"""
Gene comparison module for scaffold-contig gene integrity analysis.

Provides classes for:
- Analyzing gene integrity between scaffolds and contigs
- Detecting split, truncated, and missing genes
- Generating gene comparison reports
- Comparing gene annotations between assemblies (unique/shared genes)

Classes:
    GeneIntegrityAnalyzer: Analyze gene coverage and integrity
    GeneComparisonResult: Data class for gene comparison results
    AssemblyGeneComparator: Compare genes between scaffold and contig GFFs

GitHub: https://github.com/Aaron-Thiel/GenomeViz
License: MIT
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
from collections import defaultdict
import numpy as np

from .assembly_comparison import OverlapRegion, MultiAssemblyGFFParser


@dataclass
class GeneComparisonResult:
    """
    Result of comparing a gene between scaffold and contigs.

    Attributes:
        gene: Original gene dictionary
        source: 'scaffold' or 'contig'
        status: 'complete', 'split', 'truncated', or 'missing'
        coverage_pct: Percentage of gene covered
        num_covering_sequences: Number of contigs/scaffolds covering this gene
        covering_sequences: List of sequence names that cover this gene
        gaps: List of gap regions within the gene
        identity: Average identity in covered regions
    """
    gene: dict
    source: str  # 'scaffold' or 'contig'
    status: str  # 'complete', 'split', 'truncated', 'missing'
    coverage_pct: float
    num_covering_sequences: int
    covering_sequences: List[dict] = field(default_factory=list)
    gaps: List[Tuple[int, int]] = field(default_factory=list)
    identity: float = 0.0
    details: str = ""


class GeneIntegrityAnalyzer:
    """
    Analyze gene integrity between scaffolds and contigs.

    Determines whether genes are complete, split across multiple sequences,
    truncated, or missing entirely.
    """

    def __init__(self, overlap_regions: List[OverlapRegion],
                 gff_parser: MultiAssemblyGFFParser,
                 min_coverage_complete: float = 95.0,
                 min_coverage_partial: float = 50.0):
        """
        Initialize gene integrity analyzer.

        Args:
            overlap_regions: List of OverlapRegion objects
            gff_parser: MultiAssemblyGFFParser with gene annotations
            min_coverage_complete: Minimum coverage % to consider gene complete
            min_coverage_partial: Minimum coverage % to not consider gene missing
        """
        self.overlap_regions = overlap_regions
        self.gff_parser = gff_parser
        self.min_coverage_complete = min_coverage_complete
        self.min_coverage_partial = min_coverage_partial

        # Build coordinate maps for efficient lookup
        self._build_coverage_maps()

    def _build_coverage_maps(self):
        """Build coverage maps from overlap regions."""
        # Map: scaffold_name -> list of (ref_start, ref_end, contig_info)
        self.scaffold_to_contig_coverage = defaultdict(list)
        # Map: contig_name -> list of (ref_start, ref_end, scaffold_info)
        self.contig_to_scaffold_coverage = defaultdict(list)

        for region in self.overlap_regions:
            scaffold_name = region.scaffold_name

            for contig in region.overlapping_contigs:
                # Coverage of scaffold by this contig
                self.scaffold_to_contig_coverage[scaffold_name].append({
                    'ref_start': contig['overlap_start'],
                    'ref_end': contig['overlap_end'],
                    'contig_name': contig['contig_name'],
                    'identity': contig['identity'],
                    'strand': contig['strand']
                })

                # Coverage of contig by this scaffold
                self.contig_to_scaffold_coverage[contig['contig_name']].append({
                    'ref_start': contig['overlap_start'],
                    'ref_end': contig['overlap_end'],
                    'scaffold_name': scaffold_name,
                    'identity': contig['identity'],
                    'strand': contig['strand']
                })

    def analyze_scaffold_genes(self, scaffold_name: str) -> List[GeneComparisonResult]:
        """
        Analyze all genes in a scaffold for coverage by contigs.

        Args:
            scaffold_name: Name of scaffold to analyze

        Returns:
            List of GeneComparisonResult for each gene
        """
        results = []

        # Get genes for this scaffold
        scaffold_genes = self.gff_parser.scaffold_genes.get(scaffold_name, [])

        if not scaffold_genes:
            return results

        # Get contig coverage for this scaffold
        contig_coverages = self.scaffold_to_contig_coverage.get(scaffold_name, [])

        for gene in scaffold_genes:
            result = self._analyze_gene_coverage(
                gene, contig_coverages, 'scaffold'
            )
            results.append(result)

        return results

    def analyze_contig_genes(self, contig_name: str) -> List[GeneComparisonResult]:
        """
        Analyze all genes in a contig for coverage by scaffolds.

        Args:
            contig_name: Name of contig to analyze

        Returns:
            List of GeneComparisonResult for each gene
        """
        results = []

        # Get genes for this contig
        contig_genes = self.gff_parser.contig_genes.get(contig_name, [])

        if not contig_genes:
            return results

        # Get scaffold coverage for this contig
        scaffold_coverages = self.contig_to_scaffold_coverage.get(contig_name, [])

        for gene in contig_genes:
            result = self._analyze_gene_coverage(
                gene, scaffold_coverages, 'contig'
            )
            results.append(result)

        return results

    def analyze_reference_genes_in_range(self, ref_start: int, ref_end: int,
                                          reference_seqid: str,
                                          overlap_regions: List[OverlapRegion]) -> List[GeneComparisonResult]:
        """
        Analyze reference genes in a coordinate range based on contig coverage.

        This method is used when scaffold GFF is not available. It uses reference
        genes (from the required --gff) to show gene-level information in the
        comparison view.

        Args:
            ref_start: Start position on reference
            ref_end: End position on reference
            reference_seqid: Reference sequence ID
            overlap_regions: OverlapRegion objects for this scaffold

        Returns:
            List of GeneComparisonResult for each reference gene in range
        """
        results = []

        # Get reference genes in the specified range
        reference_genes = self.gff_parser.get_genes_in_range(
            'reference', reference_seqid, ref_start, ref_end
        )

        if not reference_genes:
            return results

        # For each reference gene, analyze contig coverage
        for gene in reference_genes:
            gene_start = gene['start']
            gene_end = gene['end']
            gene_length = gene_end - gene_start

            if gene_length <= 0:
                continue

            # Find all contigs that overlap this gene region
            covering_contigs = []
            covered_positions = np.zeros(gene_length, dtype=bool)
            identity_sum = 0.0
            identity_count = 0

            for region in overlap_regions:
                for contig in region.overlapping_contigs:
                    contig_ref_start = contig['contig_ref_start']
                    contig_ref_end = contig['contig_ref_end']

                    # Check overlap with gene in reference coordinates
                    overlap_start = max(gene_start, contig_ref_start)
                    overlap_end = min(gene_end, contig_ref_end)

                    if overlap_end > overlap_start:
                        # Map to gene-local coordinates
                        local_start = max(0, overlap_start - gene_start)
                        local_end = min(gene_length, overlap_end - gene_start)

                        covered_positions[local_start:local_end] = True

                        covering_contigs.append({
                            'name': contig['contig_name'],
                            'overlap_start': overlap_start,
                            'overlap_end': overlap_end,
                            'overlap_bp': overlap_end - overlap_start,
                            'identity': contig['identity'],
                            'strand': contig['strand'],
                            'gene_local_start': local_start,
                            'gene_local_end': local_end
                        })

                        identity_sum += contig['identity'] * (local_end - local_start)
                        identity_count += (local_end - local_start)

            # Calculate coverage statistics
            covered_bp = np.sum(covered_positions)
            coverage_pct = (covered_bp / gene_length) * 100 if gene_length > 0 else 0

            # Calculate average identity
            avg_identity = identity_sum / identity_count if identity_count > 0 else 0

            # Find gaps in coverage
            gaps = self._find_gaps(covered_positions, gene_start)

            # Determine unique covering contigs
            unique_contigs = list(set(c['name'] for c in covering_contigs))

            # Determine status
            status, details = self._determine_status(
                coverage_pct, len(unique_contigs), gaps, gene_length
            )

            # Create gene copy - keep original reference coordinates
            # since we're using reference genes as fallback and they should
            # be displayed at their actual reference positions
            gene_for_viz = gene.copy()
            # Don't convert coordinates - keep them in reference coordinate space

            results.append(GeneComparisonResult(
                gene=gene_for_viz,
                source='reference',
                status=status,
                coverage_pct=coverage_pct,
                num_covering_sequences=len(unique_contigs),
                covering_sequences=covering_contigs,
                gaps=gaps,
                identity=avg_identity,
                details=details
            ))

        return results

    def analyze_scaffold_genes_via_reference(self, scaffold_name: str,
                                              overlap_regions: List[OverlapRegion],
                                              scaffold_length: int = None) -> List[GeneComparisonResult]:
        """
        Analyze scaffold genes using reference coordinates to find contig coverage.

        This method uses the reference coordinate system to determine which
        contigs cover which parts of scaffold genes.

        When a scaffold has been reverse-complemented (indicated by _RC suffix),
        this method transforms the GFF coordinates to match the RC'd scaffold
        and then maps them to reference coordinates.

        Args:
            scaffold_name: Name of scaffold (may have _RC suffix if reverse-complemented)
            overlap_regions: OverlapRegion objects for this scaffold
            scaffold_length: Length of the scaffold sequence (needed for RC coordinate transform)

        Returns:
            List of GeneComparisonResult for each gene
        """
        results = []
        is_rc = scaffold_name.endswith('_RC')
        original_name = scaffold_name[:-3] if is_rc else scaffold_name

        # Get scaffold genes using original name (before _RC suffix was added)
        scaffold_genes = self.gff_parser.scaffold_genes.get(scaffold_name, [])
        
        if not scaffold_genes and is_rc:
            # Try original name without _RC suffix
            scaffold_genes = self.gff_parser.scaffold_genes.get(original_name, [])
            
            if scaffold_genes:
                print(f"        Found {len(scaffold_genes)} genes under original name '{original_name}'")
                
                # We need the scaffold length to transform coordinates
                # Try to get it from overlap_regions or estimate from alignment span
                if scaffold_length is None:
                    # Estimate from the maximum scaffold_local_end in overlap regions
                    max_local = 0
                    for region in overlap_regions:
                        if region.scaffold_local_end > max_local:
                            max_local = region.scaffold_local_end
                    scaffold_length = max_local
                    print(f"        Estimated scaffold length: {scaffold_length:,} bp")
                else:
                    print(f"        Using provided scaffold length: {scaffold_length:,} bp")
                
                # Debug: print overlap regions info
                print(f"        Scaffold has {len(overlap_regions)} alignment regions:")
                for i, region in enumerate(overlap_regions[:3]):  # Show first 3
                    print(f"          Region {i+1}: scaffold_local {region.scaffold_local_start:,}-{region.scaffold_local_end:,} "
                          f"-> ref {region.scaffold_ref_start:,}-{region.scaffold_ref_end:,}")
                if len(overlap_regions) > 3:
                    print(f"          ... and {len(overlap_regions) - 3} more regions")
                
                # Transform gene coordinates for reverse-complemented scaffold
                # When a sequence is RC'd: new_pos = length - old_pos + 1
                print(f"        Transforming {len(scaffold_genes)} gene coordinates for RC'd scaffold...")
                transformed_genes = []
                transform_errors = 0
                for gene in scaffold_genes:
                    # Transform coordinates
                    new_start = scaffold_length - gene['end'] + 1
                    new_end = scaffold_length - gene['start'] + 1
                    
                    # Sanity check: transformed gene should have same length
                    orig_length = gene['end'] - gene['start']
                    new_length = new_end - new_start
                    if orig_length != new_length:
                        transform_errors += 1
                        if transform_errors <= 3:
                            print(f"          ERROR: Length mismatch! orig={orig_length}, new={new_length}")
                    
                    # Flip strand
                    new_strand = '-' if gene.get('strand', '+') == '+' else '+'
                    
                    transformed_gene = gene.copy()
                    transformed_gene['start'] = new_start
                    transformed_gene['end'] = new_end
                    transformed_gene['strand'] = new_strand
                    transformed_gene['original_start'] = gene['start']
                    transformed_gene['original_end'] = gene['end']
                    transformed_gene['was_rc_transformed'] = True
                    transformed_genes.append(transformed_gene)
                
                if transform_errors > 0:
                    print(f"        WARNING: {transform_errors} genes had length mismatch after transformation!")
                
                scaffold_genes = sorted(transformed_genes, key=lambda x: x['start'])
                
                # Debug: show first few transformed genes
                print(f"        First 3 transformed genes:")
                for gene in scaffold_genes[:3]:
                    print(f"          {gene['name'][:40]}: orig {gene['original_start']:,}-{gene['original_end']:,} "
                          f"-> transformed {gene['start']:,}-{gene['end']:,}")

        if not scaffold_genes:
            return results

        # Track statistics for debugging
        total_genes = len(scaffold_genes)
        mapped_genes = 0
        skipped_unmapped = 0
        skipped_inconsistent = 0

        # For each gene, find which contigs cover it via reference coordinates
        for gene in scaffold_genes:
            gene_start = gene['start']
            gene_end = gene['end']
            gene_length = gene_end - gene_start

            if gene_length <= 0:
                continue

            # Find all contigs that overlap this gene region
            covering_contigs = []
            covered_positions = np.zeros(gene_length, dtype=bool)
            identity_sum = 0.0
            identity_count = 0
            
            # Track the reference position of this gene for visualization
            gene_ref_start = None
            gene_ref_end = None
            best_region_overlap = 0  # Track which region has best overlap with gene

            for region in overlap_regions:
                # Check if this gene falls within this scaffold alignment region
                # Gene is in scaffold-local coordinates, region has scaffold_local_start/end
                if gene_end < region.scaffold_local_start or gene_start > region.scaffold_local_end:
                    continue  # Gene not in this alignment region
                
                # Calculate how much of the gene falls within this region
                overlap_with_region_start = max(gene_start, region.scaffold_local_start)
                overlap_with_region_end = min(gene_end, region.scaffold_local_end)
                overlap_with_region_bp = overlap_with_region_end - overlap_with_region_start
                
                # Map gene coordinates to reference coordinates using the alignment
                # scaffold_local -> reference mapping
                scaffold_to_ref_offset = region.scaffold_ref_start - region.scaffold_local_start
                
                gene_ref_start_this = gene_start + scaffold_to_ref_offset
                gene_ref_end_this = gene_end + scaffold_to_ref_offset
                
                # Debug: check for suspiciously large spans
                ref_span = gene_ref_end_this - gene_ref_start_this
                if abs(ref_span - gene_length) > 1000:
                    print(f"        DEBUG: Gene '{gene.get('name', 'unknown')[:30]}' mapping issue:")
                    print(f"          gene_start={gene_start:,}, gene_end={gene_end:,}, gene_length={gene_length}")
                    print(f"          region: scaffold_local {region.scaffold_local_start:,}-{region.scaffold_local_end:,}")
                    print(f"          region: scaffold_ref {region.scaffold_ref_start:,}-{region.scaffold_ref_end:,}")
                    print(f"          offset={scaffold_to_ref_offset:,}")
                    print(f"          gene_ref={gene_ref_start_this:,}-{gene_ref_end_this:,}, span={ref_span:,}")
                
                # Use the region with the MOST overlap with the gene
                # (Don't accumulate min/max across regions - that causes huge coordinate spans)
                if gene_ref_start is None or overlap_with_region_bp > best_region_overlap:
                    gene_ref_start = gene_ref_start_this
                    gene_ref_end = gene_ref_end_this
                    best_region_overlap = overlap_with_region_bp

                for contig in region.overlapping_contigs:
                    contig_ref_start = contig['contig_ref_start']
                    contig_ref_end = contig['contig_ref_end']

                    # Check overlap with gene in reference coordinates
                    overlap_start = max(gene_ref_start_this, contig_ref_start)
                    overlap_end = min(gene_ref_end_this, contig_ref_end)

                    if overlap_end > overlap_start:
                        # Map back to gene-local coordinates
                        local_start = max(0, overlap_start - gene_ref_start_this)
                        local_end = min(gene_length, overlap_end - gene_ref_start_this)

                        covered_positions[local_start:local_end] = True

                        covering_contigs.append({
                            'name': contig['contig_name'],
                            'overlap_start': overlap_start,
                            'overlap_end': overlap_end,
                            'overlap_bp': overlap_end - overlap_start,
                            'identity': contig['identity'],
                            'strand': contig['strand'],
                            'gene_local_start': local_start,
                            'gene_local_end': local_end
                        })

                        identity_sum += contig['identity'] * (local_end - local_start)
                        identity_count += (local_end - local_start)

            # Skip genes that couldn't be mapped to reference coordinates
            if gene_ref_start is None:
                skipped_unmapped += 1
                continue  # Gene doesn't fall within any aligned region
            
            # Sanity check: gene span on reference should match gene length (approximately)
            ref_span = gene_ref_end - gene_ref_start
            if abs(ref_span - gene_length) > 100:  # Allow small differences due to indels
                print(f"        Warning: Gene '{gene.get('name', 'unknown')}' has inconsistent span: "
                      f"gene_length={gene_length}, ref_span={ref_span}. Skipping.")
                skipped_inconsistent += 1
                continue
            
            mapped_genes += 1

            # Calculate coverage statistics
            covered_bp = np.sum(covered_positions)
            coverage_pct = (covered_bp / gene_length) * 100 if gene_length > 0 else 0

            # Calculate average identity
            avg_identity = identity_sum / identity_count if identity_count > 0 else 0

            # Find gaps in coverage
            gaps = self._find_gaps(covered_positions, gene_start)

            # Determine unique covering contigs
            unique_contigs = list(set(c['name'] for c in covering_contigs))

            # Determine status
            status, details = self._determine_status(
                coverage_pct, len(unique_contigs), gaps, gene_length
            )

            # Create gene copy with reference coordinates for visualization
            gene_for_viz = gene.copy()
            if gene_ref_start is not None:
                gene_for_viz['ref_start'] = gene_ref_start
                gene_for_viz['ref_end'] = gene_ref_end

            results.append(GeneComparisonResult(
                gene=gene_for_viz,
                source='scaffold',
                status=status,
                coverage_pct=coverage_pct,
                num_covering_sequences=len(unique_contigs),
                covering_sequences=covering_contigs,
                gaps=gaps,
                identity=avg_identity,
                details=details
            ))

        # Print summary
        if skipped_unmapped > 0 or skipped_inconsistent > 0:
            print(f"        Gene mapping summary: {mapped_genes}/{total_genes} mapped, "
                  f"{skipped_unmapped} outside aligned regions, {skipped_inconsistent} inconsistent spans")

        return results

    def _analyze_gene_coverage(self, gene: dict, coverages: List[dict],
                               source: str) -> GeneComparisonResult:
        """
        Analyze coverage of a single gene.

        Args:
            gene: Gene dictionary with start, end, name, etc.
            coverages: List of coverage regions from other assembly
            source: 'scaffold' or 'contig'

        Returns:
            GeneComparisonResult with coverage analysis
        """
        gene_start = gene['start']
        gene_end = gene['end']
        gene_length = gene_end - gene_start

        # Track covered positions
        covered_positions = np.zeros(gene_length, dtype=bool)
        covering_sequences = []
        identity_sum = 0.0
        identity_count = 0

        for cov in coverages:
            cov_start = cov['ref_start']
            cov_end = cov['ref_end']

            # Check overlap
            overlap_start = max(gene_start, cov_start)
            overlap_end = min(gene_end, cov_end)

            if overlap_end > overlap_start:
                # Map to gene-local coordinates
                local_start = overlap_start - gene_start
                local_end = overlap_end - gene_start

                covered_positions[local_start:local_end] = True

                seq_name = cov.get('contig_name', cov.get('scaffold_name', 'unknown'))
                covering_sequences.append({
                    'name': seq_name,
                    'overlap_start': overlap_start,
                    'overlap_end': overlap_end,
                    'overlap_bp': overlap_end - overlap_start,
                    'identity': cov.get('identity', 0),
                    'strand': cov.get('strand', '+'),
                    'gene_local_start': local_start,
                    'gene_local_end': local_end
                })

                identity_sum += cov.get('identity', 0) * (local_end - local_start)
                identity_count += (local_end - local_start)

        # Calculate statistics
        covered_bp = np.sum(covered_positions)
        coverage_pct = (covered_bp / gene_length) * 100 if gene_length > 0 else 0
        avg_identity = identity_sum / identity_count if identity_count > 0 else 0

        # Find gaps
        gaps = self._find_gaps(covered_positions, gene_start)

        # Unique covering sequences
        unique_seqs = list(set(c['name'] for c in covering_sequences))

        # Determine status
        status, details = self._determine_status(
            coverage_pct, len(unique_seqs), gaps, gene_length
        )

        return GeneComparisonResult(
            gene=gene,
            source=source,
            status=status,
            coverage_pct=coverage_pct,
            num_covering_sequences=len(unique_seqs),
            covering_sequences=covering_sequences,
            gaps=gaps,
            identity=avg_identity,
            details=details
        )

    def _find_gaps(self, covered_positions: np.ndarray,
                   gene_start: int) -> List[Tuple[int, int]]:
        """Find gaps in coverage."""
        gaps = []
        in_gap = False
        gap_start = 0

        for i, covered in enumerate(covered_positions):
            if not covered and not in_gap:
                # Start of gap
                in_gap = True
                gap_start = i
            elif covered and in_gap:
                # End of gap
                gaps.append((gene_start + gap_start, gene_start + i))
                in_gap = False

        # Handle gap at end
        if in_gap:
            gaps.append((gene_start + gap_start, gene_start + len(covered_positions)))

        return gaps

    def _determine_status(self, coverage_pct: float, num_sequences: int,
                          gaps: List[Tuple[int, int]], gene_length: int) -> Tuple[str, str]:
        """
        Determine gene status based on coverage analysis.

        Returns:
            Tuple of (status, details_string)
        """
        if coverage_pct >= self.min_coverage_complete:
            if num_sequences == 1:
                return 'complete', f"Fully covered by single sequence ({coverage_pct:.1f}%)"
            else:
                return 'split', f"Covered by {num_sequences} sequences ({coverage_pct:.1f}%)"

        elif coverage_pct >= self.min_coverage_partial:
            if num_sequences == 0:
                return 'truncated', f"Partial coverage ({coverage_pct:.1f}%), no aligned sequences"
            elif num_sequences == 1:
                gap_info = f", {len(gaps)} gap(s)" if gaps else ""
                return 'truncated', f"Partial coverage by 1 sequence ({coverage_pct:.1f}%){gap_info}"
            else:
                gap_info = f", {len(gaps)} gap(s)" if gaps else ""
                return 'split', f"Split across {num_sequences} sequences ({coverage_pct:.1f}%){gap_info}"

        else:
            if coverage_pct > 0:
                return 'missing', f"Very low coverage ({coverage_pct:.1f}%)"
            else:
                return 'missing', "No coverage"

    def compute_summary(self, scaffold_results: Dict[str, List[GeneComparisonResult]]) -> dict:
        """
        Compute summary statistics across all scaffolds.

        Args:
            scaffold_results: Dict mapping scaffold_name -> list of GeneComparisonResult

        Returns:
            Summary statistics dictionary
        """
        total_genes = 0
        status_counts = {'complete': 0, 'split': 0, 'truncated': 0, 'missing': 0}
        coverage_values = []
        identity_values = []

        for scaffold_name, results in scaffold_results.items():
            for result in results:
                total_genes += 1
                status_counts[result.status] += 1
                coverage_values.append(result.coverage_pct)
                if result.identity > 0:
                    identity_values.append(result.identity)

        return {
            'total_genes': total_genes,
            'status_counts': status_counts,
            'status_percentages': {
                status: (count / total_genes * 100) if total_genes > 0 else 0
                for status, count in status_counts.items()
            },
            'avg_coverage': np.mean(coverage_values) if coverage_values else 0,
            'avg_identity': np.mean(identity_values) if identity_values else 0,
            'median_coverage': np.median(coverage_values) if coverage_values else 0
        }


@dataclass
class UniqueGene:
    """Represents a gene unique to one assembly."""
    gene_id: str
    seqid: str
    start: int
    end: int
    strand: str
    product: str
    name: str
    feature_type: str
    length: int
    source: str  # 'scaffold' or 'contig'
    ref_start: Optional[int] = None
    ref_end: Optional[int] = None
    ref_seqid: Optional[str] = None
    alignment_file: Optional[str] = None  # Path to gene_alignment HTML if exists


class AssemblyGeneComparator:
    """
    Compare gene annotations between scaffold and contig GFF files.

    Identifies genes that are unique to each assembly, helping to
    understand what genes are gained or lost during scaffolding.
    """

    def __init__(self, gff_parser: MultiAssemblyGFFParser,
                 scaffold_alignments: Dict = None,
                 contig_alignments: Dict = None,
                 assembly_comparator=None,
                 rc_scaffold_info: Dict[str, int] = None):
        """
        Initialize assembly gene comparator.

        Args:
            gff_parser: MultiAssemblyGFFParser with parsed GFFs
            scaffold_alignments: Scaffold alignment data for coordinate mapping
            contig_alignments: Contig alignment data for coordinate mapping
            assembly_comparator: Optional AssemblyComparator for accurate offset calculation
            rc_scaffold_info: Dict mapping RC'd scaffold names (with _RC suffix) to their lengths.
                              Used to transform gene coordinates for RC'd scaffolds.
        """
        self.gff_parser = gff_parser
        self.scaffold_alignments = scaffold_alignments or {}
        self.contig_alignments = contig_alignments or {}
        self.assembly_comparator = assembly_comparator
        self.rc_scaffold_info = rc_scaffold_info or {}

        # Build offset maps for each sequence
        # Use overlap regions if comparator is available (matches visualization exactly)
        # Otherwise fall back to alignment-based offsets
        if assembly_comparator:
            self._scaffold_offsets = self._build_offset_map_from_comparator(is_scaffold=True)
            self._contig_offsets = self._build_offset_map_from_comparator(is_scaffold=False)
        else:
            self._scaffold_offsets = self._build_offset_map(self.scaffold_alignments)
            self._contig_offsets = self._build_offset_map(self.contig_alignments)

        # Results
        self.unique_to_scaffolds: List[UniqueGene] = []
        self.unique_to_contigs: List[UniqueGene] = []
        self.shared_genes: List[dict] = []

    def _build_offset_map_from_comparator(self, is_scaffold: bool) -> Dict[str, List]:
        """
        Build offset map using the same overlap regions as the visualization.

        This ensures ref_start/ref_end values in the TSV match exactly where
        genes appear in the comparison visualization.

        Args:
            is_scaffold: True to build map for scaffolds, False for contigs

        Returns:
            Dict mapping seqid to list of overlap regions (for per-gene offset calculation)
        """
        offset_map = {}

        if is_scaffold:
            # For scaffolds, store ALL overlap regions (not just min)
            # This allows per-gene offset calculation like the visualization does
            for scaffold_name in self.assembly_comparator.get_all_scaffolds():
                overlap_regions = self.assembly_comparator.find_overlapping_contigs(scaffold_name)
                if overlap_regions:
                    # Store the regions themselves for per-gene mapping
                    offset_map[scaffold_name] = overlap_regions
                    if scaffold_name.endswith('_RC'):
                        offset_map[scaffold_name[:-3]] = overlap_regions
        else:
            # For contigs, we need to find the min ref_start from contig alignments
            # that overlap with scaffolds
            contig_min_refs = {}
            for scaffold_name in self.assembly_comparator.get_all_scaffolds():
                overlap_regions = self.assembly_comparator.find_overlapping_contigs(scaffold_name)
                for region in overlap_regions:
                    ref_seqid = region.reference_seqid
                    for contig in region.overlapping_contigs:
                        contig_name = contig['contig_name']
                        contig_ref_start = contig['contig_ref_start']

                        if contig_name not in contig_min_refs or contig_ref_start < contig_min_refs[contig_name][0]:
                            contig_min_refs[contig_name] = (contig_ref_start, 0, ref_seqid)

            for contig_name, info in contig_min_refs.items():
                offset_map[contig_name] = info
                if contig_name.endswith('_RC'):
                    offset_map[contig_name[:-3]] = info

        return offset_map

    def _build_offset_map(self, alignments: Dict) -> Dict[str, Tuple[int, int, str]]:
        """
        Build a map of sequence -> (min_ref_start, corresponding_query_start, ref_seqid).

        Fallback method when AssemblyComparator is not provided.
        Only considers primary alignments to match visualization behavior.

        Returns:
            Dict mapping seqid (with/without _RC) to (min_ref_start, query_start, ref_seqid)
        """
        offset_map = {}

        for ref_seqid, aln_list in alignments.items():
            for aln in aln_list:
                # Only use primary alignments like the visualization does
                if not aln.get('is_primary', True):
                    continue

                query_name = aln.get('query_name', '')
                ref_start = aln.get('ref_start', 0)
                query_start = aln.get('query_start', 0)

                # Track by both original name and with _RC suffix
                for name in [query_name, query_name.replace('_RC', '')]:
                    if name and (name not in offset_map or ref_start < offset_map[name][0]):
                        offset_map[name] = (ref_start, query_start, ref_seqid)

        return offset_map

    def _get_gene_signature(self, gene: dict) -> set:
        """Get set of identifiers for matching genes."""
        signatures = set()

        product = gene.get('product', '')
        if product and 'hypothetical' not in product.lower():
            # Normalize product name
            signatures.add(product.lower().strip())

        name = gene.get('name', '')
        if name:
            signatures.add(name.lower().strip())

        return signatures

    def _map_gene_to_reference(self, gene: dict, seqid: str,
                                alignments: Dict,
                                offset_map: Dict,
                                is_scaffold: bool = False) -> Tuple[Optional[int], Optional[int], Optional[str], int, int, str]:
        """
        Map gene local coordinates to reference coordinates.

        For scaffolds with overlap regions available, uses per-region offset calculation
        to match exactly how the visualization positions genes.

        For RC'd scaffolds, first transforms the gene coordinates to match
        the RC'd scaffold, then applies the appropriate offset.

        Args:
            gene: Gene dictionary with start/end positions
            seqid: Sequence ID (scaffold or contig name)
            alignments: Alignment data (kept for compatibility, not used)
            offset_map: Precomputed offset map from _build_offset_map
            is_scaffold: If True, check for RC transformation and use per-region offsets

        Returns:
            Tuple of (ref_start, ref_end, ref_seqid, local_start, local_end, effective_seqid)
            where local_start/local_end may be RC-transformed and effective_seqid has _RC suffix if applicable
        """
        gene_start = gene['start']
        gene_end = gene['end']
        local_start = gene_start
        local_end = gene_end
        effective_seqid = seqid

        # Check if this scaffold was RC'd and transform coordinates
        if is_scaffold:
            # Check if there's an RC'd version of this scaffold
            rc_name = f"{seqid}_RC"
            if rc_name in self.rc_scaffold_info:
                scaffold_length = self.rc_scaffold_info[rc_name]
                # Transform coordinates: new_pos = length - old_pos + 1
                local_start = scaffold_length - gene_end + 1
                local_end = scaffold_length - gene_start + 1
                # Use the RC name for output and offset lookup
                effective_seqid = rc_name

        # Look up offset info for this sequence
        offset_info = offset_map.get(effective_seqid) or offset_map.get(f"{seqid}_RC") or offset_map.get(seqid.replace('_RC', ''))

        if offset_info is None:
            return None, None, None, local_start, local_end, effective_seqid

        # Check if offset_info is a list of OverlapRegions (from comparator) or a tuple (from alignments)
        if isinstance(offset_info, list):
            # Per-region offset calculation (matches visualization exactly)
            # Find which overlap region contains this gene
            gene_ref_start = None
            gene_ref_end = None
            ref_seqid = None
            best_overlap = 0

            for region in offset_info:
                # Check if gene falls within this region's scaffold local coordinates
                if local_end < region.scaffold_local_start or local_start > region.scaffold_local_end:
                    continue  # Gene not in this region

                # Calculate overlap with this region
                overlap_start = max(local_start, region.scaffold_local_start)
                overlap_end = min(local_end, region.scaffold_local_end)
                overlap_bp = overlap_end - overlap_start

                if overlap_bp > best_overlap:
                    # Use this region's offset (same calculation as visualization)
                    offset = region.scaffold_ref_start - region.scaffold_local_start
                    gene_ref_start = local_start + offset
                    gene_ref_end = local_end + offset
                    ref_seqid = region.reference_seqid
                    best_overlap = overlap_bp

            if gene_ref_start is None:
                # Gene doesn't fall within any aligned region
                return None, None, None, local_start, local_end, effective_seqid

            return gene_ref_start, gene_ref_end, ref_seqid, local_start, local_end, effective_seqid

        else:
            # Tuple format from alignment-based offset (fallback)
            min_ref_start, _, ref_seqid = offset_info
            ref_start = local_start + min_ref_start
            ref_end = local_end + min_ref_start
            return ref_start, ref_end, ref_seqid, local_start, local_end, effective_seqid

    def compare(self) -> Tuple[List[UniqueGene], List[UniqueGene]]:
        """
        Compare genes between scaffold and contig GFFs.

        Returns:
            Tuple of (unique_to_scaffolds, unique_to_contigs)
        """
        # Build signature sets
        scaffold_signatures = set()
        scaffold_genes_by_sig = {}

        for seqid, genes in self.gff_parser.scaffold_genes.items():
            for gene in genes:
                sigs = self._get_gene_signature(gene)
                for sig in sigs:
                    scaffold_signatures.add(sig)
                    scaffold_genes_by_sig[sig] = (gene, seqid)

        contig_signatures = set()
        contig_genes_by_sig = {}

        for seqid, genes in self.gff_parser.contig_genes.items():
            for gene in genes:
                sigs = self._get_gene_signature(gene)
                for sig in sigs:
                    contig_signatures.add(sig)
                    contig_genes_by_sig[sig] = (gene, seqid)

        # Find unique genes
        seen_scaffold = set()
        for seqid, genes in self.gff_parser.scaffold_genes.items():
            for gene in genes:
                gene_sigs = self._get_gene_signature(gene)

                if not gene_sigs:
                    continue  # Skip hypothetical proteins without other identifiers

                is_unique = not any(sig in contig_signatures for sig in gene_sigs)

                if is_unique:
                    key = (seqid, gene['start'], gene['end'])
                    if key not in seen_scaffold:
                        seen_scaffold.add(key)

                        # Map to reference coordinates (with RC transformation if needed)
                        ref_start, ref_end, ref_seqid, local_start, local_end, effective_seqid = self._map_gene_to_reference(
                            gene, seqid, self.scaffold_alignments, self._scaffold_offsets,
                            is_scaffold=True
                        )

                        self.unique_to_scaffolds.append(UniqueGene(
                            gene_id=gene.get('name', f"gene_{gene['start']}_{gene['end']}"),
                            seqid=effective_seqid,  # Use RC'd name if applicable
                            start=local_start,  # Use RC-transformed coordinates
                            end=local_end,
                            strand=gene.get('strand', '+'),
                            product=gene.get('product', ''),
                            name=gene.get('name', ''),
                            feature_type=gene.get('type', 'CDS'),
                            length=gene['end'] - gene['start'] + 1,
                            source='scaffold',
                            ref_start=ref_start,
                            ref_end=ref_end,
                            ref_seqid=ref_seqid
                        ))

        seen_contig = set()
        for seqid, genes in self.gff_parser.contig_genes.items():
            for gene in genes:
                gene_sigs = self._get_gene_signature(gene)

                if not gene_sigs:
                    continue

                is_unique = not any(sig in scaffold_signatures for sig in gene_sigs)

                if is_unique:
                    key = (seqid, gene['start'], gene['end'])
                    if key not in seen_contig:
                        seen_contig.add(key)

                        # Map to reference coordinates
                        ref_start, ref_end, ref_seqid, local_start, local_end, effective_seqid = self._map_gene_to_reference(
                            gene, seqid, self.contig_alignments, self._contig_offsets,
                            is_scaffold=False
                        )

                        self.unique_to_contigs.append(UniqueGene(
                            gene_id=gene.get('name', f"gene_{gene['start']}_{gene['end']}"),
                            seqid=effective_seqid,
                            start=local_start,
                            end=local_end,
                            strand=gene.get('strand', '+'),
                            product=gene.get('product', ''),
                            name=gene.get('name', ''),
                            feature_type=gene.get('type', 'CDS'),
                            length=gene['end'] - gene['start'] + 1,
                            source='contig',
                            ref_start=ref_start,
                            ref_end=ref_end,
                            ref_seqid=ref_seqid
                        ))

        # Sort by reference position where available
        self.unique_to_scaffolds.sort(key=lambda g: (g.ref_start or 999999999, g.seqid, g.start))
        self.unique_to_contigs.sort(key=lambda g: (g.ref_start or 999999999, g.seqid, g.start))

        return self.unique_to_scaffolds, self.unique_to_contigs

    def link_to_alignments(self, gene_alignments_dir: str, scaffold_name: str = None):
        """
        Link unique genes to existing gene_alignment HTML files.

        Args:
            gene_alignments_dir: Path to gene_alignments directory
            scaffold_name: Optional scaffold name to filter results
        """
        from pathlib import Path
        import os
        import re

        align_dir = Path(gene_alignments_dir)
        if not align_dir.exists():
            return

        # Get all alignment files
        alignment_files = list(align_dir.glob('*_alignment.html'))

        # Create lookup by gene name/product
        file_lookup = {}
        for f in alignment_files:
            # Extract gene name from filename
            # Format: GeneName_Start-End_alignment.html
            fname = f.stem.replace('_alignment', '')
            # Normalize for matching
            normalized = re.sub(r'_\d+-\d+$', '', fname).lower().replace('_', ' ')
            file_lookup[normalized] = str(f)

        # Link scaffold unique genes
        for gene in self.unique_to_scaffolds:
            gene_key = (gene.name or gene.product).lower().replace('_', ' ')
            if gene_key in file_lookup:
                gene.alignment_file = file_lookup[gene_key]

        # Link contig unique genes
        for gene in self.unique_to_contigs:
            gene_key = (gene.name or gene.product).lower().replace('_', ' ')
            if gene_key in file_lookup:
                gene.alignment_file = file_lookup[gene_key]

    def generate_html_report(self, output_path: str, gene_alignments_base: str = None):
        """
        Generate HTML report for gene comparison with links to alignments.

        Args:
            output_path: Path for output HTML file
            gene_alignments_base: Base path to gene_alignments directory (for relative links)
        """
        from pathlib import Path

        html_content = '''<!DOCTYPE html>
<html>
<head>
    <title>Assembly Gene Comparison</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; background: #f5f5f5; }
        .container { max-width: 1400px; margin: 0 auto; }
        h1 { color: #333; border-bottom: 2px solid #2196F3; padding-bottom: 10px; }
        h2 { color: #555; margin-top: 30px; }
        .summary { background: white; padding: 20px; border-radius: 8px; margin: 20px 0;
                   box-shadow: 0 2px 4px rgba(0,0,0,0.1); display: flex; justify-content: space-around; }
        .stat { text-align: center; padding: 15px; }
        .stat-value { font-size: 2.5em; font-weight: bold; }
        .stat-label { color: #666; margin-top: 5px; }
        .new-genes .stat-value { color: #4CAF50; }
        .lost-genes .stat-value { color: #f44336; }
        table { width: 100%; border-collapse: collapse; background: white; margin: 20px 0;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1); }
        th { background: #2196F3; color: white; padding: 12px; text-align: left; position: sticky; top: 0; }
        td { padding: 10px; border-bottom: 1px solid #eee; }
        tr:hover { background: #e3f2fd; }
        .position { font-family: monospace; color: #666; }
        .mapped { color: #4CAF50; }
        .unmapped { color: #999; }
        a { color: #1976D2; text-decoration: none; }
        a:hover { text-decoration: underline; }
        .view-btn { background: #2196F3; color: white; padding: 4px 8px; border-radius: 4px;
                    font-size: 0.85em; }
        .view-btn:hover { background: #1976D2; }
        .filter-bar { background: white; padding: 15px; border-radius: 8px; margin: 20px 0; }
        input[type="text"] { padding: 8px; width: 300px; border: 1px solid #ddd; border-radius: 4px; }
        .section-new { border-left: 4px solid #4CAF50; padding-left: 15px; }
        .section-lost { border-left: 4px solid #f44336; padding-left: 15px; }
    </style>
    <script>
        function filterTable(tableId, searchTerm) {
            const table = document.getElementById(tableId);
            const rows = table.getElementsByTagName('tr');
            const term = searchTerm.toLowerCase();
            for (let i = 1; i < rows.length; i++) {
                const text = rows[i].textContent.toLowerCase();
                rows[i].style.display = text.includes(term) ? '' : 'none';
            }
        }
    </script>
</head>
<body>
<div class="container">
    <h1>Assembly Gene Comparison: Scaffolds vs Contigs</h1>

    <div class="summary">
        <div class="stat new-genes">
            <div class="stat-value">''' + str(len(self.unique_to_scaffolds)) + '''</div>
            <div class="stat-label">NEW genes in scaffolds</div>
        </div>
        <div class="stat lost-genes">
            <div class="stat-value">''' + str(len(self.unique_to_contigs)) + '''</div>
            <div class="stat-label">Genes ONLY in contigs</div>
        </div>
    </div>

    <div class="section-new">
    <h2>New Genes in Scaffolds</h2>
    <p>These genes appear in the scaffold assembly but are NOT found in the original contigs.</p>
    <div class="filter-bar">
        <input type="text" placeholder="Filter genes..." onkeyup="filterTable('new-genes', this.value)">
    </div>
    <table id="new-genes">
        <tr>
            <th>#</th>
            <th>Scaffold</th>
            <th>Local Position</th>
            <th>Ref Position</th>
            <th>Strand</th>
            <th>Length</th>
            <th>Gene</th>
            <th>Product</th>
            <th>View</th>
        </tr>
'''

        for i, gene in enumerate(self.unique_to_scaffolds, 1):
            ref_pos = f'<span class="mapped">{gene.ref_seqid}:{gene.ref_start:,}-{gene.ref_end:,}</span>' if gene.ref_start else '<span class="unmapped">unmapped</span>'

            view_link = ""
            if gene.alignment_file:
                rel_path = Path(gene.alignment_file).name
                if gene_alignments_base:
                    rel_path = f"{gene_alignments_base}/{rel_path}"
                view_link = f'<a href="{rel_path}" class="view-btn">View Alignment</a>'

            html_content += f'''        <tr>
            <td>{i}</td>
            <td>{gene.seqid}</td>
            <td class="position">{gene.start:,}-{gene.end:,}</td>
            <td class="position">{ref_pos}</td>
            <td>{gene.strand}</td>
            <td>{gene.length:,}</td>
            <td>{gene.name or '-'}</td>
            <td>{gene.product or 'hypothetical protein'}</td>
            <td>{view_link}</td>
        </tr>
'''

        html_content += '''    </table>
    </div>

    <div class="section-lost">
    <h2>Genes Only in Contigs (Lost in Scaffolding)</h2>
    <p>These genes are present in contigs but NOT found in the scaffolded assembly.</p>
    <div class="filter-bar">
        <input type="text" placeholder="Filter genes..." onkeyup="filterTable('lost-genes', this.value)">
    </div>
    <table id="lost-genes">
        <tr>
            <th>#</th>
            <th>Contig</th>
            <th>Local Position</th>
            <th>Ref Position</th>
            <th>Strand</th>
            <th>Length</th>
            <th>Gene</th>
            <th>Product</th>
            <th>View</th>
        </tr>
'''

        for i, gene in enumerate(self.unique_to_contigs, 1):
            ref_pos = f'<span class="mapped">{gene.ref_seqid}:{gene.ref_start:,}-{gene.ref_end:,}</span>' if gene.ref_start else '<span class="unmapped">unmapped</span>'

            view_link = ""
            if gene.alignment_file:
                rel_path = Path(gene.alignment_file).name
                if gene_alignments_base:
                    rel_path = f"{gene_alignments_base}/{rel_path}"
                view_link = f'<a href="{rel_path}" class="view-btn">View Alignment</a>'

            html_content += f'''        <tr>
            <td>{i}</td>
            <td>{gene.seqid}</td>
            <td class="position">{gene.start:,}-{gene.end:,}</td>
            <td class="position">{ref_pos}</td>
            <td>{gene.strand}</td>
            <td>{gene.length:,}</td>
            <td>{gene.name or '-'}</td>
            <td>{gene.product or 'hypothetical protein'}</td>
            <td>{view_link}</td>
        </tr>
'''

        html_content += '''    </table>
    </div>
</div>
</body>
</html>'''

        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(html_content)

    def write_tsv(self, output_dir: str):
        """
        Write comparison results to TSV files.

        Args:
            output_dir: Directory for output files
        """
        from pathlib import Path

        out_dir = Path(output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        # Header
        header = "source\tseqid\tlocal_start\tlocal_end\tstrand\tlength\tref_seqid\tref_start\tref_end\tgene_id\tname\tproduct\talignment_file\n"

        # New genes in scaffolds
        with open(out_dir / 'new_genes_in_scaffolds.tsv', 'w', encoding='utf-8') as f:
            f.write(header)
            for gene in self.unique_to_scaffolds:
                f.write(f"{gene.source}\t{gene.seqid}\t{gene.start}\t{gene.end}\t{gene.strand}\t"
                       f"{gene.length}\t{gene.ref_seqid or ''}\t{gene.ref_start or ''}\t{gene.ref_end or ''}\t"
                       f"{gene.gene_id}\t{gene.name}\t{gene.product}\t{gene.alignment_file or ''}\n")

        # Genes only in contigs
        with open(out_dir / 'genes_only_in_contigs.tsv', 'w', encoding='utf-8') as f:
            f.write(header)
            for gene in self.unique_to_contigs:
                f.write(f"{gene.source}\t{gene.seqid}\t{gene.start}\t{gene.end}\t{gene.strand}\t"
                       f"{gene.length}\t{gene.ref_seqid or ''}\t{gene.ref_start or ''}\t{gene.ref_end or ''}\t"
                       f"{gene.gene_id}\t{gene.name}\t{gene.product}\t{gene.alignment_file or ''}\n")

    def get_summary(self) -> dict:
        """Get summary statistics."""
        return {
            'total_scaffold_genes': sum(len(genes) for genes in self.gff_parser.scaffold_genes.values()),
            'total_contig_genes': sum(len(genes) for genes in self.gff_parser.contig_genes.values()),
            'unique_to_scaffolds': len(self.unique_to_scaffolds),
            'unique_to_contigs': len(self.unique_to_contigs),
            'scaffold_mapped': sum(1 for g in self.unique_to_scaffolds if g.ref_start is not None),
            'contig_mapped': sum(1 for g in self.unique_to_contigs if g.ref_start is not None)
        }
