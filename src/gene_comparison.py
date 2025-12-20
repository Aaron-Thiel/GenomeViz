"""
Gene comparison module for scaffold-contig gene integrity analysis.

Provides classes for:
- Analyzing gene integrity between scaffolds and contigs
- Detecting split, truncated, and missing genes
- Generating gene comparison reports

Classes:
    GeneIntegrityAnalyzer: Analyze gene coverage and integrity
    GeneComparisonResult: Data class for gene comparison results

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

            # Create gene copy with coordinates relative to scaffold region
            # (for visualization purposes, offset by ref_start)
            gene_for_viz = gene.copy()
            gene_for_viz['start'] = gene['start'] - ref_start  # Convert to scaffold-local
            gene_for_viz['end'] = gene['end'] - ref_start

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
                                              overlap_regions: List[OverlapRegion]) -> List[GeneComparisonResult]:
        """
        Analyze scaffold genes using reference coordinates to find contig coverage.

        This method uses the reference coordinate system to determine which
        contigs cover which parts of scaffold genes.

        Args:
            scaffold_name: Name of scaffold
            overlap_regions: OverlapRegion objects for this scaffold

        Returns:
            List of GeneComparisonResult for each gene
        """
        results = []

        # Get scaffold genes
        scaffold_genes = self.gff_parser.scaffold_genes.get(scaffold_name, [])

        # Handle _RC suffix: when OrientationDetector reverse-complements a scaffold,
        # it renames it with _RC suffix (e.g., 'contig_1' -> 'contig_1_RC').
        # But the scaffold GFF still has the original name, and more importantly,
        # the GFF coordinates don't match the RC'd sequence coordinates.
        # In this case, return empty to trigger fallback to reference genes.
        if not scaffold_genes and scaffold_name.endswith('_RC'):
            original_name = scaffold_name[:-3]
            if self.gff_parser.scaffold_genes.get(original_name, []):
                # Genes exist under original name but coords won't match RC'd scaffold
                print(f"        Note: '{scaffold_name}' was reverse-complemented.")
                print(f"        Scaffold GFF coords don't match - falling back to reference genes.")
                return results  # Return empty to trigger reference gene fallback

        if not scaffold_genes:
            return results

        # For each gene, find which contigs cover it via reference coordinates
        for gene in scaffold_genes:
            gene_start = gene['start']
            gene_end = gene['end']
            gene_length = gene_end - gene_start

            # Find all contigs that overlap this gene region
            covering_contigs = []
            covered_positions = np.zeros(gene_length, dtype=bool)
            identity_sum = 0.0
            identity_count = 0

            for region in overlap_regions:
                # Map gene coordinates to reference coordinates
                # This assumes scaffold local coords map to reference coords
                # via the alignment in the region
                scaffold_to_ref_offset = region.scaffold_ref_start - region.scaffold_local_start

                gene_ref_start = gene_start + scaffold_to_ref_offset
                gene_ref_end = gene_end + scaffold_to_ref_offset

                for contig in region.overlapping_contigs:
                    contig_ref_start = contig['contig_ref_start']
                    contig_ref_end = contig['contig_ref_end']

                    # Check overlap with gene in reference coordinates
                    overlap_start = max(gene_ref_start, contig_ref_start)
                    overlap_end = min(gene_ref_end, contig_ref_end)

                    if overlap_end > overlap_start:
                        # Map back to gene-local coordinates
                        local_start = max(0, overlap_start - gene_ref_start)
                        local_end = min(gene_length, overlap_end - gene_ref_start)

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

            results.append(GeneComparisonResult(
                gene=gene,
                source='scaffold',
                status=status,
                coverage_pct=coverage_pct,
                num_covering_sequences=len(unique_contigs),
                covering_sequences=covering_contigs,
                gaps=gaps,
                identity=avg_identity,
                details=details
            ))

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
