"""
Assembly comparison module for scaffold-contig analysis.

Provides classes for:
- Detecting overlapping regions between assemblies via reference coordinates
- Computing alignment statistics between scaffolds and contigs
- Managing multi-assembly GFF annotations

Classes:
    AssemblyComparator: Compare two assemblies using reference as anchor
    OverlapRegion: Data class for overlapping scaffold-contig regions
    MultiAssemblyGFFParser: Parse and organize GFFs for multiple assemblies

GitHub: https://github.com/Aaron-Thiel/GenomeViz
License: MIT
"""

from dataclasses import dataclass, field
from collections import defaultdict
from typing import Dict, List, Optional, Tuple
import numpy as np


@dataclass
class OverlapRegion:
    """
    Represents an overlapping region between a scaffold and contigs.

    Attributes:
        scaffold_name: Name of the scaffold
        scaffold_ref_start: Scaffold's start position on reference
        scaffold_ref_end: Scaffold's end position on reference
        scaffold_local_start: Start position within scaffold sequence
        scaffold_local_end: End position within scaffold sequence
        scaffold_strand: Strand of scaffold alignment
        scaffold_identity: Identity of scaffold to reference
        reference_seqid: Reference sequence ID where overlap occurs
        overlapping_contigs: List of contigs that overlap this region
    """
    scaffold_name: str
    scaffold_ref_start: int
    scaffold_ref_end: int
    scaffold_local_start: int
    scaffold_local_end: int
    scaffold_strand: str
    scaffold_identity: float
    reference_seqid: str
    overlapping_contigs: List[dict] = field(default_factory=list)


class AssemblyComparator:
    """
    Compare two assemblies (scaffolds vs contigs) using reference coordinates.

    Both assemblies are aligned to the reference genome, and overlapping
    regions are identified based on reference coordinate overlap.
    """

    def __init__(self, scaffold_alignments: Dict, contig_alignments: Dict,
                 min_overlap_bp: int = 100, min_overlap_pct: float = 0.0):
        """
        Initialize assembly comparator.

        Args:
            scaffold_alignments: Dict mapping ref_seqid -> list of scaffold alignments
            contig_alignments: Dict mapping ref_seqid -> list of contig alignments
            min_overlap_bp: Minimum overlap in base pairs to report
            min_overlap_pct: Minimum overlap as percentage of smaller region
        """
        self.scaffold_alignments = scaffold_alignments
        self.contig_alignments = contig_alignments
        self.min_overlap_bp = min_overlap_bp
        self.min_overlap_pct = min_overlap_pct

        # Index alignments by scaffold name for efficient lookup
        self._scaffold_index = self._build_scaffold_index()
        self._contig_index = self._build_contig_index()

    def _build_scaffold_index(self) -> Dict[str, List[dict]]:
        """Build index of alignments by scaffold name."""
        index = defaultdict(list)
        for ref_seqid, alignments in self.scaffold_alignments.items():
            for aln in alignments:
                if aln.get('is_primary', True):
                    aln_copy = aln.copy()
                    aln_copy['reference_seqid'] = ref_seqid
                    index[aln['query_name']].append(aln_copy)
        return dict(index)

    def _build_contig_index(self) -> Dict[str, List[dict]]:
        """Build index of alignments by contig name."""
        index = defaultdict(list)
        for ref_seqid, alignments in self.contig_alignments.items():
            for aln in alignments:
                if aln.get('is_primary', True):
                    aln_copy = aln.copy()
                    aln_copy['reference_seqid'] = ref_seqid
                    index[aln['query_name']].append(aln_copy)
        return dict(index)

    def find_overlapping_contigs(self, scaffold_name: str) -> List[OverlapRegion]:
        """
        Find all contigs that overlap with a scaffold on the reference.

        Args:
            scaffold_name: Name of the scaffold to analyze

        Returns:
            List of OverlapRegion objects describing each overlap
        """
        if scaffold_name not in self._scaffold_index:
            return []

        overlap_regions = []

        for scaffold_aln in self._scaffold_index[scaffold_name]:
            ref_seqid = scaffold_aln['reference_seqid']
            scaffold_ref_start = scaffold_aln['ref_start']
            scaffold_ref_end = scaffold_aln['ref_end']

            # Find overlapping contigs in this reference region
            overlapping_contigs = []

            if ref_seqid in self.contig_alignments:
                for contig_aln in self.contig_alignments[ref_seqid]:
                    if not contig_aln.get('is_primary', True):
                        continue

                    contig_ref_start = contig_aln['ref_start']
                    contig_ref_end = contig_aln['ref_end']

                    # Calculate overlap
                    overlap_start = max(scaffold_ref_start, contig_ref_start)
                    overlap_end = min(scaffold_ref_end, contig_ref_end)
                    overlap_bp = overlap_end - overlap_start

                    if overlap_bp >= self.min_overlap_bp:
                        # Calculate overlap percentage
                        scaffold_len = scaffold_ref_end - scaffold_ref_start
                        contig_len = contig_ref_end - contig_ref_start
                        smaller_len = min(scaffold_len, contig_len)
                        overlap_pct = (overlap_bp / smaller_len) * 100 if smaller_len > 0 else 0

                        if overlap_pct >= self.min_overlap_pct:
                            overlapping_contigs.append({
                                'contig_name': contig_aln['query_name'],
                                'contig_ref_start': contig_ref_start,
                                'contig_ref_end': contig_ref_end,
                                'contig_local_start': contig_aln['query_start'],
                                'contig_local_end': contig_aln['query_end'],
                                'contig_length': contig_aln.get('query_length', contig_ref_end - contig_ref_start),
                                'overlap_start': overlap_start,
                                'overlap_end': overlap_end,
                                'overlap_bp': overlap_bp,
                                'overlap_pct': overlap_pct,
                                'strand': contig_aln.get('strand', '+'),
                                'identity': contig_aln.get('identity', 0)
                            })

            if overlapping_contigs:
                overlap_regions.append(OverlapRegion(
                    scaffold_name=scaffold_name,
                    scaffold_ref_start=scaffold_ref_start,
                    scaffold_ref_end=scaffold_ref_end,
                    scaffold_local_start=scaffold_aln['query_start'],
                    scaffold_local_end=scaffold_aln['query_end'],
                    scaffold_strand=scaffold_aln.get('strand', '+'),
                    scaffold_identity=scaffold_aln.get('identity', 0),
                    reference_seqid=ref_seqid,
                    overlapping_contigs=sorted(overlapping_contigs,
                                              key=lambda x: x['contig_ref_start'])
                ))

        # Sort regions by reference position
        overlap_regions.sort(key=lambda x: x.scaffold_ref_start)
        return overlap_regions

    def get_all_scaffolds(self) -> List[str]:
        """Get list of all scaffold names."""
        return sorted(self._scaffold_index.keys())

    def get_all_contigs(self) -> List[str]:
        """Get list of all contig names."""
        return sorted(self._contig_index.keys())

    def compute_comparison_summary(self) -> dict:
        """
        Compute summary statistics for all scaffold-contig overlaps.

        Returns:
            Dictionary with summary statistics
        """
        scaffold_stats = {}
        total_overlap_regions = 0

        for scaffold_name in self.get_all_scaffolds():
            regions = self.find_overlapping_contigs(scaffold_name)

            total_contigs = set()
            total_overlap_bp = 0
            scaffold_coverage_bp = 0

            for region in regions:
                scaffold_coverage_bp += region.scaffold_ref_end - region.scaffold_ref_start
                for contig in region.overlapping_contigs:
                    total_contigs.add(contig['contig_name'])
                    total_overlap_bp += contig['overlap_bp']

            scaffold_stats[scaffold_name] = {
                'num_regions': len(regions),
                'num_overlapping_contigs': len(total_contigs),
                'total_overlap_bp': total_overlap_bp,
                'scaffold_coverage_bp': scaffold_coverage_bp,
                'contig_names': sorted(total_contigs),
                'regions': regions  # Store for later use
            }

            total_overlap_regions += len(regions)

        return {
            'num_scaffolds': len(scaffold_stats),
            'num_contigs': len(self.get_all_contigs()),
            'num_total_overlap_regions': total_overlap_regions,
            'scaffold_stats': scaffold_stats
        }


class MultiAssemblyGFFParser:
    """
    Parse and manage GFF annotations for multiple assemblies.

    Handles reference, scaffold, and contig GFF files separately,
    allowing gene-level comparison across assemblies.
    """

    def __init__(self, reference_gff: Optional[str] = None,
                 scaffold_gff: Optional[str] = None,
                 contig_gff: Optional[str] = None,
                 reference_gff_parser: Optional[object] = None):
        """
        Initialize multi-assembly GFF parser.

        Args:
            reference_gff: Path to reference GFF3 file
            scaffold_gff: Path to scaffold GFF3 file
            contig_gff: Path to contig GFF3 file
            reference_gff_parser: Optional pre-parsed GFFParser object for reference (overrides reference_gff)
        """
        self.reference_genes = {}
        self.scaffold_genes = {}
        self.contig_genes = {}

        # Use pre-parsed reference GFF parser if provided (e.g., for rotated coordinates)
        if reference_gff_parser:
            self.reference_genes = reference_gff_parser.genes_by_seq
            print(f"  Using pre-parsed reference GFF (rotated coordinates)")
        elif reference_gff:
            self.reference_genes = self._parse_gff(reference_gff, "reference")

        if scaffold_gff:
            self.scaffold_genes = self._parse_gff(scaffold_gff, "scaffold")
        if contig_gff:
            self.contig_genes = self._parse_gff(contig_gff, "contig")

    def _parse_gff(self, gff_file: str, assembly_type: str) -> Dict[str, List[dict]]:
        """Parse GFF3 file and return genes organized by sequence ID."""
        genes_by_seq = defaultdict(list)

        print(f"  Parsing {assembly_type} GFF: {gff_file}")

        with open(gff_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue

                feature_type = parts[2]
                if feature_type not in ['gene', 'CDS']:
                    continue

                seqid = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                attributes = parts[8]

                # Parse attributes
                attr_dict = {}
                for attr in attributes.split(';'):
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        attr_dict[key] = value

                gene_name = attr_dict.get('Name',
                            attr_dict.get('ID',
                            attr_dict.get('locus_tag',
                            f'gene_{start}_{end}')))

                # Get product/description if available
                product = attr_dict.get('product', attr_dict.get('description', ''))

                genes_by_seq[seqid].append({
                    'seqid': seqid,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'name': gene_name,
                    'product': product,
                    'type': feature_type,
                    'length': end - start + 1
                })

        # Sort by position
        for seqid in genes_by_seq:
            genes_by_seq[seqid].sort(key=lambda x: x['start'])

        total_genes = sum(len(genes) for genes in genes_by_seq.values())
        print(f"    Found {total_genes} genes across {len(genes_by_seq)} sequences")

        return dict(genes_by_seq)

    def get_genes_in_range(self, assembly_type: str, seqid: str,
                          start: int, end: int) -> List[dict]:
        """
        Get genes within a coordinate range.

        Args:
            assembly_type: 'reference', 'scaffold', or 'contig'
            seqid: Sequence ID to search
            start: Start position
            end: End position

        Returns:
            List of gene dictionaries within the range
        """
        if assembly_type == 'reference':
            genes = self.reference_genes.get(seqid, [])
        elif assembly_type == 'scaffold':
            genes = self.scaffold_genes.get(seqid, [])
        elif assembly_type == 'contig':
            genes = self.contig_genes.get(seqid, [])
        else:
            return []

        return [g for g in genes if g['start'] <= end and g['end'] >= start]

    def get_all_genes(self, assembly_type: str) -> Dict[str, List[dict]]:
        """Get all genes for an assembly type."""
        if assembly_type == 'reference':
            return self.reference_genes
        elif assembly_type == 'scaffold':
            return self.scaffold_genes
        elif assembly_type == 'contig':
            return self.contig_genes
        return {}

    def has_genes(self, assembly_type: str) -> bool:
        """Check if genes are available for an assembly type."""
        genes = self.get_all_genes(assembly_type)
        return len(genes) > 0
