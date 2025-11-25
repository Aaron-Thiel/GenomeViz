"""
Sequence handling and annotation parsing module

Provides classes for:
- Parsing GFF3 annotation files
- Managing reference sequence data
- Gene-level quality analysis

Classes:
    GFFParser: Parse GFF3 annotation files
    ReferenceSequence: Container for reference sequence and analysis results
    GeneLevelAnalyzer: Analyze alignment quality at gene level

GitHub: https://github.com/Aaron-Thiel/GenomeViz
License: MIT
"""

import numpy as np
from collections import defaultdict


class GFFParser:
    """
    Parse GFF3 annotation files to extract gene features.

    Attributes:
        genes_by_seq (dict): Dictionary mapping sequence IDs to lists of gene features
    """

    def __init__(self, gff_file):
        """
        Initialize parser and parse GFF3 file.

        Args:
            gff_file (str): Path to GFF3 annotation file
        """
        self.genes_by_seq = defaultdict(list)
        self.parse_gff(gff_file)

    def parse_gff(self, gff_file):
        """
        Parse GFF3 file and extract gene/CDS features.

        Args:
            gff_file (str): Path to GFF3 file
        """
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

                # Parse attributes field
                attr_dict = {}
                for attr in attributes.split(';'):
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        attr_dict[key] = value

                # Extract gene name from various possible attribute keys
                gene_name = attr_dict.get('Name',
                                         attr_dict.get('ID',
                                         attr_dict.get('locus_tag',
                                         f'gene_{start}_{end}')))

                self.genes_by_seq[seqid].append({
                    'seqid': seqid,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'name': gene_name,
                    'type': feature_type,
                    'length': end - start + 1
                })

        # Sort genes by position
        for seqid in self.genes_by_seq:
            self.genes_by_seq[seqid].sort(key=lambda x: x['start'])

        # Print summary
        total_genes = sum(len(genes) for genes in self.genes_by_seq.values())
        print(f"Parsed {total_genes} gene features from GFF3 across {len(self.genes_by_seq)} sequences")
        for seqid, genes in self.genes_by_seq.items():
            print(f"  {seqid}: {len(genes)} genes")


def find_oric_position(gff_file):
    """
    Find oriC position from GFF file.
    
    Returns:
        tuple: (position, seqid) or (None, None) if not found
    """
    oric_features = []
    
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            feature_type = parts[2]
            attributes = parts[8]
            
            # Check if this is an oriC feature
            if 'oriC' in feature_type or 'oriC' in attributes:
                seqid = parts[0]
                start = int(parts[3])
                oric_features.append((start, seqid, line.strip()))
    
    if len(oric_features) == 0:
        return None, None
    
    if len(oric_features) > 1:
        print(f"⚠️  Warning: Found {len(oric_features)} oriC features in GFF file:")
        for start, seqid, line in oric_features:
            print(f"    - {seqid}:{start}")
        print(f"    Using first oriC at position {oric_features[0][0]}")
    
    return oric_features[0][0], oric_features[0][1]


class ReferenceSequence:
    """
    Container for a single reference sequence and its analysis results.

    Attributes:
        seqid (str): Sequence identifier
        sequence (str): DNA sequence
        length (int): Sequence length in bp
        genes (list): List of gene features
        alignments (list): List of assembly alignments
        gene_stats (list): Per-gene quality statistics
        misassemblies (list): Detected misassemblies
        contig_mapping (dict): Contig-to-reference mapping
        seq_type (str): 'chromosome' or 'plasmid'
    """

    def __init__(self, seqid, sequence, genes):
        """
        Initialize reference sequence container.

        Args:
            seqid (str): Sequence identifier
            sequence (str): DNA sequence
            genes (list): List of gene features
        """
        self.seqid = seqid
        self.sequence = sequence
        self.length = len(sequence)
        self.genes = genes
        self.alignments = []
        self.gene_stats = []
        self.misassemblies = []
        self.contig_mapping = []

        # Classify as chromosome or plasmid based on size (>500kb = chromosome)
        self.seq_type = 'chromosome' if self.length > 500000 else 'plasmid'


class GeneLevelAnalyzer:
    """
    Analyze alignment quality at the gene level.

    Calculates coverage and identity statistics for each gene, providing
    a quality score and status classification.
    """

    def __init__(self, genes, alignments, reference_length):
        """
        Initialize gene-level analyzer.

        Args:
            genes (list): List of gene features
            alignments (list): List of alignments
            reference_length (int): Length of reference sequence in bp
        """
        self.genes = genes
        self.alignments = alignments
        self.reference_length = reference_length
        self.gene_stats = []

    def analyze(self):
        """
        Calculate per-gene alignment statistics.

        Returns:
            list: List of gene statistics dictionaries containing:
                - name: Gene name
                - coverage_pct: Percentage of gene covered
                - avg_identity: Average identity of covered regions
                - quality_score: Combined quality metric (0-100)
                - status: Classification (complete/incomplete/divergent/missing)
        """
        # Build coverage and identity maps
        coverage = np.zeros(self.reference_length)
        identity_map = np.zeros(self.reference_length)

        for aln in self.alignments:
            if aln['is_primary']:
                start = aln['ref_start']
                end = aln['ref_end']
                identity = aln['identity']

                coverage[start:end] += 1
                identity_map[start:end] = np.maximum(identity_map[start:end], identity)

        # Calculate statistics for each gene
        for gene in self.genes:
            start = max(0, gene['start'] - 1)
            end = min(self.reference_length, gene['end'])

            gene_length = end - start
            if gene_length == 0:
                continue

            gene_coverage = coverage[start:end]
            gene_identity = identity_map[start:end]

            covered_bases = np.sum(gene_coverage > 0)
            coverage_pct = (covered_bases / gene_length) * 100

            if covered_bases > 0:
                avg_identity = np.mean(gene_identity[gene_identity > 0])
            else:
                avg_identity = 0

            # Quality score: weighted combination (30% coverage, 70% identity)
            quality_score = (coverage_pct * 0.3 + avg_identity * 0.7)

            # Assign status based on metrics
            if coverage_pct < 50:
                status = 'missing'
            elif avg_identity < 90:
                status = 'divergent'
            elif coverage_pct < 95:
                status = 'incomplete'
            else:
                status = 'complete'

            self.gene_stats.append({
                'name': gene['name'],
                'seqid': gene['seqid'],
                'start': gene['start'],
                'end': gene['end'],
                'strand': gene['strand'],
                'length': gene_length,
                'coverage_pct': coverage_pct,
                'avg_identity': avg_identity,
                'quality_score': quality_score,
                'status': status
            })

        return self.gene_stats