#!/usr/bin/env python3
"""
GenomeViz - Circular Genome Assembly Visualization Tool

A comprehensive tool for visualizing and comparing bacterial genome assemblies against
reference sequences using interactive circular and linear plots.

Features:
- Circular visualizations with three information rings
- Interactive Plotly plots with hover information and zoom
- Static matplotlib plots for publication
- Automatic contig orientation detection
- Gene-level quality assessment
- Detection of gaps, inversions, and duplications

Author: Aaron Christopher Thiel <aaron.chris.thiel@gmail.com>
License: MIT
GitHub: https://github.com/Aaron-Thiel/GenomeViz
"""

import mappy as mp
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge
import plotly.graph_objects as go
from Bio import SeqIO
from Bio.Seq import Seq
import json
import tempfile
import os


__version__ = "1.0.0"


# ============================================================================
# GFF3 PARSING
# ============================================================================

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
                # Skip comments and headers
                if line.startswith('#'):
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                
                # Only extract gene and CDS features
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


# ============================================================================
# ORIENTATION DETECTION
# ============================================================================

class OrientationDetector:
    """
    Detect and correct orientation issues in assembly contigs.
    
    Analyzes alignment strand bias to identify contigs that should be
    reverse-complemented for optimal alignment to the reference.
    """
    
    def __init__(self, reference_fasta, assembly_fasta, preset='asm10'):
        """
        Initialize orientation detector.
        
        Args:
            reference_fasta (str): Path to reference genome
            assembly_fasta (str): Path to assembly to check
            preset (str): Minimap2 preset for alignment
        """
        self.reference_fasta = reference_fasta
        self.assembly_fasta = assembly_fasta
        self.preset = preset
        self.corrections = {}
        
    def detect_and_correct(self):
        """
        Detect inverted contigs and create orientation-corrected assembly.
        
        Returns:
            str: Path to corrected assembly (or original if no corrections needed)
        """
        print("\n" + "="*70)
        print("AUTOMATIC ORIENTATION DETECTION")
        print("="*70)
        
        # Create aligner
        aligner = mp.Aligner(self.reference_fasta, preset=self.preset)
        if not aligner:
            raise Exception("Failed to load reference")
        
        contig_stats = {}
        
        # Analyze each contig's alignment strand bias
        for name, seq, qual in mp.fastx_read(self.assembly_fasta):
            forward_bases = 0
            reverse_bases = 0
            
            for hit in aligner.map(seq):
                if hit.is_primary:
                    aligned_length = hit.r_en - hit.r_st
                    if hit.strand == 1:
                        forward_bases += aligned_length
                    else:
                        reverse_bases += aligned_length
            
            total_bases = forward_bases + reverse_bases
            if total_bases > 0:
                reverse_pct = (reverse_bases / total_bases) * 100
                contig_stats[name] = {
                    'forward_bases': forward_bases,
                    'reverse_bases': reverse_bases,
                    'reverse_pct': reverse_pct,
                    'needs_rc': reverse_pct > 50,  # RC if >50% aligns in reverse
                    'seq': seq
                }
        
        needs_correction = [name for name, stats in contig_stats.items() if stats['needs_rc']]
        
        # Print summary
        print(f"\nOrientation analysis for {len(contig_stats)} contigs:")
        print(f"  Contigs in correct orientation: {len(contig_stats) - len(needs_correction)}")
        print(f"  Contigs needing reverse complement: {len(needs_correction)}")
        
        if needs_correction:
            print("\nContigs to be reverse complemented:")
            for name in needs_correction:
                stats = contig_stats[name]
                print(f"  - {name}: {stats['reverse_pct']:.1f}% aligned in reverse")
        
        # Create corrected assembly if needed
        if needs_correction:
            corrected_file = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
            corrected_path = corrected_file.name
            
            print(f"\nCreating orientation-corrected assembly...")
            
            for name, seq, qual in mp.fastx_read(self.assembly_fasta):
                if name in contig_stats:
                    if contig_stats[name]['needs_rc']:
                        # Reverse complement this contig
                        seq_obj = Seq(seq)
                        rc_seq = str(seq_obj.reverse_complement())
                        corrected_file.write(f">{name}_RC\n")
                        for i in range(0, len(rc_seq), 80):
                            corrected_file.write(rc_seq[i:i+80] + "\n")
                        self.corrections[name] = 'reverse_complemented'
                    else:
                        # Keep as is
                        corrected_file.write(f">{name}\n")
                        for i in range(0, len(seq), 80):
                            corrected_file.write(seq[i:i+80] + "\n")
                        self.corrections[name] = 'kept_forward'
            
            corrected_file.close()
            print(f"Corrected assembly saved to temporary file: {corrected_path}")
            return corrected_path
        else:
            print("\nNo orientation corrections needed - all contigs in correct orientation!")
            return self.assembly_fasta


# ============================================================================
# REFERENCE SEQUENCE CONTAINER
# ============================================================================

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
        
        # Classify as chromosome or plasmid based on size
        self.seq_type = 'chromosome' if self.length > 500000 else 'plasmid'


# ============================================================================
# GENOME ALIGNMENT
# ============================================================================

class GenomeAligner:
    """
    Align assembly contigs to reference genome using minimap2.
    """
    
    def __init__(self, reference_fasta, assembly_fasta, preset='asm10'):
        """
        Initialize genome aligner.
        
        Args:
            reference_fasta (str): Path to reference genome
            assembly_fasta (str): Path to assembly
            preset (str): Minimap2 preset (asm5/asm10/asm20)
        """
        self.reference_fasta = reference_fasta
        self.assembly_fasta = assembly_fasta
        self.preset = preset
        self.reference_sequences = {}
        self.assembly_name = Path(assembly_fasta).stem
        
    def load_reference(self):
        """Load all sequences from reference genome."""
        print("\nLoading reference sequences...")
        for record in SeqIO.parse(self.reference_fasta, "fasta"):
            self.reference_sequences[record.id] = {
                'sequence': str(record.seq),
                'length': len(record.seq)
            }
            seq_type = 'chromosome' if len(record.seq) > 500000 else 'plasmid'
            print(f"  {record.id}: {len(record.seq):,} bp ({seq_type})")
        
    def align(self):
        """
        Perform alignment using minimap2.
        
        Returns:
            dict: Alignments grouped by reference sequence
        """
        print("\nAligning assembly to reference...")
        
        aligner = mp.Aligner(self.reference_fasta, preset=self.preset)
        if not aligner:
            raise Exception("Failed to load reference")
        
        alignments_by_ref = defaultdict(list)
        
        alignment_count = 0
        for name, seq, qual in mp.fastx_read(self.assembly_fasta):
            seq_len = len(seq)
            for hit in aligner.map(seq):
                alignment = {
                    'query_name': name,
                    'query_start': hit.q_st,
                    'query_end': hit.q_en,
                    'query_length': seq_len,
                    'ref_start': hit.r_st,
                    'ref_end': hit.r_en,
                    'ref_name': hit.ctg,
                    'strand': '+' if hit.strand == 1 else '-',
                    'mapq': hit.mapq,
                    'matches': hit.mlen,
                    'block_length': hit.blen,
                    'identity': hit.mlen / hit.blen * 100 if hit.blen > 0 else 0,
                    'is_primary': hit.is_primary
                }
                alignments_by_ref[hit.ctg].append(alignment)
                alignment_count += 1
        
        print(f"Found {alignment_count} total alignments")
        
        # Sort and summarize
        for ref_name in alignments_by_ref:
            alignments_by_ref[ref_name].sort(key=lambda x: x['ref_start'])
            primary_count = len([a for a in alignments_by_ref[ref_name] if a['is_primary']])
            print(f"  {ref_name}: {len(alignments_by_ref[ref_name])} alignments ({primary_count} primary)")
        
        return alignments_by_ref


# ============================================================================
# MISASSEMBLY DETECTION
# ============================================================================

class MisassemblyDetector:
    """
    Detect potential misassemblies (inversions, gaps, overlaps).
    """
    
    def __init__(self, alignments, min_gap=1000, min_inversion_size=500):
        """
        Initialize misassembly detector.
        
        Args:
            alignments (list): List of alignments
            min_gap (int): Minimum gap size to report (bp)
            min_inversion_size (int): Minimum inversion size to report (bp)
        """
        self.alignments = alignments
        self.min_gap = min_gap
        self.min_inversion_size = min_inversion_size
        self.misassemblies = []
        
    def detect(self):
        """
        Detect various types of misassemblies.
        
        Returns:
            list: List of detected misassemblies
        """
        primary_aligns = [a for a in self.alignments if a['is_primary']]
        
        # Detect inversions (reverse strand alignments)
        for aln in primary_aligns:
            if aln['strand'] == '-' and aln['block_length'] > self.min_inversion_size:
                self.misassemblies.append({
                    'type': 'inversion',
                    'ref_start': aln['ref_start'],
                    'ref_end': aln['ref_end'],
                    'query_name': aln['query_name'],
                    'size': aln['ref_end'] - aln['ref_start']
                })
        
        # Detect gaps and overlaps between consecutive alignments
        for i in range(len(primary_aligns) - 1):
            curr = primary_aligns[i]
            next_aln = primary_aligns[i + 1]
            
            gap = next_aln['ref_start'] - curr['ref_end']
            
            if gap > self.min_gap:
                # Large gap - missing sequence
                self.misassemblies.append({
                    'type': 'gap',
                    'ref_start': curr['ref_end'],
                    'ref_end': next_aln['ref_start'],
                    'query_name': f"{curr['query_name']}->{next_aln['query_name']}",
                    'size': gap
                })
            elif gap < -self.min_gap:
                # Large overlap
                self.misassemblies.append({
                    'type': 'overlap',
                    'ref_start': next_aln['ref_start'],
                    'ref_end': curr['ref_end'],
                    'query_name': f"{curr['query_name']}<->{next_aln['query_name']}",
                    'size': abs(gap)
                })
        
        return self.misassemblies


# ============================================================================
# GENE-LEVEL ANALYSIS
# ============================================================================

class GeneLevelAnalyzer:
    """
    Analyze alignment quality at the gene level.
    """
    
    def __init__(self, genes, alignments, reference_length):
        """
        Initialize gene-level analyzer.
        
        Args:
            genes (list): List of gene features
            alignments (list): List of alignments
            reference_length (int): Length of reference sequence
        """
        self.genes = genes
        self.alignments = alignments
        self.reference_length = reference_length
        self.gene_stats = []
        
    def analyze(self):
        """
        Calculate per-gene alignment statistics.
        
        Returns:
            list: List of gene statistics dictionaries
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
            
            # Quality score: weighted combination of coverage and identity
            quality_score = (coverage_pct * 0.3 + avg_identity * 0.7)
            
            # Assign status
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


# ============================================================================
# CONTIG MAPPING
# ============================================================================

class ContigMapper:
    """
    Create contig-to-reference mapping for visualization.
    """
    
    def __init__(self, alignments_by_ref):
        """
        Initialize contig mapper.
        
        Args:
            alignments_by_ref (dict): Alignments grouped by reference
        """
        self.alignments_by_ref = alignments_by_ref
        self.contig_mapping = defaultdict(lambda: defaultdict(list))
        
    def create_mapping(self):
        """
        Create detailed contig-to-reference mapping.
        
        Returns:
            dict: Nested dict of contig mappings
        """
        print("\nCreating contig-to-reference mapping...")
        
        for ref_name, alignments in self.alignments_by_ref.items():
            primary_aligns = [a for a in alignments if a['is_primary']]
            
            for aln in primary_aligns:
                contig_name = aln['query_name']
                
                self.contig_mapping[ref_name][contig_name].append({
                    'ref_start': aln['ref_start'],
                    'ref_end': aln['ref_end'],
                    'query_start': aln['query_start'],
                    'query_end': aln['query_end'],
                    'query_length': aln['query_length'],
                    'strand': aln['strand'],
                    'identity': aln['identity'],
                    'coverage': (aln['ref_end'] - aln['ref_start'])
                })
        
        # Print summary
        for ref_name in sorted(self.contig_mapping.keys()):
            contigs = self.contig_mapping[ref_name]
            print(f"\n  {ref_name}: {len(contigs)} contigs aligned")
            for contig_name in sorted(contigs.keys()):
                alns = contigs[contig_name]
                total_coverage = sum(a['coverage'] for a in alns)
                avg_identity = np.mean([a['identity'] for a in alns])
                print(f"    - {contig_name}: {len(alns)} alignment(s), "
                      f"{total_coverage:,} bp covered, {avg_identity:.1f}% identity")
        
        return dict(self.contig_mapping)


# ============================================================================
# INTERACTIVE PLOTLY VISUALIZATION
# ============================================================================

class InteractivePlotlyVisualizer:
    """
    Create interactive Plotly circular plots with hover information.
    """
    
    @staticmethod
    def create_interactive_plot(ref_seq, output_file):
        """
        Create interactive circular plot using Plotly.
        
        Args:
            ref_seq (ReferenceSequence): Reference sequence with analysis results
            output_file (str): Path to output HTML file
        """
        print(f"  Creating interactive plot: {output_file}")
        
        # Initialize figure
        fig = go.Figure()
        
        # Calculate angles for positions
        theta_values = np.linspace(0, 360, ref_seq.length, endpoint=False)
        
        # Build coverage and strand tracking arrays
        coverage = np.zeros(ref_seq.length)
        strand_info = {}
        
        for aln in ref_seq.alignments:
            if aln['is_primary']:
                for pos in range(aln['ref_start'], aln['ref_end']):
                    coverage[pos] += 1
                    if pos not in strand_info:
                        strand_info[pos] = []
                    strand_info[pos].append(aln['strand'])
        
        # ====================================================================
        # RING 1: GENE QUALITY (outer ring) - ADD FIRST
        # ====================================================================
        
        quality_colors = {
            'excellent': 'rgb(46, 204, 113)',   # Green
            'good': 'rgb(241, 196, 15)',         # Yellow
            'fair': 'rgb(230, 126, 34)',         # Orange
            'poor': 'rgb(231, 76, 60)'           # Red
        }
        
        for quality_level in ['excellent', 'good', 'fair', 'poor']:
            theta = []
            r = []
            hover_text = []
            
            for gene in ref_seq.gene_stats:
                quality = gene['quality_score']
                
                if quality_level == 'excellent' and quality < 95:
                    continue
                elif quality_level == 'good' and (quality < 85 or quality >= 95):
                    continue
                elif quality_level == 'fair' and (quality < 70 or quality >= 85):
                    continue
                elif quality_level == 'poor' and quality >= 70:
                    continue
                
                start_angle = (gene['start'] / ref_seq.length) * 360
                end_angle = (gene['end'] / ref_seq.length) * 360
                
                arc_angles = np.linspace(start_angle, end_angle, max(2, int((end_angle - start_angle) / 2)))
                theta.extend(arc_angles)
                r.extend([0.85] * len(arc_angles))
                
                hover = (f"Gene: {gene['name']}<br>"
                        f"Position: {gene['start']:,}-{gene['end']:,} bp<br>"
                        f"Length: {gene['length']:,} bp<br>"
                        f"Quality: {quality:.1f}/100<br>"
                        f"Coverage: {gene['coverage_pct']:.1f}%<br>"
                        f"Identity: {gene['avg_identity']:.1f}%<br>"
                        f"Status: {gene['status']}")
                hover_text.extend([hover] * len(arc_angles))
                
                theta.append(None)
                r.append(None)
                hover_text.append(None)
            
            if theta:
                fig.add_trace(go.Scatterpolar(
                    r=r,
                    theta=theta,
                    mode='lines',
                    name=f'Gene {quality_level}',
                    line=dict(color=quality_colors[quality_level], width=15),
                    hovertext=hover_text,
                    hoverinfo='text',
                    showlegend=True,
                    legendgroup="Ring 1: Gene Quality",
                    legendgrouptitle=dict(text="<b>Ring 1: Gene Quality</b>", font=dict(size=11))
                ))
        
        # ====================================================================
        # RING 2: ALIGNMENT STATUS (middle ring) - ADD SECOND
        # ====================================================================
        
        status_colors = {
            'complete': 'rgb(46, 204, 113)',    # Green
            'duplicated': 'rgb(243, 156, 18)',  # Orange
            'inverted': 'rgb(231, 76, 60)',     # Red
            'missing': 'rgb(149, 165, 166)'     # Gray
        }
        
        # Build segments
        segments = []
        current_segment = None
        
        for pos in range(ref_seq.length):
            if coverage[pos] == 0:
                status = 'missing'
            elif coverage[pos] > 1:
                status = 'duplicated'
            elif pos in strand_info and '-' in strand_info[pos]:
                status = 'inverted'
            else:
                status = 'complete'
            
            if current_segment is None:
                current_segment = {'start': pos, 'end': pos, 'status': status}
            elif current_segment['status'] == status:
                current_segment['end'] = pos
            else:
                segments.append(current_segment)
                current_segment = {'start': pos, 'end': pos, 'status': status}
        
        if current_segment is not None:
            segments.append(current_segment)
        
        # Create traces for alignment status
        for status_type in ['complete', 'duplicated', 'inverted', 'missing']:
            status_segments = [s for s in segments if s['status'] == status_type]
            
            if not status_segments:
                continue
            
            theta = []
            r = []
            hover_text = []
            
            for seg in status_segments:
                start_angle = (seg['start'] / ref_seq.length) * 360
                end_angle = ((seg['end'] + 1) / ref_seq.length) * 360
                
                # Create arc
                arc_angles = np.linspace(start_angle, end_angle, max(2, int((end_angle - start_angle) / 2)))
                theta.extend(arc_angles)
                r.extend([0.60] * len(arc_angles))
                
                # Hover text
                size = seg['end'] - seg['start'] + 1
                hover = f"Status: {status_type.capitalize()}<br>Position: {seg['start']:,}-{seg['end']:,} bp<br>Size: {size:,} bp"
                hover_text.extend([hover] * len(arc_angles))
                
                # Add None to separate segments
                theta.append(None)
                r.append(None)
                hover_text.append(None)
            
            fig.add_trace(go.Scatterpolar(
                r=r,
                theta=theta,
                mode='lines',
                name=f'{status_type.capitalize()}',
                line=dict(color=status_colors[status_type], width=15),
                hovertext=hover_text,
                hoverinfo='text',
                showlegend=True,
                legendgroup="Ring 2: Alignment Status",
                legendgrouptitle=dict(text="<b>Ring 2: Alignment Status</b>", font=dict(size=11))
            ))
        
        # ====================================================================
        # RING 3: CONTIG MAPPING (inner ring) - ADD THIRD
        # ====================================================================
        
        # Get unique contigs and assign colors
        contigs_in_order = []
        for aln in sorted(ref_seq.alignments, key=lambda x: x['ref_start']):
            if aln['is_primary'] and aln['query_name'] not in contigs_in_order:
                contigs_in_order.append(aln['query_name'])
        
        # Use color palette
        if len(contigs_in_order) <= 12:
            import plotly.express as px
            colors = px.colors.qualitative.Set3[:len(contigs_in_order)]
        else:
            colors = [f'hsl({i * 360 / len(contigs_in_order)}, 70%, 60%)' 
                     for i in range(len(contigs_in_order))]
        
        contig_colors = dict(zip(contigs_in_order, colors))
        
        for contig_name in contigs_in_order:
            theta = []
            r = []
            hover_text = []
            
            for aln in ref_seq.alignments:
                if not aln['is_primary'] or aln['query_name'] != contig_name:
                    continue
                
                start_angle = (aln['ref_start'] / ref_seq.length) * 360
                end_angle = (aln['ref_end'] / ref_seq.length) * 360
                
                arc_angles = np.linspace(start_angle, end_angle, max(2, int((end_angle - start_angle) / 2)))
                theta.extend(arc_angles)
                r.extend([0.35] * len(arc_angles))
                
                hover = (f"Contig: {contig_name}<br>"
                        f"Ref position: {aln['ref_start']:,}-{aln['ref_end']:,} bp<br>"
                        f"Contig position: {aln['query_start']:,}-{aln['query_end']:,} bp<br>"
                        f"Length: {aln['ref_end'] - aln['ref_start']:,} bp<br>"
                        f"Identity: {aln['identity']:.1f}%<br>"
                        f"Strand: {aln['strand']}")
                hover_text.extend([hover] * len(arc_angles))
                
                theta.append(None)
                r.append(None)
                hover_text.append(None)
            
            if theta:
                fig.add_trace(go.Scatterpolar(
                    r=r,
                    theta=theta,
                    mode='lines',
                    name=contig_name,
                    line=dict(color=contig_colors[contig_name], width=15),
                    hovertext=hover_text,
                    hoverinfo='text',
                    showlegend=True,
                    legendgroup="Ring 3: Assembly Contigs",
                    legendgrouptitle=dict(text="<b>Ring 3: Assembly Contigs</b>", font=dict(size=11))
                ))
        
        # Update layout
        seq_type = ref_seq.seq_type.capitalize()
        if ref_seq.length > 10000:
            size_str = f'{ref_seq.length/1000:.1f} kb'
        else:
            size_str = f'{ref_seq.length:,} bp'
        
        fig.update_layout(
            title=dict(
                text=f'{seq_type}: {ref_seq.seqid} ({size_str})<br>'
                     '<sub>Outer: Gene Quality | Middle: Alignment Status | Inner: Contigs</sub>',
                x=0.5,
                xanchor='center'
            ),
            polar=dict(
                radialaxis=dict(visible=False, range=[0, 1]),
                angularaxis=dict(
                    tickmode='array',
                    tickvals=np.linspace(0, 360, 8, endpoint=False),
                    ticktext=[f'{int(pos/1000)}' if ref_seq.length > 10000 else f'{int(pos):,}' 
                             for pos in np.linspace(0, ref_seq.length, 8, endpoint=False)],
                    direction='clockwise',
                    rotation=90
                )
            ),
            showlegend=True,
            legend=dict(
                orientation="v",
                yanchor="top",
                y=1,
                xanchor="left",
                x=1.02,
                title=dict(text="<b>Rings</b>", font=dict(size=14)),
                tracegroupgap=20,
                itemsizing='constant',
                font=dict(size=10)
            ),
            width=1400,
            height=1000
        )
        
        # Save as HTML
        fig.write_html(output_file)


# ============================================================================
# STATIC CIRCULAR VISUALIZATION (MATPLOTLIB)
# ============================================================================

class CircularVisualizer:
    """
    Create static circular visualizations using matplotlib (for publication).
    """
    
    @staticmethod
    def get_contig_colors(alignments):
        """
        Get consistent color mapping for contigs.
        
        Args:
            alignments (list): List of alignments
            
        Returns:
            dict: Mapping of contig names to colors
        """
        contigs_in_order = []
        for aln in sorted(alignments, key=lambda x: x['ref_start']):
            if aln['is_primary'] and aln['query_name'] not in contigs_in_order:
                contigs_in_order.append(aln['query_name'])
        
        # Use matplotlib colormap
        if len(contigs_in_order) <= 12:
            cmap = plt.cm.Set3
        elif len(contigs_in_order) <= 20:
            cmap = plt.cm.tab20
        else:
            cmap = plt.cm.hsv
        
        contig_colors = {}
        for idx, contig in enumerate(contigs_in_order):
            contig_colors[contig] = cmap(idx / max(len(contigs_in_order), 2))
        
        return contig_colors
    
    @staticmethod
    def create_circular_plot(ref_seq, output_file):
        """
        Create static circular plot using matplotlib.
        
        Args:
            ref_seq (ReferenceSequence): Reference sequence with analysis results
            output_file (str): Path to output PNG file
        """
        print(f"  Creating circular plot: {output_file}")
        
        fig = plt.figure(figsize=(16, 14))
        ax = fig.add_subplot(111, projection='polar')
        
        # Get contig colors
        contig_colors = CircularVisualizer.get_contig_colors(ref_seq.alignments)
        
        # ====================================================================
        # RING 1: GENE QUALITY (outer ring)
        # ====================================================================
        
        quality_track_r = 0.85
        quality_track_width = 0.15
        
        for gene in ref_seq.gene_stats:
            start_angle = (gene['start'] / ref_seq.length) * 2 * np.pi
            end_angle = (gene['end'] / ref_seq.length) * 2 * np.pi
            
            quality = gene['quality_score']
            
            if quality >= 95:
                color = '#2ecc71'
            elif quality >= 85:
                color = '#f1c40f'
            elif quality >= 70:
                color = '#e67e22'
            else:
                color = '#e74c3c'
            
            wedge = Wedge((0, 0), quality_track_r, 
                         np.degrees(start_angle), np.degrees(end_angle),
                         width=quality_track_width, facecolor=color, 
                         edgecolor='none', alpha=0.8,
                         transform=ax.transData._b)
            ax.add_artist(wedge)
        
        # ====================================================================
        # BUILD COVERAGE AND STRAND INFO FOR RING 2
        # ====================================================================
        
        coverage = np.zeros(ref_seq.length)
        strand_info = {}
        
        for aln in ref_seq.alignments:
            if aln['is_primary']:
                for pos in range(aln['ref_start'], aln['ref_end']):
                    coverage[pos] += 1
                    if pos not in strand_info:
                        strand_info[pos] = []
                    strand_info[pos].append(aln['strand'])
        
        # ====================================================================
        # RING 2: ALIGNMENT STATUS (middle ring)
        # ====================================================================
        
        status_track_r = 0.60
        status_track_width = 0.15
        
        status_colors = {
            'complete': '#2ecc71',      # Green
            'duplicated': '#f39c12',    # Orange
            'inverted': '#e74c3c',      # Red
            'missing': '#95a5a6'        # Gray
        }
        
        # Build segments
        segments = []
        current_segment = None
        
        for pos in range(ref_seq.length):
            if coverage[pos] == 0:
                status = 'missing'
            elif coverage[pos] > 1:
                status = 'duplicated'
            elif pos in strand_info and '-' in strand_info[pos]:
                status = 'inverted'
            else:
                status = 'complete'
            
            if current_segment is None:
                current_segment = {'start': pos, 'end': pos, 'status': status}
            elif current_segment['status'] == status:
                current_segment['end'] = pos
            else:
                segments.append(current_segment)
                current_segment = {'start': pos, 'end': pos, 'status': status}
        
        if current_segment is not None:
            segments.append(current_segment)
        
        # Draw segments as wedges
        for segment in segments:
            start_angle = (segment['start'] / ref_seq.length) * 2 * np.pi
            end_angle = ((segment['end'] + 1) / ref_seq.length) * 2 * np.pi
            color = status_colors[segment['status']]
            
            wedge = Wedge((0, 0), status_track_r,
                         np.degrees(start_angle), np.degrees(end_angle),
                         width=status_track_width, facecolor=color,
                         edgecolor='none', alpha=0.8,
                         transform=ax.transData._b)
            ax.add_artist(wedge)
        
        # ====================================================================
        # RING 3: CONTIG MAPPING (inner ring)
        # ====================================================================
        
        contig_track_r = 0.35
        contig_track_width = 0.15
        
        for aln in ref_seq.alignments:
            if aln['is_primary']:
                contig_name = aln['query_name']
                color = contig_colors.get(contig_name, '#888888')
                
                start_angle = (aln['ref_start'] / ref_seq.length) * 2 * np.pi
                end_angle = (aln['ref_end'] / ref_seq.length) * 2 * np.pi
                
                wedge = Wedge((0, 0), contig_track_r,
                            np.degrees(start_angle), np.degrees(end_angle),
                            width=contig_track_width, facecolor=color,
                            edgecolor='black', linewidth=0.5, alpha=0.7,
                            transform=ax.transData._b)
                ax.add_artist(wedge)
        
        # ====================================================================
        # CONFIGURE PLOT
        # ====================================================================
        
        ax.set_ylim(0, 1.0)
        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)
        
        # Axis labels
        n_ticks = 8
        tick_angles = np.linspace(0, 2*np.pi, n_ticks, endpoint=False)
        tick_positions_bp = np.linspace(0, ref_seq.length, n_ticks, endpoint=False)
        
        if ref_seq.length > 10000:
            tick_labels = [f'{int(pos/1000)}' for pos in tick_positions_bp]
            unit_label = 'kb'
        else:
            tick_labels = [f'{int(pos):,}' for pos in tick_positions_bp]
            unit_label = 'bp'
        
        ax.set_xticks(tick_angles)
        ax.set_xticklabels(tick_labels, fontsize=10)
        ax.set_yticks([])
        ax.grid(True, alpha=0.3, linewidth=0.5)
        
        # Title
        seq_type = ref_seq.seq_type.capitalize()
        if ref_seq.length > 10000:
            size_str = f'{ref_seq.length/1000:.1f} kb'
        else:
            size_str = f'{ref_seq.length:,} bp'
            
        plt.title(f'{seq_type}: {ref_seq.seqid}\n({size_str})\n\n'
                 f'Ring 1 (Outer): Gene Quality | Ring 2 (Middle): Alignment Status | Ring 3 (Inner): Contig Mapping',
                 fontsize=13, fontweight='bold', pad=20)
        
        ax.text(0.5, -0.15, f'Position ({unit_label})', 
               transform=ax.transAxes, ha='center', fontsize=11)
        
        # ====================================================================
        # LEGENDS
        # ====================================================================
        
        # Legend 1: Gene Quality
        quality_legend_elements = [
            plt.Line2D([0], [0], marker='s', color='w', markerfacecolor='#2ecc71', 
                      markersize=10, label='Excellent (≥95)', markeredgecolor='black', markeredgewidth=0.5),
            plt.Line2D([0], [0], marker='s', color='w', markerfacecolor='#f1c40f', 
                      markersize=10, label='Good (85-95)', markeredgecolor='black', markeredgewidth=0.5),
            plt.Line2D([0], [0], marker='s', color='w', markerfacecolor='#e67e22', 
                      markersize=10, label='Fair (70-85)', markeredgecolor='black', markeredgewidth=0.5),
            plt.Line2D([0], [0], marker='s', color='w', markerfacecolor='#e74c3c', 
                      markersize=10, label='Poor (<70)', markeredgecolor='black', markeredgewidth=0.5)
        ]
        legend1 = ax.legend(handles=quality_legend_elements, loc='upper left', 
                           bbox_to_anchor=(1.15, 1.0), title='Ring 1: Gene Quality',
                           fontsize=9, framealpha=0.9)
        ax.add_artist(legend1)
        
        # Legend 2: Alignment Status
        status_legend_elements = [
            plt.Line2D([0], [0], marker='s', color='w', markerfacecolor=status_colors['complete'], 
                      markersize=10, label='Complete (1x)', markeredgecolor='black', markeredgewidth=0.5),
            plt.Line2D([0], [0], marker='s', color='w', markerfacecolor=status_colors['duplicated'], 
                      markersize=10, label='Duplicated (>1x)', markeredgecolor='black', markeredgewidth=0.5),
            plt.Line2D([0], [0], marker='s', color='w', markerfacecolor=status_colors['inverted'], 
                      markersize=10, label='Inverted (reverse)', markeredgecolor='black', markeredgewidth=0.5),
            plt.Line2D([0], [0], marker='s', color='w', markerfacecolor=status_colors['missing'], 
                      markersize=10, label='Missing (gap)', markeredgecolor='black', markeredgewidth=0.5)
        ]
        legend2 = ax.legend(handles=status_legend_elements, loc='upper left', 
                           bbox_to_anchor=(1.15, 0.65), title='Ring 2: Alignment Status',
                           fontsize=9, framealpha=0.9)
        ax.add_artist(legend2)
        
        # Legend 3: Contig Mapping
        if len(contig_colors) <= 15:
            contig_legend_elements = [plt.Line2D([0], [0], marker='s', color='w',
                                         markerfacecolor=color, markersize=10, 
                                         label=name, markeredgecolor='black', markeredgewidth=0.5)
                             for name, color in sorted(contig_colors.items())]
            ax.legend(handles=contig_legend_elements, loc='upper left', 
                     bbox_to_anchor=(1.15, 0.30), title='Ring 3: Assembly Contigs',
                     fontsize=9, framealpha=0.9)
        
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()


# ============================================================================
# LINEAR VISUALIZATION
# ============================================================================

class LinearVisualizer:
    """
    Create linear visualizations for detailed analysis.
    """
    
    @staticmethod
    def create_linear_plot(ref_seq, output_file):
        """
        Create linear plot showing alignment coverage and quality.
        
        Args:
            ref_seq (ReferenceSequence): Reference sequence with analysis results
            output_file (str): Path to output PNG file
        """
        print(f"  Creating linear plot: {output_file}")
        
        from matplotlib.patches import Rectangle
        from matplotlib.cm import ScalarMappable
        from matplotlib.colors import Normalize
        
        fig, axes = plt.subplots(5, 1, figsize=(16, 14), 
                                gridspec_kw={'height_ratios': [2, 1, 1, 1, 1]})
        
        seq_type_label = ref_seq.seq_type.capitalize()
        fig.suptitle(f'{seq_type_label}: {ref_seq.seqid}', 
                    fontsize=14, fontweight='bold')
        
        # 1. Gene-level quality heatmap
        ax1 = axes[0]
        ax1.set_title('Gene-Level Alignment Quality (0-100 score)', fontsize=11, pad=10)
        
        for gene in ref_seq.gene_stats:
            quality = gene['quality_score']
            
            if quality >= 95:
                color = '#2ecc71'
            elif quality >= 85:
                color = '#f1c40f'
            elif quality >= 70:
                color = '#e67e22'
            else:
                color = '#e74c3c'
            
            rect = Rectangle((gene['start'], 0), gene['end'] - gene['start'], 1, 
                           facecolor=color, edgecolor='none', alpha=0.8)
            ax1.add_patch(rect)
        
        ax1.set_xlim(0, ref_seq.length)
        ax1.set_ylim(0, 1)
        ax1.set_xlabel('Position (bp)', fontsize=10)
        ax1.set_ylabel('Quality', fontsize=10)
        ax1.set_yticks([])
        
        cmap = plt.cm.RdYlGn
        norm = Normalize(vmin=0, vmax=100)
        sm = ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax1, orientation='vertical', pad=0.01, aspect=10)
        cbar.set_label('Quality Score', fontsize=9)
        
        # 2. Contig mapping - use same colors as circular plot
        ax2 = axes[1]
        ax2.set_title('Assembly Contig Mapping', fontsize=11, pad=10)
        
        # Get consistent colors using the same method as circular plot
        contig_colors = CircularVisualizer.get_contig_colors(ref_seq.alignments)
        
        for aln in ref_seq.alignments:
            if aln['is_primary']:
                y_pos = 0.5
                color = contig_colors.get(aln['query_name'], '#888888')
                
                rect = Rectangle((aln['ref_start'], y_pos - 0.3), 
                               aln['ref_end'] - aln['ref_start'], 0.6,
                               facecolor=color, edgecolor='black', linewidth=0.5, alpha=0.7)
                ax2.add_patch(rect)
                
                # Add contig label if space allows
                width = aln['ref_end'] - aln['ref_start']
                if width > ref_seq.length * 0.05:
                    ax2.text(aln['ref_start'] + width/2, 0.5, 
                            aln['query_name'], ha='center', va='center',
                            fontsize=7, fontweight='bold')
        
        ax2.set_xlim(0, ref_seq.length)
        ax2.set_ylim(0, 1)
        ax2.set_ylabel('Contigs', fontsize=10)
        ax2.set_yticks([])
        
        # 3. Coverage
        ax3 = axes[2]
        ax3.set_title('Alignment Coverage', fontsize=11, pad=10)
        
        coverage = np.zeros(ref_seq.length)
        for aln in ref_seq.alignments:
            if aln['is_primary']:
                coverage[aln['ref_start']:aln['ref_end']] += 1
        
        positions = np.arange(ref_seq.length)
        ax3.fill_between(positions, coverage, alpha=0.6, color='#3498db')
        ax3.set_xlim(0, ref_seq.length)
        ax3.set_ylabel('Coverage', fontsize=10)
        ax3.grid(True, alpha=0.3)
        
        # 4. Identity
        ax4 = axes[3]
        ax4.set_title('Alignment Identity (%)', fontsize=11, pad=10)
        
        for aln in ref_seq.alignments:
            if aln['is_primary']:
                ax4.plot([aln['ref_start'], aln['ref_end']], 
                        [aln['identity'], aln['identity']], 
                        linewidth=2, color='#9b59b6', alpha=0.6)
        
        ax4.set_xlim(0, ref_seq.length)
        ax4.set_ylim(70, 101)
        ax4.set_ylabel('Identity (%)', fontsize=10)
        ax4.axhline(y=95, color='green', linestyle='--', alpha=0.5, linewidth=1)
        ax4.axhline(y=90, color='orange', linestyle='--', alpha=0.5, linewidth=1)
        ax4.grid(True, alpha=0.3)
        
        # 5. Misassemblies
        ax5 = axes[4]
        ax5.set_title('Detected Misassemblies', fontsize=11, pad=10)
        
        y_pos = {'inversion': 0.7, 'gap': 0.4, 'overlap': 0.1}
        colors = {'inversion': '#e74c3c', 'gap': '#f39c12', 'overlap': '#9b59b6'}
        
        for mis in ref_seq.misassemblies:
            mis_type = mis['type']
            rect = Rectangle((mis['ref_start'], y_pos[mis_type] - 0.05), 
                           mis['ref_end'] - mis['ref_start'], 0.1,
                           facecolor=colors[mis_type], edgecolor='black', 
                           linewidth=0.5, alpha=0.7)
            ax5.add_patch(rect)
        
        ax5.set_xlim(0, ref_seq.length)
        ax5.set_ylim(0, 1)
        ax5.set_xlabel('Position (bp)', fontsize=10)
        ax5.set_yticks([0.1, 0.4, 0.7])
        ax5.set_yticklabels(['Overlap', 'Gap', 'Inversion'], fontsize=9)
        
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()


# ============================================================================
# MAIN PROGRAM
# ============================================================================

def main():
    """Main program entry point."""
    
    parser = argparse.ArgumentParser(
        description='GenomeViz - Circular Genome Assembly Visualization',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  %(prog)s --reference ref.fasta --assembly asm.fasta --gff genes.gff --output results/
  
  # With custom thresholds
  %(prog)s --reference ref.fasta --assembly asm.fasta --gff genes.gff --output results/ \\
           --min-gap 2000 --min-inversion 1000
  
  # Skip auto-orientation
  %(prog)s --reference ref.fasta --assembly asm.fasta --gff genes.gff --output results/ \\
           --no-auto-orient

For more information: https://github.com/Aaron-Thiel/GenomeViz
        """
    )
    
    # Required arguments
    parser.add_argument('--reference', required=True,
                       help='Reference genome (FASTA format)')
    parser.add_argument('--assembly', required=True,
                       help='Assembly to compare (FASTA format)')
    parser.add_argument('--gff', required=True,
                       help='Gene annotations (GFF3 format)')
    parser.add_argument('--output', required=True,
                       help='Output directory')
    
    # Optional arguments
    parser.add_argument('--preset', default='asm10', choices=['asm5', 'asm10', 'asm20'],
                       help='Minimap2 preset for alignment (default: asm10)')
    parser.add_argument('--min-gap', type=int, default=1000,
                       help='Minimum gap size to report (default: 1000 bp)')
    parser.add_argument('--min-inversion', type=int, default=500,
                       help='Minimum inversion size to report (default: 500 bp)')
    parser.add_argument('--no-auto-orient', action='store_true',
                       help='Skip automatic orientation detection')
    parser.add_argument('--no-circular', action='store_true',
                       help='Skip circular plot generation')
    parser.add_argument('--no-linear', action='store_true',
                       help='Skip linear plot generation')
    parser.add_argument('--no-interactive', action='store_true',
                       help='Skip interactive Plotly plot generation')
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    
    args = parser.parse_args()
    
    # Validate input files
    for file_path in [args.reference, args.assembly, args.gff]:
        if not Path(file_path).exists():
            print(f"Error: File not found: {file_path}")
            return 1
    
    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"\nOutput directory: {output_dir.absolute()}")
    
    print("=" * 70)
    print(f"GenomeViz v{__version__}")
    print("=" * 70)
    
    # Step 1: Automatic orientation detection
    assembly_to_use = args.assembly
    temp_file_to_cleanup = None
    
    if not args.no_auto_orient:
        detector = OrientationDetector(args.reference, args.assembly, preset=args.preset)
        corrected_assembly = detector.detect_and_correct()
        
        if corrected_assembly != args.assembly:
            assembly_to_use = corrected_assembly
            temp_file_to_cleanup = corrected_assembly
            print("\n✓ Using orientation-corrected assembly for analysis")
        else:
            print("\n✓ Using original assembly (no corrections needed)")
    else:
        print("\nSkipping automatic orientation detection (--no-auto-orient)")
    
    # Step 2: Parse GFF3
    print("\n" + "="*70)
    print("[1/6] PARSING GENE ANNOTATIONS")
    print("="*70)
    gff_parser = GFFParser(args.gff)
    
    # Step 3: Load reference
    print("\n" + "="*70)
    print("[2/6] LOADING REFERENCE")
    print("="*70)
    aligner = GenomeAligner(args.reference, assembly_to_use, preset=args.preset)
    aligner.load_reference()
    
    # Step 4: Align
    print("\n" + "="*70)
    print("[3/6] ALIGNING ASSEMBLY")
    print("="*70)
    alignments_by_ref = aligner.align()
    
    # Step 5: Create contig mapping
    print("\n" + "="*70)
    print("[4/6] CONTIG MAPPING")
    print("="*70)
    mapper = ContigMapper(alignments_by_ref)
    contig_mapping = mapper.create_mapping()
    
    # Save contig mapping to JSON
    mapping_file = output_dir / 'contig_mapping.json'
    with open(mapping_file, 'w') as f:
        serializable_mapping = {}
        for ref_name, contigs in contig_mapping.items():
            serializable_mapping[ref_name] = {}
            for contig_name, alns in contigs.items():
                serializable_mapping[ref_name][contig_name] = alns
        json.dump(serializable_mapping, f, indent=2)
    print(f"\nContig mapping saved to: {mapping_file}")
    
    # Step 6: Analyze sequences
    print("\n" + "="*70)
    print("[5/6] ANALYZING SEQUENCES")
    print("="*70)
    
    ref_sequences = {}
    
    for seqid in aligner.reference_sequences:
        print(f"\nAnalyzing {seqid}...")
        
        ref_data = aligner.reference_sequences[seqid]
        genes = gff_parser.genes_by_seq.get(seqid, [])
        alignments = alignments_by_ref.get(seqid, [])
        
        ref_seq = ReferenceSequence(seqid, ref_data['sequence'], genes)
        ref_seq.alignments = alignments
        ref_seq.contig_mapping = contig_mapping.get(seqid, {})
        
        if alignments:
            detector = MisassemblyDetector(alignments, 
                                          min_gap=args.min_gap,
                                          min_inversion_size=args.min_inversion)
            ref_seq.misassemblies = detector.detect()
            print(f"  Misassemblies: {len(ref_seq.misassemblies)}")
        else:
            print(f"  WARNING: No alignments found!")
        
        if genes:
            analyzer = GeneLevelAnalyzer(genes, alignments, ref_seq.length)
            ref_seq.gene_stats = analyzer.analyze()
            
            avg_quality = np.mean([g['quality_score'] for g in ref_seq.gene_stats])
            complete = len([g for g in ref_seq.gene_stats if g['status'] == 'complete'])
            print(f"  Genes: {len(genes)}, Avg Quality: {avg_quality:.1f}, Complete: {complete}")
        
        ref_sequences[seqid] = ref_seq
    
    # Step 7: Create visualizations
    print("\n" + "="*70)
    print("[6/6] CREATING VISUALIZATIONS")
    print("="*70)
    
    if not args.no_circular:
        print("\nCreating static circular plots...")
        for seqid, ref_seq in ref_sequences.items():
            output_file = output_dir / f'{seqid}_circular.png'
            CircularVisualizer.create_circular_plot(ref_seq, output_file)
    
    if not args.no_interactive:
        print("\nCreating interactive circular plots...")
        for seqid, ref_seq in ref_sequences.items():
            output_file = output_dir / f'{seqid}_interactive.html'
            InteractivePlotlyVisualizer.create_interactive_plot(ref_seq, output_file)
    
    if not args.no_linear:
        print("\nCreating linear plots...")
        for seqid, ref_seq in ref_sequences.items():
            output_file = output_dir / f'{seqid}_linear.png'
            LinearVisualizer.create_linear_plot(ref_seq, output_file)
    
    # Save statistics
    print("\nSaving gene statistics...")
    for seqid, ref_seq in ref_sequences.items():
        if ref_seq.gene_stats:
            df = pd.DataFrame(ref_seq.gene_stats)
            csv_output = output_dir / f'{seqid}_gene_stats.csv'
            df.to_csv(csv_output, index=False)
            print(f"  {csv_output}")
    
    # Create summary report
    summary_file = output_dir / 'summary_report.txt'
    with open(summary_file, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("GENOME ASSEMBLY COMPARISON SUMMARY\n")
        f.write("=" * 70 + "\n\n")
        
        f.write(f"Assembly: {Path(args.assembly).name}\n")
        f.write(f"Reference: {Path(args.reference).name}\n\n")
        
        for seqid, ref_seq in ref_sequences.items():
            f.write(f"\n{ref_seq.seq_type.upper()}: {seqid}\n")
            f.write("-" * 70 + "\n")
            f.write(f"Length: {ref_seq.length:,} bp\n")
            f.write(f"Alignments: {len(ref_seq.alignments)}\n")
            f.write(f"Contigs aligned: {len(ref_seq.contig_mapping)}\n")
            f.write(f"Genes: {len(ref_seq.genes)}\n")
            
            if ref_seq.gene_stats:
                avg_quality = np.mean([g['quality_score'] for g in ref_seq.gene_stats])
                f.write(f"Average gene quality: {avg_quality:.2f}/100\n\n")
                
                status_counts = pd.DataFrame(ref_seq.gene_stats)['status'].value_counts()
                f.write("Gene status breakdown:\n")
                for status, count in status_counts.items():
                    f.write(f"  {status.capitalize()}: {count}\n")
            
            if ref_seq.misassemblies:
                f.write(f"\nMisassemblies: {len(ref_seq.misassemblies)}\n")
                for mis_type in ['inversion', 'gap', 'overlap']:
                    count = len([m for m in ref_seq.misassemblies if m['type'] == mis_type])
                    if count > 0:
                        f.write(f"  {mis_type.capitalize()}s: {count}\n")
            
            f.write("\nContigs mapping to this sequence:\n")
            for contig_name, alns in sorted(ref_seq.contig_mapping.items()):
                total_cov = sum(a['coverage'] for a in alns)
                avg_id = np.mean([a['identity'] for a in alns])
                f.write(f"  {contig_name}: {total_cov:,} bp, {avg_id:.1f}% identity\n")
            
            f.write("\n")
    
    print(f"\nSummary report saved to: {summary_file}")
    
    # Print final summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    
    for seqid, ref_seq in ref_sequences.items():
        print(f"\n{ref_seq.seq_type.upper()}: {seqid}")
        print(f"  Length: {ref_seq.length:,} bp")
        print(f"  Contigs: {len(ref_seq.contig_mapping)}")
        print(f"  Genes: {len(ref_seq.genes)}")
        
        if ref_seq.gene_stats:
            avg_quality = np.mean([g['quality_score'] for g in ref_seq.gene_stats])
            print(f"  Avg Quality: {avg_quality:.1f}/100")
        
        if ref_seq.misassemblies:
            print(f"  Misassemblies: {len(ref_seq.misassemblies)}")
    
    print("\n" + "=" * 70)
    print(f"All outputs saved to: {output_dir.absolute()}")
    print("=" * 70)
    print("\nGenerated files:")
    for f in sorted(output_dir.iterdir()):
        print(f"  - {f.name}")
    print("\n✓ Done!")
    
    # Cleanup
    if temp_file_to_cleanup and os.path.exists(temp_file_to_cleanup):
        os.unlink(temp_file_to_cleanup)
    
    return 0


if __name__ == '__main__':
    exit(main())