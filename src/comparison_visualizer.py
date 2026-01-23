"""
Comparison visualization module for scaffold vs contig views.

Creates visualizations with the following ring structure (matching original circular plots):
- Ring 1 (outer): Gene quality - SCAFFOLD genes colored by how well CONTIGS cover them
- Ring 2 (middle): Alignment status - complete/duplicated/inverted/missing between SCAFFOLD and CONTIGS
- Ring 3 (inner): Contig mapping - which contigs overlap this scaffold region

The scaffold is the "subject" (like reference in original plots).
The contigs are the "query" (like assembly in original plots).
Only the reference region where the scaffold exists is shown.

GitHub: https://github.com/Aaron-Thiel/GenomeViz
License: MIT
"""

import json
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from pathlib import Path
from typing import List, Dict, Optional

from .assembly_comparison import OverlapRegion, MultiAssemblyGFFParser
from .gene_comparison import GeneComparisonResult
from .interactive_circular_visualizer import InteractiveCircularVisualizer


class ScaffoldSequence:
    """
    A scaffold-centric data structure mirroring ReferenceSequence.

    This allows us to use the same visualization logic as existing plots,
    but with scaffold as the subject and contigs as the query.
    """

    def __init__(self, scaffold_name: str, ref_start: int, ref_end: int,
                 overlap_regions: List[OverlapRegion],
                 gene_results: Optional[List[GeneComparisonResult]] = None,
                 reference_length: int = None,
                 genes_are_reference_coords: bool = False):
        """
        Initialize scaffold sequence for visualization.

        Args:
            scaffold_name: Name of the scaffold
            ref_start: Start position on reference
            ref_end: End position on reference
            overlap_regions: OverlapRegion objects for this scaffold
            gene_results: Gene comparison results (scaffold genes vs contig coverage)
            reference_length: Full length of reference genome (for context visualization)
            genes_are_reference_coords: If True, gene coordinates are already in reference
                                        coordinates (from reference GFF fallback).
                                        If False, genes are in scaffold-local coordinates.
        """
        self.seqid = scaffold_name
        self.seq_type = 'scaffold'

        # The scaffold's position on the reference
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.scaffold_length = ref_end - ref_start
        
        # Track if genes are in reference coordinates (don't add ref_start)
        self.genes_are_reference_coords = genes_are_reference_coords

        # Full reference length for context visualization
        # If not provided, use scaffold region as the visualization length
        self.reference_length = reference_length if reference_length else (ref_end - ref_start)
        self.length = self.reference_length  # For visualization purposes

        # Store the actual scaffold alignment regions (not just min/max)
        # This is crucial for correct gap detection - we only want to show
        # gaps WITHIN scaffold regions, not between separate alignment blocks
        self.scaffold_regions = self._extract_scaffold_regions(overlap_regions)
        print(f"      Scaffold aligns to {len(self.scaffold_regions)} region(s) on reference")
        for i, (start, end) in enumerate(self.scaffold_regions):
            print(f"        Region {i+1}: {start:,} - {end:,} bp ({end-start:,} bp)")

        # Convert contig overlaps to alignment format (like assembly alignments)
        self.alignments = self._build_contig_alignments(overlap_regions)

        # Convert gene results to gene_stats format
        if gene_results:
            self.gene_stats = self._build_gene_stats(gene_results)
            print(f"      Built {len(self.gene_stats)} gene stats for visualization")
        else:
            self.gene_stats = []
            print(f"      Warning: No gene results provided - Gene Quality ring will be empty")

        # Build misassemblies (gaps between contigs, overlaps, inversions)
        self.misassemblies = self._detect_misassemblies(overlap_regions)

    def _extract_scaffold_regions(self, overlap_regions: List[OverlapRegion]) -> List[tuple]:
        """
        Extract the actual scaffold alignment regions on the reference.
        
        Returns a list of (start, end) tuples for each separate alignment block.
        These are used to determine which parts of the reference this scaffold
        actually covers - gaps between these regions should NOT be marked as "missing".
        """
        regions = []
        for region in overlap_regions:
            regions.append((region.scaffold_ref_start, region.scaffold_ref_end))
        
        # Sort by start position and merge overlapping regions
        regions.sort(key=lambda x: x[0])
        
        merged = []
        for start, end in regions:
            if merged and start <= merged[-1][1] + 1:
                # Overlapping or adjacent - merge
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            else:
                merged.append((start, end))
        
        return merged

    def is_within_scaffold_region(self, pos: int) -> bool:
        """Check if a position falls within any scaffold alignment region."""
        for start, end in self.scaffold_regions:
            if start <= pos <= end:
                return True
        return False

    def _build_contig_alignments(self, overlap_regions: List[OverlapRegion]) -> List[dict]:
        """Convert contig overlaps to alignment format using absolute reference coordinates."""
        alignments = []

        for region in overlap_regions:
            for contig in region.overlapping_contigs:
                # Use absolute reference coordinates for context visualization
                alignments.append({
                    'query_name': contig['contig_name'],
                    'ref_start': contig['overlap_start'],  # Absolute reference position
                    'ref_end': contig['overlap_end'],
                    'query_start': contig['contig_local_start'],
                    'query_end': contig['contig_local_end'],
                    'strand': contig['strand'],
                    'identity': contig['identity'],
                    'is_primary': True
                })

        return sorted(alignments, key=lambda x: x['ref_start'])

    def _build_gene_stats(self, gene_results: List[GeneComparisonResult]) -> List[dict]:
        """Convert gene comparison results to gene_stats format with absolute reference coordinates."""
        gene_stats = []

        for result in gene_results:
            gene = result.gene

            # Map gene status to quality score
            if result.status == 'complete':
                quality_score = 95 + (result.identity - 95) if result.identity > 95 else 95
            elif result.status == 'split':
                quality_score = 85
            elif result.status == 'truncated':
                quality_score = 70
            else:  # missing
                quality_score = 30

            # Handle coordinate transformation based on gene source:
            # Priority 1: Use ref_start/ref_end if set by analyze_scaffold_genes_via_reference
            # Priority 2: Use genes_are_reference_coords flag
            # Priority 3: Add ref_start offset for scaffold-local coords
            if 'ref_start' in gene and gene['ref_start'] is not None:
                # Gene has explicit reference coordinates (from scaffold gene mapping)
                gene_start = gene['ref_start']
                gene_end = gene['ref_end']
            elif self.genes_are_reference_coords:
                # Reference genes already have absolute coordinates
                gene_start = gene['start']
                gene_end = gene['end']
            else:
                # Scaffold genes need conversion to absolute reference position
                gene_start = gene['start'] + self.ref_start
                gene_end = gene['end'] + self.ref_start

            gene_stats.append({
                'name': gene['name'],
                'start': gene_start,
                'end': gene_end,
                'scaffold_local_start': gene['start'],  # Keep original for gene alignment
                'scaffold_local_end': gene['end'],
                'length': gene['length'],
                'strand': gene.get('strand', '+'),
                'quality_score': min(100, max(0, quality_score)),
                'coverage_pct': result.coverage_pct,
                'avg_identity': result.identity,
                'status': result.status,
                'product': gene.get('product', '')
            })

        return gene_stats

    def _detect_misassemblies(self, overlap_regions: List[OverlapRegion]) -> List[dict]:
        """
        Detect gaps, overlaps, and inversions between contigs using absolute reference coordinates.
        
        IMPORTANT: Only reports gaps/overlaps that fall WITHIN scaffold alignment regions.
        Gaps between separate scaffold alignment blocks are NOT reported as they are
        expected (covered by other scaffolds).
        """
        misassemblies = []

        # Collect all contig regions sorted by position (using absolute coordinates)
        all_contigs = []
        for region in overlap_regions:
            for contig in region.overlapping_contigs:
                all_contigs.append({
                    'start': contig['overlap_start'],  # Absolute reference position
                    'end': contig['overlap_end'],
                    'strand': contig['strand'],
                    'name': contig['contig_name']
                })

        all_contigs.sort(key=lambda x: x['start'])

        # Detect gaps and overlaps between consecutive contigs
        # But ONLY if the gap falls within a scaffold alignment region
        for i in range(len(all_contigs) - 1):
            curr = all_contigs[i]
            next_c = all_contigs[i + 1]

            gap = next_c['start'] - curr['end']
            gap_start = curr['end']
            gap_end = next_c['start']

            # Check if this gap falls within a scaffold region
            # A gap is only valid if both ends are within the same scaffold region
            gap_within_scaffold = False
            for region_start, region_end in self.scaffold_regions:
                if gap_start >= region_start and gap_end <= region_end:
                    gap_within_scaffold = True
                    break

            if gap > 100 and gap_within_scaffold:  # Gap within scaffold region
                misassemblies.append({
                    'type': 'gap',
                    'ref_start': curr['end'],
                    'ref_end': next_c['start'],
                    'size': gap
                })
            elif gap < -100:  # Overlap (always report)
                misassemblies.append({
                    'type': 'overlap',
                    'ref_start': next_c['start'],
                    'ref_end': curr['end'],
                    'size': abs(gap)
                })

        # Detect inversions (reverse strand alignments)
        for contig in all_contigs:
            if contig['strand'] == '-':
                misassemblies.append({
                    'type': 'inversion',
                    'ref_start': contig['start'],
                    'ref_end': contig['end'],
                    'size': contig['end'] - contig['start']
                })

        return misassemblies


class ComparisonCircularVisualizer:
    """
    Create circular plots for scaffold-contig comparison.

    Shows the scaffold with contig coverage visualization (matching original circular plots):
    - Ring 1 (outermost): Gene quality - scaffold genes colored by contig coverage
    - Ring 2 (middle): Alignment status - how well contigs cover the scaffold
    - Ring 3 (inner): Contig mapping - which contigs cover which parts
    """

    @staticmethod
    def create_interactive_circular_plot(scaffold_seq: ScaffoldSequence, output_file: Path,
                                          reference_seqid: str = None):
        """
        Create interactive circular plot using Plotly showing scaffold in reference context.

        The plot shows the full reference genome circle with the scaffold region highlighted.
        """
        print(f"  Creating interactive circular plot: {output_file}")

        fig = go.Figure()

        ref_length = scaffold_seq.length  # This is now the full reference length
        scaffold_start = scaffold_seq.ref_start
        scaffold_end = scaffold_seq.ref_end

        # Build coverage and strand tracking arrays (only within scaffold region)
        coverage = np.zeros(ref_length)
        strand_info = {}

        for aln in scaffold_seq.alignments:
            if aln['is_primary']:
                start = max(0, aln['ref_start'])
                end = min(ref_length, aln['ref_end'])
                for pos in range(start, end):
                    coverage[pos] += 1
                    if pos not in strand_info:
                        strand_info[pos] = []
                    strand_info[pos].append(aln['strand'])

        # ====================================================================
        # RING 1: GENE QUALITY (outer ring) - Scaffold genes colored by contig coverage
        # ====================================================================

        quality_colors = {
            'excellent': 'rgb(46, 204, 113)',   # Green
            'good': 'rgb(241, 196, 15)',         # Yellow
            'fair': 'rgb(230, 126, 34)',         # Orange
            'poor': 'rgb(231, 76, 60)'           # Red
        }

        # Check if we have gene stats to display
        if scaffold_seq.gene_stats:
            for quality_level in ['excellent', 'good', 'fair', 'poor']:
                theta = []
                r = []
                hover_text = []

                for gene in scaffold_seq.gene_stats:
                    quality = gene['quality_score']

                    if quality_level == 'excellent' and quality < 95:
                        continue
                    elif quality_level == 'good' and (quality < 85 or quality >= 95):
                        continue
                    elif quality_level == 'fair' and (quality < 70 or quality >= 85):
                        continue
                    elif quality_level == 'poor' and quality >= 70:
                        continue

                    # gene['start'] and gene['end'] are already absolute reference positions
                    start_angle = (gene['start'] / scaffold_seq.length) * 360
                    end_angle = (gene['end'] / scaffold_seq.length) * 360

                    arc_angles = np.linspace(start_angle, end_angle, max(2, int((end_angle - start_angle) / 2)))
                    theta.extend(arc_angles)
                    r.extend([0.85] * len(arc_angles))  # Outermost ring

                    hover = (f"Gene: {gene['name']}<br>"
                            f"Reference position: {gene['start']:,}-{gene['end']:,} bp<br>"
                            f"Length: {gene['length']:,} bp<br>"
                            f"Quality: {quality:.1f}/100<br>"
                            f"Contig coverage: {gene['coverage_pct']:.1f}%<br>"
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
        # RING 2: ALIGNMENT STATUS (middle ring) - How contigs cover the scaffold
        # Shows complete/duplicated/inverted/missing status
        # IMPORTANT: Only show status for positions within actual scaffold alignment regions
        # Gaps between separate scaffold alignment blocks are NOT marked as "missing"
        # ====================================================================

        status_colors = {
            'complete': 'rgb(46, 204, 113)',    # Green
            'duplicated': 'rgb(243, 156, 18)',  # Orange
            'inverted': 'rgb(231, 76, 60)',     # Red
            'missing': 'rgb(149, 165, 166)'     # Gray
        }

        # Build segments ONLY within actual scaffold alignment regions
        # This prevents marking regions between separate alignment blocks as "missing"
        segments = []
        
        for region_start, region_end in scaffold_seq.scaffold_regions:
            current_segment = None
            
            for pos in range(region_start, region_end):
                if pos >= len(coverage):
                    continue
                    
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

        # Create traces for alignment status (Ring 2)
        for status_type in ['complete', 'duplicated', 'inverted', 'missing']:
            status_segments = [s for s in segments if s['status'] == status_type]

            if not status_segments:
                continue

            theta = []
            r = []
            hover_text = []

            for seg in status_segments:
                start_angle = (seg['start'] / scaffold_seq.length) * 360
                end_angle = ((seg['end'] + 1) / scaffold_seq.length) * 360

                arc_angles = np.linspace(start_angle, end_angle, max(2, int((end_angle - start_angle) / 2)))
                theta.extend(arc_angles)
                r.extend([0.60] * len(arc_angles))  # Middle ring

                size = seg['end'] - seg['start'] + 1
                hover = f"Status: {status_type.capitalize()}<br>Reference position: {seg['start']:,}-{seg['end']:,} bp<br>Size: {size:,} bp"
                hover_text.extend([hover] * len(arc_angles))

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
        # RING 3: CONTIG MAPPING (inner ring)
        # ====================================================================

        # Get unique contigs and assign colors
        contigs_in_order = []
        for aln in sorted(scaffold_seq.alignments, key=lambda x: x['ref_start']):
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

            for aln in scaffold_seq.alignments:
                if not aln['is_primary'] or aln['query_name'] != contig_name:
                    continue

                start_angle = (aln['ref_start'] / scaffold_seq.length) * 360
                end_angle = (aln['ref_end'] / scaffold_seq.length) * 360

                arc_angles = np.linspace(start_angle, end_angle, max(2, int((end_angle - start_angle) / 2)))
                theta.extend(arc_angles)
                r.extend([0.35] * len(arc_angles))  # Ring 3 position (inner)

                hover = (f"Contig: {contig_name}<br>"
                        f"Reference position: {aln['ref_start']:,}-{aln['ref_end']:,} bp<br>"
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
                    legendgroup="Ring 3: Contigs",
                    legendgrouptitle=dict(text="<b>Ring 3: Overlapping Contigs</b>", font=dict(size=11))
                ))

        # ====================================================================
        # UPDATE LAYOUT
        # ====================================================================

        if scaffold_seq.length > 10000:
            size_str = f'{scaffold_seq.length/1000:.1f} kb'
        else:
            size_str = f'{scaffold_seq.length:,} bp'

        ref_info = f' (Ref: {reference_seqid})' if reference_seqid else ''

        fig.update_layout(
            title=dict(
                text=f'Scaffold: {scaffold_seq.seqid} ({size_str}){ref_info}',
                x=0.5,
                xanchor='center',
                font=dict(family='Arial, sans-serif')
            ),
            polar=dict(
                radialaxis=dict(visible=False, range=[0, 1]),
                angularaxis=dict(
                    tickmode='array',
                    tickvals=np.linspace(0, 360, 8, endpoint=False),
                    ticktext=[f'{int(pos/1000)}' if scaffold_seq.length > 10000 else f'{int(pos):,}'
                             for pos in np.linspace(0, scaffold_seq.length, 8, endpoint=False)],
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
                title=dict(text="<b>Rings</b>", font=dict(size=14, family='Arial, sans-serif')),
                tracegroupgap=20,
                itemsizing='constant',
                font=dict(size=10, family='Arial, sans-serif')
            ),
            font=dict(family='Arial, sans-serif'),
            width=1400,
            height=1000
        )

        fig.write_html(output_file)

        # Add help button
        InteractiveCircularVisualizer.add_help_button(output_file, plot_type='circular', gene_clicking_enabled=False)

    @staticmethod
    def create_static_circular_plot(scaffold_seq: ScaffoldSequence, output_file: Path,
                                     reference_seqid: str = None):
        """Create static circular plot using matplotlib showing scaffold in reference context."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.patches import Wedge

        print(f"  Creating static circular plot: {output_file}")

        fig = plt.figure(figsize=(16, 14))
        ax = fig.add_subplot(111, projection='polar')

        ref_length = scaffold_seq.length  # Full reference length
        scaffold_start = scaffold_seq.ref_start
        scaffold_end = scaffold_seq.ref_end

        # Build coverage array (only within scaffold region for efficiency)
        coverage = np.zeros(ref_length)
        strand_info = {}

        for aln in scaffold_seq.alignments:
            if aln['is_primary']:
                start = max(0, aln['ref_start'])
                end = min(ref_length, aln['ref_end'])
                for pos in range(start, end):
                    coverage[pos] += 1
                    if pos not in strand_info:
                        strand_info[pos] = []
                    strand_info[pos].append(aln['strand'])

        # Get contig colors
        contigs_in_order = []
        for aln in sorted(scaffold_seq.alignments, key=lambda x: x['ref_start']):
            if aln['is_primary'] and aln['query_name'] not in contigs_in_order:
                contigs_in_order.append(aln['query_name'])

        if len(contigs_in_order) <= 12:
            cmap = plt.cm.Set3
        elif len(contigs_in_order) <= 20:
            cmap = plt.cm.tab20
        else:
            cmap = plt.cm.hsv

        contig_colors = {c: cmap(i / max(len(contigs_in_order), 2))
                        for i, c in enumerate(contigs_in_order)}

        # Ring parameters (3 rings: gene quality, alignment status, contig mapping)
        # Matching original circular_visualizer.py ring positions
        quality_track_r = 0.85
        quality_track_width = 0.15
        status_track_r = 0.60
        status_track_width = 0.15
        contig_track_r = 0.35
        contig_track_width = 0.15

        # ====================================================================
        # RING 1: GENE QUALITY (outermost)
        # ====================================================================

        for gene in scaffold_seq.gene_stats:
            # Calculate angles - Wedge uses degrees with 0 at East (3 o'clock), counterclockwise
            # To match Plotly (and have 0 at top going clockwise), transform:
            # new_angle = 90 - (position/length * 360)
            # This puts 0 at North and goes clockwise
            start_angle_raw = (gene['start'] / ref_length) * 360
            end_angle_raw = (gene['end'] / ref_length) * 360
            
            # Transform to 0-at-North, clockwise
            start_angle = 90 - start_angle_raw
            end_angle = 90 - end_angle_raw
            
            # Wedge needs start < end, so swap if needed (since we flipped direction)
            if start_angle > end_angle:
                start_angle, end_angle = end_angle, start_angle

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
                         start_angle, end_angle,
                         width=quality_track_width, facecolor=color,
                         edgecolor='none', alpha=0.8,
                         transform=ax.transData._b)
            ax.add_artist(wedge)

        # ====================================================================
        # RING 2: ALIGNMENT STATUS (middle)
        # IMPORTANT: Only show status for positions within actual scaffold alignment regions
        # ====================================================================

        status_colors = {
            'complete': '#2ecc71',
            'duplicated': '#f39c12',
            'inverted': '#e74c3c',
            'missing': '#95a5a6'
        }

        segments = []

        # Only process positions within actual scaffold alignment regions
        for region_start, region_end in scaffold_seq.scaffold_regions:
            current_segment = None
            
            for pos in range(region_start, region_end):
                if pos >= len(coverage):
                    continue
                    
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

        for segment in segments:
            # Transform to 0-at-North, clockwise (same as Ring 1)
            start_angle_raw = (segment['start'] / ref_length) * 360
            end_angle_raw = ((segment['end'] + 1) / ref_length) * 360
            start_angle = 90 - start_angle_raw
            end_angle = 90 - end_angle_raw
            if start_angle > end_angle:
                start_angle, end_angle = end_angle, start_angle
            
            color = status_colors[segment['status']]

            wedge = Wedge((0, 0), status_track_r,
                         start_angle, end_angle,
                         width=status_track_width, facecolor=color,
                         edgecolor='none', alpha=0.8,
                         transform=ax.transData._b)
            ax.add_artist(wedge)

        # ====================================================================
        # RING 3: CONTIG MAPPING (inner)
        # ====================================================================

        for aln in scaffold_seq.alignments:
            if aln['is_primary']:
                contig_name = aln['query_name']
                color = contig_colors.get(contig_name, '#888888')

                # Transform to 0-at-North, clockwise (same as Ring 1 and 2)
                start_angle_raw = (aln['ref_start'] / ref_length) * 360
                end_angle_raw = (aln['ref_end'] / ref_length) * 360
                start_angle = 90 - start_angle_raw
                end_angle = 90 - end_angle_raw
                if start_angle > end_angle:
                    start_angle, end_angle = end_angle, start_angle

                wedge = Wedge((0, 0), contig_track_r,
                            start_angle, end_angle,
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

        n_ticks = 8
        tick_angles = np.linspace(0, 2*np.pi, n_ticks, endpoint=False)
        tick_positions_bp = np.linspace(0, scaffold_seq.length, n_ticks, endpoint=False)

        if scaffold_seq.length > 10000:
            tick_labels = [f'{int(pos/1000)}' for pos in tick_positions_bp]
            unit_label = 'kb'
        else:
            tick_labels = [f'{int(pos):,}' for pos in tick_positions_bp]
            unit_label = 'bp'

        ax.set_xticks(tick_angles)
        ax.set_xticklabels(tick_labels, fontsize=10)
        ax.set_yticks([])
        ax.grid(True, alpha=0.3, linewidth=0.5)

        if scaffold_seq.length > 10000:
            size_str = f'{scaffold_seq.length/1000:.1f} kb'
        else:
            size_str = f'{scaffold_seq.length:,} bp'

        ref_info = f'\nReference: {reference_seqid}' if reference_seqid else ''

        plt.title(f'Scaffold: {scaffold_seq.seqid}\n({size_str}){ref_info}\n\n'
                 f'Ring 1 (Outer): Gene Quality | Ring 2 (Middle): Alignment Status | Ring 3 (Inner): Contig Mapping',
                 fontsize=13, fontweight='bold', pad=20)

        ax.text(0.5, -0.15, f'Position ({unit_label})',
               transform=ax.transAxes, ha='center', fontsize=11)

        # ====================================================================
        # LEGENDS
        # ====================================================================

        # Legend 1: Gene Quality (outermost ring)
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

        # Legend 2: Alignment Status (middle ring)
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

        # Legend 3: Contig Mapping (inner ring)
        if len(contig_colors) <= 15:
            contig_legend_elements = [plt.Line2D([0], [0], marker='s', color='w',
                                         markerfacecolor=color, markersize=10,
                                         label=name, markeredgecolor='black', markeredgewidth=0.5)
                             for name, color in sorted(contig_colors.items())]
            ax.legend(handles=contig_legend_elements, loc='upper left',
                     bbox_to_anchor=(1.15, 0.30), title='Ring 3: Contigs',
                     fontsize=9, framealpha=0.9)

        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()


class ComparisonLinearVisualizer:
    """
    Create linear plots for scaffold-contig comparison.

    Uses the EXACT same structure as LinearVisualizer:
    - Track 1: Gene quality - scaffold genes colored by contig coverage
    - Track 2: Contig mapping - which contigs cover the scaffold
    - Track 3: Coverage depth - how many contigs cover each position
    - Track 4: Alignment identity - identity of contig alignments
    - Track 5: Misassemblies - gaps, overlaps, inversions between contigs
    """

    @staticmethod
    def create_linear_plot(scaffold_seq: ScaffoldSequence, output_file: Path,
                           reference_seqid: str = None):
        """Create static linear plot using matplotlib."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.patches import Rectangle
        from matplotlib.cm import ScalarMappable
        from matplotlib.colors import Normalize

        print(f"  Creating linear plot: {output_file}")

        fig = plt.figure(figsize=(16, 14))

        ref_info = f' (Ref: {reference_seqid})' if reference_seqid else ''
        fig.suptitle(f'Scaffold: {scaffold_seq.seqid}{ref_info}',
                    fontsize=14, fontweight='bold', y=0.98)

        # Define fixed positions for all plot areas
        plot_left = 0.08
        plot_width = 0.82

        track_heights = [0.22, 0.11, 0.11, 0.11, 0.11]
        track_spacing = 0.055

        bottom_positions = []
        current_bottom = 0.06
        for i in range(4, -1, -1):
            bottom_positions.insert(0, current_bottom)
            current_bottom += track_heights[i] + track_spacing

        ax1 = fig.add_axes([plot_left, bottom_positions[0], plot_width, track_heights[0]])
        ax2 = fig.add_axes([plot_left, bottom_positions[1], plot_width, track_heights[1]])
        ax3 = fig.add_axes([plot_left, bottom_positions[2], plot_width, track_heights[2]])
        ax4 = fig.add_axes([plot_left, bottom_positions[3], plot_width, track_heights[3]])
        ax5 = fig.add_axes([plot_left, bottom_positions[4], plot_width, track_heights[4]])

        # Get contig colors
        contigs_in_order = []
        for aln in sorted(scaffold_seq.alignments, key=lambda x: x['ref_start']):
            if aln['is_primary'] and aln['query_name'] not in contigs_in_order:
                contigs_in_order.append(aln['query_name'])

        if len(contigs_in_order) <= 12:
            cmap = plt.cm.Set3
        elif len(contigs_in_order) <= 20:
            cmap = plt.cm.tab20
        else:
            cmap = plt.cm.hsv

        contig_colors = {c: cmap(i / max(len(contigs_in_order), 2))
                        for i, c in enumerate(contigs_in_order)}

        # ====================================================================
        # TRACK 1: GENE QUALITY
        # ====================================================================

        ax1.set_title('Scaffold Gene Quality (by Contig Coverage)', fontsize=11, pad=10)

        for gene in scaffold_seq.gene_stats:
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

        ax1.set_xlim(0, scaffold_seq.length)
        ax1.set_ylim(0, 1)
        ax1.set_xlabel('Position (bp)', fontsize=10)
        ax1.set_ylabel('Quality', fontsize=10)
        ax1.set_yticks([])

        # Add colorbar
        cbar_left = plot_left + plot_width + 0.01
        cbar_width = 0.02
        cax = fig.add_axes([cbar_left, bottom_positions[0], cbar_width, track_heights[0]])

        cmap_cb = plt.cm.RdYlGn
        norm = Normalize(vmin=0, vmax=100)
        sm = ScalarMappable(cmap=cmap_cb, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, cax=cax, orientation='vertical')
        cbar.set_label('Quality Score', fontsize=9)

        # ====================================================================
        # TRACK 2: CONTIG MAPPING
        # ====================================================================

        ax2.set_title('Contig Mapping', fontsize=11, pad=10)

        for aln in scaffold_seq.alignments:
            if aln['is_primary']:
                y_pos = 0.5
                color = contig_colors.get(aln['query_name'], '#888888')

                rect = Rectangle((aln['ref_start'], y_pos - 0.3),
                               aln['ref_end'] - aln['ref_start'], 0.6,
                               facecolor=color, edgecolor='black', linewidth=0.5, alpha=0.7)
                ax2.add_patch(rect)

                width = aln['ref_end'] - aln['ref_start']
                if width > scaffold_seq.length * 0.05:
                    ax2.text(aln['ref_start'] + width/2, 0.5,
                            aln['query_name'], ha='center', va='center',
                            fontsize=7, fontweight='bold')

        ax2.set_xlim(0, scaffold_seq.length)
        ax2.set_ylim(0, 1)
        ax2.set_ylabel('Contigs', fontsize=10)
        ax2.set_yticks([])

        # ====================================================================
        # TRACK 3: COVERAGE
        # ====================================================================

        ax3.set_title('Contig Coverage Depth', fontsize=11, pad=10)

        coverage = np.zeros(scaffold_seq.length)
        for aln in scaffold_seq.alignments:
            if aln['is_primary']:
                start = max(0, aln['ref_start'])
                end = min(scaffold_seq.length, aln['ref_end'])
                coverage[start:end] += 1

        positions = np.arange(scaffold_seq.length)
        ax3.fill_between(positions, coverage, alpha=0.6, color='#3498db')
        ax3.set_xlim(0, scaffold_seq.length)
        ax3.set_ylabel('Coverage', fontsize=10)
        ax3.grid(True, alpha=0.3)

        # ====================================================================
        # TRACK 4: IDENTITY
        # ====================================================================

        ax4.set_title('Contig Alignment Identity (%)', fontsize=11, pad=10)

        for aln in scaffold_seq.alignments:
            if aln['is_primary']:
                ax4.plot([aln['ref_start'], aln['ref_end']],
                        [aln['identity'], aln['identity']],
                        linewidth=2, color='#9b59b6', alpha=0.6)

        ax4.set_xlim(0, scaffold_seq.length)
        ax4.set_ylim(70, 101)
        ax4.set_ylabel('Identity (%)', fontsize=10)
        ax4.axhline(y=95, color='green', linestyle='--', alpha=0.5, linewidth=1)
        ax4.axhline(y=90, color='orange', linestyle='--', alpha=0.5, linewidth=1)
        ax4.grid(True, alpha=0.3)

        # ====================================================================
        # TRACK 5: MISASSEMBLIES
        # ====================================================================

        ax5.set_title('Detected Issues (Gaps, Overlaps, Inversions)', fontsize=11, pad=10)

        y_pos = {'inversion': 0.7, 'gap': 0.4, 'overlap': 0.1}
        colors_mis = {'inversion': '#e74c3c', 'gap': '#f39c12', 'overlap': '#9b59b6'}

        for mis in scaffold_seq.misassemblies:
            mis_type = mis['type']
            rect = Rectangle((mis['ref_start'], y_pos[mis_type] - 0.05),
                           mis['ref_end'] - mis['ref_start'], 0.1,
                           facecolor=colors_mis[mis_type], edgecolor='black',
                           linewidth=0.5, alpha=0.7)
            ax5.add_patch(rect)

        ax5.set_xlim(0, scaffold_seq.length)
        ax5.set_ylim(0, 1)
        ax5.set_xlabel('Position (bp)', fontsize=10)
        ax5.set_yticks([0.1, 0.4, 0.7])
        ax5.set_yticklabels(['Overlap', 'Gap', 'Inversion'], fontsize=9)

        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()


class ComparisonInteractiveLinearVisualizer:
    """
    Create interactive linear plots for scaffold-contig comparison using Plotly.
    
    Uses the same structure as the standard InteractiveLinearVisualizer:
    - Track 1: Gene quality - scaffold/reference genes colored by contig coverage
    - Track 2: Contig mapping - which contigs cover the scaffold
    - Track 3: Coverage depth - how many contigs cover each position
    - Track 4: Alignment identity - identity of contig alignments
    - Track 5: Misassemblies - gaps, overlaps, inversions between contigs
    """

    @staticmethod
    def create_interactive_linear_plot(scaffold_seq: ScaffoldSequence, output_file: Path,
                                        reference_seqid: str = None):
        """Create interactive linear plot using Plotly."""
        import json
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
        
        print(f"  Creating interactive linear plot: {output_file}")
        
        # Create subplot figure with 5 tracks
        fig = make_subplots(
            rows=5, cols=1,
            shared_xaxes=True,
            vertical_spacing=0.02,
            row_heights=[0.25, 0.20, 0.20, 0.15, 0.20],
            subplot_titles=(
                'Gene Quality (by Contig Coverage)',
                'Contig Mapping',
                'Coverage Depth',
                'Alignment Identity',
                'Misassemblies'
            )
        )
        
        # ====================================================================
        # TRACK 1: GENE QUALITY
        # ====================================================================
        
        quality_colors = {
            'excellent': 'rgb(46, 204, 113)',
            'good': 'rgb(241, 196, 15)',
            'fair': 'rgb(230, 126, 34)',
            'poor': 'rgb(231, 76, 60)'
        }
        
        for quality_level in ['excellent', 'good', 'fair', 'poor']:
            x_positions = []
            widths = []
            heights = []
            hover_text = []
            
            for gene in scaffold_seq.gene_stats:
                quality = gene['quality_score']
                
                if quality_level == 'excellent' and quality < 95:
                    continue
                elif quality_level == 'good' and (quality < 85 or quality >= 95):
                    continue
                elif quality_level == 'fair' and (quality < 70 or quality >= 85):
                    continue
                elif quality_level == 'poor' and quality >= 70:
                    continue
                
                x_positions.append((gene['start'] + gene['end']) / 2)
                widths.append(gene['end'] - gene['start'])
                heights.append(quality)
                
                hover = (f"<b>{gene['name']}</b><br>"
                        f"Position: {gene['start']:,}-{gene['end']:,} bp<br>"
                        f"Length: {gene['length']:,} bp<br>"
                        f"Quality: {quality:.1f}/100<br>"
                        f"Coverage: {gene['coverage_pct']:.1f}%<br>"
                        f"Identity: {gene['avg_identity']:.1f}%<br>"
                        f"Status: {gene['status']}")
                hover_text.append(hover)
            
            if x_positions:
                fig.add_trace(go.Bar(
                    x=x_positions,
                    y=heights,
                    width=widths,
                    name=f'{quality_level.capitalize()}',
                    marker=dict(
                        color=quality_colors[quality_level],
                        line=dict(color='rgba(0,0,0,0.3)', width=0)
                    ),
                    hovertext=hover_text,
                    hoverinfo='text',
                    legendgroup='Track 1: Gene Quality',
                    legendgrouptitle=dict(text="<b>Track 1: Gene Quality</b>", font=dict(size=11)),
                    showlegend=True
                ), row=1, col=1)
        
        # ====================================================================
        # TRACK 2: CONTIG MAPPING
        # ====================================================================
        
        # Get unique contigs and assign colors
        contigs_in_order = []
        for aln in sorted(scaffold_seq.alignments, key=lambda x: x['ref_start']):
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
        
        # Group by contig
        contig_data = {}
        for aln in scaffold_seq.alignments:
            if aln['is_primary']:
                contig_name = aln['query_name']
                if contig_name not in contig_data:
                    contig_data[contig_name] = {'x': [], 'widths': [], 'hover': []}
                
                center_x = (aln['ref_start'] + aln['ref_end']) / 2
                width = aln['ref_end'] - aln['ref_start']
                
                contig_data[contig_name]['x'].append(center_x)
                contig_data[contig_name]['widths'].append(width)
                
                hover = (f"<b>{contig_name}</b><br>"
                        f"Ref: {aln['ref_start']:,}-{aln['ref_end']:,} bp<br>"
                        f"Contig: {aln['query_start']:,}-{aln['query_end']:,} bp<br>"
                        f"Length: {aln['ref_end'] - aln['ref_start']:,} bp<br>"
                        f"Identity: {aln['identity']:.2f}%<br>"
                        f"Strand: {aln['strand']}")
                contig_data[contig_name]['hover'].append(hover)
        
        # Add bar traces for each contig
        for contig_name, data in contig_data.items():
            color = contig_colors.get(contig_name, 'rgb(136, 136, 136)')
            fig.add_trace(go.Bar(
                x=data['x'],
                y=[1] * len(data['x']),
                width=data['widths'],
                name=contig_name,
                marker=dict(
                    color=color,
                    line=dict(color='rgba(0,0,0,0.4)', width=1)
                ),
                hovertext=data['hover'],
                hoverinfo='text',
                legendgroup='Track 2: Contig Mapping',
                legendgrouptitle=dict(text="<b>Track 2: Contig Mapping</b>", font=dict(size=11)),
                showlegend=True
            ), row=2, col=1)
        
        # ====================================================================
        # TRACK 3: COVERAGE
        # ====================================================================
        
        coverage = np.zeros(scaffold_seq.length)
        for aln in scaffold_seq.alignments:
            if aln['is_primary']:
                start = max(0, aln['ref_start'])
                end = min(scaffold_seq.length, aln['ref_end'])
                coverage[start:end] += 1
        
        # Downsample if genome is large for performance
        window_size = max(1, scaffold_seq.length // 10000)
        if window_size > 1:
            downsampled_coverage = []
            downsampled_positions = []
            for i in range(0, scaffold_seq.length, window_size):
                window_end = min(i + window_size, scaffold_seq.length)
                downsampled_coverage.append(np.mean(coverage[i:window_end]))
                downsampled_positions.append(i + window_size // 2)
        else:
            downsampled_coverage = coverage
            downsampled_positions = np.arange(scaffold_seq.length)
        
        fig.add_trace(go.Scatter(
            x=downsampled_positions,
            y=downsampled_coverage,
            mode='lines',
            name='Coverage',
            line=dict(color='rgb(52, 152, 219)', width=1.5),
            fill='tozeroy',
            fillcolor='rgba(52, 152, 219, 0.3)',
            hovertemplate='Position: %{x:,} bp<br>Coverage: %{y:.1f}x<extra></extra>',
            legendgroup='Track 3: Coverage Depth',
            legendgrouptitle=dict(text="<b>Track 3: Coverage Depth</b>", font=dict(size=11)),
            showlegend=True
        ), row=3, col=1)
        
        # ====================================================================
        # TRACK 4: IDENTITY
        # ====================================================================
        
        identity_x = []
        identity_widths = []
        identity_heights = []
        identity_hover = []
        
        for aln in scaffold_seq.alignments:
            if aln['is_primary']:
                center_x = (aln['ref_start'] + aln['ref_end']) / 2
                width = aln['ref_end'] - aln['ref_start']
                
                identity_x.append(center_x)
                identity_widths.append(width)
                identity_heights.append(aln['identity'])
                
                hover = (f"Position: {aln['ref_start']:,}-{aln['ref_end']:,} bp<br>"
                        f"Length: {width:,} bp<br>"
                        f"Identity: {aln['identity']:.2f}%<br>"
                        f"Contig: {aln['query_name']}")
                identity_hover.append(hover)
        
        if identity_x:
            fig.add_trace(go.Bar(
                x=identity_x,
                y=identity_heights,
                width=identity_widths,
                name='Identity',
                marker=dict(
                    color='rgb(155, 89, 182)',
                    line=dict(color='rgba(0,0,0,0.3)', width=0.5)
                ),
                hovertext=identity_hover,
                hoverinfo='text',
                legendgroup='Track 4: Alignment Identity',
                legendgrouptitle=dict(text="<b>Track 4: Alignment Identity</b>", font=dict(size=11)),
                showlegend=True
            ), row=4, col=1)
        
        # Add reference lines
        fig.add_hline(y=95, line=dict(color='green', dash='dash', width=1),
                     opacity=0.5, row=4, col=1)
        fig.add_hline(y=90, line=dict(color='orange', dash='dash', width=1),
                     opacity=0.5, row=4, col=1)
        
        # ====================================================================
        # TRACK 5: MISASSEMBLIES
        # ====================================================================
        
        y_pos_map = {'inversion': 0.85, 'gap': 0.6, 'overlap': 0.35}
        color_map = {
            'inversion': 'rgb(231, 76, 60)',
            'gap': 'rgb(243, 156, 18)',
            'overlap': 'rgb(155, 89, 182)'
        }
        
        # Add misassembly markers
        for mis_type in ['inversion', 'gap', 'overlap']:
            misassemblies = [m for m in scaffold_seq.misassemblies if m['type'] == mis_type]
            
            if misassemblies:
                x_positions = []
                widths = []
                hover_text = []
                
                for mis in misassemblies:
                    center_x = (mis['ref_start'] + mis['ref_end']) / 2
                    width = mis['ref_end'] - mis['ref_start']
                    
                    x_positions.append(center_x)
                    widths.append(width)
                    
                    hover = (f"<b>{mis_type.capitalize()}</b><br>"
                            f"Position: {mis['ref_start']:,}-{mis['ref_end']:,} bp<br>"
                            f"Length: {mis['size']:,} bp")
                    hover_text.append(hover)
                
                fig.add_trace(go.Bar(
                    x=x_positions,
                    y=[0.15] * len(x_positions),
                    width=widths,
                    base=[y_pos_map[mis_type]] * len(x_positions),
                    name=mis_type.capitalize(),
                    marker=dict(
                        color=color_map[mis_type],
                        line=dict(color='rgba(0,0,0,0.3)', width=1)
                    ),
                    hovertext=hover_text,
                    hoverinfo='text',
                    legendgroup='Track 5: Misassemblies',
                    legendgrouptitle=dict(text="<b>Track 5: Misassemblies</b>", font=dict(size=11)),
                    showlegend=True,
                    orientation='v'
                ), row=5, col=1)
        
        # Add contig indicator bars at bottom of track
        for contig_name, data in contig_data.items():
            color = contig_colors.get(contig_name, 'rgb(136, 136, 136)')
            
            for i, x in enumerate(data['x']):
                fig.add_trace(go.Bar(
                    x=[x],
                    y=[0.1],
                    width=[data['widths'][i]],
                    base=[0],
                    name=contig_name,
                    marker=dict(
                        color=color,
                        opacity=0.7,
                        line=dict(color='rgba(0,0,0,0.3)', width=0.5)
                    ),
                    hovertext=[f"<b>Contig: {contig_name}</b>"],
                    hoverinfo='text',
                    legendgroup='Contig Indicators',
                    showlegend=False,
                    orientation='v'
                ), row=5, col=1)
        
        # ====================================================================
        # CONFIGURE LAYOUT
        # ====================================================================
        
        ref_info = f' (Ref: {reference_seqid})' if reference_seqid else ''
        if scaffold_seq.length > 10000:
            size_str = f'{scaffold_seq.length/1000:.1f} kb'
        else:
            size_str = f'{scaffold_seq.length:,} bp'
        
        fig.update_layout(
            title=dict(
                text=f'Scaffold: {scaffold_seq.seqid}{ref_info} ({size_str})',
                x=0.5,
                xanchor='center',
                font=dict(size=16, family='Arial, sans-serif')
            ),
            hovermode='x unified',
            height=1200,
            width=1600,
            showlegend=True,
            legend=dict(
                orientation="v",
                yanchor="top",
                y=1,
                xanchor="left",
                x=1.01,
                font=dict(size=10, family='Arial, sans-serif')
            ),
            font=dict(family='Arial, sans-serif')
        )
        
        # Update axes
        fig.update_xaxes(title_text="Position (bp)", row=5, col=1, rangeslider=dict(visible=True))
        fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='lightgray')
        
        fig.update_yaxes(title_text="Quality", row=1, col=1, range=[0, 100])
        fig.update_yaxes(title_text="Contigs", row=2, col=1, range=[0, 1.1], showticklabels=False)
        fig.update_yaxes(title_text="Coverage", row=3, col=1)
        fig.update_yaxes(title_text="Identity (%)", row=4, col=1, range=[85, 100])
        fig.update_yaxes(title_text="Type", row=5, col=1, range=[0, 1.05],
                        tickvals=[0.05, 0.42, 0.67, 0.92],
                        ticktext=['Contig', 'Overlap', 'Gap', 'Inversion'])
        
        # Configure zoom and pan
        fig.update_xaxes(
            autorange=True,
            range=[0, scaffold_seq.length],
            constrain="domain"
        )
        
        # Save HTML with config
        config = {
            'scrollZoom': True,
            'displayModeBar': True,
            'displaylogo': False,
            'modeBarButtonsToAdd': ['drawline', 'drawopenpath', 'eraseshape'],
            'toImageButtonOptions': {
                'format': 'png',
                'filename': f'{scaffold_seq.seqid}_interactive_linear',
                'height': 1200,
                'width': 1600,
                'scale': 2
            }
        }
        
        fig.write_html(output_file, config=config)
        
        # Add help button
        ComparisonInteractiveLinearVisualizer._add_help_button(output_file)
    
    @staticmethod
    def _add_help_button(html_file):
        """Add help button to the interactive linear plot."""
        
        # Read the HTML file
        with open(html_file, 'r', encoding='utf-8') as f:
            html_content = f.read()
        
        help_html = '''
        <style>
        .help-button {
            position: fixed;
            top: 10px;
            right: 10px;
            z-index: 1000;
            background: #3498db;
            color: white;
            border: none;
            border-radius: 50%;
            width: 30px;
            height: 30px;
            cursor: pointer;
            font-size: 18px;
            box-shadow: 0 2px 5px rgba(0,0,0,0.2);
        }
        .help-button:hover {
            background: #2980b9;
        }
        .help-popup {
            display: none;
            position: fixed;
            top: 50px;
            right: 10px;
            z-index: 999;
            background: white;
            padding: 15px;
            border-radius: 8px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.2);
            max-width: 350px;
            font-family: Arial, sans-serif;
            font-size: 13px;
        }
        .help-popup h3 {
            margin-top: 0;
            color: #2c3e50;
        }
        .help-popup ul {
            padding-left: 20px;
        }
        .help-popup li {
            margin: 8px 0;
        }
        </style>
        <button class="help-button" onclick="toggleHelp()">?</button>
        <div class="help-popup" id="helpPopup">
            <h3>Navigation Guide</h3>
            <ul>
                <li><b>Zoom:</b> Use scroll wheel or drag the range slider at the bottom</li>
                <li><b>Pan:</b> Click and drag on the plot</li>
                <li><b>Reset:</b> Double-click to reset view</li>
                <li><b>Hover:</b> Move mouse over elements for details</li>
            </ul>
            <h4>Track Information</h4>
            <ul>
                <li><b>Track 1:</b> Gene quality (green=excellent, yellow=good, orange=fair, red=poor)</li>
                <li><b>Track 2:</b> Contig positions on reference</li>
                <li><b>Track 3:</b> Coverage depth</li>
                <li><b>Track 4:</b> Alignment identity (%)</li>
                <li><b>Track 5:</b> Issues (gaps, inversions, overlaps)</li>
            </ul>
        </div>
        <script>
        function toggleHelp() {
            var popup = document.getElementById('helpPopup');
            popup.style.display = popup.style.display === 'none' ? 'block' : 'none';
        }
        </script>
        '''
        
        # Insert before closing body tag
        html_content = html_content.replace('</body>', help_html + '\n</body>')
        
        with open(html_file, 'w', encoding='utf-8') as f:
            f.write(html_content)


class ComparisonIndexGenerator:
    """Generate index pages linking to all scaffold comparison views."""

    def __init__(self, comparison_summary: dict, output_dir: Path,
                 gene_summary: Optional[dict] = None):
        self.summary = comparison_summary
        self.output_dir = output_dir
        self.gene_summary = gene_summary

    def generate_index_html(self) -> Path:
        """Generate main index HTML page."""
        gene_section = ""
        if self.gene_summary:
            gs = self.gene_summary
            gene_section = f'''
        <h2>Gene Integrity Summary</h2>
        <div class="summary-grid">
            <div class="stat-card">
                <div class="stat-value">{gs.get('total_genes', 0)}</div>
                <div class="stat-label">Total Genes</div>
            </div>
            <div class="stat-card complete">
                <div class="stat-value">{gs['status_counts'].get('complete', 0)}</div>
                <div class="stat-label">Complete</div>
            </div>
            <div class="stat-card split">
                <div class="stat-value">{gs['status_counts'].get('split', 0)}</div>
                <div class="stat-label">Split</div>
            </div>
            <div class="stat-card truncated">
                <div class="stat-value">{gs['status_counts'].get('truncated', 0)}</div>
                <div class="stat-label">Truncated</div>
            </div>
            <div class="stat-card missing">
                <div class="stat-value">{gs['status_counts'].get('missing', 0)}</div>
                <div class="stat-label">Missing</div>
            </div>
        </div>
'''

        html_content = f'''<!DOCTYPE html>
<html>
<head>
    <title>Scaffold-Contig Comparison - GenomeViz</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 0; padding: 20px; background: #f5f5f5; }}
        .container {{ max-width: 1400px; margin: 0 auto; background: white; padding: 30px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
        h1 {{ color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 15px; }}
        h2 {{ color: #34495e; margin-top: 30px; }}
        .summary-grid {{ display: flex; flex-wrap: wrap; gap: 15px; margin: 20px 0; }}
        .stat-card {{ background: #ecf0f1; padding: 20px; border-radius: 8px; text-align: center; min-width: 120px; }}
        .stat-card .stat-value {{ font-size: 28px; font-weight: bold; color: #3498db; }}
        .stat-card .stat-label {{ font-size: 12px; color: #7f8c8d; margin-top: 5px; }}
        .stat-card.complete .stat-value {{ color: #2ecc71; }}
        .stat-card.split .stat-value {{ color: #f1c40f; }}
        .stat-card.truncated .stat-value {{ color: #e67e22; }}
        .stat-card.missing .stat-value {{ color: #e74c3c; }}
        table {{ width: 100%; border-collapse: collapse; margin: 20px 0; }}
        th, td {{ padding: 12px 15px; text-align: left; border-bottom: 1px solid #eee; }}
        th {{ background: #3498db; color: white; }}
        tr:hover {{ background: #f8f9fa; }}
        .btn {{ display: inline-block; padding: 6px 14px; text-decoration: none; border-radius: 4px; font-size: 13px; margin: 2px; }}
        .btn-circular {{ background: #9b59b6; color: white; }}
        .btn-circular:hover {{ background: #8e44ad; }}
        .btn-linear {{ background: #3498db; color: white; }}
        .btn-linear:hover {{ background: #2980b9; }}
        .btn-png {{ background: #27ae60; color: white; }}
        .btn-png:hover {{ background: #219a52; }}
        footer {{ margin-top: 40px; padding-top: 20px; border-top: 1px solid #eee; color: #7f8c8d; font-size: 12px; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>Scaffold-Contig Comparison Results</h1>

        <div class="summary-grid">
            <div class="stat-card">
                <div class="stat-value">{self.summary['num_scaffolds']}</div>
                <div class="stat-label">Scaffolds</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{self.summary['num_contigs']}</div>
                <div class="stat-label">Contigs</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{self.summary['num_total_overlap_regions']}</div>
                <div class="stat-label">Overlap Regions</div>
            </div>
        </div>

        {gene_section}

        <h2>Scaffold Comparisons</h2>
        <table>
            <thead>
                <tr>
                    <th>Scaffold</th>
                    <th>Overlapping Contigs</th>
                    <th>Total Overlap</th>
                    <th>Views</th>
                </tr>
            </thead>
            <tbody>
'''

        for scaffold_name, stats in self.summary['scaffold_stats'].items():
            safe_name = "".join(c if c.isalnum() or c in '-_' else '_' for c in scaffold_name)

            html_content += f'''
                <tr>
                    <td><strong>{scaffold_name}</strong></td>
                    <td>{stats['num_overlapping_contigs']}</td>
                    <td>{stats['total_overlap_bp']:,} bp</td>
                    <td>
                        <a href="{safe_name}/{safe_name}_circular.html" class="btn btn-circular">Circular</a>
                        <a href="{safe_name}/{safe_name}_linear.png" class="btn btn-linear">Linear</a>
                        <a href="{safe_name}/{safe_name}_circular.png" class="btn btn-png">PNG</a>
                    </td>
                </tr>
'''

        html_content += '''
            </tbody>
        </table>

        <footer>
            Generated by GenomeViz - Scaffold-Contig Comparison Mode<br>
            <a href="https://github.com/Aaron-Thiel/GenomeViz">https://github.com/Aaron-Thiel/GenomeViz</a>
        </footer>
    </div>
</body>
</html>
'''

        index_file = self.output_dir / 'index.html'
        with open(index_file, 'w', encoding='utf-8') as f:
            f.write(html_content)

        return index_file


class ComparisonGeneAlignmentVisualizer:
    """
    Create detailed gene alignment visualizations for scaffold-contig comparison.

    Shows scaffold genes with contig coverage, highlighting:
    - Which contigs cover which parts of the gene
    - Gaps between contig coverage
    - Parts of the gene only present in scaffold (not covered by contigs)
    """

    @staticmethod
    def create_gene_alignment_html(gene_result, scaffold_seq: ScaffoldSequence,
                                    output_dir: Path) -> str:
        """
        Create a detailed alignment visualization for a scaffold gene.

        Args:
            gene_result: GeneComparisonResult with coverage analysis
            scaffold_seq: ScaffoldSequence containing alignment data
            output_dir: Directory to save the gene alignment HTML

        Returns:
            Path to the generated HTML file
        """
        gene = gene_result.gene
        safe_gene_name = "".join(c if c.isalnum() or c in ('-', '_') else '_' for c in gene['name'])
        unique_name = f"{safe_gene_name}_{gene['start']}-{gene['end']}"
        output_file = output_dir / f"{unique_name}_alignment.html"

        gene_start = gene['start']
        gene_end = gene['end']
        gene_length = gene_end - gene_start

        # Build coverage map showing which contigs cover which positions
        coverage_map = []  # List of (contig_name, start, end) for each position
        position_contigs = [[] for _ in range(gene_length)]

        for cov_seq in gene_result.covering_sequences:
            local_start = cov_seq.get('gene_local_start', 0)
            local_end = cov_seq.get('gene_local_end', gene_length)
            contig_name = cov_seq.get('name', 'unknown')

            for pos in range(max(0, local_start), min(gene_length, local_end)):
                position_contigs[pos].append({
                    'name': contig_name,
                    'identity': cov_seq.get('identity', 0),
                    'strand': cov_seq.get('strand', '+')
                })

        # Generate visual representation
        html_content = ComparisonGeneAlignmentVisualizer._generate_html(
            gene, gene_result, position_contigs, scaffold_seq
        )

        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html_content)

        # Add help button
        InteractiveCircularVisualizer.add_help_button(output_file, plot_type='gene', gene_clicking_enabled=False)

        return str(output_file)

    @staticmethod
    def _generate_html(gene: dict, gene_result, position_contigs: List,
                       scaffold_seq: ScaffoldSequence) -> str:
        """Generate HTML content for gene alignment visualization."""
        gene_length = gene['end'] - gene['start']

        # Get unique contigs covering this gene
        all_contigs = set()
        for pos_contigs in position_contigs:
            for c in pos_contigs:
                all_contigs.add(c['name'])

        # Assign colors to contigs
        contig_colors = {}
        color_palette = [
            '#3498db', '#e74c3c', '#2ecc71', '#9b59b6', '#f39c12',
            '#1abc9c', '#e67e22', '#34495e', '#16a085', '#c0392b'
        ]
        for i, contig in enumerate(sorted(all_contigs)):
            contig_colors[contig] = color_palette[i % len(color_palette)]

        # Build coverage visualization
        coverage_html = []
        block_size = 100  # Show in blocks of 100bp

        for block_start in range(0, gene_length, block_size):
            block_end = min(block_start + block_size, gene_length)
            block_html = f'<div class="block"><div class="block-header">Position {gene["start"] + block_start:,} - {gene["start"] + block_end:,}</div>'

            # Create a track for each contig
            contig_tracks = {c: [] for c in all_contigs}
            gap_track = []

            for pos in range(block_start, block_end):
                pos_coverage = position_contigs[pos]

                if not pos_coverage:
                    gap_track.append(('gap', pos))
                else:
                    gap_track.append(('covered', pos))
                    for c in pos_coverage:
                        contig_tracks[c['name']].append((pos, c['identity'], c['strand']))

            # Render gap track (scaffold regions not covered by contigs)
            gap_cells = []
            for status, pos in gap_track:
                if status == 'gap':
                    gap_cells.append('<span class="gap" title="Not covered by contigs">-</span>')
                else:
                    gap_cells.append('<span class="covered">|</span>')

            block_html += f'<div class="track scaffold-track"><span class="track-label">Scaffold:</span>{"".join(gap_cells)}</div>'

            # Render contig tracks
            for contig_name in sorted(all_contigs):
                contig_cells = []
                track_data = {item[0]: (item[1], item[2]) for item in contig_tracks[contig_name]}

                for pos in range(block_start, block_end):
                    if pos in track_data:
                        identity, strand = track_data[pos]
                        color = contig_colors[contig_name]
                        strand_char = '>' if strand == '+' else '<'
                        contig_cells.append(
                            f'<span class="contig-base" style="background-color: {color};" '
                            f'title="{contig_name}: {identity:.1f}% identity, strand: {strand}">{strand_char}</span>'
                        )
                    else:
                        contig_cells.append('<span class="no-coverage">.</span>')

                block_html += f'<div class="track contig-track"><span class="track-label" style="color: {contig_colors[contig_name]};">{contig_name[:15]}:</span>{"".join(contig_cells)}</div>'

            block_html += '</div>'
            coverage_html.append(block_html)

        # Build summary statistics
        covered_positions = sum(1 for pc in position_contigs if pc)
        gap_positions = gene_length - covered_positions

        # Find gap regions
        gaps = gene_result.gaps
        gap_summary = ""
        if gaps:
            gap_items = [f"<li>{start:,} - {end:,} ({end-start:,} bp)</li>" for start, end in gaps]
            gap_summary = f"<h3>Gap Regions (Not Covered by Contigs)</h3><ul>{''.join(gap_items)}</ul>"

        html = f'''<!DOCTYPE html>
<html>
<head>
    <title>{gene['name']} - Scaffold Gene Alignment</title>
    <style>
        body {{ font-family: 'Consolas', 'Monaco', monospace; margin: 20px; background: #f5f5f5; }}
        .container {{ max-width: 1400px; margin: 0 auto; background: white; padding: 30px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
        h1 {{ color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 15px; }}
        h2 {{ color: #34495e; margin-top: 30px; }}
        h3 {{ color: #7f8c8d; }}
        .stats {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin: 20px 0; }}
        .stat {{ background: #ecf0f1; padding: 15px; border-radius: 8px; text-align: center; }}
        .stat-value {{ font-size: 24px; font-weight: bold; color: #3498db; }}
        .stat-label {{ font-size: 12px; color: #7f8c8d; margin-top: 5px; }}
        .stat.good .stat-value {{ color: #2ecc71; }}
        .stat.warning .stat-value {{ color: #f39c12; }}
        .stat.bad .stat-value {{ color: #e74c3c; }}
        .block {{ margin: 15px 0; padding: 10px; background: #fafafa; border-radius: 4px; overflow-x: auto; }}
        .block-header {{ font-size: 12px; color: #7f8c8d; margin-bottom: 5px; }}
        .track {{ white-space: nowrap; margin: 2px 0; font-size: 11px; }}
        .track-label {{ display: inline-block; width: 120px; font-weight: bold; text-align: right; padding-right: 10px; }}
        .scaffold-track {{ background: #f0f0f0; padding: 2px 0; }}
        .gap {{ color: #e74c3c; font-weight: bold; }}
        .covered {{ color: #2ecc71; }}
        .contig-base {{ padding: 0 1px; color: white; font-weight: bold; }}
        .no-coverage {{ color: #ccc; }}
        .legend {{ margin: 20px 0; padding: 15px; background: #f9f9f9; border-radius: 8px; }}
        .legend-item {{ display: inline-block; margin: 5px 15px 5px 0; }}
        .legend-color {{ display: inline-block; width: 20px; height: 14px; vertical-align: middle; margin-right: 5px; border-radius: 2px; }}
        ul {{ line-height: 1.8; }}
        footer {{ margin-top: 40px; padding-top: 20px; border-top: 1px solid #eee; color: #7f8c8d; font-size: 12px; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>Gene: {gene['name']}</h1>
        <p><strong>Scaffold:</strong> {scaffold_seq.seqid} | <strong>Position:</strong> {gene['start']:,} - {gene['end']:,} bp | <strong>Length:</strong> {gene_length:,} bp</p>

        <div class="stats">
            <div class="stat {'good' if gene_result.coverage_pct >= 95 else 'warning' if gene_result.coverage_pct >= 70 else 'bad'}">
                <div class="stat-value">{gene_result.coverage_pct:.1f}%</div>
                <div class="stat-label">Contig Coverage</div>
            </div>
            <div class="stat">
                <div class="stat-value">{gene_result.identity:.1f}%</div>
                <div class="stat-label">Average Identity</div>
            </div>
            <div class="stat">
                <div class="stat-value">{gene_result.num_covering_sequences}</div>
                <div class="stat-label">Covering Contigs</div>
            </div>
            <div class="stat {'bad' if gap_positions > 0 else 'good'}">
                <div class="stat-value">{gap_positions:,}</div>
                <div class="stat-label">Gap Positions (bp)</div>
            </div>
        </div>

        <p><strong>Status:</strong> <span style="color: {'#2ecc71' if gene_result.status == 'complete' else '#f39c12' if gene_result.status == 'split' else '#e67e22' if gene_result.status == 'truncated' else '#e74c3c'};">{gene_result.status.upper()}</span> - {gene_result.details}</p>

        {gap_summary}

        <h2>Contig Legend</h2>
        <div class="legend">
            {''.join(f'<div class="legend-item"><span class="legend-color" style="background-color: {color};"></span>{name}</div>' for name, color in sorted(contig_colors.items()))}
            <div class="legend-item"><span class="legend-color" style="background-color: #e74c3c;"></span>Gap (not covered)</div>
        </div>

        <h2>Coverage Alignment</h2>
        <p><em>Scaffold track shows "-" for positions not covered by any contig, "|" for covered positions.<br>
        Contig tracks show ">" for forward strand, "<" for reverse strand coverage.</em></p>

        {''.join(coverage_html)}

        <footer>
            Generated by GenomeViz - Scaffold-Contig Comparison Mode<br>
            <a href="https://github.com/Aaron-Thiel/GenomeViz">https://github.com/Aaron-Thiel/GenomeViz</a>
        </footer>
    </div>
</body>
</html>'''

        return html

    @staticmethod
    def add_circular_click_handler(html_file: Path, scaffold_seq: ScaffoldSequence,
                                    gene_results: List):
        """
        Add JavaScript to circular plot HTML to make genes clickable.

        Args:
            html_file: Path to the HTML file
            scaffold_seq: ScaffoldSequence with gene information
            gene_results: List of GeneComparisonResult objects
        """
        if not gene_results:
            return

        # Build gene mapping with unique filenames
        gene_list = []
        for result in gene_results:
            gene = result.gene
            safe_gene_name = "".join(c if c.isalnum() or c in ('-', '_') else '_' for c in gene['name'])
            unique_name = f"{safe_gene_name}_{gene['start']}-{gene['end']}"
            link_path = f"gene_alignments/{unique_name}_alignment.html"

            gene_list.append({
                'name': gene['name'],
                'start': gene['start'],
                'end': gene['end'],
                'link': link_path
            })

        # Read the HTML file
        with open(html_file, 'r', encoding='utf-8') as f:
            html_content = f.read()

        # Create JavaScript code for click handling in polar coordinates
        js_code = '''
        <script>
        // Gene positions and links
        const geneList = ''' + json.dumps(gene_list, indent=2) + ''';
        const refLength = ''' + str(scaffold_seq.length) + ''';

        // Add click handler to the plot
        document.addEventListener('DOMContentLoaded', function() {
            const plotDiv = document.getElementsByClassName('plotly-graph-div')[0];

            if (plotDiv) {
                plotDiv.on('plotly_click', function(data) {
                    if (data.points && data.points.length > 0) {
                        const point = data.points[0];

                        // Check if click is on Ring 1 (gene quality)
                        if (point.data.name && point.data.name.includes('Gene')) {
                            // Get theta (angle) from click
                            let theta = point.theta;

                            // Convert angle to genome position
                            while (theta < 0) theta += 360;
                            while (theta >= 360) theta -= 360;

                            const position = (theta / 360) * refLength;

                            // Find which gene was clicked
                            for (let i = 0; i < geneList.length; i++) {
                                const gene = geneList[i];
                                if (position >= gene.start && position <= gene.end) {
                                    // Open gene alignment in new tab
                                    window.open(gene.link, '_blank');
                                    console.log('Opening alignment for:', gene.name, 'at position', Math.round(position));
                                    break;
                                }
                            }
                        }
                    }
                });
            }
        });
        </script>
        '''

        # Insert JavaScript before closing body tag
        html_content = html_content.replace('</body>', js_code + '\n</body>')

        # Write back to file
        with open(html_file, 'w', encoding='utf-8') as f:
            f.write(html_content)
