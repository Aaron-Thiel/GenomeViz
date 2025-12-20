# Implementation Plan: Multi-Assembly Comparison Mode

## Overview

Add the ability to compare scaffolds and contigs directly using the reference genome as a coordinate anchor. Both assemblies are aligned to the reference once, and overlapping regions are identified via reference coordinates.

## Critical Design Decisions

### 1. GFF3 Files for All Assemblies

**Recommendation**: Accept separate GFF3 files for scaffolds and contigs.

```
--gff            Reference GFF3 (required, existing)
--scaffold-gff   Scaffold GFF3 (optional)
--contig-gff     Contig GFF3 (optional)
```

**Visualization approach**:
- In comparison views, show genes from the scaffold/contig GFFs at their native coordinates
- This allows users to see "does contig gene X align with scaffold gene Y?"
- When GFFs are not provided for scaffold/contig, those tracks simply won't have gene annotations

### 2. Coordinate Mapping Strategy

```
REFERENCE (anchor)
    |
    +-- Scaffold alignments --> {scaffold_name: [(ref_start, ref_end, scaffold_start, scaffold_end), ...]}
    |
    +-- Contig alignments   --> {contig_name: [(ref_start, ref_end, contig_start, contig_end), ...]}

For each scaffold region [ref_start, ref_end]:
    Find all contigs where their ref_start..ref_end overlaps
```

### 3. Output Structure

```
output/
├── (existing outputs unchanged)
├── assembly_comparison/
│   ├── index.html                           # Main navigation page
│   ├── comparison_summary.json              # Overlap statistics
│   ├── scaffold_A/
│   │   ├── scaffold_A_vs_contigs_linear.html
│   │   ├── scaffold_A_vs_contigs_circular.html
│   │   ├── scaffold_A_vs_contigs_circular.png
│   │   └── scaffold_A_gene_comparison.csv   # Gene-level comparison if GFFs provided
│   ├── scaffold_B/
│   │   └── ...
│   └── overview_matrix.html                 # Matrix showing all scaffold-contig overlaps
```

---

## New CLI Arguments

```python
# In genomeViz.py argument parser, add:

# Multi-assembly comparison mode
parser.add_argument('--scaffolds',
    help='Scaffold assembly (FASTA format) for scaffold-contig comparison')
parser.add_argument('--contigs',
    help='Contig assembly (FASTA format) for scaffold-contig comparison')
parser.add_argument('--scaffold-gff',
    help='Gene annotations for scaffolds (GFF3 format)')
parser.add_argument('--contig-gff',
    help='Gene annotations for contigs (GFF3 format)')
parser.add_argument('--no-comparison-circular', action='store_true',
    help='Skip circular plots in assembly comparison')
parser.add_argument('--no-comparison-linear', action='store_true',
    help='Skip linear plots in assembly comparison')
```

**Backward Compatibility**: When `--scaffolds` and `--contigs` are not provided, the tool behaves exactly as before. The new comparison mode is entirely opt-in.

---

## New Modules

### 1. `src/assembly_comparison.py`

Core logic for comparing two assemblies via reference coordinates.

```python
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
"""

from dataclasses import dataclass
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
        overlapping_contigs: List of contigs that overlap this region
        reference_seqid: Reference sequence ID where overlap occurs
    """
    scaffold_name: str
    scaffold_ref_start: int
    scaffold_ref_end: int
    scaffold_local_start: int
    scaffold_local_end: int
    reference_seqid: str
    overlapping_contigs: List[dict]  # Each dict contains contig alignment info


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
                        overlap_pct = (overlap_bp / smaller_len) * 100

                        if overlap_pct >= self.min_overlap_pct:
                            overlapping_contigs.append({
                                'contig_name': contig_aln['query_name'],
                                'contig_ref_start': contig_ref_start,
                                'contig_ref_end': contig_ref_end,
                                'contig_local_start': contig_aln['query_start'],
                                'contig_local_end': contig_aln['query_end'],
                                'overlap_start': overlap_start,
                                'overlap_end': overlap_end,
                                'overlap_bp': overlap_bp,
                                'overlap_pct': overlap_pct,
                                'strand': contig_aln['strand'],
                                'identity': contig_aln.get('identity', 0)
                            })

            if overlapping_contigs:
                overlap_regions.append(OverlapRegion(
                    scaffold_name=scaffold_name,
                    scaffold_ref_start=scaffold_ref_start,
                    scaffold_ref_end=scaffold_ref_end,
                    scaffold_local_start=scaffold_aln['query_start'],
                    scaffold_local_end=scaffold_aln['query_end'],
                    reference_seqid=ref_seqid,
                    overlapping_contigs=sorted(overlapping_contigs,
                                              key=lambda x: x['contig_ref_start'])
                ))

        return overlap_regions

    def get_all_scaffolds(self) -> List[str]:
        """Get list of all scaffold names."""
        return list(self._scaffold_index.keys())

    def compute_comparison_summary(self) -> dict:
        """
        Compute summary statistics for all scaffold-contig overlaps.

        Returns:
            Dictionary with summary statistics
        """
        all_overlaps = []
        scaffold_stats = {}

        for scaffold_name in self.get_all_scaffolds():
            regions = self.find_overlapping_contigs(scaffold_name)

            total_contigs = set()
            total_overlap_bp = 0

            for region in regions:
                for contig in region.overlapping_contigs:
                    total_contigs.add(contig['contig_name'])
                    total_overlap_bp += contig['overlap_bp']

            scaffold_stats[scaffold_name] = {
                'num_regions': len(regions),
                'num_overlapping_contigs': len(total_contigs),
                'total_overlap_bp': total_overlap_bp,
                'contig_names': list(total_contigs)
            }

            all_overlaps.extend(regions)

        return {
            'num_scaffolds': len(scaffold_stats),
            'num_total_overlaps': len(all_overlaps),
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
                 contig_gff: Optional[str] = None):
        """
        Initialize multi-assembly GFF parser.

        Args:
            reference_gff: Path to reference GFF3 file
            scaffold_gff: Path to scaffold GFF3 file
            contig_gff: Path to contig GFF3 file
        """
        self.reference_genes = {}
        self.scaffold_genes = {}
        self.contig_genes = {}

        if reference_gff:
            self.reference_genes = self._parse_gff(reference_gff)
        if scaffold_gff:
            self.scaffold_genes = self._parse_gff(scaffold_gff)
        if contig_gff:
            self.contig_genes = self._parse_gff(contig_gff)

    def _parse_gff(self, gff_file: str) -> Dict[str, List[dict]]:
        """Parse GFF3 file and return genes organized by sequence ID."""
        genes_by_seq = defaultdict(list)

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

                genes_by_seq[seqid].append({
                    'seqid': seqid,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'name': gene_name,
                    'type': feature_type,
                    'length': end - start + 1
                })

        # Sort by position
        for seqid in genes_by_seq:
            genes_by_seq[seqid].sort(key=lambda x: x['start'])

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
```

---

### 2. `src/comparison_visualizer.py`

Visualizations for scaffold-contig comparisons.

```python
"""
Comparison visualization module for scaffold vs contig views.

Creates:
- Interactive linear plots showing scaffold aligned with overlapping contigs
- Interactive circular plots showing scaffold-contig relationships
- Static publication-quality versions

Classes:
    ComparisonLinearVisualizer: Linear multi-track comparison views
    ComparisonCircularVisualizer: Circular comparison views
    ComparisonIndexGenerator: Generate index/navigation pages
"""

import json
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from pathlib import Path
from typing import List, Dict, Optional

from .assembly_comparison import OverlapRegion, MultiAssemblyGFFParser
from .circular_visualizer import CircularVisualizer


class ComparisonLinearVisualizer:
    """
    Create interactive linear plots comparing scaffold to overlapping contigs.

    Track layout:
    1. Reference region (genes if GFF provided)
    2. Scaffold track (with genes if scaffold GFF provided)
    3-N. Contig tracks (one per overlapping contig, with genes if contig GFF provided)
    N+1. Alignment identity between scaffold and contigs (via reference coordinates)
    """

    def __init__(self, overlap_regions: List[OverlapRegion],
                 reference_sequence: dict,
                 scaffold_sequences: dict,
                 contig_sequences: dict,
                 gff_parser: Optional[MultiAssemblyGFFParser] = None):
        """
        Initialize comparison linear visualizer.

        Args:
            overlap_regions: List of OverlapRegion objects for one scaffold
            reference_sequence: Dict with reference sequence data
            scaffold_sequences: Dict mapping scaffold names to sequence data
            contig_sequences: Dict mapping contig names to sequence data
            gff_parser: MultiAssemblyGFFParser for gene annotations
        """
        self.overlap_regions = overlap_regions
        self.reference_sequence = reference_sequence
        self.scaffold_sequences = scaffold_sequences
        self.contig_sequences = contig_sequences
        self.gff_parser = gff_parser

        # Determine scaffold name from regions
        self.scaffold_name = overlap_regions[0].scaffold_name if overlap_regions else None

    def create_comparison_plot(self, output_file: Path) -> None:
        """
        Create interactive linear comparison plot.

        Args:
            output_file: Path to output HTML file
        """
        if not self.overlap_regions:
            return

        # Collect all unique contigs
        all_contigs = set()
        for region in self.overlap_regions:
            for contig in region.overlapping_contigs:
                all_contigs.add(contig['contig_name'])
        all_contigs = sorted(all_contigs)

        # Calculate number of tracks needed
        # Track 1: Reference region
        # Track 2: Scaffold
        # Tracks 3-N: Contigs (one per contig)
        # Track N+1: Overlap/identity comparison
        num_contig_tracks = len(all_contigs)
        total_tracks = 3 + num_contig_tracks  # ref + scaffold + contigs + comparison

        # Calculate row heights
        base_height = 0.15
        remaining = 1.0 - (base_height * 3)  # Reserve for ref, scaffold, comparison
        contig_height = remaining / max(num_contig_tracks, 1)

        row_heights = [base_height]  # Reference
        row_heights.append(base_height)  # Scaffold
        row_heights.extend([contig_height] * num_contig_tracks)  # Contigs
        row_heights.append(base_height)  # Comparison

        # Normalize to sum to 1
        total = sum(row_heights)
        row_heights = [h/total for h in row_heights]

        # Create subplot titles
        titles = ['Reference Region', f'Scaffold: {self.scaffold_name}']
        titles.extend([f'Contig: {c}' for c in all_contigs])
        titles.append('Scaffold-Contig Identity')

        fig = make_subplots(
            rows=total_tracks, cols=1,
            shared_xaxes=True,
            vertical_spacing=0.02,
            row_heights=row_heights,
            subplot_titles=titles
        )

        # Determine x-axis range (reference coordinates)
        all_ref_starts = [r.scaffold_ref_start for r in self.overlap_regions]
        all_ref_ends = [r.scaffold_ref_end for r in self.overlap_regions]
        x_min = min(all_ref_starts)
        x_max = max(all_ref_ends)

        # Add 5% padding
        padding = (x_max - x_min) * 0.05
        x_min = max(0, x_min - padding)
        x_max = x_max + padding

        # ================================================================
        # TRACK 1: REFERENCE REGION
        # ================================================================
        self._add_reference_track(fig, row=1, x_min=x_min, x_max=x_max)

        # ================================================================
        # TRACK 2: SCAFFOLD
        # ================================================================
        self._add_scaffold_track(fig, row=2)

        # ================================================================
        # TRACKS 3-N: CONTIGS
        # ================================================================
        for idx, contig_name in enumerate(all_contigs):
            self._add_contig_track(fig, contig_name, row=3+idx)

        # ================================================================
        # FINAL TRACK: COMPARISON/IDENTITY
        # ================================================================
        self._add_comparison_track(fig, row=total_tracks)

        # ================================================================
        # LAYOUT
        # ================================================================
        ref_seqid = self.overlap_regions[0].reference_seqid
        fig.update_layout(
            title=dict(
                text=f'Scaffold-Contig Comparison: {self.scaffold_name}<br>'
                     f'<sub>Reference region: {ref_seqid}:{int(x_min):,}-{int(x_max):,}</sub>',
                x=0.5,
                xanchor='center'
            ),
            height=200 + (150 * total_tracks),
            width=1600,
            showlegend=True,
            legend=dict(orientation='v', x=1.02, y=1)
        )

        # Update x-axis
        fig.update_xaxes(
            title_text='Reference Position (bp)',
            row=total_tracks, col=1,
            range=[x_min, x_max],
            rangeslider=dict(visible=True)
        )

        # Save
        config = {
            'scrollZoom': True,
            'displayModeBar': True,
            'displaylogo': False
        }
        fig.write_html(output_file, config=config)

    def _add_reference_track(self, fig, row: int, x_min: float, x_max: float):
        """Add reference region track with genes."""
        # Add background bar for reference region
        fig.add_trace(go.Bar(
            x=[(x_min + x_max) / 2],
            y=[1],
            width=[x_max - x_min],
            marker=dict(color='rgba(200, 200, 200, 0.3)'),
            name='Reference Region',
            showlegend=False,
            hoverinfo='skip'
        ), row=row, col=1)

        # Add genes if GFF available
        if self.gff_parser and self.overlap_regions:
            ref_seqid = self.overlap_regions[0].reference_seqid
            genes = self.gff_parser.get_genes_in_range(
                'reference', ref_seqid, int(x_min), int(x_max)
            )

            for gene in genes:
                color = 'rgb(52, 152, 219)'  # Blue for reference genes
                center = (gene['start'] + gene['end']) / 2
                width = gene['end'] - gene['start']

                fig.add_trace(go.Bar(
                    x=[center],
                    y=[0.8],
                    width=[width],
                    marker=dict(color=color, line=dict(color='black', width=0.5)),
                    name=gene['name'],
                    hovertext=f"<b>{gene['name']}</b><br>"
                             f"Position: {gene['start']:,}-{gene['end']:,}<br>"
                             f"Strand: {gene['strand']}",
                    hoverinfo='text',
                    showlegend=False
                ), row=row, col=1)

        fig.update_yaxes(range=[0, 1.1], showticklabels=False, row=row, col=1)

    def _add_scaffold_track(self, fig, row: int):
        """Add scaffold track with alignment blocks."""
        colors = ['rgb(46, 204, 113)', 'rgb(39, 174, 96)']  # Greens

        for idx, region in enumerate(self.overlap_regions):
            color = colors[idx % len(colors)]
            center = (region.scaffold_ref_start + region.scaffold_ref_end) / 2
            width = region.scaffold_ref_end - region.scaffold_ref_start

            fig.add_trace(go.Bar(
                x=[center],
                y=[1],
                width=[width],
                marker=dict(color=color, line=dict(color='black', width=1)),
                name=f'{region.scaffold_name}',
                hovertext=f"<b>{region.scaffold_name}</b><br>"
                         f"Ref: {region.scaffold_ref_start:,}-{region.scaffold_ref_end:,}<br>"
                         f"Local: {region.scaffold_local_start:,}-{region.scaffold_local_end:,}",
                hoverinfo='text',
                legendgroup='scaffold',
                showlegend=(idx == 0)
            ), row=row, col=1)

        # Add scaffold genes if available
        if self.gff_parser and self.scaffold_name:
            # Note: Scaffold genes are at scaffold coordinates, need to map to reference
            # For simplicity, we show them at their scaffold positions if scaffold == reference coord mapping exists
            # This is a simplification - full implementation would need coordinate transformation
            pass

        fig.update_yaxes(range=[0, 1.1], showticklabels=False, row=row, col=1)

    def _add_contig_track(self, fig, contig_name: str, row: int):
        """Add track for a single contig."""
        contig_colors = ['rgb(155, 89, 182)', 'rgb(142, 68, 173)']  # Purples

        trace_idx = 0
        for region in self.overlap_regions:
            for contig in region.overlapping_contigs:
                if contig['contig_name'] == contig_name:
                    color = contig_colors[trace_idx % len(contig_colors)]
                    center = (contig['contig_ref_start'] + contig['contig_ref_end']) / 2
                    width = contig['contig_ref_end'] - contig['contig_ref_start']

                    fig.add_trace(go.Bar(
                        x=[center],
                        y=[1],
                        width=[width],
                        marker=dict(color=color, line=dict(color='black', width=1)),
                        name=contig_name,
                        hovertext=f"<b>{contig_name}</b><br>"
                                 f"Ref: {contig['contig_ref_start']:,}-{contig['contig_ref_end']:,}<br>"
                                 f"Identity: {contig['identity']:.1f}%<br>"
                                 f"Strand: {contig['strand']}<br>"
                                 f"Overlap: {contig['overlap_bp']:,} bp ({contig['overlap_pct']:.1f}%)",
                        hoverinfo='text',
                        legendgroup=f'contig_{contig_name}',
                        showlegend=(trace_idx == 0)
                    ), row=row, col=1)
                    trace_idx += 1

        fig.update_yaxes(range=[0, 1.1], showticklabels=False, row=row, col=1)

    def _add_comparison_track(self, fig, row: int):
        """Add scaffold-contig comparison track showing alignment relationships."""
        # Show overlap regions as colored bars with identity
        for region in self.overlap_regions:
            for contig in region.overlapping_contigs:
                identity = contig['identity']

                # Color based on identity
                if identity >= 99:
                    color = 'rgb(46, 204, 113)'  # Green
                elif identity >= 95:
                    color = 'rgb(241, 196, 15)'  # Yellow
                elif identity >= 90:
                    color = 'rgb(230, 126, 34)'  # Orange
                else:
                    color = 'rgb(231, 76, 60)'   # Red

                center = (contig['overlap_start'] + contig['overlap_end']) / 2
                width = contig['overlap_bp']

                fig.add_trace(go.Bar(
                    x=[center],
                    y=[identity],
                    width=[width],
                    marker=dict(color=color, line=dict(color='rgba(0,0,0,0.3)', width=0.5)),
                    name=f"Overlap: {contig['contig_name']}",
                    hovertext=f"<b>Overlap Region</b><br>"
                             f"Scaffold: {region.scaffold_name}<br>"
                             f"Contig: {contig['contig_name']}<br>"
                             f"Position: {contig['overlap_start']:,}-{contig['overlap_end']:,}<br>"
                             f"Length: {contig['overlap_bp']:,} bp<br>"
                             f"Identity: {identity:.2f}%",
                    hoverinfo='text',
                    showlegend=False
                ), row=row, col=1)

        fig.update_yaxes(
            title_text='Identity %',
            range=[85, 100],
            row=row, col=1
        )

        # Add reference lines
        fig.add_hline(y=99, line=dict(color='green', dash='dash', width=1),
                     opacity=0.5, row=row, col=1)
        fig.add_hline(y=95, line=dict(color='orange', dash='dash', width=1),
                     opacity=0.5, row=row, col=1)


class ComparisonCircularVisualizer:
    """
    Create circular plots showing scaffold-contig relationships.

    Ring layout (outer to inner):
    1. Reference region with genes
    2. Scaffold alignment
    3. Contig alignments (stacked)
    4. Overlap indicators
    """

    def __init__(self, overlap_regions: List[OverlapRegion],
                 gff_parser: Optional[MultiAssemblyGFFParser] = None):
        """
        Initialize comparison circular visualizer.

        Args:
            overlap_regions: List of OverlapRegion objects for one scaffold
            gff_parser: MultiAssemblyGFFParser for gene annotations
        """
        self.overlap_regions = overlap_regions
        self.gff_parser = gff_parser
        self.scaffold_name = overlap_regions[0].scaffold_name if overlap_regions else None

    def create_circular_plot(self, output_file: Path, interactive: bool = True):
        """
        Create circular comparison plot.

        Args:
            output_file: Path to output file (.html for interactive, .png for static)
            interactive: Whether to create interactive Plotly version
        """
        if not self.overlap_regions:
            return

        if interactive:
            self._create_interactive_circular(output_file)
        else:
            self._create_static_circular(output_file)

    def _create_interactive_circular(self, output_file: Path):
        """Create interactive circular plot using Plotly."""
        fig = go.Figure()

        # Determine angular range based on reference coordinates
        all_ref_starts = [r.scaffold_ref_start for r in self.overlap_regions]
        all_ref_ends = [r.scaffold_ref_end for r in self.overlap_regions]
        ref_min = min(all_ref_starts)
        ref_max = max(all_ref_ends)
        ref_length = ref_max - ref_min

        # Ring radii
        r_ref = 0.95       # Outer: Reference
        r_scaffold = 0.80  # Middle: Scaffold
        r_contig = 0.65    # Inner: Contigs
        r_overlap = 0.50   # Innermost: Overlap indicators

        ring_width = 0.08

        def pos_to_angle(pos):
            """Convert reference position to angle (radians)."""
            return ((pos - ref_min) / ref_length) * 2 * np.pi - np.pi/2

        def add_arc(r_inner, r_outer, start_pos, end_pos, color, name, hover):
            """Add an arc segment to the figure."""
            theta_start = pos_to_angle(start_pos)
            theta_end = pos_to_angle(end_pos)

            # Create arc points
            n_points = max(20, int((theta_end - theta_start) * 50))
            theta = np.linspace(theta_start, theta_end, n_points)

            # Outer arc
            x_outer = r_outer * np.cos(theta)
            y_outer = r_outer * np.sin(theta)

            # Inner arc (reversed)
            x_inner = r_inner * np.cos(theta[::-1])
            y_inner = r_inner * np.sin(theta[::-1])

            # Combine to form closed shape
            x = np.concatenate([x_outer, x_inner, [x_outer[0]]])
            y = np.concatenate([y_outer, y_inner, [y_outer[0]]])

            fig.add_trace(go.Scatter(
                x=x, y=y,
                fill='toself',
                fillcolor=color,
                line=dict(color='rgba(0,0,0,0.3)', width=0.5),
                name=name,
                hovertext=hover,
                hoverinfo='text',
                mode='lines'
            ))

        # ================================================================
        # RING 1: REFERENCE
        # ================================================================
        add_arc(r_ref - ring_width, r_ref, ref_min, ref_max,
               'rgba(200, 200, 200, 0.5)', 'Reference Region',
               f'Reference: {ref_min:,}-{ref_max:,}')

        # ================================================================
        # RING 2: SCAFFOLD
        # ================================================================
        for region in self.overlap_regions:
            add_arc(r_scaffold - ring_width, r_scaffold,
                   region.scaffold_ref_start, region.scaffold_ref_end,
                   'rgb(46, 204, 113)', region.scaffold_name,
                   f'<b>{region.scaffold_name}</b><br>'
                   f'Ref: {region.scaffold_ref_start:,}-{region.scaffold_ref_end:,}')

        # ================================================================
        # RING 3: CONTIGS
        # ================================================================
        contig_colors = ['rgb(155, 89, 182)', 'rgb(52, 152, 219)',
                        'rgb(230, 126, 34)', 'rgb(231, 76, 60)']
        color_idx = 0

        for region in self.overlap_regions:
            for contig in region.overlapping_contigs:
                color = contig_colors[color_idx % len(contig_colors)]
                add_arc(r_contig - ring_width, r_contig,
                       contig['contig_ref_start'], contig['contig_ref_end'],
                       color, contig['contig_name'],
                       f"<b>{contig['contig_name']}</b><br>"
                       f"Ref: {contig['contig_ref_start']:,}-{contig['contig_ref_end']:,}<br>"
                       f"Identity: {contig['identity']:.1f}%")
                color_idx += 1

        # ================================================================
        # RING 4: OVERLAP INDICATORS
        # ================================================================
        for region in self.overlap_regions:
            for contig in region.overlapping_contigs:
                identity = contig['identity']
                if identity >= 99:
                    color = 'rgb(46, 204, 113)'
                elif identity >= 95:
                    color = 'rgb(241, 196, 15)'
                else:
                    color = 'rgb(230, 126, 34)'

                add_arc(r_overlap - ring_width/2, r_overlap,
                       contig['overlap_start'], contig['overlap_end'],
                       color, 'Overlap',
                       f"Overlap: {contig['overlap_bp']:,} bp<br>"
                       f"Identity: {identity:.1f}%")

        # ================================================================
        # LAYOUT
        # ================================================================
        ref_seqid = self.overlap_regions[0].reference_seqid
        fig.update_layout(
            title=dict(
                text=f'Scaffold-Contig Comparison: {self.scaffold_name}<br>'
                     f'<sub>Reference: {ref_seqid}</sub>',
                x=0.5
            ),
            showlegend=True,
            width=900,
            height=900,
            xaxis=dict(scaleanchor='y', scaleratio=1, showgrid=False,
                      zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)
        )

        # Add ring labels
        fig.add_annotation(x=0, y=r_ref+0.05, text='Reference', showarrow=False)
        fig.add_annotation(x=0, y=r_scaffold+0.05, text='Scaffold', showarrow=False)
        fig.add_annotation(x=0, y=r_contig+0.05, text='Contigs', showarrow=False)

        config = {'scrollZoom': True, 'displayModeBar': True, 'displaylogo': False}
        fig.write_html(output_file, config=config)

    def _create_static_circular(self, output_file: Path):
        """Create static circular plot using matplotlib."""
        import matplotlib.pyplot as plt
        from matplotlib.patches import Wedge
        from matplotlib.collections import PatchCollection

        fig, ax = plt.subplots(figsize=(10, 10))

        # Similar logic to interactive but using matplotlib
        # Implementation follows same pattern as _create_interactive_circular
        # but uses Wedge patches instead of Plotly traces

        # [Detailed matplotlib implementation would go here]
        # For brevity, saving basic version

        ax.set_xlim(-1.2, 1.2)
        ax.set_ylim(-1.2, 1.2)
        ax.set_aspect('equal')
        ax.axis('off')

        ax.set_title(f'Scaffold-Contig Comparison: {self.scaffold_name}')

        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()


class ComparisonIndexGenerator:
    """
    Generate index pages linking to all scaffold comparison views.
    """

    def __init__(self, comparison_summary: dict, output_dir: Path):
        """
        Initialize index generator.

        Args:
            comparison_summary: Summary from AssemblyComparator.compute_comparison_summary()
            output_dir: Base output directory for comparison files
        """
        self.summary = comparison_summary
        self.output_dir = output_dir

    def generate_index_html(self) -> Path:
        """
        Generate main index HTML page.

        Returns:
            Path to generated index.html
        """
        html_content = '''<!DOCTYPE html>
<html>
<head>
    <title>Assembly Comparison - GenomeViz</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; background: #f5f5f5; }
        .container { max-width: 1200px; margin: 0 auto; background: white; padding: 30px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }
        h1 { color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }
        h2 { color: #34495e; margin-top: 30px; }
        .summary { background: #ecf0f1; padding: 20px; border-radius: 5px; margin: 20px 0; }
        .summary-stat { display: inline-block; margin-right: 30px; }
        .summary-stat .value { font-size: 24px; font-weight: bold; color: #3498db; }
        .summary-stat .label { font-size: 12px; color: #7f8c8d; }
        table { width: 100%; border-collapse: collapse; margin: 20px 0; }
        th, td { padding: 12px; text-align: left; border-bottom: 1px solid #ddd; }
        th { background: #3498db; color: white; }
        tr:hover { background: #f5f5f5; }
        .link-btn { background: #3498db; color: white; padding: 5px 15px; text-decoration: none; border-radius: 3px; margin-right: 5px; }
        .link-btn:hover { background: #2980b9; }
        .link-btn.circular { background: #9b59b6; }
        .link-btn.circular:hover { background: #8e44ad; }
    </style>
</head>
<body>
    <div class="container">
        <h1>Assembly Comparison Results</h1>

        <div class="summary">
            <div class="summary-stat">
                <div class="value">''' + str(self.summary['num_scaffolds']) + '''</div>
                <div class="label">Scaffolds Analyzed</div>
            </div>
            <div class="summary-stat">
                <div class="value">''' + str(self.summary['num_total_overlaps']) + '''</div>
                <div class="label">Total Overlap Regions</div>
            </div>
        </div>

        <h2>Scaffold Comparisons</h2>
        <table>
            <thead>
                <tr>
                    <th>Scaffold</th>
                    <th>Overlapping Contigs</th>
                    <th>Total Overlap (bp)</th>
                    <th>Views</th>
                </tr>
            </thead>
            <tbody>
'''

        for scaffold_name, stats in self.summary['scaffold_stats'].items():
            safe_name = scaffold_name.replace(' ', '_')
            html_content += f'''
                <tr>
                    <td><strong>{scaffold_name}</strong></td>
                    <td>{stats['num_overlapping_contigs']}</td>
                    <td>{stats['total_overlap_bp']:,}</td>
                    <td>
                        <a href="{safe_name}/{safe_name}_vs_contigs_linear.html" class="link-btn">Linear</a>
                        <a href="{safe_name}/{safe_name}_vs_contigs_circular.html" class="link-btn circular">Circular</a>
                    </td>
                </tr>
'''

        html_content += '''
            </tbody>
        </table>

        <h2>Contig Overlap Details</h2>
        <table>
            <thead>
                <tr>
                    <th>Scaffold</th>
                    <th>Overlapping Contigs</th>
                </tr>
            </thead>
            <tbody>
'''

        for scaffold_name, stats in self.summary['scaffold_stats'].items():
            contigs = ', '.join(stats['contig_names'][:5])
            if len(stats['contig_names']) > 5:
                contigs += f' (+{len(stats["contig_names"]) - 5} more)'

            html_content += f'''
                <tr>
                    <td>{scaffold_name}</td>
                    <td>{contigs}</td>
                </tr>
'''

        html_content += '''
            </tbody>
        </table>

        <p style="color: #7f8c8d; font-size: 12px; margin-top: 30px;">
            Generated by GenomeViz - Assembly Comparison Mode
        </p>
    </div>
</body>
</html>
'''

        index_file = self.output_dir / 'index.html'
        with open(index_file, 'w') as f:
            f.write(html_content)

        return index_file
```

---

## Main Pipeline Changes

### Updates to `genomeViz.py`

Add a new step after existing pipeline (Step 11):

```python
# After Step 10, add:

# Step 11: Assembly comparison (if enabled)
if args.scaffolds and args.contigs:
    print("\n" + "="*70)
    print("[11/11] ASSEMBLY COMPARISON")
    print("="*70)

    from src.assembly_comparison import AssemblyComparator, MultiAssemblyGFFParser
    from src.comparison_visualizer import (
        ComparisonLinearVisualizer,
        ComparisonCircularVisualizer,
        ComparisonIndexGenerator
    )

    # Step 11a: Align scaffolds to reference
    print("\n  [11a] Aligning scaffolds to reference...")
    scaffold_aligner = GenomeAligner(reference_to_use, args.scaffolds, preset=args.preset)
    scaffold_aligner.load_reference()
    scaffold_alignments = scaffold_aligner.align()

    # Step 11b: Align contigs to reference
    print("\n  [11b] Aligning contigs to reference...")
    contig_aligner = GenomeAligner(reference_to_use, args.contigs, preset=args.preset)
    contig_aligner.load_reference()
    contig_alignments = contig_aligner.align()

    # Step 11c: Parse additional GFF files
    print("\n  [11c] Parsing assembly GFF files...")
    multi_gff = MultiAssemblyGFFParser(
        reference_gff=args.gff,
        scaffold_gff=args.scaffold_gff if hasattr(args, 'scaffold_gff') else None,
        contig_gff=args.contig_gff if hasattr(args, 'contig_gff') else None
    )

    # Step 11d: Find overlapping regions
    print("\n  [11d] Finding scaffold-contig overlaps...")
    comparator = AssemblyComparator(scaffold_alignments, contig_alignments)
    comparison_summary = comparator.compute_comparison_summary()

    print(f"  Found {comparison_summary['num_scaffolds']} scaffolds with overlapping contigs")

    # Step 11e: Create comparison output directory
    comparison_dir = output_dir / 'assembly_comparison'
    comparison_dir.mkdir(exist_ok=True)

    # Step 11f: Generate visualizations for each scaffold
    print("\n  [11e] Generating comparison visualizations...")

    for scaffold_name in comparator.get_all_scaffolds():
        overlap_regions = comparator.find_overlapping_contigs(scaffold_name)

        if not overlap_regions:
            continue

        # Create scaffold output directory
        safe_name = scaffold_name.replace(' ', '_')
        scaffold_dir = comparison_dir / safe_name
        scaffold_dir.mkdir(exist_ok=True)

        print(f"    Processing {scaffold_name}...")

        # Create linear visualization
        if not args.no_comparison_linear:
            linear_viz = ComparisonLinearVisualizer(
                overlap_regions,
                reference_sequence=aligner.reference_sequences,
                scaffold_sequences={},  # Load if needed
                contig_sequences={},    # Load if needed
                gff_parser=multi_gff
            )
            linear_viz.create_comparison_plot(
                scaffold_dir / f'{safe_name}_vs_contigs_linear.html'
            )

        # Create circular visualization
        if not args.no_comparison_circular:
            circular_viz = ComparisonCircularVisualizer(
                overlap_regions,
                gff_parser=multi_gff
            )
            circular_viz.create_circular_plot(
                scaffold_dir / f'{safe_name}_vs_contigs_circular.html',
                interactive=True
            )
            # Also create static version
            circular_viz.create_circular_plot(
                scaffold_dir / f'{safe_name}_vs_contigs_circular.png',
                interactive=False
            )

    # Step 11f: Generate index page
    print("\n  [11f] Generating comparison index...")
    index_gen = ComparisonIndexGenerator(comparison_summary, comparison_dir)
    index_file = index_gen.generate_index_html()
    print(f"  Index page: {index_file}")

    # Save comparison summary JSON
    import json
    summary_file = comparison_dir / 'comparison_summary.json'
    with open(summary_file, 'w') as f:
        json.dump(comparison_summary, f, indent=2)
```

---

## Testing Strategy

### Unit Tests

```python
# tests/test_assembly_comparison.py

def test_overlap_detection():
    """Test that overlapping regions are correctly identified."""
    scaffold_alns = {
        'chr1': [
            {'query_name': 'scaffold_1', 'ref_start': 0, 'ref_end': 10000,
             'query_start': 0, 'query_end': 10000, 'is_primary': True}
        ]
    }
    contig_alns = {
        'chr1': [
            {'query_name': 'contig_1', 'ref_start': 5000, 'ref_end': 15000,
             'query_start': 0, 'query_end': 10000, 'is_primary': True,
             'strand': '+', 'identity': 99.5}
        ]
    }

    comparator = AssemblyComparator(scaffold_alns, contig_alns)
    regions = comparator.find_overlapping_contigs('scaffold_1')

    assert len(regions) == 1
    assert regions[0].overlapping_contigs[0]['contig_name'] == 'contig_1'
    assert regions[0].overlapping_contigs[0]['overlap_bp'] == 5000


def test_no_overlap():
    """Test that non-overlapping regions are not reported."""
    scaffold_alns = {
        'chr1': [
            {'query_name': 'scaffold_1', 'ref_start': 0, 'ref_end': 5000,
             'query_start': 0, 'query_end': 5000, 'is_primary': True}
        ]
    }
    contig_alns = {
        'chr1': [
            {'query_name': 'contig_1', 'ref_start': 10000, 'ref_end': 15000,
             'query_start': 0, 'query_end': 5000, 'is_primary': True,
             'strand': '+', 'identity': 99.5}
        ]
    }

    comparator = AssemblyComparator(scaffold_alns, contig_alns)
    regions = comparator.find_overlapping_contigs('scaffold_1')

    assert len(regions) == 0 or len(regions[0].overlapping_contigs) == 0
```

### Integration Tests

```bash
# Test with example data
python genomeViz.py \
    --reference examples/reference.fasta \
    --scaffolds examples/scaffolds.fasta \
    --contigs examples/contigs.fasta \
    --gff examples/reference.gff \
    --scaffold-gff examples/scaffolds.gff \
    --contig-gff examples/contigs.gff \
    --output test_comparison/
```

---

## Implementation Order

1. **Phase 1**: Core comparison logic
   - [ ] Create `src/assembly_comparison.py` with `AssemblyComparator` class
   - [ ] Implement coordinate overlap detection
   - [ ] Add unit tests

2. **Phase 2**: Linear visualization
   - [ ] Create `ComparisonLinearVisualizer` class
   - [ ] Implement multi-track layout
   - [ ] Add gene annotation support

3. **Phase 3**: Circular visualization
   - [ ] Create `ComparisonCircularVisualizer` class
   - [ ] Implement interactive Plotly version
   - [ ] Implement static matplotlib version

4. **Phase 4**: Index and navigation
   - [ ] Create `ComparisonIndexGenerator` class
   - [ ] Generate summary statistics
   - [ ] Create navigable index page

5. **Phase 5**: Integration
   - [ ] Add CLI arguments to `genomeViz.py`
   - [ ] Integrate comparison pipeline
   - [ ] Add GFF parsing for scaffolds/contigs

6. **Phase 6**: Testing and documentation
   - [ ] Integration tests with example data
   - [ ] Update README and USAGE.md
   - [ ] Update CHANGELOG.md

---

## Example Usage

```bash
# Basic comparison (scaffolds vs contigs, using reference as anchor)
python genomeViz.py \
    --reference ref.fasta \
    --scaffolds scaffolds.fasta \
    --contigs contigs.fasta \
    --gff ref.gff \
    --output results/

# With gene annotations for all assemblies
python genomeViz.py \
    --reference ref.fasta \
    --scaffolds scaffolds.fasta \
    --contigs contigs.fasta \
    --gff ref.gff \
    --scaffold-gff scaffolds.gff \
    --contig-gff contigs.gff \
    --output results/

# Skip circular comparison plots
python genomeViz.py \
    --reference ref.fasta \
    --scaffolds scaffolds.fasta \
    --contigs contigs.fasta \
    --gff ref.gff \
    --output results/ \
    --no-comparison-circular

# Standard mode (unchanged behavior)
python genomeViz.py \
    --reference ref.fasta \
    --assembly assembly.fasta \
    --gff ref.gff \
    --output results/
```
