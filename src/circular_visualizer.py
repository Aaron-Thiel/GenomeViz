"""
CircularVisualizer - Static circular genome visualization

Creates publication-quality static circular plots using matplotlib showing:
- Gene quality (outer ring)
- Alignment status (middle ring)
- Contig mapping (inner ring)

GitHub: https://github.com/Aaron-Thiel/GenomeViz
License: MIT
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge


class CircularVisualizer:
    """
    Create static circular visualizations using matplotlib for publication.
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
