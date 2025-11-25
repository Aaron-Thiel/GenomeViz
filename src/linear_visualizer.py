"""
LinearVisualizer - Static linear genome visualization

Creates publication-quality static linear plots using matplotlib showing:
- Gene-level quality track
- Contig mapping track
- Coverage depth track
- Alignment identity track
- Misassemblies track

GitHub: https://github.com/Aaron-Thiel/GenomeViz
License: MIT
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

# Import CircularVisualizer for consistent contig coloring
from .circular_visualizer import CircularVisualizer


class LinearVisualizer:
    """
    Create linear visualizations for detailed genome analysis.
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

        # Create figure with tight_layout disabled so we can control positions manually
        fig = plt.figure(figsize=(16, 14))
        
        seq_type_label = ref_seq.seq_type.capitalize()
        fig.suptitle(f'{seq_type_label}: {ref_seq.seqid}',
                    fontsize=14, fontweight='bold', y=0.98)

        # Define fixed positions for all plot areas [left, bottom, width, height]
        # These positions ensure all tracks have identical width
        plot_left = 0.08
        plot_width = 0.82  # Leave room for colorbar on the right
        
        # Heights and positions for 5 tracks
        track_heights = [0.22, 0.11, 0.11, 0.11, 0.11]  # Proportional heights
        track_spacing = 0.055  # Space between tracks
        
        # Calculate bottom positions (from bottom to top)
        bottom_positions = []
        current_bottom = 0.06
        for i in range(4, -1, -1):  # Start from bottom track (track 5) to top (track 1)
            bottom_positions.insert(0, current_bottom)
            current_bottom += track_heights[i] + track_spacing
        
        # Create axes with exact positions
        ax1 = fig.add_axes([plot_left, bottom_positions[0], plot_width, track_heights[0]])
        ax2 = fig.add_axes([plot_left, bottom_positions[1], plot_width, track_heights[1]])
        ax3 = fig.add_axes([plot_left, bottom_positions[2], plot_width, track_heights[2]])
        ax4 = fig.add_axes([plot_left, bottom_positions[3], plot_width, track_heights[3]])
        ax5 = fig.add_axes([plot_left, bottom_positions[4], plot_width, track_heights[4]])

        # ====================================================================
        # TRACK 1: GENE QUALITY
        # ====================================================================

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

        # Add colorbar to the right of the plot area
        cbar_left = plot_left + plot_width + 0.01
        cbar_width = 0.02
        cax = fig.add_axes([cbar_left, bottom_positions[0], cbar_width, track_heights[0]])
        
        cmap = plt.cm.RdYlGn
        norm = Normalize(vmin=0, vmax=100)
        sm = ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, cax=cax, orientation='vertical')
        cbar.set_label('Quality Score', fontsize=9)

        # ====================================================================
        # TRACK 2: CONTIG MAPPING
        # ====================================================================

        ax2.set_title('Assembly Contig Mapping', fontsize=11, pad=10)

        # Get consistent colors using CircularVisualizer method
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

        # ====================================================================
        # TRACK 3: COVERAGE
        # ====================================================================

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

        # ====================================================================
        # TRACK 4: IDENTITY
        # ====================================================================

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

        # ====================================================================
        # TRACK 5: MISASSEMBLIES
        # ====================================================================

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

        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()