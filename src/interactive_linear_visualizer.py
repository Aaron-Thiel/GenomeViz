"""
InteractiveLinearVisualizer - Interactive linear genome visualization

Creates interactive HTML plots using Plotly with:
- Gene quality track (clickable to view gene alignments)
- Contig mapping track
- Coverage depth track
- Alignment identity track
- Misassemblies track with contig indicators
- Dynamic border visibility based on zoom level
- Help button with navigation guide

GitHub: https://github.com/Aaron-Thiel/GenomeViz
License: MIT
"""

import json
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Import other visualizer classes
from .circular_visualizer import CircularVisualizer
from .interactive_circular_visualizer import InteractiveCircularVisualizer
from .gene_alignment_visualizer import GeneAlignmentVisualizer


class InteractiveLinearVisualizer:
    """
    Create interactive linear plots with multi-level zooming and clickable genes.
    """

    def __init__(self, ref_seq, aligner, gene_clicking_enabled=True):
        """
        Initialize interactive linear visualizer.

        Args:
            ref_seq (ReferenceSequence): Reference sequence with analysis results
            aligner (GenomeAligner): Aligner object to access assembly sequences
            gene_clicking_enabled (bool): Whether gene clicking feature is available
        """
        self.ref_seq = ref_seq
        self.aligner = aligner
        self.gene_clicking_enabled = gene_clicking_enabled
        self.gene_visualizer = GeneAlignmentVisualizer(ref_seq, aligner.assembly_fasta)

    def create_interactive_linear_plot(self, output_file):
        """
        Create interactive linear plot with multi-level zoom.

        Args:
            output_file (str): Path to output HTML file
        """
        print(f"  Creating interactive linear plot: {output_file}")
        print(f"    Extracting nucleotide-level differences...")

        # Create subplot figure with 5 tracks
        fig = make_subplots(
            rows=5, cols=1,
            shared_xaxes=True,
            vertical_spacing=0.02,
            row_heights=[0.25, 0.20, 0.20, 0.15, 0.20],
            subplot_titles=(
                'Gene Quality',
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

        # Build gene annotations for click handler
        gene_annotations = []
        for gene in self.ref_seq.gene_stats:
            safe_gene_name = "".join(c if c.isalnum() or c in ('-', '_') else '_' for c in gene['name'])
            unique_name = f"{safe_gene_name}_{gene['start']}-{gene['end']}"
            gene_annotations.append({
                'gene': gene,
                'safe_name': unique_name
            })

        for quality_level in ['excellent', 'good', 'fair', 'poor']:
            x_positions = []
            widths = []
            heights = []
            hover_text = []

            for gene in self.ref_seq.gene_stats:
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

                hover = (f"<b>{gene['name']}</b> (click to view alignment)<br>"
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

        # Get contig colors
        contig_colors = CircularVisualizer.get_contig_colors(self.ref_seq.alignments)

        # Convert matplotlib colors to RGB strings
        contig_colors_rgb = {}
        for contig_name, color in contig_colors.items():
            if isinstance(color, tuple):
                r, g, b = int(color[0]*255), int(color[1]*255), int(color[2]*255)
                contig_colors_rgb[contig_name] = f'rgb({r}, {g}, {b})'
            else:
                contig_colors_rgb[contig_name] = color

        # Group by contig
        contig_data = {}
        for aln in self.ref_seq.alignments:
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
            color = contig_colors_rgb.get(contig_name, 'rgb(136, 136, 136)')
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

        coverage = np.zeros(self.ref_seq.length)
        for aln in self.ref_seq.alignments:
            if aln['is_primary']:
                coverage[aln['ref_start']:aln['ref_end']] += 1

        # Downsample if genome is large for performance
        window_size = max(1, self.ref_seq.length // 10000)
        if window_size > 1:
            downsampled_coverage = []
            downsampled_positions = []
            for i in range(0, self.ref_seq.length, window_size):
                window_end = min(i + window_size, self.ref_seq.length)
                downsampled_coverage.append(np.mean(coverage[i:window_end]))
                downsampled_positions.append(i + window_size // 2)
        else:
            downsampled_coverage = coverage
            downsampled_positions = np.arange(self.ref_seq.length)

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

        for aln in self.ref_seq.alignments:
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
            misassemblies = [m for m in self.ref_seq.misassemblies if m['type'] == mis_type]

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
                            f"Length: {mis['size']:,} bp<br>"
                            f"Contig: {mis['query_name']}")
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
        for aln in self.ref_seq.alignments:
            if aln['is_primary']:
                contig_name = aln['query_name']
                color = contig_colors_rgb.get(contig_name, 'rgb(136, 136, 136)')

                center_x = (aln['ref_start'] + aln['ref_end']) / 2
                width = aln['ref_end'] - aln['ref_start']

                hover = (f"<b>Contig: {contig_name}</b><br>"
                        f"Position: {aln['ref_start']:,}-{aln['ref_end']:,} bp<br>"
                        f"Strand: {aln['strand']}")

                fig.add_trace(go.Bar(
                    x=[center_x],
                    y=[0.1],
                    width=[width],
                    base=[0],
                    name=contig_name,
                    marker=dict(
                        color=color,
                        opacity=0.7,
                        line=dict(color='rgba(0,0,0,0.3)', width=0.5)
                    ),
                    hovertext=[hover],
                    hoverinfo='text',
                    legendgroup='Contig Indicators',
                    showlegend=False,
                    orientation='v'
                ), row=5, col=1)

        # ====================================================================
        # CONFIGURE LAYOUT
        # ====================================================================

        seq_type = self.ref_seq.seq_type.capitalize()
        if self.ref_seq.length > 10000:
            size_str = f'{self.ref_seq.length/1000:.1f} kb'
        else:
            size_str = f'{self.ref_seq.length:,} bp'

        fig.update_layout(
            title=dict(
                text=f'{seq_type}: {self.ref_seq.seqid} ({size_str})',
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
            range=[0, self.ref_seq.length],
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
                'filename': f'{self.ref_seq.seqid}_interactive_linear',
                'height': 1200,
                'width': 1600,
                'scale': 2
            }
        }

        # Write HTML with custom click handler
        fig.write_html(output_file, config=config)

        # Add custom JavaScript and help button
        self._add_click_handler_to_html(output_file, gene_annotations)
        InteractiveCircularVisualizer.add_help_button(
            output_file,
            plot_type='linear',
            gene_clicking_enabled=self.gene_clicking_enabled
        )

    def _add_click_handler_to_html(self, html_file, gene_annotations):
        """Add JavaScript to HTML file to make gene blocks clickable."""

        # Build gene mapping
        gene_list = []
        for gene_info in gene_annotations:
            gene = gene_info['gene']
            safe_name = gene_info['safe_name']
            link_path = f"gene_alignments/{safe_name}_alignment.html"

            gene_list.append({
                'name': gene['name'],
                'start': gene['start'],
                'end': gene['end'],
                'link': link_path
            })

        # Read the HTML file
        with open(html_file, 'r', encoding='utf-8') as f:
            html_content = f.read()

        # Create JavaScript code for click handling and dynamic borders
        js_code = '''
        <script>
        // Gene positions and links
        const geneList = ''' + json.dumps(gene_list, indent=2) + ''';

        // Add click handler and dynamic border adjustment
        document.addEventListener('DOMContentLoaded', function() {
            const plotDiv = document.getElementsByClassName('plotly-graph-div')[0];

            if (plotDiv) {
                // Function to update border visibility based on zoom level
                function updateBorders() {
                    const xaxis = plotDiv.layout.xaxis;
                    if (xaxis && xaxis.range) {
                        const range = xaxis.range[1] - xaxis.range[0];
                        const borderWidth = range <= 200000 ? 1 : 0;

                        // Update gene bar traces (first 4 traces are gene quality)
                        const update = {
                            'marker.line.width': borderWidth
                        };

                        Plotly.restyle(plotDiv, update, [0, 1, 2, 3]);
                    }
                }

                // Update borders on zoom/pan
                plotDiv.on('plotly_relayout', function(eventData) {
                    updateBorders();
                });

                // Initial border update
                setTimeout(updateBorders, 100);

                // Click handler for genes
                plotDiv.on('plotly_click', function(data) {
                    if (data.points && data.points.length > 0) {
                        const point = data.points[0];
                        const xValue = point.x;

                        // Check if click is on track 1 (gene quality)
                        if (point.data.name && (point.data.name.includes('xcellent') ||
                                                point.data.name.includes('ood') ||
                                                point.data.name.includes('air') ||
                                                point.data.name.includes('oor'))) {

                            // Find which gene was clicked
                            for (let i = 0; i < geneList.length; i++) {
                                const gene = geneList[i];
                                if (xValue >= gene.start && xValue <= gene.end) {
                                    // Open gene alignment in new tab
                                    window.open(gene.link, '_blank');
                                    console.log('Opening alignment for:', gene.name, 'at', xValue);
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

    def generate_gene_alignments(self, gene_align_dir):
        """
        Generate individual gene alignment HTML files.

        Args:
            gene_align_dir (Path): Directory to save gene alignment files
        """
        print(f"  Generating gene alignments for {self.ref_seq.seqid}...")

        for gene in self.ref_seq.gene_stats:
            try:
                self.gene_visualizer.create_gene_alignment_html(gene, gene_align_dir)
            except Exception as e:
                print(f"    Warning: Could not create alignment for {gene['name']}: {e}")

        print(f"  Created {len(self.ref_seq.gene_stats)} gene alignment files in {gene_align_dir.name}/")
