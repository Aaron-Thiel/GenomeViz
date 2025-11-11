"""
InteractiveCircularVisualizer - Interactive circular genome visualization

Creates interactive HTML plots using Plotly with:
- Gene quality (outer ring) - clickable to view gene alignments
- Alignment status (middle ring)
- Contig mapping (inner ring)
- Help button with navigation guide

GitHub: https://github.com/Aaron-Thiel/GenomeViz
License: MIT
"""

import json
import numpy as np
import plotly.graph_objects as go


class InteractiveCircularVisualizer:
    """
    Create interactive Plotly circular plots with hover information and clickable genes.
    """

    @staticmethod
    def create_interactive_circular_plot(ref_seq, output_file):
        """
        Create interactive circular plot using Plotly.

        Args:
            ref_seq (ReferenceSequence): Reference sequence with analysis results
            output_file (str): Path to output HTML file
        """
        print(f"  Creating interactive circular plot: {output_file}")

        # Initialize figure
        fig = go.Figure()

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
        # RING 1: GENE QUALITY (outer ring)
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
        # RING 2: ALIGNMENT STATUS (middle ring)
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

                arc_angles = np.linspace(start_angle, end_angle, max(2, int((end_angle - start_angle) / 2)))
                theta.extend(arc_angles)
                r.extend([0.60] * len(arc_angles))

                size = seg['end'] - seg['start'] + 1
                hover = f"Status: {status_type.capitalize()}<br>Position: {seg['start']:,}-{seg['end']:,} bp<br>Size: {size:,} bp"
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

        # ====================================================================
        # UPDATE LAYOUT
        # ====================================================================

        seq_type = ref_seq.seq_type.capitalize()
        if ref_seq.length > 10000:
            size_str = f'{ref_seq.length/1000:.1f} kb'
        else:
            size_str = f'{ref_seq.length:,} bp'

        fig.update_layout(
            title=dict(
                text=f'{seq_type}: {ref_seq.seqid} ({size_str})',
                x=0.5,
                xanchor='center',
                font=dict(family='Arial, sans-serif')
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
                title=dict(text="<b>Rings</b>", font=dict(size=14, family='Arial, sans-serif')),
                tracegroupgap=20,
                itemsizing='constant',
                font=dict(size=10, family='Arial, sans-serif')
            ),
            font=dict(family='Arial, sans-serif'),
            width=1400,
            height=1000
        )

        # Save as HTML
        fig.write_html(output_file)

    @staticmethod
    def add_help_button(html_file, plot_type='linear', gene_clicking_enabled=True):
        """
        Add a help/info button and modal to HTML file.

        Args:
            html_file (str): Path to the HTML file
            plot_type (str): Type of plot ('linear', 'circular', or 'gene')
            gene_clicking_enabled (bool): Whether gene clicking feature is available
        """
        # Read the HTML file
        with open(html_file, 'r', encoding='utf-8') as f:
            html_content = f.read()

        # Define help content based on plot type
        if plot_type == 'linear':
            help_title = "Linear Plot Navigation Guide"

            # Build interactive features list
            if gene_clicking_enabled:
                click_instruction = '<li><strong>Click on genes</strong> in Track 1 to view detailed alignments</li>'
            else:
                click_instruction = '<li><strong>Note:</strong> Gene clicking feature is disabled (--no-gene-alignments flag was used)</li>'

            help_content = f'''
                <h3>Understanding the Plot</h3>
                <ul>
                    <li><strong>Track 1 (Gene Quality):</strong> Shows genes colored by assembly quality
                        <ul>
                            <li>Green (Excellent): ≥95% identity</li>
                            <li>Yellow (Good): 85-95% identity</li>
                            <li>Orange (Fair): 70-85% identity</li>
                            <li>Red (Poor): &lt;70% identity</li>
                        </ul>
                    </li>
                    <li><strong>Track 2 (Contig Mapping):</strong> Shows which assembly contigs map to the reference</li>
                    <li><strong>Track 3 (Coverage Depth):</strong> Shows how many times each position is covered</li>
                    <li><strong>Track 4 (Alignment Identity):</strong> Shows sequence identity percentage across the genome</li>
                    <li><strong>Track 5 (Misassemblies):</strong> Highlights regions with assembly errors
                        <ul>
                            <li>Bottom layer: Contig positions (colored bars showing which contig)</li>
                            <li>Upper layers: Overlaps, Gaps, and Inversions</li>
                        </ul>
                    </li>
                </ul>

                <h3>Interactive Features</h3>
                <ul>
                    {click_instruction}
                    <li><strong>Hover</strong> over any element to see detailed information</li>
                    <li><strong>Zoom:</strong> Click and drag to select a region, double-click to reset</li>
                    <li><strong>Pan:</strong> Use the pan tool in the toolbar</li>
                </ul>
            '''
        elif plot_type == 'circular':
            help_title = "Circular Plot Navigation Guide"

            # Build interactive features list
            if gene_clicking_enabled:
                click_instruction = '<li><strong>Click on genes</strong> in Ring 1 to view detailed alignments</li>'
            else:
                click_instruction = '<li><strong>Note:</strong> Gene clicking feature is disabled (--no-gene-alignments flag was used)</li>'

            help_content = f'''
                <h3>Understanding the Rings</h3>
                <ul>
                    <li><strong>Ring 1 (Outer - Gene Quality):</strong> Shows genes colored by assembly quality
                        <ul>
                            <li>Green (Excellent): ≥95% identity</li>
                            <li>Yellow (Good): 85-95% identity</li>
                            <li>Orange (Fair): 70-85% identity</li>
                            <li>Red (Poor): &lt;70% identity</li>
                        </ul>
                    </li>
                    <li><strong>Ring 2 (Middle - Alignment Status):</strong> Shows assembly coverage status
                        <ul>
                            <li>Green (Complete): Single coverage</li>
                            <li>Orange (Duplicated): Multiple coverage</li>
                            <li>Red (Inverted): Reverse strand alignment</li>
                            <li>Gray (Missing): No coverage</li>
                        </ul>
                    </li>
                    <li><strong>Ring 3 (Inner - Contig Mapping):</strong> Shows which assembly contigs map to each region</li>
                </ul>

                <h3>Interactive Features</h3>
                <ul>
                    {click_instruction}
                    <li><strong>Hover</strong> over any ring to see detailed information</li>
                </ul>
            '''
        else:  # gene alignment
            help_title = "Gene Alignment Guide"
            help_content = '''
                <h3>Understanding the Alignment</h3>
                <ul>
                    <li><strong>Reference:</strong> The top sequence from the reference genome</li>
                    <li><strong>Assembly:</strong> The bottom sequence from your assembly</li>
                    <li><strong>Match indicators (|):</strong> Show positions where sequences match</li>
                    <li><strong>Gaps (-):</strong> Indicate insertions/deletions or missing coverage</li>
                </ul>

                <h3>Alignment Method</h3>
                <ul>
                    <li><strong>Blue box:</strong> CIGAR-based alignment (accurate gap placement from minimap2)</li>
                    <li><strong>Yellow box:</strong> Fallback alignment (no CIGAR data available)</li>
                </ul>

                <h3>Coverage Statistics</h3>
                <ul>
                    <li><strong>Identity:</strong> Percentage of matching bases (from minimap2)</li>
                    <li><strong>Coverage:</strong> How many positions are covered by assembly alignments</li>
                    <li><strong>Gaps in assembly:</strong> Regions not covered by any alignment</li>
                </ul>

                <h3>Navigation</h3>
                <ul>
                    <li><strong>Scroll horizontally</strong> to view the entire gene sequence</li>
                </ul>
            '''

        # Create HTML/CSS/JavaScript for help button and modal
        help_code = '''
        <style>
        /* Help button styles */
        .help-button {
            position: fixed;
            bottom: 20px;
            right: 20px;
            width: 50px;
            height: 50px;
            border-radius: 50%;
            background-color: #007bff;
            color: white;
            border: none;
            font-size: 24px;
            font-weight: bold;
            cursor: pointer;
            box-shadow: 0 4px 6px rgba(0,0,0,0.3);
            z-index: 9999;
            transition: all 0.3s ease;
        }

        .help-button:hover {
            background-color: #0056b3;
            transform: scale(1.1);
        }

        /* Modal styles */
        .help-modal {
            display: none;
            position: fixed;
            z-index: 10000;
            left: 0;
            top: 0;
            width: 100%;
            height: 100%;
            overflow: auto;
            background-color: rgba(0,0,0,0.5);
        }

        .help-modal-content {
            background-color: #fefefe;
            margin: 5% auto;
            padding: 30px;
            border: 1px solid #888;
            border-radius: 8px;
            width: 80%;
            max-width: 700px;
            box-shadow: 0 4px 20px rgba(0,0,0,0.3);
            font-family: 'Courier New', monospace;
        }

        .help-modal-header {
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 20px;
            padding-bottom: 15px;
            border-bottom: 2px solid #007bff;
        }

        .help-modal-header h2 {
            margin: 0;
            color: #007bff;
        }

        .help-close {
            color: #aaa;
            font-size: 32px;
            font-weight: bold;
            cursor: pointer;
            transition: color 0.3s ease;
        }

        .help-close:hover,
        .help-close:focus {
            color: #000;
        }

        .help-modal-content h3 {
            color: #333;
            margin-top: 20px;
            margin-bottom: 10px;
        }

        .help-modal-content ul {
            line-height: 1.8;
            color: #555;
        }

        .help-modal-content li {
            margin-bottom: 8px;
        }

        .help-modal-content strong {
            color: #007bff;
        }
        </style>

        <!-- Help Button -->
        <button class="help-button" onclick="openHelpModal()" title="Help & Navigation Guide">?</button>

        <!-- Help Modal -->
        <div id="helpModal" class="help-modal">
            <div class="help-modal-content">
                <div class="help-modal-header">
                    <h2>''' + help_title + '''</h2>
                    <span class="help-close" onclick="closeHelpModal()">&times;</span>
                </div>
                <div class="help-modal-body">
                    ''' + help_content + '''
                </div>
            </div>
        </div>

        <script>
        function openHelpModal() {
            document.getElementById('helpModal').style.display = 'block';
        }

        function closeHelpModal() {
            document.getElementById('helpModal').style.display = 'none';
        }

        // Close modal when clicking outside of it
        window.onclick = function(event) {
            const modal = document.getElementById('helpModal');
            if (event.target == modal) {
                modal.style.display = 'none';
            }
        }

        // Close modal with Escape key
        document.addEventListener('keydown', function(event) {
            if (event.key === 'Escape') {
                closeHelpModal();
            }
        });
        </script>
        '''

        # Insert help button and modal before closing body tag
        html_content = html_content.replace('</body>', help_code + '\n</body>')

        # Write back to file
        with open(html_file, 'w', encoding='utf-8') as f:
            f.write(html_content)

    @staticmethod
    def add_circular_click_handler(html_file, ref_seq):
        """
        Add JavaScript to circular plot HTML to make genes clickable.

        Args:
            html_file (str): Path to the HTML file
            ref_seq (ReferenceSequence): Reference sequence with gene information
        """
        # Build gene mapping with unique filenames
        gene_list = []
        for gene in ref_seq.gene_stats:
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
        const refLength = ''' + str(ref_seq.length) + ''';

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
