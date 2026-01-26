"""
Dashboard generators for GenomeViz output.

Creates HTML dashboard pages for:
- Individual alignment modes (contig vs reference, scaffold vs reference)
- Root-level sample overview linking all analysis modes

GitHub: https://github.com/Aaron-Thiel/GenomeViz
License: MIT
"""

from pathlib import Path
from typing import Dict, List, Optional
import json


class AlignmentDashboardGenerator:
    """
    Generate index.html dashboards for contig_alignment and scaffold_alignment folders.

    These dashboards show alignment statistics and link to visualizations.
    """

    def __init__(self, output_dir: Path, ref_sequences: dict, assembly_name: str,
                 reference_name: str, mode: str = 'contig'):
        """
        Initialize the dashboard generator.

        Args:
            output_dir: Directory where index.html will be created
            ref_sequences: Dictionary of ReferenceSequence objects with analysis results
            assembly_name: Name of the assembly file
            reference_name: Name of the reference file
            mode: 'contig' or 'scaffold' - affects display labels
        """
        self.output_dir = Path(output_dir)
        self.ref_sequences = ref_sequences
        self.assembly_name = assembly_name
        self.reference_name = reference_name
        self.mode = mode
        self.mode_label = "Contig" if mode == 'contig' else "Scaffold"

    def generate_index_html(self) -> Path:
        """Generate the index.html dashboard."""

        # Collect stats from all reference sequences
        total_genes = 0
        complete_genes = 0
        incomplete_genes = 0
        divergent_genes = 0
        missing_genes = 0
        total_alignments = 0
        total_contigs = 0
        inversions = 0
        gaps = 0
        overlaps = 0

        chromosome_sections = []

        for seqid, ref_seq in self.ref_sequences.items():
            total_alignments += len(ref_seq.alignments)
            total_contigs += len(ref_seq.contig_mapping)

            # Gene stats
            if ref_seq.gene_stats:
                total_genes += len(ref_seq.gene_stats)
                complete_genes += len([g for g in ref_seq.gene_stats if g['status'] == 'complete'])
                incomplete_genes += len([g for g in ref_seq.gene_stats if g['status'] == 'incomplete'])
                divergent_genes += len([g for g in ref_seq.gene_stats if g['status'] == 'divergent'])
                missing_genes += len([g for g in ref_seq.gene_stats if g['status'] == 'missing'])

                avg_coverage = sum(g['coverage_pct'] for g in ref_seq.gene_stats) / len(ref_seq.gene_stats)
                avg_identity = sum(g['avg_identity'] for g in ref_seq.gene_stats) / len(ref_seq.gene_stats)
                avg_quality = sum(g['quality_score'] for g in ref_seq.gene_stats) / len(ref_seq.gene_stats)
            else:
                avg_coverage = avg_identity = avg_quality = 0

            # Misassemblies
            if ref_seq.misassemblies:
                inversions += len([m for m in ref_seq.misassemblies if m['type'] == 'inversion'])
                gaps += len([m for m in ref_seq.misassemblies if m['type'] == 'gap'])
                overlaps += len([m for m in ref_seq.misassemblies if m['type'] == 'overlap'])

            # Check for visualization files
            seqid_dir = self.output_dir / seqid
            has_circular = (seqid_dir / f'{seqid}_interactive_circular.html').exists()
            has_linear = (seqid_dir / f'{seqid}_interactive_linear.html').exists()
            has_gene_stats = (seqid_dir / f'{seqid}_gene_stats.csv').exists()

            viz_links = []
            if has_circular:
                viz_links.append(f'<a href="{seqid}/{seqid}_interactive_circular.html" class="viz-link circular">Circular</a>')
            if has_linear:
                viz_links.append(f'<a href="{seqid}/{seqid}_interactive_linear.html" class="viz-link linear">Linear</a>')
            if has_gene_stats:
                viz_links.append(f'<a href="{seqid}/{seqid}_gene_stats.csv" class="viz-link csv">Gene Stats</a>')

            chromosome_sections.append(f'''
            <div class="chromosome-card">
                <h3>{seqid}</h3>
                <div class="chromosome-stats">
                    <div class="stat-row">
                        <span class="stat-label">Length:</span>
                        <span class="stat-value">{ref_seq.length:,} bp</span>
                    </div>
                    <div class="stat-row">
                        <span class="stat-label">Alignments:</span>
                        <span class="stat-value">{len(ref_seq.alignments)}</span>
                    </div>
                    <div class="stat-row">
                        <span class="stat-label">{self.mode_label}s Aligned:</span>
                        <span class="stat-value">{len(ref_seq.contig_mapping)}</span>
                    </div>
                    <div class="stat-row">
                        <span class="stat-label">Genes:</span>
                        <span class="stat-value">{len(ref_seq.genes)}</span>
                    </div>
                    <div class="stat-row">
                        <span class="stat-label">Avg Coverage:</span>
                        <span class="stat-value">{avg_coverage:.1f}%</span>
                    </div>
                    <div class="stat-row">
                        <span class="stat-label">Avg Identity:</span>
                        <span class="stat-value">{avg_identity:.1f}%</span>
                    </div>
                </div>
                <div class="viz-links">
                    {' '.join(viz_links) if viz_links else '<span class="no-viz">No visualizations available</span>'}
                </div>
            </div>
            ''')

        html_content = f'''<!DOCTYPE html>
<html>
<head>
    <title>{self.mode_label} vs Reference Analysis - GenomeViz</title>
    <style>
        :root {{
            --primary-color: #2563eb;
            --success-color: #16a34a;
            --warning-color: #d97706;
            --danger-color: #dc2626;
            --info-color: #0891b2;
            --purple-color: #7c3aed;
            --bg-color: #f8fafc;
            --card-bg: #ffffff;
            --text-color: #1e293b;
            --text-muted: #64748b;
            --border-color: #e2e8f0;
        }}

        * {{ box-sizing: border-box; margin: 0; padding: 0; }}

        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: var(--bg-color);
            color: var(--text-color);
            line-height: 1.6;
            padding: 2rem;
        }}

        .container {{ max-width: 1200px; margin: 0 auto; }}

        header {{
            text-align: center;
            margin-bottom: 2rem;
            padding-bottom: 1rem;
            border-bottom: 2px solid var(--border-color);
        }}

        h1 {{ color: var(--primary-color); margin-bottom: 0.5rem; }}
        .subtitle {{ color: var(--text-muted); font-size: 1.1rem; }}

        .back-link {{
            display: inline-block;
            margin-bottom: 1rem;
            color: var(--primary-color);
            text-decoration: none;
        }}
        .back-link:hover {{ text-decoration: underline; }}

        .files-info {{
            background: var(--card-bg);
            padding: 1rem;
            border-radius: 8px;
            margin-bottom: 2rem;
            border: 1px solid var(--border-color);
        }}
        .files-info p {{ margin: 0.25rem 0; }}
        .files-info strong {{ color: var(--text-muted); }}

        .summary-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
            gap: 1rem;
            margin-bottom: 2rem;
        }}

        .stat-card {{
            background: var(--card-bg);
            padding: 1.5rem;
            border-radius: 8px;
            text-align: center;
            border: 1px solid var(--border-color);
        }}
        .stat-card.complete {{ border-left: 4px solid var(--success-color); }}
        .stat-card.incomplete {{ border-left: 4px solid var(--warning-color); }}
        .stat-card.missing {{ border-left: 4px solid var(--danger-color); }}

        .stat-value {{
            font-size: 2rem;
            font-weight: bold;
            color: var(--primary-color);
        }}
        .stat-card.complete .stat-value {{ color: var(--success-color); }}
        .stat-card.incomplete .stat-value {{ color: var(--warning-color); }}
        .stat-card.missing .stat-value {{ color: var(--danger-color); }}

        .stat-label {{
            color: var(--text-muted);
            font-size: 0.9rem;
            margin-top: 0.25rem;
        }}

        h2 {{
            color: var(--text-color);
            margin: 2rem 0 1rem 0;
            padding-bottom: 0.5rem;
            border-bottom: 1px solid var(--border-color);
        }}

        .chromosome-card {{
            background: var(--card-bg);
            padding: 1.5rem;
            border-radius: 8px;
            margin-bottom: 1rem;
            border: 1px solid var(--border-color);
        }}
        .chromosome-card h3 {{
            color: var(--primary-color);
            margin-bottom: 1rem;
        }}

        .chromosome-stats {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 0.5rem;
            margin-bottom: 1rem;
        }}

        .stat-row {{
            display: flex;
            justify-content: space-between;
            padding: 0.25rem 0;
        }}
        .stat-row .stat-label {{ color: var(--text-muted); }}
        .stat-row .stat-value {{ font-weight: 500; }}

        .viz-links {{
            display: flex;
            gap: 0.5rem;
            flex-wrap: wrap;
            padding-top: 1rem;
            border-top: 1px solid var(--border-color);
        }}

        .viz-link {{
            display: inline-block;
            padding: 0.375rem 0.75rem;
            color: white;
            text-decoration: none;
            border-radius: 4px;
            font-size: 0.8rem;
            font-weight: 500;
        }}
        .viz-link.circular {{ background: var(--purple-color); }}
        .viz-link.linear {{ background: var(--primary-color); }}
        .viz-link.csv {{ background: var(--success-color); }}
        .viz-link:hover {{ opacity: 0.9; }}

        .no-viz {{ color: var(--text-muted); font-style: italic; }}

        .misassembly-section {{
            background: var(--card-bg);
            padding: 1.5rem;
            border-radius: 8px;
            border: 1px solid var(--border-color);
        }}

        .misassembly-grid {{
            display: grid;
            grid-template-columns: repeat(3, 1fr);
            gap: 1rem;
            margin-top: 1rem;
        }}

        .misassembly-item {{
            text-align: center;
            padding: 1rem;
            background: var(--bg-color);
            border-radius: 4px;
        }}
        .misassembly-item .count {{
            font-size: 1.5rem;
            font-weight: bold;
        }}
        .misassembly-item .label {{
            color: var(--text-muted);
            font-size: 0.9rem;
        }}

        footer {{
            text-align: center;
            margin-top: 3rem;
            padding-top: 1rem;
            border-top: 1px solid var(--border-color);
            color: var(--text-muted);
            font-size: 0.9rem;
        }}
        footer a {{ color: var(--primary-color); }}
    </style>
</head>
<body>
    <div class="container">
        <a href="../index.html" class="back-link">← Back to Sample Overview</a>

        <header>
            <h1>{self.mode_label} vs Reference Analysis</h1>
            <p class="subtitle">Alignment quality assessment</p>
        </header>

        <div class="files-info">
            <p><strong>Assembly:</strong> {self.assembly_name}</p>
            <p><strong>Reference:</strong> {self.reference_name}</p>
        </div>

        <h2>Gene Quality Summary</h2>
        <div class="summary-grid">
            <div class="stat-card">
                <div class="stat-value">{total_genes}</div>
                <div class="stat-label">Total Genes</div>
            </div>
            <div class="stat-card complete">
                <div class="stat-value">{complete_genes}</div>
                <div class="stat-label">Complete</div>
            </div>
            <div class="stat-card incomplete">
                <div class="stat-value">{incomplete_genes}</div>
                <div class="stat-label">Incomplete</div>
            </div>
            <div class="stat-card missing">
                <div class="stat-value">{missing_genes}</div>
                <div class="stat-label">Missing</div>
            </div>
        </div>

        <div class="misassembly-section">
            <h3>Misassemblies Detected</h3>
            <div class="misassembly-grid">
                <div class="misassembly-item">
                    <div class="count">{inversions}</div>
                    <div class="label">Inversions</div>
                </div>
                <div class="misassembly-item">
                    <div class="count">{gaps}</div>
                    <div class="label">Gaps</div>
                </div>
                <div class="misassembly-item">
                    <div class="count">{overlaps}</div>
                    <div class="label">Overlaps</div>
                </div>
            </div>
        </div>

        <h2>Per-Chromosome Details</h2>
        {''.join(chromosome_sections)}

        <footer>
            Generated by <a href="https://github.com/Aaron-Thiel/GenomeViz">GenomeViz</a>
        </footer>
    </div>
</body>
</html>
'''

        index_file = self.output_dir / 'index.html'
        with open(index_file, 'w', encoding='utf-8') as f:
            f.write(html_content)

        return index_file


class SampleDashboardGenerator:
    """
    Generate root-level index.html that provides an overview of all analysis modes.

    Links to:
    - contig_alignment/index.html
    - scaffold_alignment/index.html
    - assembly_comparison/index.html
    """

    def __init__(self, output_dir: Path, sample_name: str,
                 contig_summary: Optional[dict] = None,
                 scaffold_summary: Optional[dict] = None,
                 comparison_summary: Optional[dict] = None):
        """
        Initialize the sample dashboard generator.

        Args:
            output_dir: Root output directory for the sample
            sample_name: Name of the sample
            contig_summary: Summary data from contig analysis (optional)
            scaffold_summary: Summary data from scaffold analysis (optional)
            comparison_summary: Summary data from comparison analysis (optional)
        """
        self.output_dir = Path(output_dir)
        self.sample_name = sample_name
        self.contig_summary = contig_summary
        self.scaffold_summary = scaffold_summary
        self.comparison_summary = comparison_summary

    def generate_index_html(self) -> Path:
        """Generate the root index.html dashboard."""

        # Check which analyses are available
        contig_dir = self.output_dir / 'contig_alignment'
        scaffold_dir = self.output_dir / 'scaffold_alignment'
        comparison_dir = self.output_dir / 'assembly_comparison'

        has_contig = contig_dir.exists() and (contig_dir / 'index.html').exists()
        has_scaffold = scaffold_dir.exists() and (scaffold_dir / 'index.html').exists()
        has_comparison = comparison_dir.exists() and (comparison_dir / 'index.html').exists()

        # Build analysis cards
        analysis_cards = []

        if has_contig or (contig_dir / 'summary_report.txt').exists():
            contig_stats = self._parse_summary_file(contig_dir / 'summary_report.txt')
            analysis_cards.append(self._build_analysis_card(
                'Contig vs Reference',
                'contig_alignment',
                'Alignment of individual contigs against the reference genome',
                contig_stats,
                has_contig
            ))

        if has_scaffold or (scaffold_dir / 'summary_report.txt').exists():
            scaffold_stats = self._parse_summary_file(scaffold_dir / 'summary_report.txt')
            analysis_cards.append(self._build_analysis_card(
                'Scaffold vs Reference',
                'scaffold_alignment',
                'Alignment of assembled scaffolds against the reference genome',
                scaffold_stats,
                has_scaffold
            ))

        if has_comparison or comparison_dir.exists():
            comparison_stats = self._parse_comparison_summary(comparison_dir / 'comparison_summary.json')
            analysis_cards.append(self._build_comparison_card(
                'Scaffold vs Contig Comparison',
                'assembly_comparison',
                'Comparison of scaffold assembly against original contigs',
                comparison_stats,
                has_comparison
            ))

        html_content = f'''<!DOCTYPE html>
<html>
<head>
    <title>{self.sample_name} - GenomeViz Analysis</title>
    <style>
        :root {{
            --primary-color: #2563eb;
            --success-color: #16a34a;
            --warning-color: #d97706;
            --danger-color: #dc2626;
            --info-color: #0891b2;
            --purple-color: #7c3aed;
            --bg-color: #f8fafc;
            --card-bg: #ffffff;
            --text-color: #1e293b;
            --text-muted: #64748b;
            --border-color: #e2e8f0;
        }}

        * {{ box-sizing: border-box; margin: 0; padding: 0; }}

        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: var(--bg-color);
            color: var(--text-color);
            line-height: 1.6;
            padding: 2rem;
        }}

        .container {{ max-width: 1200px; margin: 0 auto; }}

        header {{
            text-align: center;
            margin-bottom: 2rem;
            padding-bottom: 1rem;
            border-bottom: 2px solid var(--border-color);
        }}

        h1 {{
            color: var(--primary-color);
            margin-bottom: 0.5rem;
            font-size: 2.5rem;
        }}
        .subtitle {{ color: var(--text-muted); font-size: 1.1rem; }}

        .sample-info {{
            background: var(--card-bg);
            padding: 1.5rem;
            border-radius: 8px;
            margin-bottom: 2rem;
            border: 1px solid var(--border-color);
            text-align: center;
        }}
        .sample-info h2 {{
            color: var(--text-color);
            margin-bottom: 0.5rem;
        }}

        .analysis-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(350px, 1fr));
            gap: 1.5rem;
            margin-bottom: 2rem;
        }}

        .analysis-card {{
            background: var(--card-bg);
            border-radius: 12px;
            border: 1px solid var(--border-color);
            overflow: hidden;
            transition: box-shadow 0.2s, transform 0.2s;
        }}
        .analysis-card:hover {{
            box-shadow: 0 4px 12px rgba(0,0,0,0.1);
            transform: translateY(-2px);
        }}

        .analysis-card.contig {{ border-top: 4px solid var(--info-color); }}
        .analysis-card.scaffold {{ border-top: 4px solid var(--primary-color); }}
        .analysis-card.comparison {{ border-top: 4px solid var(--purple-color); }}

        .card-header {{
            padding: 1.5rem;
            border-bottom: 1px solid var(--border-color);
        }}
        .card-header h3 {{
            margin-bottom: 0.5rem;
            font-size: 1.25rem;
        }}
        .card-header p {{
            color: var(--text-muted);
            font-size: 0.9rem;
        }}

        .card-body {{
            padding: 1.5rem;
        }}

        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(2, 1fr);
            gap: 1rem;
            margin-bottom: 1rem;
        }}

        .stat-item {{
            text-align: center;
            padding: 0.75rem;
            background: var(--bg-color);
            border-radius: 6px;
        }}
        .stat-item .value {{
            font-size: 1.5rem;
            font-weight: bold;
            color: var(--primary-color);
        }}
        .stat-item.success .value {{ color: var(--success-color); }}
        .stat-item.warning .value {{ color: var(--warning-color); }}
        .stat-item.danger .value {{ color: var(--danger-color); }}
        .stat-item .label {{
            font-size: 0.8rem;
            color: var(--text-muted);
        }}

        .card-footer {{
            padding: 1rem 1.5rem;
            border-top: 1px solid var(--border-color);
            background: var(--bg-color);
        }}

        .view-btn {{
            display: inline-block;
            width: 100%;
            padding: 0.75rem 1.5rem;
            background: var(--primary-color);
            color: white;
            text-decoration: none;
            border-radius: 6px;
            text-align: center;
            font-weight: 500;
            transition: opacity 0.2s;
        }}
        .view-btn:hover {{ opacity: 0.9; }}
        .view-btn.disabled {{
            background: var(--text-muted);
            cursor: not-allowed;
        }}

        .no-data {{
            color: var(--text-muted);
            font-style: italic;
            text-align: center;
            padding: 2rem;
        }}

        footer {{
            text-align: center;
            margin-top: 3rem;
            padding-top: 1rem;
            border-top: 1px solid var(--border-color);
            color: var(--text-muted);
            font-size: 0.9rem;
        }}
        footer a {{ color: var(--primary-color); }}
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>GenomeViz Analysis</h1>
            <p class="subtitle">Genome assembly quality assessment dashboard</p>
        </header>

        <div class="sample-info">
            <h2>{self.sample_name}</h2>
            <p>Analysis results overview</p>
        </div>

        <div class="analysis-grid">
            {''.join(analysis_cards) if analysis_cards else '<div class="no-data">No analysis results found</div>'}
        </div>

        <footer>
            Generated by <a href="https://github.com/Aaron-Thiel/GenomeViz">GenomeViz</a>
        </footer>
    </div>
</body>
</html>
'''

        index_file = self.output_dir / 'index.html'
        with open(index_file, 'w', encoding='utf-8') as f:
            f.write(html_content)

        return index_file

    def _parse_summary_file(self, summary_file: Path) -> dict:
        """Parse summary_report.txt to extract key statistics."""
        stats = {
            'complete': 0,
            'incomplete': 0,
            'missing': 0,
            'coverage': 0,
            'identity': 0,
            'gaps': 0
        }

        if not summary_file.exists():
            return stats

        try:
            with open(summary_file, 'r') as f:
                content = f.read()

            # Parse gene quality
            for line in content.split('\n'):
                line = line.strip()
                if line.startswith('Complete:'):
                    stats['complete'] = int(line.split(':')[1].strip())
                elif line.startswith('Incomplete:'):
                    stats['incomplete'] = int(line.split(':')[1].strip())
                elif line.startswith('Missing:'):
                    stats['missing'] = int(line.split(':')[1].strip())
                elif line.startswith('Coverage:'):
                    stats['coverage'] = float(line.split(':')[1].strip().replace('%', ''))
                elif line.startswith('Identity:'):
                    stats['identity'] = float(line.split(':')[1].strip().replace('%', ''))
                elif line.startswith('Gaps:'):
                    stats['gaps'] = int(line.split(':')[1].strip())
        except Exception:
            pass

        return stats

    def _parse_comparison_summary(self, json_file: Path) -> dict:
        """Parse comparison_summary.json to extract key statistics."""
        stats = {
            'num_scaffolds': 0,
            'num_contigs': 0,
            'overlap_regions': 0
        }

        if not json_file.exists():
            return stats

        try:
            with open(json_file, 'r') as f:
                data = json.load(f)
            stats['num_scaffolds'] = data.get('num_scaffolds', 0)
            stats['num_contigs'] = data.get('num_contigs', 0)
            stats['overlap_regions'] = data.get('num_total_overlap_regions', 0)
        except Exception:
            pass

        return stats

    def _build_analysis_card(self, title: str, folder: str, description: str,
                             stats: dict, has_dashboard: bool) -> str:
        """Build HTML for an analysis card (contig or scaffold mode)."""
        card_class = 'contig' if 'Contig' in title else 'scaffold'

        return f'''
        <div class="analysis-card {card_class}">
            <div class="card-header">
                <h3>{title}</h3>
                <p>{description}</p>
            </div>
            <div class="card-body">
                <div class="stats-grid">
                    <div class="stat-item success">
                        <div class="value">{stats['complete']}</div>
                        <div class="label">Complete Genes</div>
                    </div>
                    <div class="stat-item warning">
                        <div class="value">{stats['incomplete']}</div>
                        <div class="label">Incomplete</div>
                    </div>
                    <div class="stat-item danger">
                        <div class="value">{stats['missing']}</div>
                        <div class="label">Missing</div>
                    </div>
                    <div class="stat-item">
                        <div class="value">{stats['gaps']}</div>
                        <div class="label">Gaps</div>
                    </div>
                </div>
                <div class="stats-grid">
                    <div class="stat-item">
                        <div class="value">{stats['coverage']:.1f}%</div>
                        <div class="label">Avg Coverage</div>
                    </div>
                    <div class="stat-item">
                        <div class="value">{stats['identity']:.1f}%</div>
                        <div class="label">Avg Identity</div>
                    </div>
                </div>
            </div>
            <div class="card-footer">
                <a href="{folder}/index.html" class="view-btn{'' if has_dashboard else ' disabled'}">
                    {'View Details' if has_dashboard else 'Dashboard Not Available'}
                </a>
            </div>
        </div>
        '''

    def _build_comparison_card(self, title: str, folder: str, description: str,
                               stats: dict, has_dashboard: bool) -> str:
        """Build HTML for the comparison analysis card."""
        return f'''
        <div class="analysis-card comparison">
            <div class="card-header">
                <h3>{title}</h3>
                <p>{description}</p>
            </div>
            <div class="card-body">
                <div class="stats-grid">
                    <div class="stat-item">
                        <div class="value">{stats['num_scaffolds']}</div>
                        <div class="label">Scaffolds</div>
                    </div>
                    <div class="stat-item">
                        <div class="value">{stats['num_contigs']}</div>
                        <div class="label">Contigs</div>
                    </div>
                    <div class="stat-item success">
                        <div class="value">{stats['overlap_regions']}</div>
                        <div class="label">Overlap Regions</div>
                    </div>
                </div>
            </div>
            <div class="card-footer">
                <a href="{folder}/index.html" class="view-btn{'' if has_dashboard else ' disabled'}">
                    {'View Details' if has_dashboard else 'Dashboard Not Available'}
                </a>
            </div>
        </div>
        '''
