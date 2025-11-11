"""
GeneAlignmentVisualizer - Detailed gene-level alignment visualization

Creates detailed HTML pages for individual genes showing:
- CIGAR-based alignment with accurate gap placement
- Nucleotide and amino acid sequences
- Coverage statistics and quality metrics
- Visual highlighting of matches, mismatches, and gaps

GitHub: https://github.com/Aaron-Thiel/GenomeViz
License: MIT
"""

from pathlib import Path
from Bio.Seq import Seq
import mappy as mp

# Import InteractiveCircularVisualizer for help button functionality
from .interactive_circular_visualizer import InteractiveCircularVisualizer


class GeneAlignmentVisualizer:
    """
    Create detailed gene-level alignment visualizations with MSA-style output.
    """

    def __init__(self, ref_seq, assembly_fasta):
        """
        Initialize gene alignment visualizer.

        Args:
            ref_seq (ReferenceSequence): Reference sequence with analysis results
            assembly_fasta (str): Path to assembly FASTA file
        """
        self.ref_seq = ref_seq
        self.assembly_fasta = assembly_fasta
        self.assembly_sequences = self._load_assembly_sequences()

    def _load_assembly_sequences(self):
        """Load all assembly sequences for nucleotide-level comparison."""
        sequences = {}
        for name, seq, qual in mp.fastx_read(self.assembly_fasta):
            sequences[name] = seq
        return sequences

    def create_gene_alignment_html(self, gene, output_dir):
        """
        Create a detailed alignment visualization for a single gene.

        Args:
            gene (dict): Gene information including start, end, name
            output_dir (Path): Directory to save the gene alignment HTML

        Returns:
            str: Path to the generated HTML file
        """
        # Sanitize gene name for filename and make it unique with position
        safe_gene_name = "".join(c if c.isalnum() or c in ('-', '_') else '_' for c in gene['name'])
        unique_name = f"{safe_gene_name}_{gene['start']}-{gene['end']}"
        output_file = Path(output_dir) / f"{unique_name}_alignment.html"

        # Get gene region from reference
        gene_start = gene['start'] - 1  # Convert to 0-based
        gene_end = gene['end']
        ref_seq = self.ref_seq.sequence[gene_start:gene_end]

        # Find overlapping alignments
        overlapping_alns = []
        for aln in self.ref_seq.alignments:
            if aln['is_primary']:
                # Check if alignment overlaps with gene
                if not (aln['ref_end'] <= gene_start or aln['ref_start'] >= gene_end):
                    overlapping_alns.append(aln)

        if not overlapping_alns:
            # No assembly coverage for this gene
            html_content = f"""
            <html>
            <head><title>{gene['name']} - No Coverage</title></head>
            <body style="font-family: Arial, sans-serif; margin: 20px;">
                <h1>Gene: {gene['name']}</h1>
                <p><b>Position:</b> {gene['start']:,} - {gene['end']:,} bp</p>
                <p><b>Length:</b> {len(ref_seq):,} bp</p>
                <p style="color: red;"><b>No assembly coverage for this gene region.</b></p>
            </body>
            </html>
            """
            with open(output_file, 'w') as f:
                f.write(html_content)
            return str(output_file)

        # Build complete gene alignment showing all positions
        gene_length = gene_end - gene_start

        # Initialize with full reference sequence
        full_ref_aligned = list(ref_seq)
        full_asm_aligned = ['-'] * gene_length
        coverage_map = [False] * gene_length

        # Process all overlapping alignments
        has_cigar = False
        for aln in overlapping_alns:
            query_name = aln['query_name']
            if query_name not in self.assembly_sequences:
                continue

            # Get sequences for CIGAR-based alignment
            query_seq_full = self.assembly_sequences[query_name]
            query_seq = query_seq_full[aln['query_start']:aln['query_end']]
            ref_seq_aln = self.ref_seq.sequence[aln['ref_start']:aln['ref_end']]

            # Calculate overlap with gene
            overlap_start = max(gene_start, aln['ref_start'])
            overlap_end = min(gene_end, aln['ref_end'])

            if overlap_start >= overlap_end:
                continue

            # Use CIGAR to create properly aligned sequences
            cigar = aln.get('cigar', None)
            if cigar:
                has_cigar = True
                aligned_ref, aligned_query, match_str = self._parse_cigar_alignment(
                    ref_seq_aln, query_seq, cigar,
                    aln['ref_start'], aln['query_start'], aln['strand']
                )

                # Map aligned sequences to gene coordinates
                gene_offset_start = overlap_start - gene_start
                aln_offset_start = overlap_start - aln['ref_start']

                # Find position in aligned sequences
                ref_pos = 0
                aln_idx = 0
                for i, base in enumerate(aligned_ref):
                    if base != '-':
                        if ref_pos == aln_offset_start:
                            aln_idx = i
                            break
                        ref_pos += 1

                # Copy aligned bases to full gene alignment
                gene_pos = gene_offset_start
                for i in range(aln_idx, len(aligned_ref)):
                    if gene_pos >= gene_length:
                        break

                    if aligned_ref[i] == '-':
                        # Insertion in assembly
                        full_ref_aligned.insert(gene_pos, '-')
                        full_asm_aligned.insert(gene_pos, aligned_query[i])
                        coverage_map.insert(gene_pos, True)
                        gene_pos += 1
                    else:
                        # Match or mismatch
                        full_asm_aligned[gene_pos] = aligned_query[i]
                        coverage_map[gene_pos] = True
                        gene_pos += 1
                        ref_pos += 1
                        if ref_pos >= (overlap_end - aln['ref_start']):
                            break
            else:
                # Fallback to simple extraction if no CIGAR
                offset_in_aln = overlap_start - aln['ref_start']
                length = overlap_end - overlap_start
                query_start = aln['query_start'] + offset_in_aln
                query_end = query_start + length

                asm_seq = query_seq_full[query_start:query_end]
                if aln['strand'] == '-':
                    asm_seq = str(Seq(asm_seq).reverse_complement())

                # Insert into full alignment
                gene_offset = overlap_start - gene_start
                for i, base in enumerate(asm_seq):
                    if gene_offset + i < gene_length:
                        full_asm_aligned[gene_offset + i] = base
                        coverage_map[gene_offset + i] = True

        # Convert to strings
        full_ref_str = ''.join(full_ref_aligned)
        full_asm_str = ''.join(full_asm_aligned)

        # Calculate coverage statistics
        covered_positions = sum(coverage_map)
        coverage_pct = (covered_positions / len(coverage_map) * 100) if coverage_map else 0

        assembly_seqs = [{
            'sequence': full_asm_str,
            'aligned_ref': full_ref_str,
            'ref_start': gene_start,
            'ref_end': gene_end,
            'has_cigar': has_cigar,
            'coverage_map': coverage_map,
            'covered_positions': covered_positions,
            'coverage_pct': coverage_pct
        }] if overlapping_alns else []

        # Translate sequences if length is divisible by 3
        show_translation = len(ref_seq) % 3 == 0 and len(ref_seq) >= 3

        if show_translation:
            try:
                ref_aa = str(Seq(ref_seq).translate())
            except:
                ref_aa = None
                show_translation = False

        # Create HTML visualization
        html_content = self._generate_gene_alignment_html(
            gene, ref_seq, assembly_seqs, ref_aa if show_translation else None
        )

        with open(output_file, 'w') as f:
            f.write(html_content)

        # Add help button to gene alignment file
        InteractiveCircularVisualizer.add_help_button(output_file, plot_type='gene')

        return str(output_file)

    def _generate_gene_alignment_html(self, gene, ref_seq, assembly_seqs, ref_aa=None):
        """Generate HTML content for gene alignment visualization."""

        # Use identity from gene stats
        identity = gene.get('avg_identity', 0)

        # Get assembly sequence for display
        if assembly_seqs:
            asm_seq = assembly_seqs[0]['sequence']
        else:
            asm_seq = ""

        # Build HTML
        html = f"""
        <html>
        <head>
            <title>{gene['name']} - Gene Alignment</title>
            <style>
                body {{
                    font-family: 'Courier New', monospace;
                    margin: 20px;
                    background-color: #f5f5f5;
                }}
                .header {{
                    font-family: Arial, sans-serif;
                    background-color: white;
                    padding: 20px;
                    border-radius: 5px;
                    margin-bottom: 20px;
                    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
                }}
                .alignment-container {{
                    background-color: white;
                    padding: 20px;
                    border-radius: 5px;
                    overflow-x: auto;
                    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
                }}
                .seq-line {{
                    margin: 5px 0;
                    white-space: pre;
                    line-height: 1.6;
                }}
                .label {{
                    display: inline-block;
                    width: 150px;
                    font-weight: bold;
                }}
                .match {{
                    background-color: #d4edda;
                }}
                .mismatch {{
                    background-color: #f8d7da;
                    color: #721c24;
                    font-weight: bold;
                }}
                .gap {{
                    background-color: #fff3cd;
                }}
                .missing {{
                    background-color: #f8d7da;
                    color: #721c24;
                }}
                .stats {{
                    margin: 10px 0;
                }}
            </style>
        </head>
        <body>
            <div class="header">
                <h1>Gene Alignment: {gene['name']}</h1>
                <div class="stats">
                    <p><b>Position:</b> {gene['start']:,} - {gene['end']:,} bp</p>
                    <p><b>Length:</b> {len(ref_seq):,} bp ({len(ref_seq)//3} codons)</p>
                    <p><b>Strand:</b> {gene['strand']}</p>
                    <p><b>Identity (minimap2):</b> {identity:.2f}%</p>
                    <p><b>Coverage:</b> {gene.get('coverage_pct', 0):.1f}%</p>
                    <p><b>Quality Score:</b> {gene.get('quality_score', 0):.1f}/100</p>
                </div>
                <div style="background-color: #e7f3ff; padding: 10px; border-left: 4px solid #2196F3; margin-top: 10px;">
                    <p style="margin: 0; font-size: 13px;"><b>About Identity Calculation:</b></p>
                    <p style="margin: 5px 0 0 0; font-size: 12px;">
                        The <b>minimap2 identity</b> shown above is calculated from the alignment's CIGAR string, which properly accounts for insertions, deletions, and gaps.
                        The alignment below uses the CIGAR string to insert gaps in the correct positions, so the visual alignment should match the identity score.
                    </p>
                </div>
            </div>
        """

        if assembly_seqs:
            html += '<div class="alignment-container">'

            # Use CIGAR-aligned sequences if available
            if assembly_seqs[0].get('has_cigar', False):
                aligned_ref_seq = assembly_seqs[0]['aligned_ref']
                aligned_asm_seq = assembly_seqs[0]['sequence']
                coverage_pct_display = assembly_seqs[0].get('coverage_pct', 0)
                covered_pos = assembly_seqs[0].get('covered_positions', 0)
                total_pos = len(assembly_seqs[0].get('coverage_map', []))

                html += '<div style="background-color: #d1ecf1; padding: 10px; border-left: 4px solid #0c5460; margin-bottom: 15px;">'
                html += '<p style="margin: 0; font-size: 13px; color: #0c5460;"><b>CIGAR-based Alignment - Complete Gene View</b></p>'
                html += f'<p style="margin: 5px 0 0 0; font-size: 12px; color: #0c5460;">Showing full gene ({len(ref_seq):,} bp). '
                html += f'Coverage: {covered_pos:,}/{total_pos:,} positions ({coverage_pct_display:.1f}%). '
                html += 'Gaps (-) in assembly indicate missing coverage.</p>'
                html += '</div>'

                # Show amino acid alignment if available
                if ref_aa:
                    html += '<h2>Amino Acid Alignment</h2>'
                    html += '<p style="color: #666; font-size: 14px;">Scroll horizontally | Green = match, Red = mismatch, Yellow = gap/indel</p>'
                    try:
                        asm_aa_aligned = str(Seq(aligned_asm_seq.replace('-', 'N')).translate()).replace('X', '-')
                        ref_aa_aligned = str(Seq(aligned_ref_seq.replace('-', 'N')).translate()).replace('X', '-')
                        html += self._format_alignment_html(ref_aa_aligned, asm_aa_aligned, "Reference AA", "Assembly AA")
                    except Exception as e:
                        html += f'<p style="color: orange;">Could not translate aligned sequences: {e}</p>'

                # Show nucleotide alignment
                html += '<h2>Nucleotide Alignment</h2>'
                html += '<p style="color: #666; font-size: 14px;">Scroll horizontally | Green = match, Red = mismatch, Yellow = gap/indel</p>'
                html += self._format_alignment_html(aligned_ref_seq, aligned_asm_seq, "Reference", "Assembly")

            else:
                # Fallback to simple alignment
                html += '<div style="background-color: #fff3cd; padding: 10px; border-left: 4px solid #856404; margin-bottom: 15px;">'
                html += '<p style="margin: 0; font-size: 13px; color: #856404;"><b>Simple Character Alignment</b></p>'
                html += '<p style="margin: 5px 0 0 0; font-size: 12px; color: #856404;">CIGAR data not available. Showing raw sequences without gap insertion.</p>'
                html += '</div>'

                # Show amino acid alignment if available
                if ref_aa:
                    html += '<h2>Amino Acid Alignment</h2>'
                    html += '<p style="color: #666; font-size: 14px;">Scroll horizontally | Green = match, Red = mismatch</p>'
                    try:
                        asm_aa = str(Seq(asm_seq).translate())
                        html += self._format_alignment_html(ref_aa, asm_aa, "Reference AA", "Assembly AA")
                    except:
                        pass

                # Show nucleotide alignment
                html += '<h2>Nucleotide Alignment</h2>'
                html += '<p style="color: #666; font-size: 14px;">Scroll horizontally | Green = match, Red = mismatch</p>'
                html += self._format_alignment_html(ref_seq, asm_seq, "Reference", "Assembly")

            html += '</div>'
        else:
            html += '<div class="alignment-container"><p style="color: red;">No assembly coverage for this gene.</p></div>'

        html += """
        </body>
        </html>
        """

        return html

    def _parse_cigar_alignment(self, ref_seq, query_seq, cigar, ref_start, query_start, strand):
        """
        Parse CIGAR string to create properly aligned sequences with gaps.

        Args:
            ref_seq (str): Reference sequence
            query_seq (str): Query sequence
            cigar (list): CIGAR operations as list of tuples [(length, op), ...]
            ref_start (int): Start position in reference
            query_start (int): Start position in query
            strand (str): '+' or '-'

        Returns:
            tuple: (aligned_ref, aligned_query, match_string)
        """
        # CIGAR operations: 0=M, 1=I, 2=D, 3=N, 4=S, 5=H, 7==, 8=X
        if cigar is None or len(cigar) == 0:
            return ref_seq, query_seq, ''

        # Reverse complement query if needed
        if strand == '-':
            query_seq = str(Seq(query_seq).reverse_complement())

        aligned_ref = []
        aligned_query = []
        match_string = []

        ref_pos = 0
        query_pos = 0

        for length, op in cigar:
            if op == 0 or op == 7 or op == 8:  # M, =, X (match/mismatch)
                for i in range(length):
                    if ref_pos < len(ref_seq) and query_pos < len(query_seq):
                        r = ref_seq[ref_pos]
                        q = query_seq[query_pos]
                        aligned_ref.append(r)
                        aligned_query.append(q)
                        match_string.append('|' if r == q else ' ')
                        ref_pos += 1
                        query_pos += 1

            elif op == 1:  # I (insertion to reference)
                for i in range(length):
                    if query_pos < len(query_seq):
                        aligned_ref.append('-')
                        aligned_query.append(query_seq[query_pos])
                        match_string.append(' ')
                        query_pos += 1

            elif op == 2:  # D (deletion from reference)
                for i in range(length):
                    if ref_pos < len(ref_seq):
                        aligned_ref.append(ref_seq[ref_pos])
                        aligned_query.append('-')
                        match_string.append(' ')
                        ref_pos += 1

            elif op == 4:  # S (soft clip) - skip in query
                query_pos += length

            elif op == 5:  # H (hard clip) - already removed
                pass

            elif op == 3:  # N (skipped region in reference)
                ref_pos += length

        return ''.join(aligned_ref), ''.join(aligned_query), ''.join(match_string)

    def _format_alignment_html(self, seq1, seq2, label1, label2, chars_per_line=100):
        """Format two sequences as horizontally scrollable alignment with coloring."""
        min_len = min(len(seq1), len(seq2))

        # Build colored sequences
        colored1 = ""
        colored2 = ""
        match_line = ""

        # Calculate identity
        matches = 0
        for i, (c1, c2) in enumerate(zip(seq1[:min_len], seq2[:min_len])):
            if c1 == c2:
                colored1 += f'<span class="match">{c1}</span>'
                colored2 += f'<span class="match">{c2}</span>'
                match_line += '|'
                matches += 1
            else:
                colored1 += f'<span class="mismatch">{c1}</span>'
                colored2 += f'<span class="mismatch">{c2}</span>'
                match_line += ' '

        # Handle sequence length differences
        if len(seq1) > min_len:
            for c in seq1[min_len:]:
                colored1 += f'<span class="gap">{c}</span>'
            colored2 += '<span class="gap">' + '-' * (len(seq1) - min_len) + '</span>'
            match_line += ' ' * (len(seq1) - min_len)

        if len(seq2) > min_len:
            colored1 += '<span class="gap">' + '-' * (len(seq2) - min_len) + '</span>'
            for c in seq2[min_len:]:
                colored2 += f'<span class="gap">{c}</span>'
            match_line += ' ' * (len(seq2) - min_len)

        # Calculate identity
        total_len = max(len(seq1), len(seq2))
        simple_identity = (matches / min_len * 100) if min_len > 0 else 0

        html = f'''
        <div style="overflow-x: auto; max-width: 100%; border: 1px solid #dee2e6; padding: 10px; background-color: #f8f9fa;">
            <div class="seq-line"><span class="label">{label1}:</span> {colored1}</div>
            <div class="seq-line"><span class="label">{"Match:":>15}</span> {match_line}</div>
            <div class="seq-line"><span class="label">{label2}:</span> {colored2}</div>
            <div style="margin-top: 10px; font-size: 12px; color: #666;">
                <b>Alignment stats:</b> {matches} matches, {min_len - matches} mismatches/gaps out of {max(len(seq1), len(seq2))} positions
                ({simple_identity:.2f}% identity excluding gaps)
            </div>
        </div>
        '''
        return html
