"""
New Gene Analyzer - Identify and classify new genes in scaffolds compared to contigs.

Integrated module for GenomeViz that:
- Finds genes in scaffolds not present in contigs (by protein hash)
- Classifies them based on contig-scaffold alignment overlap
- Generates CSV reports and visualization files
- Uses rotated genome coordinates from GenomeViz pipeline

GitHub: https://github.com/Aaron-Thiel/GenomeViz
License: MIT
"""

import os
import hashlib
import pandas as pd
import mappy as mp
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq


class NewGeneAnalyzer:
    """
    Identify and classify new genes in scaffolds compared to contigs.
    
    A gene is "new" if its protein sequence (by MD5 hash) exists in scaffolds
    but not in contigs. Classification depends on alignment overlap patterns.
    """

    def __init__(self, scaffolds_fasta, contigs_fasta, scaffolds_gff, contigs_gff,
                 output_dir, origin_offset=0, preset='asm5', scaffold_alignments=None):
        """
        Initialize new gene analyzer.

        Args:
            scaffolds_fasta (str): Path to scaffold sequences
            contigs_fasta (str): Path to contig sequences
            scaffolds_gff (str): Path to scaffold GFF3 annotations
            contigs_gff (str): Path to contig GFF3 annotations
            output_dir (str/Path): Base output directory
            origin_offset (int): Coordinate offset from genome rotation (default: 0)
            preset (str): Minimap2 preset for contig-to-scaffold alignment
            scaffold_alignments (dict): Scaffold-to-reference alignments for coordinate mapping
        """
        self.scaffolds_fasta = scaffolds_fasta
        self.contigs_fasta = contigs_fasta
        self.scaffolds_gff = scaffolds_gff
        self.contigs_gff = contigs_gff
        self.output_dir = Path(output_dir)
        self.origin_offset = origin_offset
        self.preset = preset
        self.scaffold_alignments = scaffold_alignments or {}

        # Create output subfolders
        self.new_genes_dir = self.output_dir / 'new_genes'
        self.visualizations_dir = self.new_genes_dir / 'visualizations'
        self.new_genes_dir.mkdir(parents=True, exist_ok=True)
        self.visualizations_dir.mkdir(parents=True, exist_ok=True)

        # Build scaffold alignment index for coordinate mapping
        self._build_scaffold_alignment_index()

    def _build_scaffold_alignment_index(self):
        """
        Build an index mapping scaffold names to their reference alignments.

        This is used to convert scaffold-local gene coordinates to reference coordinates.
        Uses only primary alignments for consistency with visualizations.
        """
        from collections import defaultdict

        self.scaffold_to_ref_map = defaultdict(list)

        for ref_seqid, alignments in self.scaffold_alignments.items():
            for aln in alignments:
                # Only use primary alignments
                if not aln.get('is_primary', True):
                    continue

                scaffold_name = aln.get('query_name', '')
                if scaffold_name:
                    self.scaffold_to_ref_map[scaffold_name].append({
                        'ref_seqid': ref_seqid,
                        'ref_start': aln.get('ref_start', 0),
                        'ref_end': aln.get('ref_end', 0),
                        'scaffold_local_start': aln.get('query_start', 0),
                        'scaffold_local_end': aln.get('query_end', 0),
                        'strand': aln.get('strand', '+')
                    })

                    # Also index by name without _RC suffix
                    if scaffold_name.endswith('_RC'):
                        original_name = scaffold_name[:-3]
                        self.scaffold_to_ref_map[original_name].append({
                            'ref_seqid': ref_seqid,
                            'ref_start': aln.get('ref_start', 0),
                            'ref_end': aln.get('ref_end', 0),
                            'scaffold_local_start': aln.get('query_start', 0),
                            'scaffold_local_end': aln.get('query_end', 0),
                            'strand': aln.get('strand', '+')
                        })

    def _map_gene_to_reference(self, gene_start, gene_end, scaffold_name):
        """
        Map gene coordinates from scaffold-local to reference genome coordinates.

        Args:
            gene_start: Gene start position on scaffold (local coordinates)
            gene_end: Gene end position on scaffold (local coordinates)
            scaffold_name: Name of the scaffold

        Returns:
            Tuple of (ref_start, ref_end, ref_seqid) or (None, None, None) if unmapped
        """
        alignments = self.scaffold_to_ref_map.get(scaffold_name, [])

        if not alignments:
            # Try with _RC suffix
            alignments = self.scaffold_to_ref_map.get(f"{scaffold_name}_RC", [])

        if not alignments:
            return None, None, None

        # Find the alignment region that contains (or best overlaps) this gene
        best_alignment = None
        best_overlap = 0

        for aln in alignments:
            scaffold_local_start = aln['scaffold_local_start']
            scaffold_local_end = aln['scaffold_local_end']

            # Calculate overlap between gene and this alignment region
            overlap_start = max(gene_start, scaffold_local_start)
            overlap_end = min(gene_end, scaffold_local_end)
            overlap_bp = max(0, overlap_end - overlap_start)

            if overlap_bp > best_overlap:
                best_overlap = overlap_bp
                best_alignment = aln

        if not best_alignment or best_overlap == 0:
            return None, None, None

        # Map gene coordinates to reference using the best alignment
        aln = best_alignment

        # Calculate offset: reference_pos = gene_pos + offset
        offset = aln['ref_start'] - aln['scaffold_local_start']

        ref_start = gene_start + offset
        ref_end = gene_end + offset
        ref_seqid = aln['ref_seqid']

        return ref_start, ref_end, ref_seqid

    def get_gene_hashes(self, faa_file):
        """Returns a set of protein MD5 hashes from a .faa file."""
        hashes = set()
        if os.path.exists(faa_file):
            for record in SeqIO.parse(faa_file, "fasta"):
                seq_str = str(record.seq).strip().upper()
                if seq_str:
                    hashes.add(hashlib.md5(seq_str.encode('utf-8')).hexdigest())
        return hashes

    def get_seq_data(self, fna_file):
        """Returns a dict of {id: sequence_string} from a fasta file."""
        seqs = {}
        if os.path.exists(fna_file):
            for record in SeqIO.parse(fna_file, "fasta"):
                seqs[record.id] = str(record.seq).upper()
        return seqs

    def parse_scaffold_genes(self, gff_file, faa_file):
        """Parses scaffold GFF3 and FAA to get gene details."""
        genes = []
        seq_map = {}
        
        if os.path.exists(faa_file):
            for record in SeqIO.parse(faa_file, "fasta"):
                seq_id = record.id.split()[0]
                seq_map[seq_id] = str(record.seq).strip().upper()
        
        if not os.path.exists(gff_file):
            return genes
            
        with open(gff_file, 'r') as f:
            for line in f:
                if line.startswith("#") or line.strip() == "":
                    continue
                parts = line.strip().split('\t')
                if len(parts) != 9:
                    continue
                if parts[2] == "CDS":
                    attrs = {x.split('=')[0]: x.split('=')[1] for x in parts[8].split(';') if '=' in x}
                    gene_id = attrs.get('ID', '')
                    prot_seq = seq_map.get(gene_id, "")
                    genes.append({
                        'id': gene_id,
                        'scaffold': parts[0],
                        'start': int(parts[3]),
                        'end': int(parts[4]),
                        'strand': parts[6],
                        'product': attrs.get('product', 'hypothetical protein'),
                        'hash': hashlib.md5(prot_seq.encode('utf-8')).hexdigest() if prot_seq else None,
                        'prot_seq': prot_seq
                    })
        return genes

    def parse_cigar_to_aligned_blocks(self, cigar_str, ref_start, query_start):
        """
        Parse CIGAR string to identify actual aligned blocks in reference coordinates.
        
        This properly handles gaps, soft-clips, and insertions in the alignment.
        Returns list of (ref_block_start, ref_block_end, query_block_start, query_block_end) tuples
        representing contiguous aligned regions.
        
        Args:
            cigar_str: CIGAR string from minimap2 (e.g., "100M5D95M")
            ref_start: Starting position on reference
            query_start: Starting position on query
            
        Returns:
            List of tuples: (ref_start, ref_end, query_start, query_end) for each aligned block
        """
        import re
        
        # If no CIGAR string, return empty list - will be handled by caller
        if not cigar_str:
            return []
        
        blocks = []
        ref_pos = ref_start
        query_pos = query_start
        block_start_ref = ref_start
        block_start_query = query_start
        in_aligned_block = False
        
        # Parse CIGAR string: extract operation and length
        for match in re.finditer(r'(\d+)([MIDNSHPX=])', cigar_str):
            length = int(match.group(1))
            op = match.group(2)
            
            # Operations: M=match/mismatch, I=insertion, D=deletion, N=skip, S=soft-clip, H=hard-clip, P=pad, X=mismatch, ==match
            if op in ['M', '=', 'X']:  # Aligned region (match, mismatch, or sequence match)
                if not in_aligned_block:
                    block_start_ref = ref_pos
                    block_start_query = query_pos
                    in_aligned_block = True
                ref_pos += length
                query_pos += length
                
            elif op == 'D' or op == 'N':  # Deletion or skip in query (gap in query)
                if in_aligned_block:
                    # End current block
                    blocks.append((block_start_ref, ref_pos, block_start_query, query_pos))
                    in_aligned_block = False
                ref_pos += length
                
            elif op == 'I':  # Insertion in query (gap in reference)
                if in_aligned_block:
                    # Insertions are still part of the aligned block but don't advance ref_pos
                    query_pos += length
                else:
                    query_pos += length
                
            elif op in ['S', 'H']:  # Soft-clip or hard-clip
                if in_aligned_block:
                    # End current block
                    blocks.append((block_start_ref, ref_pos, block_start_query, query_pos))
                    in_aligned_block = False
                # For soft-clips, advance query_pos; for hard-clips, don't
                if op == 'S':
                    query_pos += length
        
        # Close any remaining open block
        if in_aligned_block:
            blocks.append((block_start_ref, ref_pos, block_start_query, query_pos))
        
        return blocks

    def map_contigs_to_scaffolds(self):
        """
        Align contigs to scaffolds using mappy (minimap2 Python binding).
        
        Properly handles CIGAR strings to identify actual aligned regions,
        ensuring that soft-clipped or gapped regions are not counted as coverage.
        
        Returns:
            dict: Mappings grouped by scaffold ID, with list of aligned blocks per contig
        """
        print("\n  Aligning contigs to scaffolds...")
        print(f"    Using preset: {self.preset}")
        mappings = defaultdict(list)
        
        try:
            # Load contigs as REFERENCE, map scaffolds as QUERIES
            # This way we get scaffold gene coordinates mapped back to contig space
            aligner = mp.Aligner(self.contigs_fasta, preset=self.preset)
            if not aligner:
                print("  ERROR: Failed to load contigs for alignment")
                return mappings
            
            print(f"    Aligner initialized successfully")
            
            scaffold_count = 0
            hit_count = 0
            block_count = 0
            for name, seq, qual in mp.fastx_read(self.scaffolds_fasta):
                scaffold_count += 1
                hits_for_scaffold = list(aligner.map(seq))
                if scaffold_count <= 5:
                    print(f"    Scaffold {scaffold_count} ({name}): {len(hits_for_scaffold)} hits (seq_len={len(seq)})")
                
                for hit in hits_for_scaffold:
                    hit_count += 1
                    # Extract alignment statistics
                    matches = hit.mlen
                    block_len = hit.blen
                    identity = matches / block_len if block_len > 0 else 0

                    # Calculate alignment fraction
                    s_len = len(seq)  # scaffold length (query)
                    c_len = hit.ctg_len  # contig length (reference)
                    s_start = hit.q_st
                    s_end = hit.q_en
                    alignment_fraction = (s_end - s_start) / s_len if s_len > 0 else 0
                    
                    # Parse CIGAR to get actual aligned blocks (handles gaps, soft-clips, etc.)
                    # Note: If CIGAR is not available, use the full alignment span as a single block
                    if hit.cigar_str:
                        cigar_blocks = self.parse_cigar_to_aligned_blocks(hit.cigar_str, hit.r_st, hit.q_st)
                    else:
                        # No CIGAR available - treat entire alignment as one block
                        cigar_blocks = [(hit.r_st, hit.r_en, hit.q_st, hit.q_en)]
                    
                    if len(cigar_blocks) == 0:
                        # If still empty, force at least one block
                        cigar_blocks = [(hit.r_st, hit.r_en, hit.q_st, hit.q_en)]
                    
                    if scaffold_count <= 3 and hit_count <= 5:
                        print(f"      Hit {hit_count}: {name} -> {hit.ctg}[{hit.r_st}:{hit.r_en}] cigar={hit.cigar_str} blocks={len(cigar_blocks)}")
                    
                    # Store each aligned block separately for accurate overlap detection
                    for ref_block_start, ref_block_end, query_block_start, query_block_end in cigar_blocks:
                        block_count += 1
                        # Tuple: (scaffold_start, scaffold_end, contig_id, contig_start, contig_end, contig_len, strand, identity, alignment_fraction)
                        # Note: ref is CONTIG (reference in mappy), query is SCAFFOLD (query in mappy)
                        mappings[name].append((
                            query_block_start, query_block_end, hit.ctg, ref_block_start, ref_block_end,
                            c_len, '+' if hit.strand == 1 else '-', identity, alignment_fraction
                        ))
                        if scaffold_count <= 3 and hit_count <= 5:
                            print(f"        -> Block {block_count}: scaffold[{query_block_start}:{query_block_end}] contig[{ref_block_start}:{ref_block_end}]")
            
            print(f"    Processed {scaffold_count} scaffolds, found {hit_count} hits with {block_count} aligned blocks")
            total_align_blocks = sum(len(v) for v in mappings.values())
            print(f"  Found {total_align_blocks} total aligned blocks across {len(mappings)} scaffolds")
            if len(mappings) > 0:
                print(f"    Scaffolds with alignments: {sorted(mappings.keys())[:5]}...")
                for scaffold_name in sorted(mappings.keys())[:2]:
                    print(f"      {scaffold_name}: {len(mappings[scaffold_name])} blocks")
        except Exception as e:
            print(f"  ERROR during alignment: {e}")
            import traceback
            traceback.print_exc()
        
        return mappings

    def get_translated_context(self, scaffold_seq, start_bp, end_bp, strand):
        """Translates a specific window of the scaffold."""
        if start_bp >= end_bp or start_bp < 0:
            return ""
        
        sub_seq = Seq(scaffold_seq[start_bp:end_bp])
        if strand == '-':
            sub_seq = sub_seq.reverse_complement()
        
        translated = str(sub_seq[:(len(sub_seq)//3)*3].translate())
        return translated

    def generate_visualization(self, gene, overlapping_contigs, full_scaffold_seq, contig_seqs):
        """
        Generates both nucleotide and amino acid alignment visualizations with detailed positional context.
        Shows aligned parts (UPPER) and soft-clips (lower) in a stacked layout.

        Returns:
            tuple: (visualization_string, mapping_stats_string)
        """
        prot_seq = gene['prot_seq']
        if not prot_seq:
            return "No protein sequence available", ""

        scaffold_len = len(full_scaffold_seq)
        g_start_0, g_end_0, strand = gene['start'] - 1, gene['end'], gene['strand']

        context_bp = 30
        win_start = max(0, g_start_0 - context_bp)
        win_end = min(scaffold_len, g_end_0 + context_bp)

        # Get both nucleotide and amino acid context windows
        win_nt_seq = Seq(full_scaffold_seq[win_start:win_end])
        if strand == '-':
            win_nt_seq = win_nt_seq.reverse_complement()
        win_nt = str(win_nt_seq)

        win_aa = self.get_translated_context(full_scaffold_seq, win_start, win_end, strand)
        win_len_aa = len(win_aa)
        win_len_nt = len(win_nt)

        if win_len_aa == 0:
            return "Could not generate visualization", ""

        vis_lines_aa = []
        vis_lines_nt = []
        mapping_stats = []
        contig_details = []  # Store detailed info for header
        
        for s_start, s_end, c_id, c_start, c_end, c_len, c_strand, identity, align_frac in sorted(overlapping_contigs, key=lambda x: x[0]):
            line_chars_aa = [" "] * win_len_aa
            line_chars_nt = [" "] * win_len_nt

            mapping_stats.append(f"{c_id}: align: {align_frac:.2f}, id: {identity:.2f}")

            # Generate detailed info for header
            contig_info = (f"  {c_id}: Length={c_len:,}bp, Aligned contig[{c_start:,}-{c_end:,}] to scaffold[{s_start:,}-{s_end:,}], "
                          f"Strand={c_strand}, Identity={identity:.1%}")
            contig_details.append(contig_info)

            # Find overlap between window and alignment
            intersect_start = max(win_start, s_start)
            intersect_end = min(win_end, s_end)

            # Calculate how much of the contig alignment extends beyond the visible window
            # For LEFT: How many bp of the alignment are BEFORE the window starts?
            if s_start < win_start:
                left_extend_aa = (win_start - s_start) // 3
                left_extend_nt = win_start - s_start
            else:
                left_extend_aa = 0
                left_extend_nt = 0

            # For RIGHT: How many bp of the alignment continue AFTER the window ends?
            if s_end > win_end:
                right_extend_aa = (s_end - win_end) // 3
                right_extend_nt = s_end - win_end
            else:
                right_extend_aa = 0
                right_extend_nt = 0

            # Track if alignment starts/ends within the visible window
            starts_in_window = s_start >= win_start
            ends_in_window = s_end <= win_end

            if intersect_start < intersect_end:
                # Get the actual contig sequence for the aligned region
                # Try with original ID first, then without _RC suffix as fallback
                full_c_seq = contig_seqs.get(c_id, "")

                if not full_c_seq and c_id.endswith('_RC'):
                    # Fallback: get original contig and RC it to match what was aligned
                    base_id = c_id[:-3]
                    full_c_seq = contig_seqs.get(base_id, "")
                    if full_c_seq:
                        # RC the entire sequence to simulate the oriented _RC version
                        full_c_seq = str(Seq(full_c_seq).reverse_complement())

                # Track if we're using the scaffold-based fallback (already display-ordered)
                used_scaffold_fallback = False

                if full_c_seq:
                    # Calculate which part of the contig corresponds to the visible intersection
                    # intersect_start/end are in scaffold coordinates
                    # We need to map back to contig coordinates

                    # Offset from alignment start (in scaffold coords) to intersection start
                    offset_from_align_start = intersect_start - s_start
                    visible_length = intersect_end - intersect_start

                    # Map to contig coordinates - depends on alignment strand
                    if c_strand == '-':
                        # For minus strand: scaffold forward = contig backward
                        # c_end corresponds to s_start, c_start corresponds to s_end
                        contig_visible_end = c_end - offset_from_align_start
                        contig_visible_start = contig_visible_end - visible_length
                    else:
                        # For plus strand: scaffold forward = contig forward
                        contig_visible_start = c_start + offset_from_align_start
                        contig_visible_end = contig_visible_start + visible_length

                    # Ensure coordinates are within bounds
                    contig_visible_start = max(0, contig_visible_start)
                    contig_visible_end = min(len(full_c_seq), contig_visible_end)

                    # Extract the contig sequence for this region
                    contig_region_seq = Seq(full_c_seq[contig_visible_start:contig_visible_end])

                    # Apply strand transformation if needed
                    if c_strand == '-':
                        contig_region_seq = contig_region_seq.reverse_complement()

                    contig_nt = str(contig_region_seq).upper()

                    # For amino acid: use the scaffold's AA at the corresponding positions
                    # This ensures correct reading frame alignment since nucleotides match
                    # Calculate the AA range that corresponds to this region
                    nt_start_in_window = intersect_start - win_start
                    nt_end_in_window = intersect_end - win_start

                    # For minus strand genes, win_aa is already in display order (5'->3' gene direction)
                    # We need to convert scaffold coordinates to display coordinates
                    if strand == '-':
                        # Flip coordinates for display space
                        nt_start_display = win_len_nt - nt_end_in_window
                        nt_end_display = win_len_nt - nt_start_in_window
                    else:
                        nt_start_display = nt_start_in_window
                        nt_end_display = nt_end_in_window

                    # Find the codon boundaries that overlap with this NT region
                    # First full codon that starts at or after nt_start
                    first_codon_idx = (nt_start_display + 2) // 3  # Round up to next codon
                    # Last codon that ends at or before nt_end
                    last_codon_idx = nt_end_display // 3

                    if first_codon_idx < last_codon_idx and last_codon_idx <= len(win_aa):
                        contig_aa = win_aa[first_codon_idx:last_codon_idx]
                    else:
                        contig_aa = ""
                else:
                    # Fallback to scaffold sequence if contig not found
                    # Since win_nt is already in display order, we extract directly
                    used_scaffold_fallback = True
                    nt_start_in_window = intersect_start - win_start
                    nt_end_in_window = intersect_end - win_start

                    # For minus strand, use display coordinates for extraction from win_nt/win_aa
                    if strand == '-':
                        nt_start_display = win_len_nt - nt_end_in_window
                        nt_end_display = win_len_nt - nt_start_in_window
                        contig_nt = win_nt[nt_start_display:nt_end_display].upper()
                        first_codon_idx = (nt_start_display + 2) // 3
                        last_codon_idx = nt_end_display // 3
                    else:
                        contig_nt = win_nt[nt_start_in_window:nt_end_in_window].upper()
                        first_codon_idx = (nt_start_in_window + 2) // 3
                        last_codon_idx = nt_end_in_window // 3
                    contig_aa = win_aa[first_codon_idx:last_codon_idx] if first_codon_idx < last_codon_idx else ""

                # Position in the display window
                # For minus strand, the scaffold display is RC'd so we need to:
                # 1. Flip positions (original end -> display start)
                # 2. RC the contig sequence to match the RC'd scaffold display (only if from actual contig)
                if strand == '-':
                    # In RC'd display: intersection maps to flipped position
                    rel_start_nt = win_len_nt - (intersect_end - win_start)
                    # RC the contig sequence to match the RC'd scaffold display
                    # But NOT if we used the scaffold fallback (already in display order)
                    if not used_scaffold_fallback:
                        contig_nt = str(Seq(contig_nt).reverse_complement())
                    # Note: contig_aa is already extracted in display order (see above),
                    # so no reversal is needed here
                else:
                    rel_start_nt = intersect_start - win_start

                rel_start_aa = (rel_start_nt + 2) // 3  # First full codon index

                # Fill in nucleotide sequence
                for i, nt in enumerate(contig_nt):
                    idx = rel_start_nt + i
                    if 0 <= idx < win_len_nt:
                        line_chars_nt[idx] = nt

                # Fill in amino acid sequence
                for i, aa in enumerate(contig_aa):
                    idx = rel_start_aa + i
                    if 0 <= idx < win_len_aa:
                        line_chars_aa[idx] = aa

            # Add unaligned contig regions (overhang beyond alignment)
            # Calculate contig overhang amounts
            left_overhang_bp = c_start  # bp of contig before aligned region
            right_overhang_bp = c_len - c_end  # bp of contig after aligned region

            # Get contig sequence (try with _RC suffix fallback)
            full_c_seq = contig_seqs.get(c_id, "")
            if not full_c_seq and c_id.endswith('_RC'):
                full_c_seq = contig_seqs.get(c_id[:-3], "")
                if full_c_seq:
                    # RC the sequence to match the oriented _RC version
                    full_c_seq = str(Seq(full_c_seq).reverse_complement())

            # Track how much contig extends beyond the visible window
            left_beyond_window_nt = 0
            right_beyond_window_nt = 0

            if full_c_seq:
                # Left overhang: contig sequence before the aligned region
                if left_overhang_bp > 0:
                    clip_seq = Seq(full_c_seq[:c_start])
                    if c_strand == '-':
                        clip_seq = clip_seq.reverse_complement()
                    clip_nt = str(clip_seq).lower()
                    clip_aa = str(clip_seq[:(len(clip_seq)//3)*3].translate()).lower()

                    # Calculate how much fits in the window vs extends beyond
                    space_in_window_left = s_start - win_start  # scaffold space available to the left

                    if left_overhang_bp > space_in_window_left:
                        # Some contig extends beyond the window
                        left_beyond_window_nt = left_overhang_bp - space_in_window_left
                        # Only show the part that fits in window
                        clip_nt_visible = clip_nt[left_beyond_window_nt:]
                        clip_aa_visible = clip_aa[left_beyond_window_nt // 3:]
                    else:
                        clip_nt_visible = clip_nt
                        clip_aa_visible = clip_aa

                    # Fill in visible overhang (right-aligned to alignment start)
                    for i, nt in enumerate(reversed(clip_nt_visible)):
                        idx = (s_start - win_start) - 1 - i
                        if 0 <= idx < win_len_nt:
                            line_chars_nt[idx] = nt

                    for i, aa in enumerate(reversed(clip_aa_visible)):
                        idx = (s_start - win_start) // 3 - 1 - i
                        if 0 <= idx < win_len_aa:
                            line_chars_aa[idx] = aa

                # Right overhang: contig sequence after the aligned region
                if right_overhang_bp > 0:
                    clip_seq = Seq(full_c_seq[c_end:])
                    if c_strand == '-':
                        clip_seq = clip_seq.reverse_complement()
                    clip_nt = str(clip_seq).lower()
                    clip_aa = str(clip_seq[:(len(clip_seq)//3)*3].translate()).lower()

                    # Calculate how much fits in the window vs extends beyond
                    space_in_window_right = win_end - s_end  # scaffold space available to the right

                    if right_overhang_bp > space_in_window_right:
                        # Some contig extends beyond the window
                        right_beyond_window_nt = right_overhang_bp - space_in_window_right
                        # Only show the part that fits in window
                        clip_nt_visible = clip_nt[:space_in_window_right]
                        clip_aa_visible = clip_aa[:space_in_window_right // 3]
                    else:
                        clip_nt_visible = clip_nt
                        clip_aa_visible = clip_aa

                    # Fill in visible overhang (left-aligned from alignment end)
                    for i, nt in enumerate(clip_nt_visible):
                        idx = (s_end - win_start) + i
                        if 0 <= idx < win_len_nt:
                            line_chars_nt[idx] = nt

                    for i, aa in enumerate(clip_aa_visible):
                        idx = (s_end - win_start) // 3 + i
                        if 0 <= idx < win_len_aa:
                            line_chars_aa[idx] = aa

            # Add < and > markers to show alignment boundaries within the window
            # < marks start of aligned region, > marks end of aligned region
            if starts_in_window:
                # Insert < at the start of the aligned region
                start_pos_nt = s_start - win_start
                start_pos_aa = start_pos_nt // 3
                if 0 <= start_pos_nt < win_len_nt:
                    line_chars_nt[start_pos_nt] = '<'
                if 0 <= start_pos_aa < win_len_aa:
                    line_chars_aa[start_pos_aa] = '<'

            if ends_in_window:
                # Insert > at the end of the aligned region
                end_pos_nt = s_end - win_start - 1
                end_pos_aa = (s_end - win_start - 1) // 3
                if 0 <= end_pos_nt < win_len_nt:
                    line_chars_nt[end_pos_nt] = '>'
                if 0 <= end_pos_aa < win_len_aa:
                    line_chars_aa[end_pos_aa] = '>'

            # Format the lines with proper prefixes and suffixes
            total_left_extend_nt = left_extend_nt + left_beyond_window_nt
            total_left_extend_aa = left_extend_aa + (left_beyond_window_nt // 3)

            prefix_nt = f"({total_left_extend_nt}...)".ljust(12) if total_left_extend_nt > 0 else " " * 12
            prefix_aa = f"({total_left_extend_aa}...)".ljust(12) if total_left_extend_aa > 0 else " " * 12

            seq_aa = ''.join(line_chars_aa)
            seq_nt = ''.join(line_chars_nt)

            # Suffix shows: alignment extension + contig overhang beyond window
            total_right_extend_nt = right_extend_nt + right_beyond_window_nt
            total_right_extend_aa = right_extend_aa + (right_beyond_window_nt // 3)

            suffix_nt = f" (...{total_right_extend_nt})" if total_right_extend_nt > 0 else ""
            suffix_aa = f" (...{total_right_extend_aa})" if total_right_extend_aa > 0 else ""

            vis_lines_aa.append(f"{c_id:<25} {prefix_aa}{seq_aa}{suffix_aa}")
            vis_lines_nt.append(f"{c_id:<25} {prefix_nt}{seq_nt}{suffix_nt}")
        
        # Add scaffold sequence lines (both AA and NT)
        s_chars_aa = list(win_aa)
        s_chars_nt = list(win_nt)

        # Mark gene boundaries in AA
        g_rel_s_aa = (g_start_0 - win_start) // 3
        g_rel_e_aa = (g_end_0 - win_start) // 3

        if 0 <= g_rel_s_aa < win_len_aa:
            s_chars_aa[g_rel_s_aa] = "|"
        if 0 <= (g_rel_e_aa - 1) < win_len_aa:
            s_chars_aa[g_rel_e_aa - 1] = "|"

        # Mark gene boundaries in NT
        g_rel_s_nt = g_start_0 - win_start
        g_rel_e_nt = g_end_0 - win_start

        if 0 <= g_rel_s_nt < win_len_nt:
            s_chars_nt[g_rel_s_nt] = "|"
        if 0 <= (g_rel_e_nt - 1) < win_len_nt:
            s_chars_nt[g_rel_e_nt - 1] = "|"

        # Scaffold context - use fixed-width prefix for alignment
        left_context = win_start
        right_context = scaffold_len - win_end

        # Fixed-width prefix (12 characters) for alignment with contig lines
        if left_context > 0:
            s_prefix_aa = f"({left_context}...)".ljust(12)
            s_prefix_nt = f"({left_context}...)".ljust(12)
        else:
            s_prefix_aa = " " * 12
            s_prefix_nt = " " * 12

        s_suffix_aa = f" (...{right_context})" if right_context > 0 else ""
        s_suffix_nt = f" (...{right_context})" if right_context > 0 else ""

        vis_lines_aa.append(f"{'Scaffold (Context)':<25} {s_prefix_aa}{''.join(s_chars_aa)}{s_suffix_aa}")
        vis_lines_nt.append(f"{'Scaffold (Context)':<25} {s_prefix_nt}{''.join(s_chars_nt)}{s_suffix_nt}")

        # Build header with contig details
        header = []
        if contig_details:
            header.append("Contig alignment details:")
            header.extend(contig_details)
            header.append("")

        # Combine both visualizations in stacked format
        final_vis = []

        # Header with contig details
        if header:
            final_vis.extend(header)

        # Amino acid visualization
        final_vis.append("Amino Acid Alignment:")
        final_vis.append("-" * 100)
        final_vis.extend(vis_lines_aa)
        final_vis.append("")

        # Nucleotide visualization
        final_vis.append("Nucleotide Alignment:")
        final_vis.append("-" * 100)
        final_vis.extend(vis_lines_nt)

        return "\n".join(final_vis), "; ".join(mapping_stats)

    def classify_new_gene(self, gene, mappings, full_scaffold_seq, contig_seqs, contig_genes=None):
        """
        Classify a new gene based on alignment overlap patterns with CIGAR-aware blocks.
        
        Now properly handles gaps and soft-clips in alignments. Only regions that are
        actually aligned (not soft-clipped or in gaps) are considered for overlap.
        
        Maps scaffold coordinates back to contig coordinates to find which contig genes
        may have been combined or re-predicted to form this new scaffold gene.
        
        Classifications:
        - "Other Source": No contig overlap or unaligned region (gene created by assembly)
        - "Contextual Re-prediction": Single contig overlaps the gene region
        - "Combined Contigs (N joined)": Multiple contigs overlap the gene region
        
        Returns:
            tuple: (classification, joined_contigs_str, visualization, stats)
        """
        scaffold = gene['scaffold']
        g_s, g_e = gene['start'] - 1, gene['end']
        
        # Get alignments for this scaffold
        scaffold_alignments = mappings.get(scaffold, [])
        
        # Find overlapping alignments (now based on CIGAR-parsed aligned blocks)
        ovl = [m for m in scaffold_alignments if max(g_s, m[0]) < min(g_e, m[1])]
        
        # Collect information about overlapping contigs
        contig_gene_sources = {}  # Map: contig_name -> list of genes in the overlapping region
        
        if ovl and contig_genes:
            # For each overlapping alignment block, map coordinates back to contig space
            for aln in ovl:
                ref_block_start, ref_block_end, contig_id, query_block_start, query_block_end, c_len, strand, identity, align_frac = aln

                # Map the overlapping region [g_s:g_e] from scaffold to contig coordinates
                overlap_start = max(g_s, ref_block_start)
                overlap_end = min(g_e, ref_block_end)

                if overlap_start < overlap_end:
                    # Calculate offset: how much into the alignment block is the overlap start?
                    offset_in_block = overlap_start - ref_block_start

                    # Map to contig coordinates
                    contig_overlap_start = query_block_start + offset_in_block
                    contig_overlap_end = contig_overlap_start + (overlap_end - overlap_start)

                    if contig_id not in contig_gene_sources:
                        contig_gene_sources[contig_id] = []

                    # Check which contig genes overlap with this region
                    for contig_gene in contig_genes:
                        if contig_gene['scaffold'] == contig_id:
                            c_g_s, c_g_e = contig_gene['start'] - 1, contig_gene['end']
                            if max(c_g_s, contig_overlap_start) < min(c_g_e, contig_overlap_end):
                                contig_gene_sources[contig_id].append({
                                    'id': contig_gene['id'],
                                    'product': contig_gene['product'],
                                    'start': contig_gene['start'],
                                    'end': contig_gene['end']
                                })
        
        visual, stats = self.generate_visualization(gene, ovl, full_scaffold_seq, contig_seqs)
        
        # Determine classification based on what we found
        all_contributing_contigs = set(contig_gene_sources.keys())
        
        if not all_contributing_contigs:
            # No contigs overlap this region
            return "Other Source", "", visual, ""
        elif len(all_contributing_contigs) == 1:
            # Single contig overlaps
            return "Contextual Re-prediction", list(all_contributing_contigs)[0], visual, stats
        else:
            # Multiple contigs overlap
            return f"Combined Contigs ({len(all_contributing_contigs)} joined)", ", ".join(sorted(all_contributing_contigs)), visual, stats

    def analyze(self):
        """
        Main analysis workflow.
        
        Returns:
            dict: Results dataframe and visualization data
        """
        print("\n" + "="*70)
        print("NEW GENE ANALYZER")
        print("="*70)
        
        # Step 1: Get protein hashes
        print("\n  Step 1: Loading protein sequences...")
        contig_hashes = self.get_gene_hashes(self.contigs_gff.replace('.gff3', '.faa').replace('.gff', '.faa'))
        print(f"    Contig proteins: {len(contig_hashes)} unique hashes")
        
        # Step 2: Parse scaffold genes
        print("\n  Step 2: Parsing scaffold genes...")
        scaffold_genes = self.parse_scaffold_genes(self.scaffolds_gff, 
                                                   self.scaffolds_gff.replace('.gff3', '.faa').replace('.gff', '.faa'))
        print(f"    Scaffold genes: {len(scaffold_genes)} total CDS")
        
        # Step 3: Align contigs to scaffolds
        scaffold_to_contig_mappings = self.map_contigs_to_scaffolds()
        
        # Step 4: Load sequences
        print("\n  Step 4: Loading sequences...")
        scaffold_seqs = self.get_seq_data(self.scaffolds_fasta)
        contig_seqs = self.get_seq_data(self.contigs_fasta)
        print(f"    Scaffold sequences: {len(scaffold_seqs)}")
        print(f"    Contig sequences: {len(contig_seqs)}")
        
        # Step 4b: Parse contig genes (to map back overlaps)
        print("\n  Step 4b: Parsing contig genes...")
        contig_genes = self.parse_scaffold_genes(self.contigs_gff,
                                                 self.contigs_gff.replace('.gff3', '.faa').replace('.gff', '.faa'))
        print(f"    Contig genes: {len(contig_genes)} total CDS")
        
        # Step 5: Find and classify new genes
        print("\n  Step 5: Finding and classifying new genes...")
        results = []
        vis_lines = []
        
        new_gene_count = 0
        for gene in scaffold_genes:
            if gene['hash'] and gene['hash'] not in contig_hashes:
                new_gene_count += 1
                
                # Classify the gene (now with contig gene information)
                classification, joined, visualization, stats = self.classify_new_gene(
                    gene, scaffold_to_contig_mappings,
                    scaffold_seqs.get(gene['scaffold'], ""), contig_seqs,
                    contig_genes=contig_genes
                )
                
                # Map to reference genome coordinates (global coordinate system)
                ref_start, ref_end, ref_seqid = self._map_gene_to_reference(
                    gene['start'], gene['end'], gene['scaffold']
                )

                results.append({
                    "GeneID": gene['id'],
                    "Product": gene['product'],
                    "Scaffold": gene['scaffold'],
                    "Scaffold_Start": gene['start'],
                    "Scaffold_End": gene['end'],
                    "Ref_SeqID": ref_seqid if ref_seqid else '',
                    "Ref_Start": ref_start if ref_start else '',
                    "Ref_End": ref_end if ref_end else '',
                    "Classification": classification,
                    "Joined": joined,
                    "Mapping_Stats": stats
                })
                
                # Format reference coordinates for display
                ref_location = f"{ref_seqid}:{ref_start:,}-{ref_end:,}" if ref_seqid else "unmapped"

                vis_lines.append(
                    f"Gene: {gene['id']} ({gene['product']})\n"
                    f"Scaffold locus: {gene['scaffold']}:{gene['start']}-{gene['end']}\n"
                    f"Reference position: {ref_location}\n"
                    f"Class: {classification}\n"
                    f"Joined: {joined}\n"
                    f"Stats: {stats}\n"
                    f"{'-'*100}\n{visualization}\n{'='*100}\n"
                )
        
        print(f"    Found {new_gene_count} new genes")
        
        # Step 6: Save results
        print("\n  Step 6: Saving results...")
        
        if results:
            # Save CSV files
            df = pd.DataFrame(results)
            
            # Main report
            report_file = self.new_genes_dir / 'gene_comparison_report.csv'
            df.to_csv(report_file, index=False)
            print(f"    Report: {report_file}")
            
            # Combined contigs
            df_combined = df[df['Classification'].str.contains("Combined Contigs", case=False, na=False)]
            if not df_combined.empty:
                combined_file = self.new_genes_dir / 'combined_contigs.csv'
                df_combined.to_csv(combined_file, index=False)
                print(f"    Combined contigs: {combined_file}")
            
            # Contextual re-predictions
            df_reprediction = df[df['Classification'].str.contains("Contextual Re-prediction", case=False, na=False)]
            if not df_reprediction.empty:
                repred_file = self.new_genes_dir / 'contextual_repredictions.csv'
                df_reprediction.to_csv(repred_file, index=False)
                print(f"    Contextual re-predictions: {repred_file}")
            
            # Save visualization file
            if vis_lines:
                legend = """
NEW GENE VISUALIZATIONS
=======================================================================================

This file shows detailed alignments of new genes found in scaffolds, displaying both
amino acid and nucleotide sequences in a stacked format for easy comparison.

LEGEND:
-------
  UPPERCASE letters  = Aligned region (primary alignment between scaffold and contig)
  lowercase letters  = Unaligned contig sequence (extends beyond aligned region)
  <                  = Start of contig alignment within the visible window
  >                  = End of contig alignment within the visible window
  |                  = Gene boundary markers (start and end positions)
  (N...)             = N bp of contig sequence extend to the LEFT (before visible window)
  (...N)             = N bp of contig sequence continue to the RIGHT (after visible window)

HEADER INFORMATION:
-------------------
Each gene includes detailed contig alignment information in the header:
  - Length: Total contig sequence length
  - Aligned: Positions aligned [contig coordinates] to [scaffold coordinates]
  - Strand: Alignment orientation (+/-)
  - Identity: Sequence identity percentage in aligned region

SCAFFOLD CONTEXT:
-----------------
The scaffold line shows:
  - (N...) prefix: N bp of scaffold extend before the visible window
  - |...|: Gene start and end positions marked with vertical bars
  - (...N) suffix: N bp of scaffold continue after the visible window

CONTIG DISPLAY:
---------------
Each contig line shows:
  - (N...) prefix: N bp of contig extend before visible window (alignment + overhang)
  - UPPERCASE: Aligned portion of the contig
  - lowercase: Unaligned contig sequence (overhang) within the visible window
  - (...N) suffix: N bp of contig extend after visible window (alignment + overhang)

FORMAT:
-------
Each gene shows TWO stacked visualizations:
  1. Amino Acid Alignment - Protein sequence view
  2. Nucleotide Alignment - DNA sequence view

All sequences for a gene are aligned vertically for easy visual comparison.

=======================================================================================

"""
                # Save TXT file
                vis_file_txt = self.visualizations_dir / 'new_genes_visualizations.txt'
                with open(vis_file_txt, 'w') as f:
                    f.write(legend)
                    f.write("\n".join(vis_lines))
                print(f"    Visualizations (TXT): {vis_file_txt}")

                # Save HTML file with proper styling for horizontal scroll
                vis_file_html = self.visualizations_dir / 'new_genes_visualizations.html'
                import html as html_module
                escaped_legend = html_module.escape(legend)
                escaped_content = html_module.escape("\n".join(vis_lines))
                html_content = f'''<!DOCTYPE html>
<html>
<head>
    <title>New Gene Visualizations - GenomeViz</title>
    <style>
        :root {{
            --bg-color: #1e1e1e;
            --text-color: #d4d4d4;
            --header-bg: #2d2d2d;
            --border-color: #404040;
            --link-color: #569cd6;
        }}
        body {{
            font-family: 'Consolas', 'Monaco', 'Courier New', monospace;
            background: var(--bg-color);
            color: var(--text-color);
            margin: 0;
            padding: 20px;
        }}
        .container {{
            max-width: 100%;
        }}
        h1 {{
            color: #569cd6;
            border-bottom: 2px solid var(--border-color);
            padding-bottom: 10px;
        }}
        .back-link {{
            display: inline-block;
            margin-bottom: 15px;
            color: var(--link-color);
            text-decoration: none;
        }}
        .back-link:hover {{
            text-decoration: underline;
        }}
        .legend {{
            background: var(--header-bg);
            padding: 15px;
            border-radius: 8px;
            margin-bottom: 20px;
            border: 1px solid var(--border-color);
        }}
        .legend pre {{
            margin: 0;
            white-space: pre;
            font-size: 13px;
        }}
        .alignments {{
            background: #252526;
            border: 1px solid var(--border-color);
            border-radius: 8px;
            overflow-x: auto;
            padding: 15px;
        }}
        .alignments pre {{
            margin: 0;
            white-space: pre;
            font-size: 13px;
            line-height: 1.4;
        }}
        footer {{
            margin-top: 30px;
            padding-top: 15px;
            border-top: 1px solid var(--border-color);
            color: #808080;
            font-size: 12px;
        }}
        footer a {{
            color: var(--link-color);
        }}
    </style>
</head>
<body>
    <div class="container">
        <a href="../../index.html" class="back-link">← Back to Comparison Overview</a>
        <h1>New Gene Visualizations</h1>

        <div class="legend">
            <pre>{escaped_legend}</pre>
        </div>

        <div class="alignments">
            <pre>{escaped_content}</pre>
        </div>

        <footer>
            Generated by <a href="https://github.com/Aaron-Thiel/GenomeViz">GenomeViz</a>
        </footer>
    </div>
</body>
</html>'''
                with open(vis_file_html, 'w', encoding='utf-8') as f:
                    f.write(html_content)
                print(f"    Visualizations (HTML): {vis_file_html}")
            
            return {'results': df, 'new_genes_count': new_gene_count, 'vis_dir': self.visualizations_dir}
        else:
            print("    No new genes found")
            return {'results': pd.DataFrame(), 'new_genes_count': 0, 'vis_dir': self.visualizations_dir}
