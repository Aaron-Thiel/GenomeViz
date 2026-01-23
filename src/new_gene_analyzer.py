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
        print(f"    Scaffolds file: {self.scaffolds_fasta}")
        print(f"    Contigs file: {self.contigs_fasta}")
        print(f"    DEBUG: Checking files exist...")
        import os
        print(f"      Scaffolds exists: {os.path.exists(self.scaffolds_fasta)}")
        print(f"      Contigs exists: {os.path.exists(self.contigs_fasta)}")
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
                    # Note: hit.q_len is NOT available in mappy, so we use len(seq) from the loop
                    s_len = len(seq)
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
                        # Tuple: (ref_start, ref_end, contig_id, query_start, query_end, s_len, strand, identity, alignment_fraction)
                        # Now: ref is CONTIG (reference), query is SCAFFOLD (query)
                        mappings[name].append((
                            query_block_start, query_block_end, hit.ctg, ref_block_start, ref_block_end,
                            s_len, '+' if hit.strand == 1 else '-', identity, alignment_fraction
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
        Generates an alignment visualization showing aligned parts (UPPER) and soft-clips (lower).
        
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
        win_aa = self.get_translated_context(full_scaffold_seq, win_start, win_end, strand)
        win_len_aa = len(win_aa)
        
        if win_len_aa == 0:
            return "Could not generate visualization", ""
        
        vis_lines = []
        mapping_stats = []
        
        for s_start, s_end, c_id, c_start, c_end, c_len, c_strand, identity, align_frac in sorted(overlapping_contigs, key=lambda x: x[0]):
            line_chars = [" "] * win_len_aa
            mapping_stats.append(f"{c_id}: align: {align_frac:.2f}, id: {identity:.2f}")
            
            # Find overlap between window and alignment
            intersect_start = max(win_start, s_start)
            intersect_end = min(win_end, s_end)
            
            if intersect_start < intersect_end:
                rel_start_aa = (intersect_start - win_start) // 3
                rel_end_aa = (intersect_end - win_start) // 3
                
                # Mark aligned region with uppercase
                for i in range(rel_start_aa, min(rel_end_aa, win_len_aa)):
                    line_chars[i] = win_aa[i]
                
                # Mark boundaries
                if s_start >= win_start and rel_start_aa < win_len_aa:
                    line_chars[rel_start_aa] = "<"
                if s_end <= win_end and (rel_end_aa - 1) < win_len_aa:
                    line_chars[max(0, rel_end_aa - 1)] = ">"
            
            # Add soft-clip regions (unaligned contig sequence)
            full_c_seq = contig_seqs.get(c_id, "")
            if full_c_seq:
                # Left soft-clip
                if c_start > 0 and s_start >= win_start:
                    clip_seq = Seq(full_c_seq[:c_start])
                    if c_strand == '-':
                        clip_seq = clip_seq.reverse_complement()
                    clip_aa = str(clip_seq[:(len(clip_seq)//3)*3].translate()).lower()
                    for i, aa in enumerate(reversed(clip_aa)):
                        idx = (s_start - win_start) // 3 - 1 - i
                        if 0 <= idx < win_len_aa:
                            line_chars[idx] = aa
                
                # Right soft-clip
                if c_end < c_len and s_end <= win_end:
                    clip_seq = Seq(full_c_seq[c_end:])
                    if c_strand == '-':
                        clip_seq = clip_seq.reverse_complement()
                    clip_aa = str(clip_seq[:(len(clip_seq)//3)*3].translate()).lower()
                    for i, aa in enumerate(clip_aa):
                        idx = (s_end - win_start) // 3 + i
                        if 0 <= idx < win_len_aa:
                            line_chars[idx] = aa
            
            vis_lines.append(f"{c_id:<20} {' ':<12} {''.join(line_chars)}")
        
        # Add scaffold sequence line
        s_chars = list(win_aa)
        g_rel_s = (g_start_0 - win_start) // 3
        g_rel_e = (g_end_0 - win_start) // 3
        
        if 0 <= g_rel_s < win_len_aa:
            s_chars[g_rel_s] = "|"
        if 0 <= (g_rel_e - 1) < win_len_aa:
            s_chars[g_rel_e - 1] = "|"
        
        s_prefix = f"({win_start}...)" if win_start > 0 else " " * 12
        s_suffix = f"(...{scaffold_len - win_end})" if win_end < scaffold_len else ""
        vis_lines.append(f"{'Scaffold (Context)':<20} {s_prefix:>12} {''.join(s_chars)} {s_suffix}")
        
        return "\n".join(vis_lines), "; ".join(mapping_stats)

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
            # Debug: Show details
            if gene['id'].startswith('MJAJBK_00'):
                print(f"    DEBUG: Processing {len(ovl)} overlapping blocks")
            
            # For each overlapping alignment block, map coordinates back to contig space
            # and find what genes are there
            for aln in ovl:
                ref_block_start, ref_block_end, contig_id, query_block_start, query_block_end, c_len, strand, identity, align_frac = aln
                
                if gene['id'].startswith('MJAJBK_00'):
                    print(f"      Block: {contig_id}[{ref_block_start}:{ref_block_end}] -> query[{query_block_start}:{query_block_end}]")
                
                # Map the overlapping region [g_s:g_e] from scaffold to contig coordinates
                # Find the intersection with this alignment block
                overlap_start = max(g_s, ref_block_start)
                overlap_end = min(g_e, ref_block_end)
                
                if overlap_start < overlap_end:
                    # Calculate offset: how much into the alignment block is the overlap start?
                    offset_in_block = overlap_start - ref_block_start
                    
                    # Map to contig coordinates
                    contig_overlap_start = query_block_start + offset_in_block
                    contig_overlap_end = contig_overlap_start + (overlap_end - overlap_start)
                    
                    if gene['id'].startswith('MJAJBK_00'):
                        print(f"        Mapped to contig: [{contig_overlap_start}:{contig_overlap_end}]")
                    
                    # Find genes in the contig that overlap this region
                    if contig_id not in contig_gene_sources:
                        contig_gene_sources[contig_id] = []
                    
                    # Check which contig genes overlap with this region
                    genes_found_in_contig = 0
                    if contig_genes:
                        for contig_gene in contig_genes:
                            if contig_gene['scaffold'] == contig_id:  # Same contig
                                c_g_s, c_g_e = contig_gene['start'] - 1, contig_gene['end']
                                # Check if contig gene overlaps the mapped region
                                if max(c_g_s, contig_overlap_start) < min(c_g_e, contig_overlap_end):
                                    genes_found_in_contig += 1
                                    contig_gene_sources[contig_id].append({
                                        'id': contig_gene['id'],
                                        'product': contig_gene['product'],
                                        'start': contig_gene['start'],
                                        'end': contig_gene['end']
                                    })
                    
                    if gene['id'].startswith('MJAJBK_00'):
                        print(f"        Found {genes_found_in_contig} genes in {contig_id} at [{contig_overlap_start}:{contig_overlap_end}]")
        
        # Debug output
        if gene['id'].startswith('MJAJBK_00'):
            print(f"    DEBUG: {gene['id']} [{g_s}:{g_e}] on {scaffold}")
            print(f"      Contig_genes available: {contig_genes is not None and len(contig_genes) > 0}")
            print(f"      Alignments on this scaffold: {len(scaffold_alignments)}")
            print(f"      Overlapping blocks: {len(ovl)}")
            if len(ovl) == 0 and len(scaffold_alignments) > 0:
                print(f"      WARNING: No overlaps found! Showing scaffold alignments:")
                for i, aln in enumerate(scaffold_alignments[:3]):
                    ref_s, ref_e = aln[0], aln[1]
                    contig_id = aln[2]
                    print(f"        Alignment {i}: scaffold[{ref_s}:{ref_e}] <- {contig_id}")
            if scaffold_alignments:
                for i, aln in enumerate(scaffold_alignments[:2]):
                    ref_s, ref_e = aln[0], aln[1]
                    contig_id = aln[2]
                    print(f"        Block {i}: scaffold[{ref_s}:{ref_e}] <- {contig_id}")
                    # Check if this block overlaps the gene
                    if max(g_s, ref_s) < min(g_e, ref_e):
                        print(f"          -> OVERLAPS gene at [{max(g_s, ref_s)}:{min(g_e, ref_e)}]")
            for contig_id, genes in contig_gene_sources.items():
                print(f"        Contig {contig_id}: {len(genes)} genes")
                for g in genes[:2]:
                    print(f"          - {g['id']}: {g['product']}")
        
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
                vis_file = self.visualizations_dir / 'new_genes_visualizations.txt'
                with open(vis_file, 'w') as f:
                    f.write(
                        "Legend: UPPER = Aligned region, lower = Unaligned contig overhang, "
                        "| = Gene position, < > = Contig alignment boundaries\n\n"
                    )
                    f.write("\n".join(vis_lines))
                print(f"    Visualizations: {vis_file}")
            
            return {'results': df, 'new_genes_count': new_gene_count, 'vis_dir': self.visualizations_dir}
        else:
            print("    No new genes found")
            return {'results': pd.DataFrame(), 'new_genes_count': 0, 'vis_dir': self.visualizations_dir}
