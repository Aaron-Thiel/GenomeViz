"""
Genome alignment and analysis module

Provides classes for:
- Automatic orientation detection and correction
- Genome alignment using minimap2
- Misassembly detection (gaps, inversions, overlaps)
- Contig-to-reference mapping

Classes:
    OrientationDetector: Detect and correct contig orientation
    GenomeAligner: Align assembly to reference genome
    MisassemblyDetector: Detect potential misassemblies
    ContigMapper: Create contig-to-reference mappings

GitHub: https://github.com/Aaron-Thiel/GenomeViz
License: MIT
"""

import mappy as mp
import numpy as np
import tempfile
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq


class OrientationDetector:
    """
    Detect and correct orientation issues in assembly contigs.

    Analyzes alignment strand bias to identify contigs that should be
    reverse-complemented for optimal alignment to the reference.
    """

    def __init__(self, reference_fasta, assembly_fasta, preset='asm10'):
        """
        Initialize orientation detector.

        Args:
            reference_fasta (str): Path to reference genome
            assembly_fasta (str): Path to assembly to check
            preset (str): Minimap2 preset for alignment
        """
        self.reference_fasta = reference_fasta
        self.assembly_fasta = assembly_fasta
        self.preset = preset
        self.corrections = {}

    def detect_and_correct(self):
        """
        Detect inverted contigs and create orientation-corrected assembly.

        Returns:
            tuple: (corrected_assembly_path, temp_file_to_cleanup or None)
        """

        aligner = mp.Aligner(self.reference_fasta, preset=self.preset)
        if not aligner:
            raise Exception("Failed to load reference")

        contig_stats = {}

        # Analyze each contig's alignment strand bias
        for name, seq, qual in mp.fastx_read(self.assembly_fasta):
            forward_bases = 0
            reverse_bases = 0

            for hit in aligner.map(seq):
                if hit.is_primary:
                    aligned_length = hit.r_en - hit.r_st
                    if hit.strand == 1:
                        forward_bases += aligned_length
                    else:
                        reverse_bases += aligned_length

            total_bases = forward_bases + reverse_bases
            if total_bases > 0:
                reverse_pct = (reverse_bases / total_bases) * 100
                contig_stats[name] = {
                    'forward_bases': forward_bases,
                    'reverse_bases': reverse_bases,
                    'reverse_pct': reverse_pct,
                    'needs_rc': reverse_pct > 50,  # RC if >50% aligns in reverse
                    'seq': seq
                }

        needs_correction = [name for name, stats in contig_stats.items() if stats['needs_rc']]

        # Print summary
        print(f"\nOrientation analysis for {len(contig_stats)} contigs:")
        print(f"  Contigs in correct orientation: {len(contig_stats) - len(needs_correction)}")
        print(f"  Contigs needing reverse complement: {len(needs_correction)}")

        if needs_correction:
            print("\nContigs to be reverse complemented:")
            for name in needs_correction:
                stats = contig_stats[name]
                print(f"  - {name}: {stats['reverse_pct']:.1f}% aligned in reverse")

        # Create corrected assembly if needed
        if needs_correction:
            corrected_file = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
            corrected_path = corrected_file.name

            print(f"\nCreating orientation-corrected assembly...")

            for name, seq, qual in mp.fastx_read(self.assembly_fasta):
                if name in contig_stats:
                    if contig_stats[name]['needs_rc']:
                        # Reverse complement this contig
                        seq_obj = Seq(seq)
                        rc_seq = str(seq_obj.reverse_complement())
                        corrected_file.write(f">{name}_RC\n")
                        for i in range(0, len(rc_seq), 80):
                            corrected_file.write(rc_seq[i:i+80] + "\n")
                        self.corrections[name] = 'reverse_complemented'
                    else:
                        # Keep as is
                        corrected_file.write(f">{name}\n")
                        for i in range(0, len(seq), 80):
                            corrected_file.write(seq[i:i+80] + "\n")
                        self.corrections[name] = 'kept_forward'

            corrected_file.close()
            print(f"Corrected assembly saved to temporary file: {corrected_path}")
            print("\n✓ Using orientation-corrected assembly for analysis")
            return corrected_path, corrected_path
        else:
            print("\nNo orientation corrections needed - all contigs in correct orientation!")
            print("✓ Using original assembly (no corrections needed)")
            return self.assembly_fasta, None


class GenomeAligner:
    """
    Align assembly contigs to reference genome using minimap2.

    Provides methods to load reference sequences and perform alignment
    with customizable minimap2 presets.
    """

    def __init__(self, reference_fasta, assembly_fasta, preset='asm10'):
        """
        Initialize genome aligner.

        Args:
            reference_fasta (str): Path to reference genome
            assembly_fasta (str): Path to assembly
            preset (str): Minimap2 preset (asm5/asm10/asm20)
        """
        self.reference_fasta = reference_fasta
        self.assembly_fasta = assembly_fasta
        self.preset = preset
        self.reference_sequences = {}
        self.assembly_name = Path(assembly_fasta).stem

    def load_reference(self):
        """Load all sequences from reference genome."""
        print("\nLoading reference sequences...")
        for record in SeqIO.parse(self.reference_fasta, "fasta"):
            self.reference_sequences[record.id] = {
                'sequence': str(record.seq),
                'length': len(record.seq)
            }
            seq_type = 'chromosome' if len(record.seq) > 500000 else 'plasmid'
            print(f"  {record.id}: {len(record.seq):,} bp ({seq_type})")

    def align(self):
        """
        Perform alignment using minimap2.

        Returns:
            dict: Alignments grouped by reference sequence ID
        """
        print("\nAligning assembly to reference...")

        aligner = mp.Aligner(self.reference_fasta, preset=self.preset)
        if not aligner:
            raise Exception("Failed to load reference")

        alignments_by_ref = defaultdict(list)

        alignment_count = 0
        for name, seq, qual in mp.fastx_read(self.assembly_fasta):
            seq_len = len(seq)
            for hit in aligner.map(seq):
                alignment = {
                    'query_name': name,
                    'query_start': hit.q_st,
                    'query_end': hit.q_en,
                    'query_length': seq_len,
                    'ref_start': hit.r_st,
                    'ref_end': hit.r_en,
                    'ref_name': hit.ctg,
                    'strand': '+' if hit.strand == 1 else '-',
                    'mapq': hit.mapq,
                    'matches': hit.mlen,
                    'block_length': hit.blen,
                    'identity': hit.mlen / hit.blen * 100 if hit.blen > 0 else 0,
                    'is_primary': hit.is_primary,
                    'cigar': hit.cigar  # Store CIGAR for proper alignment display
                }
                alignments_by_ref[hit.ctg].append(alignment)
                alignment_count += 1

        print(f"Found {alignment_count} total alignments")

        # Sort and summarize
        for ref_name in alignments_by_ref:
            alignments_by_ref[ref_name].sort(key=lambda x: x['ref_start'])
            primary_count = len([a for a in alignments_by_ref[ref_name] if a['is_primary']])
            print(f"  {ref_name}: {len(alignments_by_ref[ref_name])} alignments ({primary_count} primary)")

        return alignments_by_ref


class MisassemblyDetector:
    """
    Detect potential misassemblies (inversions, gaps, overlaps).

    Analyzes alignment patterns to identify potential assembly errors
    including strand inversions, sequence gaps, and overlaps.
    """

    def __init__(self, alignments, min_gap=1000, min_inversion_size=500):
        """
        Initialize misassembly detector.

        Args:
            alignments (list): List of alignments
            min_gap (int): Minimum gap size to report (bp)
            min_inversion_size (int): Minimum inversion size to report (bp)
        """
        self.alignments = alignments
        self.min_gap = min_gap
        self.min_inversion_size = min_inversion_size
        self.misassemblies = []

    def detect(self):
        """
        Detect various types of misassemblies.

        Returns:
            list: List of detected misassemblies with type, position, and size
        """
        primary_aligns = [a for a in self.alignments if a['is_primary']]

        # Detect inversions (reverse strand alignments)
        for aln in primary_aligns:
            if aln['strand'] == '-' and aln['block_length'] > self.min_inversion_size:
                self.misassemblies.append({
                    'type': 'inversion',
                    'ref_start': aln['ref_start'],
                    'ref_end': aln['ref_end'],
                    'query_name': aln['query_name'],
                    'size': aln['ref_end'] - aln['ref_start']
                })

        # Detect gaps and overlaps between consecutive alignments
        for i in range(len(primary_aligns) - 1):
            curr = primary_aligns[i]
            next_aln = primary_aligns[i + 1]

            gap = next_aln['ref_start'] - curr['ref_end']

            if gap > self.min_gap:
                # Large gap - missing sequence
                self.misassemblies.append({
                    'type': 'gap',
                    'ref_start': curr['ref_end'],
                    'ref_end': next_aln['ref_start'],
                    'query_name': f"{curr['query_name']}->{next_aln['query_name']}",
                    'size': gap
                })
            elif gap < -self.min_gap:
                # Large overlap
                self.misassemblies.append({
                    'type': 'overlap',
                    'ref_start': next_aln['ref_start'],
                    'ref_end': curr['ref_end'],
                    'query_name': f"{curr['query_name']}<->{next_aln['query_name']}",
                    'size': abs(gap)
                })

        return self.misassemblies


class ContigMapper:
    """
    Create contig-to-reference mapping for visualization.

    Organizes alignment data into a hierarchical structure
    for easy access during visualization generation.
    """

    def __init__(self, alignments_by_ref):
        """
        Initialize contig mapper.

        Args:
            alignments_by_ref (dict): Alignments grouped by reference
        """
        self.alignments_by_ref = alignments_by_ref
        self.contig_mapping = defaultdict(lambda: defaultdict(list))

    def create_mapping(self):
        """
        Create detailed contig-to-reference mapping.

        Returns:
            dict: Nested dict mapping reference -> contig -> alignment details
        """
        print("\nCreating contig-to-reference mapping...")

        for ref_name, alignments in self.alignments_by_ref.items():
            primary_aligns = [a for a in alignments if a['is_primary']]

            for aln in primary_aligns:
                contig_name = aln['query_name']

                self.contig_mapping[ref_name][contig_name].append({
                    'ref_start': aln['ref_start'],
                    'ref_end': aln['ref_end'],
                    'query_start': aln['query_start'],
                    'query_end': aln['query_end'],
                    'query_length': aln['query_length'],
                    'strand': aln['strand'],
                    'identity': aln['identity'],
                    'coverage': (aln['ref_end'] - aln['ref_start'])
                })

        # Print summary
        for ref_name in sorted(self.contig_mapping.keys()):
            contigs = self.contig_mapping[ref_name]
            print(f"\n  {ref_name}: {len(contigs)} contigs aligned")
            for contig_name in sorted(contigs.keys()):
                alns = contigs[contig_name]
                total_coverage = sum(a['coverage'] for a in alns)
                avg_identity = np.mean([a['identity'] for a in alns])
                print(f"    - {contig_name}: {len(alns)} alignment(s), "
                      f"{total_coverage:,} bp covered, {avg_identity:.1f}% identity")

        return dict(self.contig_mapping)
