#!/usr/bin/env python3
"""
GenomeViz - Genome Assembly Visualization Tool

A comprehensive tool for visualizing and comparing bacterial genome assemblies against
reference sequences using interactive circular and linear plots.

Features:
- Circular and linear visualizations with multiple information tracks
- Interactive Plotly plots with hover information and zoom
- Static matplotlib plots for publication
- Automatic contig orientation detection
- Gene-level quality assessment and clickable alignments
- Detection of gaps, inversions, and duplications

GitHub: https://github.com/Aaron-Thiel/GenomeViz
License: MIT
"""

import argparse
import os
import sys
from pathlib import Path

from src.alignment import GenomeAligner, ContigMapper, OrientationDetector
from src.sequence import GFFParser
from src.utils import (validate_files, create_summary_report, save_gene_statistics,
                       print_final_summary, print_directory_tree, create_output_directories,
                       analyze_sequences, save_contig_mapping, create_visualizations,
                       generate_gene_alignments)
from src import __version__


def main():
    """Main program entry point."""

    parser = argparse.ArgumentParser(
        description='GenomeViz - Circular Genome Assembly Visualization',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  %(prog)s --reference ref.fasta --assembly asm.fasta --gff genes.gff --output results/

  # With custom thresholds
  %(prog)s --reference ref.fasta --assembly asm.fasta --gff genes.gff --output results/ \\
           --min-gap 2000 --min-inversion 1000

  # Skip auto-orientation
  %(prog)s --reference ref.fasta --assembly asm.fasta --gff genes.gff --output results/ \\
           --no-auto-orient

For more information: https://github.com/Aaron-Thiel/GenomeViz
        """
    )

    # Required arguments
    parser.add_argument('--reference', required=True,
                       help='Reference genome (FASTA format)')
    parser.add_argument('--assembly', required=True,
                       help='Assembly to compare (FASTA format)')
    parser.add_argument('--gff', required=True,
                       help='Gene annotations (GFF3 format)')
    parser.add_argument('--output', required=True,
                       help='Output directory')

    # Optional arguments
    parser.add_argument('--preset', default='asm10', choices=['asm5', 'asm10', 'asm20'],
                       help='Minimap2 preset for alignment (default: asm10)')
    parser.add_argument('--min-gap', type=int, default=1000,
                       help='Minimum gap size to report (default: 1000 bp)')
    parser.add_argument('--min-inversion', type=int, default=500,
                       help='Minimum inversion size to report (default: 500 bp)')
    parser.add_argument('--no-auto-orient', action='store_true',
                       help='Skip automatic orientation detection')
    parser.add_argument('--no-circular', action='store_true',
                       help='Skip circular plot generation')
    parser.add_argument('--no-linear', action='store_true',
                       help='Skip linear plot generation')
    parser.add_argument('--no-interactive', action='store_true',
                       help='Skip interactive Plotly plot generation')
    parser.add_argument('--no-interactive-linear', action='store_true',
                       help='Skip interactive linear plot generation')
    parser.add_argument('--no-gene-alignments', action='store_true',
                       help='Skip individual gene alignment file generation (disables gene clicking feature)')
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')

    args = parser.parse_args()

    # Validate input files
    try:
        validate_files(args.reference, args.assembly, args.gff)
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return 1

    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print(f"GenomeViz v{__version__}")
    print("=" * 70)
    print(f"\nOutput directory: {output_dir.absolute()}")

    # Step 1: Automatic orientation detection
    print("\n" + "="*70)
    print("[1/10] AUTOMATIC ORIENTATION DETECTION")
    print("="*70)

    if not args.no_auto_orient:
        detector = OrientationDetector(args.reference, args.assembly, preset=args.preset)
        assembly_to_use, temp_file_to_cleanup = detector.detect_and_correct()
    else:
        assembly_to_use = args.assembly
        temp_file_to_cleanup = None

    # Step 2: Parse gene annotations
    print("\n" + "="*70)
    print("[2/10] PARSING GENE ANNOTATIONS")
    print("="*70)
    gff_parser = GFFParser(args.gff)

    # Step 3: Load reference
    print("\n" + "="*70)
    print("[3/10] LOADING REFERENCE")
    print("="*70)
    aligner = GenomeAligner(args.reference, assembly_to_use, preset=args.preset)
    aligner.load_reference()

    # Step 4: Align assembly
    print("\n" + "="*70)
    print("[4/10] ALIGNING ASSEMBLY")
    print("="*70)
    alignments_by_ref = aligner.align()

    # Step 5: Create contig mapping
    print("\n" + "="*70)
    print("[5/10] CONTIG MAPPING")
    print("="*70)
    mapper = ContigMapper(alignments_by_ref)
    contig_mapping = mapper.create_mapping()

    # Step 6: Analyze sequences
    print("\n" + "="*70)
    print("[6/10] ANALYZING SEQUENCES")
    print("="*70)

    ref_sequences = analyze_sequences(
        aligner, gff_parser, alignments_by_ref, contig_mapping,
        args.min_gap, args.min_inversion
    )

    # Create output directories and save contig mapping
    seqid_dirs = create_output_directories(output_dir, ref_sequences)
    save_contig_mapping(output_dir, contig_mapping)

    # Step 7: Create visualizations
    print("\n" + "="*70)
    print("[7/10] CREATING VISUALIZATIONS")
    print("="*70)

    interactive_circular_files, interactive_linear_files = create_visualizations(
        ref_sequences, seqid_dirs, aligner, args
    )

    # Step 8: Generate gene alignments
    print("\n" + "="*70)
    print("[8/10] GENERATING GENE ALIGNMENTS")
    print("="*70)
    
    if not args.no_gene_alignments and interactive_linear_files:
        generate_gene_alignments(
            interactive_linear_files, interactive_circular_files,
            seqid_dirs, args.no_gene_alignments
        )

    # Step 9: Save statistics
    print("\n" + "="*70)
    print("[9/10] SAVING STATISTICS")
    print("="*70)
    
    if ref_sequences:
        save_gene_statistics(output_dir, ref_sequences)

    # Step 10: Summary
    print("\n" + "="*70)
    print("[10/10] SUMMARY")
    print("="*70)

    summary_file = create_summary_report(output_dir, ref_sequences, args.assembly, args.reference)
    print(f"Summary report saved to: {summary_file}")

    print_final_summary(output_dir, ref_sequences)

    print("\n" + "=" * 70)
    print("OUTPUT DIRECTORY STRUCTURE")
    print("=" * 70)
    print()
    print_directory_tree(output_dir, max_depth=2, skip_contents_of={'gene_alignments'})
    print(f"\n✓ All outputs saved to: {output_dir.absolute()}")
    print("\n✓ Done!")

    # Cleanup
    if temp_file_to_cleanup and os.path.exists(temp_file_to_cleanup):
        os.unlink(temp_file_to_cleanup)

    return 0


if __name__ == '__main__':
    sys.exit(main())
