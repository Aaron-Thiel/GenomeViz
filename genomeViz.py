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
- Origin of replication alignment
- Gene-level quality assessment and clickable alignments
- Detection of gaps, inversions, and duplications
- Scaffold-contig comparison mode for assembly QC

GitHub: https://github.com/Aaron-Thiel/GenomeViz
License: MIT
"""

import argparse
import os
import sys
import warnings
from pathlib import Path

# Suppress specific warnings
warnings.filterwarnings('ignore', category=UserWarning, module='Bio.Seq')
os.environ['QT_QPA_PLATFORM'] = 'offscreen'  # Suppress Qt warnings

from src.alignment import GenomeAligner, ContigMapper, OrientationDetector
from src.sequence import GFFParser, find_oric_position
from src.utils import (validate_files, create_summary_report, save_gene_statistics,
                       print_final_summary, print_directory_tree, create_output_directories,
                       analyze_sequences, save_contig_mapping, create_visualizations,
                       generate_gene_alignments, rotate_sequence_and_features)
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

  # Manually set origin position
  %(prog)s --reference ref.fasta --assembly asm.fasta --gff genes.gff --output results/ \\
           --origin 150000

  # Disable origin rotation
  %(prog)s --reference ref.fasta --assembly asm.fasta --gff genes.gff --output results/ \\
           --origin 0

  # Scaffold-contig comparison mode
  %(prog)s --reference ref.fasta --scaffolds scaffolds.fasta --contigs contigs.fasta \\
           --gff ref.gff --output results/

  # With GFF annotations for all assemblies
  %(prog)s --reference ref.fasta --scaffolds scaffolds.fasta --contigs contigs.fasta \\
           --gff ref.gff --scaffold-gff scaffolds.gff --contig-gff contigs.gff --output results/

For more information: https://github.com/Aaron-Thiel/GenomeViz
        """
    )

    # Required arguments
    parser.add_argument('--reference', required=True,
                       help='Reference genome (FASTA format)')
    parser.add_argument('--assembly', default=None,
                       help='Assembly to compare (FASTA format)')
    parser.add_argument('--gff', required=True,
                       help='Gene annotations for reference (GFF3 format)')
    parser.add_argument('--output', required=True,
                       help='Output directory')

    # Scaffold-contig comparison mode
    parser.add_argument('--scaffolds',
                       help='Scaffold assembly for comparison (FASTA format)')
    parser.add_argument('--contigs',
                       help='Contig assembly for comparison (FASTA format)')
    parser.add_argument('--scaffold-gff',
                       help='Gene annotations for scaffolds (GFF3 format)')
    parser.add_argument('--contig-gff',
                       help='Gene annotations for contigs (GFF3 format)')
    parser.add_argument('--no-comparison-circular', action='store_true',
                       help='Skip circular plots in assembly comparison')
    parser.add_argument('--no-comparison-linear', action='store_true',
                       help='Skip linear plots in assembly comparison')

    # Optional arguments
    parser.add_argument('--preset', default='asm10', choices=['asm5', 'asm10', 'asm20'],
                       help='Minimap2 preset for alignment (default: asm10)')
    parser.add_argument('--min-gap', type=int, default=1000,
                       help='Minimum gap size to report (default: 1000 bp)')
    parser.add_argument('--min-inversion', type=int, default=500,
                       help='Minimum inversion size to report (default: 500 bp)')
    parser.add_argument('--origin', type=int, default=None,
                       help='Manually set origin position (bp). Use 0 to disable origin rotation. '
                            'If not specified, will auto-detect oriC from GFF.')
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

    # Validate arguments - either --assembly or both --scaffolds and --contigs are required
    if not args.assembly and not (args.scaffolds and args.contigs):
        parser.error("Either --assembly OR both --scaffolds and --contigs are required")

    # Validate input files
    try:
        files_to_check = [args.reference, args.gff]
        if args.assembly:
            files_to_check.append(args.assembly)
        if args.scaffolds:
            files_to_check.append(args.scaffolds)
        if args.contigs:
            files_to_check.append(args.contigs)
        if args.scaffold_gff:
            files_to_check.append(args.scaffold_gff)
        if args.contig_gff:
            files_to_check.append(args.contig_gff)

        for filepath in files_to_check:
            if not os.path.exists(filepath):
                raise FileNotFoundError(f"File not found: {filepath}")
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return 1

    # Determine run mode
    run_standard_mode = args.assembly is not None
    run_comparison_mode = args.scaffolds is not None and args.contigs is not None

    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print(f"GenomeViz v{__version__}")
    print("=" * 70)
    print(f"\nOutput directory: {output_dir.absolute()}")

    mode_info = []
    if run_standard_mode:
        mode_info.append("Standard (assembly vs reference)")
    if run_comparison_mode:
        mode_info.append("Comparison (scaffolds vs contigs)")
    print(f"Mode: {' + '.join(mode_info)}")

    # Track temporary files for cleanup
    temp_files_to_cleanup = []

    # Step 1: Parse gene annotations (do this first to find oriC)
    print("\n" + "="*70)
    print("[1] PARSING GENE ANNOTATIONS")
    print("="*70)
    gff_parser = GFFParser(args.gff)

    # Step 2: Sequence preparation (origin rotation for reference)
    print("\n" + "="*70)
    print("[2] SEQUENCE PREPARATION")
    print("="*70)

    print("\n  [2a] Origin Detection & Reference Rotation")
    print("  " + "-"*66)

    reference_to_use = args.reference
    origin_position = args.origin

    if origin_position is None:
        # Auto-detect oriC from GFF
        origin_position, origin_seqid = find_oric_position(args.gff)

        if origin_position:
            print(f"  ✓ Found oriC at position {origin_position:,} bp on {origin_seqid}")
        else:
            print("  ℹ️  No oriC found in GFF file - using sequence as-is")
            origin_position = 0

    if origin_position and origin_position != 0:
        # Determine which seqid to rotate
        if 'origin_seqid' not in locals():
            # Need to load reference temporarily to find largest sequence
            temp_aligner = GenomeAligner(args.reference, args.reference, preset=args.preset)
            temp_aligner.load_reference()
            origin_seqid = max(temp_aligner.reference_sequences.keys(),
                              key=lambda k: temp_aligner.reference_sequences[k]['length'])
            print(f"  ℹ️  No oriC seqid detected, rotating largest sequence: {origin_seqid}")

        reference_to_use, gff_parser = rotate_sequence_and_features(
            args.reference, gff_parser, origin_position, origin_seqid
        )
        temp_files_to_cleanup.append(reference_to_use)
        print(f"  ✓ Reference rotated to start at oriC position {origin_position:,}")
    else:
        print("  ℹ️  Origin rotation disabled (--origin 0 or not found)")

    # =========================================================================
    # STANDARD MODE: Assembly vs Reference
    # =========================================================================
    if run_standard_mode:
        print("\n" + "="*70)
        print("STANDARD MODE: Assembly vs Reference")
        print("="*70)

        # Orientation detection for assembly contigs
        print("\n  [2b] Assembly Orientation Detection")
        print("  " + "-"*66)

        if not args.no_auto_orient:
            detector = OrientationDetector(reference_to_use, args.assembly, preset=args.preset)
            assembly_to_use, temp_asm = detector.detect_and_correct()
            if temp_asm:
                temp_files_to_cleanup.append(temp_asm)
        else:
            assembly_to_use = args.assembly
            print("  ℹ️  Orientation detection skipped (--no-auto-orient)")

        # Load reference and align
        print("\n" + "-"*70)
        print("Loading reference and aligning assembly...")
        print("-"*70)
        aligner = GenomeAligner(reference_to_use, assembly_to_use, preset=args.preset)
        aligner.load_reference()
        alignments_by_ref = aligner.align()

        # Create contig mapping
        mapper = ContigMapper(alignments_by_ref)
        contig_mapping = mapper.create_mapping()

        # Analyze sequences
        print("\nAnalyzing sequences...")
        ref_sequences = analyze_sequences(
            aligner, gff_parser, alignments_by_ref, contig_mapping,
            args.min_gap, args.min_inversion
        )

        # Create output directories and save contig mapping
        seqid_dirs = create_output_directories(output_dir, ref_sequences)
        save_contig_mapping(output_dir, contig_mapping)

        # Create visualizations
        print("\nCreating visualizations...")
        interactive_circular_files, interactive_linear_files = create_visualizations(
            ref_sequences, seqid_dirs, aligner, args
        )

        # Generate gene alignments
        if not args.no_gene_alignments and interactive_linear_files:
            print("\nGenerating gene alignments...")
            generate_gene_alignments(
                interactive_linear_files, interactive_circular_files,
                seqid_dirs, args.no_gene_alignments
            )

        # Save statistics
        if ref_sequences:
            save_gene_statistics(output_dir, ref_sequences)

        # Generate summary report
        summary_file = create_summary_report(output_dir, ref_sequences, args.assembly, args.reference)
        print(f"\nSummary report: {summary_file}")
        print_final_summary(output_dir, ref_sequences)

    # =========================================================================
    # COMPARISON MODE: Scaffolds vs Contigs
    # =========================================================================
    if run_comparison_mode:
        print("\n" + "="*70)
        print("COMPARISON MODE: Scaffolds vs Contigs")
        print("="*70)

        from src.assembly_comparison import AssemblyComparator, MultiAssemblyGFFParser
        from src.comparison_visualizer import (
            ScaffoldSequence,
            ComparisonLinearVisualizer,
            ComparisonCircularVisualizer,
            ComparisonIndexGenerator,
            ComparisonGeneAlignmentVisualizer
        )
        from src.gene_comparison import GeneIntegrityAnalyzer

        # Step C1: Orientation detection for scaffolds
        print("\n  [C1] Scaffold Orientation Detection")
        print("  " + "-"*66)
        scaffolds_to_use = args.scaffolds
        if not args.no_auto_orient:
            scaffold_detector = OrientationDetector(reference_to_use, args.scaffolds, preset=args.preset)
            scaffolds_to_use, temp_scaffolds = scaffold_detector.detect_and_correct()
            if temp_scaffolds:
                temp_files_to_cleanup.append(temp_scaffolds)
        else:
            print("  ℹ️  Orientation detection skipped (--no-auto-orient)")

        # Step C2: Orientation detection for contigs
        print("\n  [C2] Contig Orientation Detection")
        print("  " + "-"*66)
        contigs_to_use = args.contigs
        if not args.no_auto_orient:
            contig_detector = OrientationDetector(reference_to_use, args.contigs, preset=args.preset)
            contigs_to_use, temp_contigs = contig_detector.detect_and_correct()
            if temp_contigs:
                temp_files_to_cleanup.append(temp_contigs)
        else:
            print("  ℹ️  Orientation detection skipped (--no-auto-orient)")

        # Step C3: Align scaffolds to reference
        print("\n  [C3] Aligning scaffolds to reference...")
        scaffold_aligner = GenomeAligner(reference_to_use, scaffolds_to_use, preset=args.preset)
        scaffold_aligner.load_reference()
        scaffold_alignments = scaffold_aligner.align()

        # Step C4: Align contigs to reference
        print("\n  [C4] Aligning contigs to reference...")
        contig_aligner = GenomeAligner(reference_to_use, contigs_to_use, preset=args.preset)
        contig_aligner.load_reference()
        contig_alignments = contig_aligner.align()

        # Step C5: Parse additional GFF files
        print("\n  [C5] Parsing assembly GFF files...")
        multi_gff = MultiAssemblyGFFParser(
            reference_gff=args.gff,
            scaffold_gff=args.scaffold_gff,
            contig_gff=args.contig_gff
        )

        # Step C6: Find overlapping regions
        print("\n  [C6] Finding scaffold-contig overlaps...")
        comparator = AssemblyComparator(scaffold_alignments, contig_alignments)
        comparison_summary = comparator.compute_comparison_summary()

        print(f"  Found {comparison_summary['num_scaffolds']} scaffolds")
        print(f"  Found {comparison_summary['num_contigs']} contigs")
        print(f"  Found {comparison_summary['num_total_overlap_regions']} overlap regions")

        # Get reference length for context visualization
        # Use the largest reference sequence length
        reference_lengths = {seqid: info['length']
                            for seqid, info in scaffold_aligner.reference_sequences.items()}
        total_reference_length = max(reference_lengths.values()) if reference_lengths else 0
        print(f"  Reference genome length: {total_reference_length:,} bp")

        # Create comparison output directory
        comparison_dir = output_dir / 'assembly_comparison'
        comparison_dir.mkdir(exist_ok=True)

        # Step C7: Generate visualizations for each scaffold
        print("\n  [C7] Generating comparison visualizations...")

        all_gene_results = {}

        for scaffold_name in comparator.get_all_scaffolds():
            overlap_regions = comparator.find_overlapping_contigs(scaffold_name)

            if not overlap_regions:
                print(f"    Skipping {scaffold_name} (no overlapping contigs)")
                continue

            # Create scaffold output directory
            safe_name = "".join(c if c.isalnum() or c in '-_' else '_' for c in scaffold_name)
            scaffold_dir = comparison_dir / safe_name
            scaffold_dir.mkdir(exist_ok=True)

            print(f"    Processing {scaffold_name} ({len(overlap_regions)} regions)...")

            # Calculate scaffold reference span
            ref_start = min(r.scaffold_ref_start for r in overlap_regions)
            ref_end = max(r.scaffold_ref_end for r in overlap_regions)
            reference_seqid = overlap_regions[0].reference_seqid

            # Analyze gene integrity
            # Priority: scaffold GFF > reference GFF (required --gff)
            gene_results = None
            if multi_gff.has_genes('scaffold'):
                # Use scaffold genes if scaffold GFF is provided
                gene_analyzer = GeneIntegrityAnalyzer(overlap_regions, multi_gff)
                gene_results = gene_analyzer.analyze_scaffold_genes_via_reference(
                    scaffold_name, overlap_regions
                )
                if gene_results:
                    all_gene_results[scaffold_name] = gene_results
                    print(f"      Using {len(gene_results)} scaffold genes")
            
            # Fall back to reference genes if scaffold genes not available
            # This happens when: no scaffold GFF, or scaffold was RC'd (coords don't match)
            if not gene_results and multi_gff.has_genes('reference'):
                gene_analyzer = GeneIntegrityAnalyzer(overlap_regions, multi_gff)
                gene_results = gene_analyzer.analyze_reference_genes_in_range(
                    ref_start, ref_end, reference_seqid, overlap_regions
                )
                if gene_results:
                    all_gene_results[scaffold_name] = gene_results
                    print(f"      Using {len(gene_results)} reference genes in region")

            # Create ScaffoldSequence object for visualization
            # Pass reference_length for full genome context visualization
            scaffold_seq = ScaffoldSequence(
                scaffold_name=scaffold_name,
                ref_start=ref_start,
                ref_end=ref_end,
                overlap_regions=overlap_regions,
                gene_results=gene_results,
                reference_length=total_reference_length
            )

            # Create linear visualization
            if not args.no_comparison_linear:
                ComparisonLinearVisualizer.create_linear_plot(
                    scaffold_seq,
                    scaffold_dir / f'{safe_name}_linear.png',
                    reference_seqid=reference_seqid
                )

            # Create circular visualizations
            circular_html_file = None
            if not args.no_comparison_circular:
                circular_html_file = scaffold_dir / f'{safe_name}_circular.html'
                # Interactive HTML
                ComparisonCircularVisualizer.create_interactive_circular_plot(
                    scaffold_seq,
                    circular_html_file,
                    reference_seqid=reference_seqid
                )
                # Static PNG
                ComparisonCircularVisualizer.create_static_circular_plot(
                    scaffold_seq,
                    scaffold_dir / f'{safe_name}_circular.png',
                    reference_seqid=reference_seqid
                )

            # Generate gene alignment files if genes are available
            if gene_results and not args.no_gene_alignments:
                gene_align_dir = scaffold_dir / 'gene_alignments'
                gene_align_dir.mkdir(exist_ok=True)

                print(f"      Generating {len(gene_results)} gene alignments...")
                for gene_result in gene_results:
                    try:
                        ComparisonGeneAlignmentVisualizer.create_gene_alignment_html(
                            gene_result, scaffold_seq, gene_align_dir
                        )
                    except Exception as e:
                        print(f"      Warning: Could not create alignment for {gene_result.gene['name']}: {e}")

                # Add click handlers to circular plot
                if circular_html_file and circular_html_file.exists():
                    ComparisonGeneAlignmentVisualizer.add_circular_click_handler(
                        circular_html_file, scaffold_seq, gene_results
                    )

        # Step C8: Compute gene summary if available
        gene_summary = None
        if all_gene_results:
            gene_analyzer = GeneIntegrityAnalyzer([], multi_gff)
            gene_summary = gene_analyzer.compute_summary(all_gene_results)
            print(f"\n  Gene Integrity Summary:")
            print(f"    Complete: {gene_summary['status_counts']['complete']}")
            print(f"    Split: {gene_summary['status_counts']['split']}")
            print(f"    Truncated: {gene_summary['status_counts']['truncated']}")
            print(f"    Missing: {gene_summary['status_counts']['missing']}")

        # Step C9: Generate index page
        print("\n  [C9] Generating comparison index...")
        index_gen = ComparisonIndexGenerator(comparison_summary, comparison_dir, gene_summary)
        index_file = index_gen.generate_index_html()
        print(f"  Index page: {index_file}")

        # Save comparison summary JSON
        import json
        summary_file = comparison_dir / 'comparison_summary.json'
        # Remove non-serializable 'regions' from summary before saving
        serializable_summary = {
            'num_scaffolds': comparison_summary['num_scaffolds'],
            'num_contigs': comparison_summary['num_contigs'],
            'num_total_overlap_regions': comparison_summary['num_total_overlap_regions'],
            'scaffold_stats': {
                name: {k: v for k, v in stats.items() if k != 'regions'}
                for name, stats in comparison_summary['scaffold_stats'].items()
            }
        }
        with open(summary_file, 'w') as f:
            json.dump(serializable_summary, f, indent=2)

    # =========================================================================
    # FINAL OUTPUT
    # =========================================================================
    print("\n" + "=" * 70)
    print("OUTPUT DIRECTORY STRUCTURE")
    print("=" * 70)
    print()
    print_directory_tree(output_dir, max_depth=2, skip_contents_of={'gene_alignments'})

    print("\n" + "=" * 70)
    print(f"✓ All outputs saved to: {output_dir.absolute()}")
    print("✓ Analysis complete!")
    print("=" * 70)

    # Cleanup temporary files
    for temp_file in temp_files_to_cleanup:
        if temp_file and os.path.exists(temp_file):
            os.unlink(temp_file)

    return 0


if __name__ == '__main__':
    sys.exit(main())