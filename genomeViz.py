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
import json
import hashlib
from pathlib import Path
from datetime import datetime

# Suppress specific warnings
warnings.filterwarnings('ignore', category=UserWarning, module='Bio.Seq')
os.environ['QT_QPA_PLATFORM'] = 'offscreen'  # Suppress Qt warnings

from src.alignment import GenomeAligner, ContigMapper, OrientationDetector
from src.sequence import GFFParser, find_oric_position
from src.utils import (validate_files, create_summary_report, save_gene_statistics,
                       print_final_summary, print_directory_tree, create_output_directories,
                       analyze_sequences, save_contig_mapping, create_visualizations,
                       generate_gene_alignments, rotate_sequence_and_features)
from src.new_gene_analyzer import NewGeneAnalyzer
from src import __version__


def get_file_hash(filepath):
    """Get MD5 hash of file contents for cache validation."""
    if not os.path.exists(filepath):
        return None
    md5_hash = hashlib.md5()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            md5_hash.update(chunk)
    return md5_hash.hexdigest()


def load_cache(output_dir):
    """Load cache metadata from output directory."""
    cache_file = Path(output_dir) / '.genomeviz_cache.json'
    if cache_file.exists():
        try:
            with open(cache_file, 'r') as f:
                return json.load(f)
        except:
            return {}
    return {}


def save_cache(output_dir, cache_data):
    """Save cache metadata to output directory."""
    cache_file = Path(output_dir) / '.genomeviz_cache.json'
    with open(cache_file, 'w') as f:
        json.dump(cache_data, f, indent=2)


def check_cache_valid(cache_data, mode_key, input_files):
    """Check if cached results are still valid based on input file hashes."""
    if mode_key not in cache_data:
        return False
    
    mode_cache = cache_data[mode_key]
    if 'input_hashes' not in mode_cache:
        return False
    
    current_hashes = {fname: get_file_hash(fpath) for fname, fpath in input_files.items()}
    cached_hashes = mode_cache.get('input_hashes', {})
    
    # Check if all input files match
    return current_hashes == cached_hashes


def auto_detect_files(input_dir):
    """
    Automatically detect standard file names in the input directory.

    Args:
        input_dir: Path to input directory

    Returns:
        dict: Detected file paths with keys: reference, gff, scaffolds,
              contigs, scaffold_gff, contig_gff, scaffold_faa, contig_faa
    """
    from pathlib import Path

    input_path = Path(input_dir)
    detected = {}

    # Reference files
    for ext in ['fna', 'fasta', 'fa']:
        ref_file = input_path / f'reference.{ext}'
        if ref_file.exists():
            detected['reference'] = str(ref_file)
            break

    for ext in ['gff3', 'gff']:
        gff_file = input_path / f'reference.{ext}'
        if gff_file.exists():
            detected['gff'] = str(gff_file)
            break

    # Scaffold files
    for ext in ['fna', 'fasta', 'fa']:
        scaffold_file = input_path / f'scaffolds.{ext}'
        if scaffold_file.exists():
            detected['scaffolds'] = str(scaffold_file)
            break

    for ext in ['gff3', 'gff']:
        scaffold_gff = input_path / f'scaffolds.{ext}'
        if scaffold_gff.exists():
            detected['scaffold_gff'] = str(scaffold_gff)
            break

    scaffold_faa = input_path / 'scaffolds.faa'
    if scaffold_faa.exists():
        detected['scaffold_faa'] = str(scaffold_faa)

    # Contig files
    for ext in ['fna', 'fasta', 'fa']:
        contig_file = input_path / f'contigs.{ext}'
        if contig_file.exists():
            detected['contigs'] = str(contig_file)
            break

    for ext in ['gff3', 'gff']:
        contig_gff = input_path / f'contigs.{ext}'
        if contig_gff.exists():
            detected['contig_gff'] = str(contig_gff)
            break

    contig_faa = input_path / 'contigs.faa'
    if contig_faa.exists():
        detected['contig_faa'] = str(contig_faa)

    return detected


def create_smart_output_dir(base_output_dir, mode_type):
    """
    Create smart subfolder structure for output.

    Args:
        base_output_dir: Base output directory
        mode_type: 'contig', 'scaffold', 'comparison', or 'all'

    Returns:
        dict: Paths for contig, scaffold, and/or comparison output
    """
    base = Path(base_output_dir)

    if mode_type == 'contig':
        output_path = base / 'contig_alignment'
        return {'contig': output_path}
    elif mode_type == 'scaffold':
        output_path = base / 'scaffold_alignment'
        return {'scaffold': output_path}
    elif mode_type == 'comparison':
        output_path = base / 'assembly_comparison'
        return {'comparison': output_path}
    else:  # all
        return {
            'contig': base / 'contig_alignment',
            'scaffold': base / 'scaffold_alignment',
            'comparison': base / 'assembly_comparison'
        }


def main():
    """Main program entry point."""

    parser = argparse.ArgumentParser(
        description='GenomeViz - Circular Genome Assembly Visualization',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Simplified usage with auto-detection (recommended)
  %(prog)s --input input_dir/ --output results/ --no-gene-alignments

  # Run all modes (contig, scaffold, and comparison)
  %(prog)s --input input_dir/ --output results/ --mode all

  # Run only scaffold mode
  %(prog)s --input input_dir/ --output results/ --mode scaffold

  # Run only comparison mode
  %(prog)s --input input_dir/ --output results/ --mode comparison

  # Manual file specification - contig mode
  %(prog)s --reference ref.fasta --contig asm.fasta --reference-gff genes.gff --output results/ --mode contig

  # Manual file specification - scaffold mode with custom thresholds
  %(prog)s --reference ref.fasta --scaffold asm.fasta --reference-gff genes.gff --output results/ --mode scaffold \\
           --min-gap 2000 --min-inversion 1000

  # Skip auto-orientation
  %(prog)s --reference ref.fasta --contig asm.fasta --reference-gff genes.gff --output results/ --mode contig \\
           --no-auto-orient

  # Manually set origin position
  %(prog)s --reference ref.fasta --scaffold asm.fasta --reference-gff genes.gff --output results/ --mode scaffold \\
           --origin 150000

  # Disable origin rotation
  %(prog)s --reference ref.fasta --contig asm.fasta --reference-gff genes.gff --output results/ --mode contig \\
           --origin 0

  # Scaffold-contig comparison mode (manual files)
  %(prog)s --reference ref.fasta --scaffold scaffolds.fasta --contig contigs.fasta \\
           --reference-gff ref.gff --output results/ --mode comparison

  # With GFF annotations for all assemblies
  %(prog)s --reference ref.fasta --scaffold scaffolds.fasta --contig contigs.fasta \\
           --reference-gff ref.gff --scaffold-gff scaffolds.gff --contig-gff contigs.gff --output results/ --mode comparison

Auto-detection expects files named:
  - reference.{fna,fasta,fa} and reference.{gff,gff3}
  - scaffolds.{fna,fasta,fa}, scaffolds.{gff,gff3}, scaffolds.faa
  - contigs.{fna,fasta,fa}, contigs.{gff,gff3}, contigs.faa

For more information: https://github.com/Aaron-Thiel/GenomeViz
        """
    )

    # Input directory for automatic file detection
    parser.add_argument('--input',
                       help='Input directory containing reference, scaffold, and contig files. '
                            'Automatically detects files named: reference.{fna,fasta,gff3}, '
                            'scaffolds.{fna,fasta,gff3,faa}, contigs.{fna,fasta,gff3,faa}. '
                            'Manual file arguments override auto-detection.')

    # Required arguments (or can be auto-detected with --input)
    parser.add_argument('--reference',
                       help='Reference genome (FASTA format)')
    parser.add_argument('--reference-gff',
                       help='Gene annotations for reference (GFF3 format)')
    parser.add_argument('--output', required=True,
                       help='Output directory')

    # Mode selection
    parser.add_argument('--mode', choices=['contig', 'scaffold', 'comparison', 'all'],
                       help='Analysis mode: "contig" (contig vs reference), '
                            '"scaffold" (scaffold vs reference), "comparison" (scaffold vs contig), '
                            'or "all" (run all three). Auto-detected if not specified.')

    # Input assembly files
    parser.add_argument('--scaffold',
                       help='Scaffold assembly (FASTA format)')
    parser.add_argument('--contig',
                       help='Contig assembly (FASTA format)')

    # Comparison mode annotations (optional, if --scaffold-gff and --contig-gff not provided, uses --reference-gff)
    parser.add_argument('--scaffold-gff',
                       help='Gene annotations for scaffolds (GFF3 format) - optional for comparison mode')
    parser.add_argument('--contig-gff',
                       help='Gene annotations for contigs (GFF3 format) - optional for comparison mode')

    # Analysis parameters
    parser.add_argument('--preset', default='asm10', choices=['asm5', 'asm10', 'asm20'],
                       help='Minimap2 preset for alignment (default: asm10)')
    parser.add_argument('--origin', type=int, default=None,
                       help='Manually set origin position (bp). Use 0 to disable origin rotation. '
                            'If not specified, will auto-detect oriC from GFF.')
    parser.add_argument('--min-gap', type=int, default=1000,
                       help='Minimum gap size to report (default: 1000 bp)')
    parser.add_argument('--min-inversion', type=int, default=500,
                       help='Minimum inversion size to report (default: 500 bp)')

    # Orientation and reference preparation
    parser.add_argument('--no-auto-orient', action='store_true',
                       help='Skip automatic orientation detection')

    # Visualization output options
    parser.add_argument('--no-static', action='store_true',
                       help='Skip static PNG plots (circular and linear)')
    parser.add_argument('--no-interactive', action='store_true',
                       help='Skip interactive HTML plots (circular and linear)')
    parser.add_argument('--no-comparison', action='store_true',
                       help='Skip comparison mode visualizations')
    parser.add_argument('--no-gene-alignments', action='store_true',
                       help='Skip individual gene alignment file generation (disables gene clicking feature)')

    # Execution control
    parser.add_argument('--force', action='store_true',
                       help='Force reprocessing of all analyses, ignoring cache')
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')

    args = parser.parse_args()

    # Auto-detect files if --input provided
    detected_files = {}
    if args.input:
        if not os.path.exists(args.input):
            parser.error(f"Input directory does not exist: {args.input}")

        detected_files = auto_detect_files(args.input)

        # Use detected files if manual arguments not provided
        if not args.reference and 'reference' in detected_files:
            args.reference = detected_files['reference']
            print(f"Auto-detected reference: {args.reference}")

        if not args.reference_gff and 'gff' in detected_files:
            args.reference_gff = detected_files['gff']
            print(f"Auto-detected reference GFF: {args.reference_gff}")

        if not args.scaffold and 'scaffolds' in detected_files:
            args.scaffold = detected_files['scaffolds']
            print(f"Auto-detected scaffold: {args.scaffold}")

        if not args.contig and 'contigs' in detected_files:
            args.contig = detected_files['contigs']
            print(f"Auto-detected contig: {args.contig}")

        if not args.scaffold_gff and 'scaffold_gff' in detected_files:
            args.scaffold_gff = detected_files['scaffold_gff']
            print(f"Auto-detected scaffold GFF: {args.scaffold_gff}")

        if not args.contig_gff and 'contig_gff' in detected_files:
            args.contig_gff = detected_files['contig_gff']
            print(f"Auto-detected contig GFF: {args.contig_gff}")

    # Auto-detect mode if not specified
    if not args.mode:
        if args.scaffold and args.contig:
            args.mode = 'all'  # If both scaffold and contig present, run all modes
            print("Auto-detected mode: all")
        elif args.scaffold:
            args.mode = 'scaffold'
            print("Auto-detected mode: scaffold")
        elif args.contig:
            args.mode = 'contig'
            print("Auto-detected mode: contig")
        else:
            parser.error("Cannot auto-detect mode. Please specify --mode or provide --scaffold and/or --contig.")

    # Validate required arguments
    if not args.reference:
        parser.error("--reference is required (or use --input with reference.{fna,fasta})")
    if not args.reference_gff:
        parser.error("--reference-gff is required (or use --input with reference.{gff,gff3})")

    # Validate mode-specific requirements
    if args.mode == 'all':
        if not (args.scaffold and args.contig):
            parser.error("--mode all requires both --scaffold and --contig")
    elif args.mode == 'scaffold':
        if not args.scaffold:
            parser.error("--mode scaffold requires --scaffold")
    elif args.mode == 'contig':
        if not args.contig:
            parser.error("--mode contig requires --contig")
    elif args.mode == 'comparison':
        if not (args.scaffold and args.contig):
            parser.error("--mode comparison requires both --scaffold and --contig")

    # Validate input files
    try:
        files_to_check = [args.reference, args.reference_gff]
        if args.scaffold:
            files_to_check.append(args.scaffold)
        if args.contig:
            files_to_check.append(args.contig)
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

    # Initialize cache variables
    contig_mode_cache_key = None
    contig_input_files = {}
    scaffold_mode_cache_key = None
    scaffold_input_files = {}
    comparison_mode_cache_key = None
    comparison_input_files = {}

    # Determine run mode
    run_scaffold_mode = args.mode in ['scaffold', 'all']
    run_contig_mode = args.mode in ['contig', 'all']
    run_comparison_mode = args.mode in ['comparison', 'all']

    # Create output directory
    base_output_dir = Path(args.output)
    base_output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create smart subfolders for each mode
    output_dirs = create_smart_output_dir(base_output_dir, args.mode)
    for mode_dir in output_dirs.values():
        mode_dir.mkdir(parents=True, exist_ok=True)
    
    # Load cache from base output directory
    cache = load_cache(base_output_dir) if not args.force else {}
    
    output_dir = Path(args.output)

    print("=" * 70)
    print(f"GenomeViz v{__version__}")
    print("=" * 70)
    print(f"\nOutput directory: {output_dir.absolute()}")

    mode_info = []
    if run_contig_mode:
        mode_info.append("Contig (contig vs reference)")
    if run_scaffold_mode:
        mode_info.append("Scaffold (scaffold vs reference)")
    if run_comparison_mode:
        mode_info.append("Comparison (scaffold vs contig)")
    print(f"Mode: {' + '.join(mode_info)}")

    # Track temporary files for cleanup
    temp_files_to_cleanup = []

    # Step 1: Parse gene annotations (do this first to find oriC)
    print("\n" + "="*70)
    print("[1] PARSING GENE ANNOTATIONS")
    print("="*70)
    gff_parser = GFFParser(args.reference_gff)

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
        origin_position, origin_seqid = find_oric_position(args.reference_gff)

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
    # CONTIG MODE: Contig vs Reference
    # =========================================================================
    if run_contig_mode:
        print("\n" + "="*70)
        print("CONTIG MODE: Contig vs Reference")
        print("="*70)

        # Set output directory for this mode
        if args.mode == 'all':
            output_dir = output_dirs['contig']
            print(f"Output directory: {output_dir.absolute()}")

        # Check cache for contig mode
        contig_mode_cache_key = 'contig_mode'
        contig_input_files = {
            'reference': args.reference,
            'contig': args.contig,
            'gff': args.reference_gff
        }

        # Initialize cache entry if it doesn't exist
        if contig_mode_cache_key not in cache:
            cache[contig_mode_cache_key] = {}

        if not args.force and check_cache_valid(cache, contig_mode_cache_key, contig_input_files):
            print("\n✓ Using cached results (unchanged input files)")
        else:
            if 'input_hashes' in cache[contig_mode_cache_key]:
                print("\n⚠ Input files changed, reprocessing...")

        # Orientation detection for contigs
        print("\n  [2b] Contig Orientation Detection")
        print("  " + "-"*66)

        if not args.no_auto_orient:
            base_output_dir.mkdir(parents=True, exist_ok=True)
            oriented_contig_path = base_output_dir / '.oriented_contig.fasta'

            if oriented_contig_path.exists() and not args.force:
                contig_to_use = str(oriented_contig_path)
                print(f"  ✓ Using cached oriented contig: {oriented_contig_path.name}")
            else:
                detector = OrientationDetector(reference_to_use, args.contig, preset=args.preset)
                contig_to_use, temp_asm = detector.detect_and_correct()

                if contig_to_use != args.contig:
                    with open(contig_to_use, 'r') as src:
                        with open(oriented_contig_path, 'w') as dst:
                            dst.write(src.read())
                    contig_to_use = str(oriented_contig_path)
                    if temp_asm:
                        temp_files_to_cleanup.append(temp_asm)
        else:
            contig_to_use = args.contig
            print("  ℹ️  Orientation detection skipped (--no-auto-orient)")

        # Load reference and align
        print("\n" + "-"*70)
        print("Loading reference and aligning contig...")
        print("-"*70)
        aligner = GenomeAligner(reference_to_use, contig_to_use, preset=args.preset)
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
        summary_file = create_summary_report(output_dir, ref_sequences, args.contig, args.reference)
        print(f"\nSummary report: {summary_file}")
        print_final_summary(output_dir, ref_sequences)

    # =========================================================================
    # SCAFFOLD MODE: Scaffold vs Reference
    # =========================================================================
    if run_scaffold_mode:
        print("\n" + "="*70)
        print("SCAFFOLD MODE: Scaffold vs Reference")
        print("="*70)

        # Set output directory for this mode
        if args.mode == 'all':
            output_dir = output_dirs['scaffold']
            print(f"Output directory: {output_dir.absolute()}")

        # Check cache for scaffold mode
        scaffold_mode_cache_key = 'scaffold_mode'
        scaffold_input_files = {
            'reference': args.reference,
            'scaffold': args.scaffold,
            'gff': args.reference_gff
        }

        # Initialize cache entry if it doesn't exist
        if scaffold_mode_cache_key not in cache:
            cache[scaffold_mode_cache_key] = {}

        if not args.force and check_cache_valid(cache, scaffold_mode_cache_key, scaffold_input_files):
            print("\n✓ Using cached results (unchanged input files)")
        else:
            if 'input_hashes' in cache[scaffold_mode_cache_key]:
                print("\n⚠ Input files changed, reprocessing...")

        # Orientation detection for scaffolds
        print("\n  [2b] Scaffold Orientation Detection")
        print("  " + "-"*66)

        if not args.no_auto_orient:
            base_output_dir.mkdir(parents=True, exist_ok=True)
            oriented_scaffold_path = base_output_dir / '.oriented_scaffold.fasta'

            if oriented_scaffold_path.exists() and not args.force:
                scaffold_to_use = str(oriented_scaffold_path)
                print(f"  ✓ Using cached oriented scaffold: {oriented_scaffold_path.name}")
            else:
                detector = OrientationDetector(reference_to_use, args.scaffold, preset=args.preset)
                scaffold_to_use, temp_asm = detector.detect_and_correct()

                if scaffold_to_use != args.scaffold:
                    with open(scaffold_to_use, 'r') as src:
                        with open(oriented_scaffold_path, 'w') as dst:
                            dst.write(src.read())
                    scaffold_to_use = str(oriented_scaffold_path)
                    if temp_asm:
                        temp_files_to_cleanup.append(temp_asm)
        else:
            scaffold_to_use = args.scaffold
            print("  ℹ️  Orientation detection skipped (--no-auto-orient)")

        # Load reference and align
        print("\n" + "-"*70)
        print("Loading reference and aligning scaffold...")
        print("-"*70)
        aligner = GenomeAligner(reference_to_use, scaffold_to_use, preset=args.preset)
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
        summary_file = create_summary_report(output_dir, ref_sequences, args.scaffold, args.reference)
        print(f"\nSummary report: {summary_file}")
        print_final_summary(output_dir, ref_sequences)

    # =========================================================================
    # COMPARISON MODE: Scaffolds vs Contigs
    # =========================================================================
    if run_comparison_mode:
        print("\n" + "="*70)
        print("COMPARISON MODE: Scaffolds vs Contigs")
        print("="*70)
        
        # Set output directory for this mode
        if args.mode == 'all':
            output_dir = output_dirs['comparison']
            print(f"Output directory: {output_dir.absolute()}")

        # Check cache for comparison mode
        comparison_mode_cache_key = 'comparison_mode'
        comparison_input_files = {
            'reference': args.reference,
            'scaffold': args.scaffold,
            'contig': args.contig,
            'gff': args.reference_gff,
            'scaffold_gff': args.scaffold_gff or '',
            'contig_gff': args.contig_gff or ''
        }
        
        # Initialize cache entry if it doesn't exist
        if comparison_mode_cache_key not in cache:
            cache[comparison_mode_cache_key] = {}
        
        if not args.force and check_cache_valid(cache, comparison_mode_cache_key, comparison_input_files):
            print("\n✓ Using cached results (unchanged input files)")
        else:
            if 'input_hashes' in cache[comparison_mode_cache_key]:
                print("\n⚠ Input files changed, reprocessing...")

        from src.assembly_comparison import AssemblyComparator, MultiAssemblyGFFParser
        from src.comparison_visualizer import (
            ScaffoldSequence,
            ComparisonLinearVisualizer,
            ComparisonInteractiveLinearVisualizer,
            ComparisonCircularVisualizer,
            ComparisonIndexGenerator,
            ComparisonGeneAlignmentVisualizer
        )
        from src.gene_comparison import GeneIntegrityAnalyzer

        # Step C1: Orientation detection for scaffolds
        print("\n  [C1] Scaffold Orientation Detection")
        print("  " + "-"*66)
        scaffolds_to_use = args.scaffold
        if not args.no_auto_orient:
            # Save oriented scaffolds to parent output directory for consistent caching across modes
            # This allows standard → both transition without reprocessing
            base_output_dir.mkdir(parents=True, exist_ok=True)
            oriented_scaffolds_path = base_output_dir / '.oriented_scaffolds.fasta'
            
            # Check if we already have a cached oriented scaffolds
            if oriented_scaffolds_path.exists() and not args.force:
                scaffolds_to_use = str(oriented_scaffolds_path)
                print(f"  ✓ Using cached oriented scaffolds: {oriented_scaffolds_path.name}")
            else:
                scaffold_detector = OrientationDetector(reference_to_use, args.scaffold, preset=args.preset)
                scaffolds_to_use, temp_scaffolds = scaffold_detector.detect_and_correct()
                
                # Save to persistent location instead of temp
                if scaffolds_to_use != args.scaffold:
                    # Copy the corrected scaffolds to output directory
                    with open(scaffolds_to_use, 'r') as src:
                        with open(oriented_scaffolds_path, 'w') as dst:
                            dst.write(src.read())
                    scaffolds_to_use = str(oriented_scaffolds_path)
                    if temp_scaffolds:
                        temp_files_to_cleanup.append(temp_scaffolds)
        else:
            print("  ℹ️  Orientation detection skipped (--no-auto-orient)")

        # Step C2: Orientation detection for contigs
        print("\n  [C2] Contig Orientation Detection")
        print("  " + "-"*66)
        contigs_to_use = args.contig
        if not args.no_auto_orient:
            # Save oriented contigs to parent output directory for consistent caching across modes
            # This allows standard → both transition without reprocessing
            base_output_dir.mkdir(parents=True, exist_ok=True)
            oriented_contigs_path = base_output_dir / '.oriented_contigs.fasta'
            
            # Check if we already have a cached oriented contigs
            if oriented_contigs_path.exists() and not args.force:
                contigs_to_use = str(oriented_contigs_path)
                print(f"  ✓ Using cached oriented contigs: {oriented_contigs_path.name}")
            else:
                contig_detector = OrientationDetector(reference_to_use, args.contig, preset=args.preset)
                contigs_to_use, temp_contigs = contig_detector.detect_and_correct()
                
                # Save to persistent location instead of temp
                if contigs_to_use != args.contig:
                    # Copy the corrected contigs to output directory
                    with open(contigs_to_use, 'r') as src:
                        with open(oriented_contigs_path, 'w') as dst:
                            dst.write(src.read())
                    contigs_to_use = str(oriented_contigs_path)
                    if temp_contigs:
                        temp_files_to_cleanup.append(temp_contigs)
        else:
            print("  ℹ️  Orientation detection skipped (--no-auto-orient)")

        # Step C3: Align scaffolds to reference
        print("\n  [C3] Aligning scaffolds to reference...")
        if origin_position and origin_position != 0:
            print(f"        Using rotated reference (origin at {origin_position:,} bp)")
        scaffold_aligner = GenomeAligner(reference_to_use, scaffolds_to_use, preset=args.preset)
        scaffold_aligner.load_reference()
        scaffold_alignments = scaffold_aligner.align()

        # Step C4: Align contigs to reference
        print("\n  [C4] Aligning contigs to reference...")
        if origin_position and origin_position != 0:
            print(f"        Using rotated reference (origin at {origin_position:,} bp)")
        contig_aligner = GenomeAligner(reference_to_use, contigs_to_use, preset=args.preset)
        contig_aligner.load_reference()
        contig_alignments = contig_aligner.align()

        # Step C5: Parse additional GFF files
        print("\n  [C5] Parsing assembly GFF files...")
        multi_gff = MultiAssemblyGFFParser(
            reference_gff_parser=gff_parser,  # Use rotated gff_parser for correct coordinates
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

        # Track origin offset for coordinate consistency across modes
        if origin_position and origin_position != 0:
            print(f"  Origin offset applied: {origin_position:,} bp (oriC at position 0)")

        # Use output_dir directly (it's already set to assembly_comparison/)
        comparison_dir = output_dir

        # Step C7: Generate visualizations for each scaffold
        print("\n  [C7] Generating comparison visualizations...")

        all_gene_results = {}
        rc_scaffold_info = {}  # Track RC'd scaffolds and their lengths for TSV coordinate consistency

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

            # Log coordinates for debugging rotation consistency
            if origin_position and origin_position != 0:
                print(f"      Scaffold maps to reference: {ref_start:,} - {ref_end:,} bp (rotated coords)")

            # Analyze gene integrity
            # Priority: scaffold GFF > reference GFF (required --gff)
            gene_results = None
            genes_are_reference_coords = False  # Track if genes are in reference coordinates
            
            # Get scaffold length from alignment data (needed for RC coordinate transformation)
            scaffold_length = None
            for ref_seqid_aln, alns in scaffold_alignments.items():
                for aln in alns:
                    if aln['query_name'] == scaffold_name or aln['query_name'] == scaffold_name.replace('_RC', ''):
                        scaffold_length = aln.get('query_length')
                        print(f"      Found scaffold_length={scaffold_length:,} from alignment of '{aln['query_name']}'")
                        break
                if scaffold_length:
                    break
            
            if scaffold_length is None:
                print(f"      Warning: Could not find scaffold_length for '{scaffold_name}'")
            
            # Track RC'd scaffolds for TSV coordinate consistency
            if scaffold_name.endswith('_RC') and scaffold_length:
                rc_scaffold_info[scaffold_name] = scaffold_length
                print(f"      Tracking RC'd scaffold: {scaffold_name} (length={scaffold_length:,})")
            
            if multi_gff.has_genes('scaffold'):
                # Use scaffold genes if scaffold GFF is provided
                gene_analyzer = GeneIntegrityAnalyzer(overlap_regions, multi_gff)
                gene_results = gene_analyzer.analyze_scaffold_genes_via_reference(
                    scaffold_name, overlap_regions, scaffold_length=scaffold_length
                )
                if gene_results:
                    all_gene_results[scaffold_name] = gene_results
                    # Scaffold genes now have ref_start/ref_end set if they map to reference
                    genes_are_reference_coords = True  # Use ref coords from gene results
                    print(f"      Using {len(gene_results)} scaffold genes")
            
            # Fall back to reference genes if scaffold genes not available
            # This happens when: no scaffold GFF, or scaffold was RC'd (coords don't match)
            if not gene_results and multi_gff.has_genes('reference'):
                gene_analyzer = GeneIntegrityAnalyzer(overlap_regions, multi_gff)
                
                # Get genes from EACH scaffold alignment region separately
                # This prevents fetching genes from regions between scaffold alignments
                gene_results = []
                for region in overlap_regions:
                    region_genes = gene_analyzer.analyze_reference_genes_in_range(
                        region.scaffold_ref_start, region.scaffold_ref_end, 
                        reference_seqid, [region]  # Only pass this region for coverage analysis
                    )
                    gene_results.extend(region_genes)
                
                if gene_results:
                    all_gene_results[scaffold_name] = gene_results
                    genes_are_reference_coords = True  # Reference genes already have absolute coords
                    print(f"      Using {len(gene_results)} reference genes in region")

            # Create ScaffoldSequence object for visualization
            # Pass reference_length for full genome context visualization
            # Pass origin_position to ensure coordinate consistency with contig/scaffold modes
            scaffold_seq = ScaffoldSequence(
                scaffold_name=scaffold_name,
                ref_start=ref_start,
                ref_end=ref_end,
                overlap_regions=overlap_regions,
                gene_results=gene_results,
                reference_length=total_reference_length,
                genes_are_reference_coords=genes_are_reference_coords,
                origin_offset=origin_position if origin_position else 0
            )

            # Create linear visualization
            if not args.no_comparison:
                # Static PNG
                if not args.no_static:
                    ComparisonLinearVisualizer.create_linear_plot(
                        scaffold_seq,
                        scaffold_dir / f'{safe_name}_linear.png',
                        reference_seqid=reference_seqid
                    )
                # Interactive HTML
                if not args.no_interactive:
                    ComparisonInteractiveLinearVisualizer.create_interactive_linear_plot(
                        scaffold_seq,
                        scaffold_dir / f'{safe_name}_interactive_linear.html',
                        reference_seqid=reference_seqid
                    )

            # Create circular visualizations
            circular_html_file = None
            if not args.no_comparison:
                # Interactive HTML
                if not args.no_interactive:
                    circular_html_file = scaffold_dir / f'{safe_name}_circular.html'
                    ComparisonCircularVisualizer.create_interactive_circular_plot(
                        scaffold_seq,
                        circular_html_file,
                        reference_seqid=reference_seqid
                    )
                # Static PNG
                if not args.no_static:
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


        # Step C11: New gene analyzer (find/classify genes new to scaffolds)
        if args.scaffold_gff and args.contig_gff:
            print("\n  [C11] Analyzing new genes in scaffolds...")
            
            # Get FAA file paths (derive from GFF paths)
            scaffold_faa = args.scaffold_gff.replace('.gff3', '.faa').replace('.gff', '.faa')
            contig_faa = args.contig_gff.replace('.gff3', '.faa').replace('.gff', '.faa')
            
            # Check if FAA files exist
            if os.path.exists(scaffold_faa) and os.path.exists(contig_faa):
                # Create analyzer with scaffold alignments for coordinate mapping
                # This maps genes to the global reference coordinate system used in visualizations
                analyzer = NewGeneAnalyzer(
                    scaffolds_fasta=scaffolds_to_use,
                    contigs_fasta=contigs_to_use,
                    scaffolds_gff=args.scaffold_gff,
                    contigs_gff=args.contig_gff,
                    output_dir=output_dir,
                    origin_offset=0,  # Not used when scaffold_alignments provided
                    preset=args.preset,
                    scaffold_alignments=scaffold_alignments  # Pass alignments for reference coordinate mapping
                )
                
                analysis_results = analyzer.analyze()
                
                if analysis_results['new_genes_count'] > 0:
                    print(f"\n  ✓ New gene analysis complete: {analysis_results['new_genes_count']} new genes found")
                    print(f"    Output directory: {analyzer.new_genes_dir}")
                else:
                    print("\n  ℹ️  No new genes found (all scaffold genes present in contigs)")
            else:
                print("\n  ℹ️  Skipping new gene analysis - FAA files not found")
                if not os.path.exists(scaffold_faa):
                    print(f"     Missing: {scaffold_faa}")
                if not os.path.exists(contig_faa):
                    print(f"     Missing: {contig_faa}")

    # =========================================================================
    # CACHE MANAGEMENT
    # =========================================================================

    # Update and save cache with input file hashes
    if run_contig_mode and contig_mode_cache_key in cache:
        cache[contig_mode_cache_key]['input_hashes'] = {
            fname: get_file_hash(fpath) for fname, fpath in contig_input_files.items()
        }
        cache[contig_mode_cache_key]['timestamp'] = datetime.now().isoformat()

    if run_scaffold_mode and scaffold_mode_cache_key in cache:
        cache[scaffold_mode_cache_key]['input_hashes'] = {
            fname: get_file_hash(fpath) for fname, fpath in scaffold_input_files.items()
        }
        cache[scaffold_mode_cache_key]['timestamp'] = datetime.now().isoformat()

    if run_comparison_mode and comparison_mode_cache_key in cache:
        cache[comparison_mode_cache_key]['input_hashes'] = {
            fname: get_file_hash(fpath) for fname, fpath in comparison_input_files.items()
        }
        cache[comparison_mode_cache_key]['timestamp'] = datetime.now().isoformat()

    # Save cache to base output directory
    save_cache(base_output_dir, cache)

    # =========================================================================
    # FINAL OUTPUT
    # =========================================================================
    print("\n" + "=" * 70)
    print("OUTPUT DIRECTORY STRUCTURE")
    print("=" * 70)
    print()
    
    # Print structure for the mode that was run
    if args.mode == 'all':
        if run_contig_mode:
            print("CONTIG ALIGNMENT MODE:")
            print("-" * 70)
            print_directory_tree(output_dirs['contig'], max_depth=2, skip_contents_of={'gene_alignments'})
            print()
        if run_scaffold_mode:
            print("SCAFFOLD ALIGNMENT MODE:")
            print("-" * 70)
            print_directory_tree(output_dirs['scaffold'], max_depth=2, skip_contents_of={'gene_alignments'})
            print()
        if run_comparison_mode:
            print("ASSEMBLY COMPARISON MODE:")
            print("-" * 70)
            print_directory_tree(output_dirs['comparison'], max_depth=2, skip_contents_of={'gene_alignments'})
    else:
        print_directory_tree(output_dir, max_depth=2, skip_contents_of={'gene_alignments'})

    print("\n" + "=" * 70)
    if args.mode == 'all':
        print(f"✓ All outputs saved to: {base_output_dir.absolute()}")
        if run_contig_mode:
            print(f"  - Contig alignment: {output_dirs['contig'].absolute()}")
        if run_scaffold_mode:
            print(f"  - Scaffold alignment: {output_dirs['scaffold'].absolute()}")
        if run_comparison_mode:
            print(f"  - Assembly comparison: {output_dirs['comparison'].absolute()}")
    else:
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