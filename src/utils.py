"""
Utility functions for GenomeViz

Provides helper functions for:
- File validation
- Summary report generation
- Output directory management
- Directory tree visualization

GitHub: https://github.com/Aaron-Thiel/GenomeViz
License: MIT
"""

import os
import numpy as np
import pandas as pd
from pathlib import Path


def validate_files(reference, assembly, gff):
    """
    Validate that input files exist and are readable.

    Args:
        reference (str): Path to reference FASTA file
        assembly (str): Path to assembly FASTA file
        gff (str): Path to GFF3 annotation file

    Raises:
        FileNotFoundError: If any required file is missing
    """
    files = {
        'Reference': reference,
        'Assembly': assembly,
        'GFF3': gff
    }

    for name, path in files.items():
        if not os.path.exists(path):
            raise FileNotFoundError(f"{name} file not found: {path}")


def create_summary_report(output_dir, ref_sequences, assembly_path, reference_path):
    """
    Create a summary report of the analysis.

    Args:
        output_dir (Path): Output directory path
        ref_sequences (dict): Dictionary of ReferenceSequence objects
        assembly_path (str): Path to assembly file
        reference_path (str): Path to reference file

    Returns:
        Path: Path to created summary report file
    """
    summary_file = output_dir / 'summary_report.txt'

    with open(summary_file, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("GENOME ASSEMBLY COMPARISON SUMMARY\n")
        f.write("=" * 70 + "\n\n")

        f.write(f"Assembly: {Path(assembly_path).name}\n")
        f.write(f"Reference: {Path(reference_path).name}\n\n")

        for seqid, ref_seq in ref_sequences.items():
            f.write(f"\n{ref_seq.seq_type.upper()}: {seqid}\n")
            f.write("-" * 70 + "\n")
            f.write(f"Length: {ref_seq.length:,} bp\n")
            f.write(f"Alignments: {len(ref_seq.alignments)}\n")
            f.write(f"Contigs aligned: {len(ref_seq.contig_mapping)}\n")
            f.write(f"Genes: {len(ref_seq.genes)}\n")

            if ref_seq.gene_stats:
                # Gene quality summary
                complete = len([g for g in ref_seq.gene_stats if g['status'] == 'complete'])
                incomplete = len([g for g in ref_seq.gene_stats if g['status'] == 'incomplete'])
                divergent = len([g for g in ref_seq.gene_stats if g['status'] == 'divergent'])
                missing = len([g for g in ref_seq.gene_stats if g['status'] == 'missing'])

                f.write(f"\nGene Quality Summary:\n")
                f.write(f"  Complete: {complete}\n")
                f.write(f"  Incomplete: {incomplete}\n")
                f.write(f"  Divergent: {divergent}\n")
                f.write(f"  Missing: {missing}\n")

                # Average statistics
                avg_coverage = sum(g['coverage_pct'] for g in ref_seq.gene_stats) / len(ref_seq.gene_stats)
                avg_identity = sum(g['avg_identity'] for g in ref_seq.gene_stats) / len(ref_seq.gene_stats)
                avg_quality = sum(g['quality_score'] for g in ref_seq.gene_stats) / len(ref_seq.gene_stats)

                f.write(f"\nAverage Gene Metrics:\n")
                f.write(f"  Coverage: {avg_coverage:.1f}%\n")
                f.write(f"  Identity: {avg_identity:.1f}%\n")
                f.write(f"  Quality Score: {avg_quality:.1f}/100\n")

            # Misassembly summary
            if ref_seq.misassemblies:
                inversions = len([m for m in ref_seq.misassemblies if m['type'] == 'inversion'])
                gaps = len([m for m in ref_seq.misassemblies if m['type'] == 'gap'])
                overlaps = len([m for m in ref_seq.misassemblies if m['type'] == 'overlap'])

                f.write(f"\nMisassemblies Detected:\n")
                f.write(f"  Inversions: {inversions}\n")
                f.write(f"  Gaps: {gaps}\n")
                f.write(f"  Overlaps: {overlaps}\n")

        f.write("\n" + "=" * 70 + "\n")
        f.write("Analysis complete!\n")
        f.write("=" * 70 + "\n")

    return summary_file


def save_gene_statistics(output_dir, ref_sequences):
    """
    Save gene statistics to CSV files.

    Args:
        output_dir (Path): Output directory path
        ref_sequences (dict): Dictionary of ReferenceSequence objects
    """
    for seqid, ref_seq in ref_sequences.items():
        if ref_seq.gene_stats:
            df = pd.DataFrame(ref_seq.gene_stats)
            seqid_dir = output_dir / seqid
            csv_output = seqid_dir / f'{seqid}_gene_stats.csv'
            df.to_csv(csv_output, index=False)
            print(f"  {csv_output}")


def print_final_summary(output_dir, ref_sequences):
    """
    Print final summary of analysis results.

    Args:
        output_dir (Path): Output directory path
        ref_sequences (dict): Dictionary of ReferenceSequence objects
    """
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    for seqid, ref_seq in ref_sequences.items():
        print(f"\n{ref_seq.seq_type.upper()}: {seqid}")
        print(f"  Length: {ref_seq.length:,} bp")
        print(f"  Contigs: {len(ref_seq.contig_mapping)}")
        print(f"  Genes: {len(ref_seq.genes)}")

        if ref_seq.gene_stats:
            avg_quality = np.mean([g['quality_score'] for g in ref_seq.gene_stats])
            print(f"  Avg Quality: {avg_quality:.1f}/100")

        if ref_seq.misassemblies:
            print(f"  Misassemblies: {len(ref_seq.misassemblies)}")


def print_directory_tree(directory, prefix="", is_last=True, max_depth=3, current_depth=0, skip_contents_of=None):
    """
    Print directory structure as a tree.

    Args:
        directory (Path): Directory to visualize
        prefix (str): Prefix for tree branches
        is_last (bool): Whether this is the last item in current level
        max_depth (int): Maximum depth to traverse
        current_depth (int): Current depth in tree
        skip_contents_of (set): Set of directory names to show but not traverse into
    """
    if skip_contents_of is None:
        skip_contents_of = set()

    if current_depth > max_depth:
        return

    directory = Path(directory)
    if not directory.exists():
        return

    # Print current directory/file
    connector = "└── " if is_last else "├── "
    if current_depth == 0:
        print(f"{directory.name}/")
    else:
        print(f"{prefix}{connector}{directory.name}/")

    # Get contents
    try:
        contents = sorted(directory.iterdir(), key=lambda x: (not x.is_dir(), x.name))
    except PermissionError:
        return

    # Filter out hidden files and __pycache__
    contents = [item for item in contents
                if not item.name.startswith('.')
                and item.name != '__pycache__']

    # Update prefix for children
    extension = "    " if is_last else "│   "
    new_prefix = prefix + extension if current_depth > 0 else ""

    # Print contents
    for i, item in enumerate(contents):
        is_last_item = (i == len(contents) - 1)

        if item.is_dir():
            # Check if we should skip the contents of this directory
            if item.name in skip_contents_of:
                # Show the directory but don't recurse into it
                connector = "└── " if is_last_item else "├── "
                print(f"{new_prefix}{connector}{item.name}/")
            else:
                # Recurse normally
                print_directory_tree(item, new_prefix, is_last_item, max_depth, current_depth + 1, skip_contents_of)
        else:
            connector = "└── " if is_last_item else "├── "
            print(f"{new_prefix}{connector}{item.name}")


def create_output_directories(output_dir, ref_sequences):
    """
    Create separate output directories for each sequence.

    Args:
        output_dir (Path): Base output directory
        ref_sequences (dict): Dictionary of ReferenceSequence objects

    Returns:
        dict: Mapping of sequence IDs to their output directories
    """
    output_dirs = {}
    for seqid in ref_sequences.keys():
        seqid_dir = output_dir / seqid
        seqid_dir.mkdir(exist_ok=True)
        output_dirs[seqid] = seqid_dir

    return output_dirs


def analyze_sequences(aligner, gff_parser, alignments_by_ref, contig_mapping, min_gap, min_inversion):
    """
    Analyze all reference sequences with alignments and gene data.

    Args:
        aligner: GenomeAligner object with reference sequences
        gff_parser: GFFParser object with gene annotations
        alignments_by_ref (dict): Alignments grouped by reference sequence
        contig_mapping (dict): Contig-to-reference mapping
        min_gap (int): Minimum gap size to report
        min_inversion (int): Minimum inversion size to report

    Returns:
        dict: Dictionary of ReferenceSequence objects with analysis results
    """
    from src.sequence import ReferenceSequence, GeneLevelAnalyzer
    from src.alignment import MisassemblyDetector

    ref_sequences = {}

    for seqid in aligner.reference_sequences:
        print(f"\nAnalyzing {seqid}...")

        ref_data = aligner.reference_sequences[seqid]
        genes = gff_parser.genes_by_seq.get(seqid, [])
        alignments = alignments_by_ref.get(seqid, [])

        ref_seq = ReferenceSequence(seqid, ref_data['sequence'], genes)
        ref_seq.alignments = alignments
        ref_seq.contig_mapping = contig_mapping.get(seqid, {})

        if alignments:
            detector = MisassemblyDetector(alignments,
                                          min_gap=min_gap,
                                          min_inversion_size=min_inversion)
            ref_seq.misassemblies = detector.detect()
            print(f"  Misassemblies: {len(ref_seq.misassemblies)}")
        else:
            print(f"  WARNING: No alignments found!")

        if genes:
            analyzer = GeneLevelAnalyzer(genes, alignments, ref_seq.length)
            ref_seq.gene_stats = analyzer.analyze()

            avg_quality = np.mean([g['quality_score'] for g in ref_seq.gene_stats])
            complete = len([g for g in ref_seq.gene_stats if g['status'] == 'complete'])
            print(f"  Genes: {len(genes)}, Avg Quality: {avg_quality:.1f}, Complete: {complete}")

        ref_sequences[seqid] = ref_seq

    return ref_sequences


def save_contig_mapping(output_dir, contig_mapping):
    """
    Save contig mapping to JSON file.

    Args:
        output_dir (Path): Output directory path
        contig_mapping (dict): Contig-to-reference mapping

    Returns:
        Path: Path to saved JSON file
    """
    import json

    mapping_file = output_dir / 'contig_mapping.json'
    with open(mapping_file, 'w') as f:
        serializable_mapping = {}
        for ref_name, contigs in contig_mapping.items():
            serializable_mapping[ref_name] = {}
            for contig_name, alns in contigs.items():
                serializable_mapping[ref_name][contig_name] = alns
        json.dump(serializable_mapping, f, indent=2)

    return mapping_file


def create_visualizations(ref_sequences, seqid_dirs, aligner, args):
    """
    Create all requested visualizations.

    Args:
        ref_sequences (dict): Dictionary of ReferenceSequence objects
        seqid_dirs (dict): Mapping of sequence IDs to output directories
        aligner: GenomeAligner object
        args: Command-line arguments

    Returns:
        tuple: (interactive_circular_files, interactive_linear_files)
    """
    interactive_circular_files = []
    interactive_linear_files = []

    # Conditionally import visualizers based on flags
    if not args.no_circular or not args.no_interactive:
        from src.circular_visualizer import CircularVisualizer
        from src.interactive_circular_visualizer import InteractiveCircularVisualizer

    if not args.no_linear or not args.no_interactive_linear:
        from src.linear_visualizer import LinearVisualizer
        from src.interactive_linear_visualizer import InteractiveLinearVisualizer

    # Create static circular plots
    if not args.no_circular:
        print("\nCreating static circular plots...")
        for seqid, ref_seq in ref_sequences.items():
            output_file = seqid_dirs[seqid] / f'{seqid}_circular.png'
            CircularVisualizer.create_circular_plot(ref_seq, output_file)

    # Create interactive circular plots
    if not args.no_interactive:
        print("\nCreating interactive circular plots...")
        for seqid, ref_seq in ref_sequences.items():
            output_file = seqid_dirs[seqid] / f'{seqid}_interactive_circular.html'
            InteractiveCircularVisualizer.create_interactive_circular_plot(ref_seq, output_file)
            InteractiveCircularVisualizer.add_help_button(
                output_file,
                plot_type='circular',
                gene_clicking_enabled=(not args.no_gene_alignments)
            )
            interactive_circular_files.append((output_file, ref_seq))

    # Create static linear plots
    if not args.no_linear:
        print("\nCreating static linear plots...")
        for seqid, ref_seq in ref_sequences.items():
            output_file = seqid_dirs[seqid] / f'{seqid}_linear.png'
            LinearVisualizer.create_linear_plot(ref_seq, output_file)

    # Create interactive linear plots
    if not args.no_interactive_linear:
        print("\nCreating interactive linear plots...")
        for seqid, ref_seq in ref_sequences.items():
            output_file = seqid_dirs[seqid] / f'{seqid}_interactive_linear.html'
            visualizer = InteractiveLinearVisualizer(
                ref_seq, aligner,
                gene_clicking_enabled=(not args.no_gene_alignments)
            )
            visualizer.create_interactive_linear_plot(output_file)
            interactive_linear_files.append((visualizer, ref_seq, seqid))

    return interactive_circular_files, interactive_linear_files


def generate_gene_alignments(interactive_linear_files, interactive_circular_files,
                            seqid_dirs, no_gene_alignments):
    """
    Generate gene alignment HTML files and add click handlers.

    Args:
        interactive_linear_files (list): List of (visualizer, ref_seq, seqid) tuples
        interactive_circular_files (list): List of (output_file, ref_seq) tuples
        seqid_dirs (dict): Mapping of sequence IDs to output directories
        no_gene_alignments (bool): Whether to skip gene alignment generation
    """
    if no_gene_alignments or not interactive_linear_files:
        return

    from src.interactive_circular_visualizer import InteractiveCircularVisualizer

    # Generate gene alignments
    for visualizer, ref_seq, seqid in interactive_linear_files:
        gene_align_dir = seqid_dirs[seqid] / 'gene_alignments'
        gene_align_dir.mkdir(exist_ok=True)
        visualizer.generate_gene_alignments(gene_align_dir)

    # Add click handlers to circular plots
    if interactive_circular_files:
        print("\nAdding click handlers to circular plots...")
        for output_file, ref_seq in interactive_circular_files:
            try:
                InteractiveCircularVisualizer.add_circular_click_handler(output_file, ref_seq)
            except Exception as e:
                print(f"  Warning: Could not add click handler to {output_file.name}: {e}")
