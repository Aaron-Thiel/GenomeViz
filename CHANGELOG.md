# Changelog

All notable changes to GenomeViz will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.3.0] - 2026-02-23

### Added
- **Scaffold-Contig Comparison Mode**
  - New `--mode comparison` to compare scaffolds against contigs using the reference as anchor
  - Detects overlapping regions between scaffold and contig assemblies via reference coordinates
  - Generates per-scaffold circular and linear comparison visualizations
  - New modules: `assembly_comparison.py`, `comparison_visualizer.py`, `gene_comparison.py`

- **Multi-Mode Analysis (`--mode all`)**
  - Run contig alignment, scaffold alignment, and comparison modes in a single invocation
  - Smart output directory structure: `contig_alignment/`, `scaffold_alignment/`, `assembly_comparison/`
  - Root-level sample dashboard linking all analysis modes

- **Auto-Detection Input Mode (`--input`)**
  - New `--input` flag for automatic file detection from a directory
  - Expects standard file names: `reference.{fna,fasta}`, `scaffolds.{fna,fasta,gff3,faa}`, `contigs.{fna,fasta,gff3,faa}`
  - Auto-detects analysis mode based on available files
  - Manual file arguments override auto-detected files

- **Result Caching**
  - MD5-based cache validation of input files
  - Skips expensive alignment/analysis when inputs unchanged
  - Visualization options can be toggled without reprocessing alignments
  - Use `--force` to override cache and reprocess everything

- **New Gene Analyzer**
  - Identifies genes in scaffolds not present in contigs (by protein hash)
  - Classifies new genes based on contig-scaffold alignment overlap
  - Generates CSV reports and interactive visualizations
  - Requires scaffold and contig GFF + FAA files

- **Gene Integrity Analysis**
  - Per-scaffold gene integrity assessment (complete, split, truncated, missing)
  - Supports both scaffold-derived and reference-derived gene annotations
  - Gene comparison reports with coverage and identity statistics

- **HTML Dashboards**
  - `AlignmentDashboardGenerator` for contig/scaffold alignment index pages
  - `SampleDashboardGenerator` for root-level sample overview
  - Interactive dashboards with alignment statistics and links to visualizations

- **Docker Support**
  - Dockerfile for containerized execution
  - Nextflow-compatible with `/app` on PATH

### Changed
- Version updated from 1.2.0 to 1.3.0
- CLI reorganized: `--assembly` replaced by `--scaffold` and `--contig` for explicit assembly type
- `--gff` renamed to `--reference-gff` for clarity; added `--scaffold-gff` and `--contig-gff`
- New `--mode` flag selects analysis mode: `contig`, `scaffold`, `comparison`, or `all`
- New `--no-static` flag replaces separate `--no-circular` and `--no-linear` flags
- `--no-interactive` now skips both circular and linear interactive plots
- Main script expanded from ~214 lines to ~1200 lines with multi-mode pipeline
- Import ordering cleaned up (stdlib first, then local modules)
- Bare `except` replaced with explicit exception types

### New Modules
- `src/assembly_comparison.py` - Scaffold-contig overlap detection via reference
- `src/comparison_visualizer.py` - Circular and linear comparison visualizations
- `src/gene_comparison.py` - Gene integrity analysis between assemblies
- `src/new_gene_analyzer.py` - Novel gene identification in scaffolds
- `src/dashboard_generator.py` - HTML dashboard generation

## [1.2.0] - 2025-11-25

### Added
- **Origin of Replication (oriC) Alignment**
  - Automatic detection of origin of replication (oriC) from GFF annotations
  - `find_oric_position()` function to locate oriC in genome annotations
  - `rotate_sequence_and_features()` function to rotate reference sequences and gene coordinates
  - Ensures visualizations consistently start at the biological origin of replication
  - Improves comparability across different genome assemblies and references
  - Works automatically when oriC annotation is present in GFF file
  - New `--origin` command-line flag to manually control origin position:
    - Auto-detects oriC by default (when flag not specified)
    - Use `--origin 0` to disable rotation
    - Use `--origin <position>` to manually set origin position in base pairs

### Changed
- Version updated from 1.1.0 to 1.2.0
- Enhanced README with better visualization examples and descriptions
- Improved .gitignore to exclude additional output directories

### Documentation
- Updated README.md with improved formatting and image placement
- Added centered visualization examples with descriptions
- Enhanced feature descriptions for better clarity

## [1.1.0]

### Added
- **Interactive Linear Plot with Multi-Level Zooming**
  - New `InteractiveLinearVisualizer` class for interactive linear genome visualization
  - Five comprehensive tracks:
    - Gene quality track with color-coded quality scores
    - Contig mapping track showing assembly alignment
    - Coverage depth track with smooth visualization
    - Alignment identity track with reference lines (90%, 95%)
    - Misassemblies track highlighting inversions, gaps, and overlaps
  - Advanced zoom functionality:
    - Genome level: Full genome overview with all tracks
    - Gene level: Zoomed view showing individual gene annotations
    - Nucleotide level: Detailed sequence differences, SNPs, and indels
  - Interactive features:
    - Smooth scroll-to-zoom and click-to-pan
    - Range slider for quick navigation across the genome
    - Detailed hover information at all zoom levels
    - Export capability to high-resolution PNG
    - Toggle tracks on/off in legend with organized legend groups
  - Performance optimizations:
    - Automatic downsampling for large genomes (>10 Mbp)
    - Efficient rendering using Plotly subplots

- **Clickable Gene Feature**
  - Click any gene in interactive plots to view detailed alignment information
  - Individual gene alignment HTML files generated automatically
  - MSA-style visualization showing nucleotide-level differences
  - CIGAR string parsing for accurate gap and mismatch display
  - Works in both circular and linear interactive plots
  - Command-line option `--no-gene-alignments` to skip generation if not needed

- **Help Button System**
  - Built-in help modal in all interactive HTML plots
  - Context-specific instructions for each plot type
  - Explains zoom, pan, legend, and export functionality
  - Documents gene clicking feature when enabled

- **Organized Output Structure**
  - Separate directory created for each reference sequence
  - All sequence-specific files grouped together
  - Gene alignment files organized in `gene_alignments/` subdirectory
  - Directory tree visualization printed to terminal
  - Cleaner, more intuitive output organization

- **Modular Code Architecture**
  - Split visualizers into separate files for better maintainability:
    - `circular_visualizer.py`: Static circular plots
    - `interactive_circular_visualizer.py`: Interactive circular plots with help and click handlers
    - `linear_visualizer.py`: Static linear plots
    - `interactive_linear_visualizer.py`: Interactive linear plots with multi-level zoom
    - `gene_alignment_visualizer.py`: Detailed gene alignment HTML generation
  - Utility functions extracted to `utils.py` for cleaner main file
  - Consistent file headers with GitHub repository link
  - Improved code organization following bioinformatics tool conventions

- **10-Step Workflow**
  - Clearer terminal output with numbered steps (1/10 through 10/10)
  - Steps: Orientation, Parse GFF, Load Reference, Align, Map Contigs, Analyze, Visualize, Gene Alignments, Statistics, Summary
  - Each step follows consistent pattern: title → check → execute
  - Progress tracking throughout analysis pipeline

- **Legend Improvements**
  - All tracks now have legend group titles for better organization
  - Contigs now visible in Track 2 legend with their assigned colors
  - Consistent legend styling across circular and linear plots
  - Easier navigation of complex multi-track visualizations

### Changed
- Version updated from 1.0.0 to 1.1.0
- Enhanced visualization workflow to include interactive linear plots by default
- Renamed `InteractivePlotlyVisualizer` to `InteractiveCircularVisualizer` for clarity
- Output files renamed: `*_interactive.html` → `*_interactive_circular.html`
- Main `genomeViz.py` file streamlined to ~214 lines (from ~3000 lines)
- Automatic orientation detection now returns tuple for better error handling
- Statistics and summary separated into distinct steps (9/10 and 10/10)
- Directory tree shows `gene_alignments/` folders without listing contents

### Performance
- Optimized coverage track rendering with intelligent downsampling
- Reduced memory footprint for large genomes through windowed processing
- Conditional imports for visualizers based on command-line flags
- More efficient file organization reduces directory clutter

## [1.0.0]

### Added
- Initial release of GenomeRings
- Three-ring circular visualization system
  - Ring 1: Gene quality assessment
  - Ring 2: Alignment status (complete/duplicated/inverted/missing)
  - Ring 3: Contig mapping
- Interactive Plotly plots with hover information and zoom
- Static matplotlib plots for publication (300 DPI)
- Linear multi-track plots for detailed analysis
- Automatic contig orientation detection
- Gene-level quality statistics export (CSV)
- Detailed contig mapping (JSON)
- Summary report generation (TXT)
- Support for multiple sequences (chromosome + plasmids)
- Configurable alignment presets (asm5/asm10/asm20)
- Adjustable detection thresholds for gaps and inversions
- Comprehensive documentation and examples

### Features
- Visualizes bacterial genome assemblies against reference
- Detects gaps, inversions, and duplications
- Color-coded quality assessment
- Works with standard bioinformatics formats (FASTA, GFF3)
- Compatible with popular annotation tools (Bakta, Prokka, PGAP)

### Dependencies
- Python 3.7+
- mappy (minimap2)
- biopython
- numpy
- pandas
- matplotlib
- plotly

## [Unreleased]

### Planned Features
- Batch comparison mode (multiple assemblies)
- Export to SVG format
- Customizable color schemes
- Integration with NCBI datasets
- Command-line auto-completion
- GUI interface option

---

## Version History

- **1.3.0**: Scaffold-contig comparison mode, multi-mode pipeline, caching, dashboards, Docker support
- **1.2.0**: Origin of replication alignment and documentation improvements
- **1.1.0**: Interactive linear plots with multi-level zoom and modular architecture
- **1.0.0**: Initial public release
