# Changelog

All notable changes to GenomeViz will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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

- **1.0.0**: Initial public release
