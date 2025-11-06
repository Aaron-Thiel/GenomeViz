# Changelog

All notable changes to GenomeRings will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2024-11-06

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
