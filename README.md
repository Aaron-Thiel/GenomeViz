# GenomeViz

<p align="center">
  <img src="https://img.shields.io/badge/version-1.0.0-blue.svg" alt="Version">
  <img src="https://img.shields.io/badge/license-MIT-green.svg" alt="License">
  <img src="https://img.shields.io/badge/python-3.7+-blue.svg" alt="Python">
</p>

**GenomeViz** is a comprehensive tool for visualizing and comparing bacterial genome assemblies against reference sequences. It creates beautiful circular and linear plots with three information rings showing gene quality, alignment status, and contig mapping.

## ✨ Features

### Three-Ring Visualization System

![Circular Plot Example](examples/output/contig_1_circular.png)

- **Ring 1 (Outer)**: Gene quality assessment
  - Color-coded by quality score (excellent/good/fair/poor)
  - Shows gene coverage and identity
  
- **Ring 2 (Middle)**: Alignment status
  - 🟢 Green: Complete coverage (1x)
  - 🟠 Orange: Duplicated regions (>1x)
  - 🔴 Red: Inverted regions (reverse strand)
  - ⚪ Gray: Missing regions (gaps)
  
- **Ring 3 (Inner)**: Contig mapping
  - Each contig has a unique color
  - Shows which assembly contigs align to reference

### Interactive Plots

- **Zoom and pan** on any region
- **Hover information** for:
  - Gene names, positions, and quality scores
  - Contig names, positions, and identities
  - Alignment status and coverage
- **Export** as high-resolution images

### Additional Features

- ✅ Automatic contig orientation detection
- ✅ Static matplotlib plots for publication
- ✅ Linear plots for detailed analysis
- ✅ Gene-level quality statistics (CSV export)
- ✅ Detection of gaps, inversions, and duplications
- ✅ Support for multiple sequences (chromosome + plasmids)

## 🚀 Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/Aaron-Thiel/GenomeViz.git
cd GenomeViz

# Install dependencies
pip install -r requirements.txt
```

### Basic Usage

```bash
python genomeViz.py \
  --reference examples/input/reference.fna \
  --assembly examples/input/sample.fna \
  --gff examples/input/reference.gff3 \
  --output results/
```

### Output Files

The tool generates several output files:

```
results/
├── contig_1_circular.png          # Static circular plot
├── contig_1_interactive.html      # Interactive circular plot (zoom/hover)
├── contig_1_linear.png            # Linear multi-track plot
├── contig_1_gene_stats.csv        # Per-gene quality statistics
├── contig_mapping.json            # Detailed contig alignment data
└── summary_report.txt             # Overall summary report
```

## 📖 Detailed Documentation

### Command-Line Options

#### Required Arguments

- `--reference`: Reference genome in FASTA format
- `--assembly`: Assembly to compare in FASTA format
- `--gff`: Gene annotations in GFF3 format
- `--output`: Output directory for results

#### Optional Arguments

- `--preset`: Minimap2 alignment preset
  - `asm5`: More sensitive (>99% similar sequences)
  - `asm10`: Default, balanced (recommended)
  - `asm20`: Faster, less sensitive
  
- `--min-gap`: Minimum gap size to report (default: 1000 bp)
- `--min-inversion`: Minimum inversion size to report (default: 500 bp)
- `--no-auto-orient`: Skip automatic orientation detection
- `--no-circular`: Skip circular plot generation
- `--no-linear`: Skip linear plot generation
- `--no-interactive`: Skip interactive Plotly plot generation

### Understanding the Visualizations

#### Circular Plots

The circular plots provide a genome-wide overview:

1. **Gene Quality (Ring 1)**:
   - **Green** (Excellent): ≥95% quality
   - **Yellow** (Good): 85-95% quality
   - **Orange** (Fair): 70-85% quality
   - **Red** (Poor): <70% quality

2. **Alignment Status (Ring 2)**:
   - **Green**: Normal 1x coverage, correctly oriented
   - **Orange**: Duplicated regions (>1x coverage)
   - **Red**: Inverted contigs (need reverse complement)
   - **Gray**: Missing sequence (gaps in assembly)

3. **Contig Mapping (Ring 3)**:
   - Each color represents a different assembly contig
   - Shows which contigs align to each region

#### Linear Plots

Linear plots show five tracks:
1. Gene-level quality
2. Coverage depth
3. Sequence identity
4. Alignment strand
5. Contig mapping

These are useful for detailed analysis of specific regions.

#### Interactive HTML Plots

The interactive plots allow you to:
- **Zoom** by scrolling or dragging
- **Pan** by clicking and dragging
- **Hover** to see detailed information
- **Toggle** individual traces on/off
- **Export** as PNG images

### Interpreting Results

#### Quality Scores

Each gene receives a quality score (0-100) based on:
- **Coverage (30% weight)**: Percentage of gene covered by assembly
- **Identity (70% weight)**: Average sequence identity in covered regions

#### Gene Status Categories

- **Complete**: ≥95% coverage AND ≥90% identity
- **Incomplete**: 50-95% coverage
- **Divergent**: <90% identity (may indicate gene variation)
- **Missing**: <50% coverage (likely absent from assembly)

#### Common Patterns

**Pattern 1: High-Quality Assembly**
```
Ring 1: Mostly green (excellent gene quality)
Ring 2: Mostly green (complete coverage)
Ring 3: Few colors (assembled into few contigs)
```

**Pattern 2: Fragmented Assembly**
```
Ring 1: Mix of colors with gaps
Ring 2: Green + gray segments (gaps present)
Ring 3: Many colors (many small contigs)
```

**Pattern 3: Misassembled Regions**
```
Ring 1: Good overall
Ring 2: Red segments (inversions)
Ring 3: Some contigs in unusual positions
```

## 🔧 Advanced Usage

### Adjusting Detection Thresholds

For small plasmids (<100 kb):
```bash
python genomeViz.py \
  --reference plasmid.fna \
  --assembly assembly.fna \
  --gff plasmid.gff3 \
  --output results/ \
  --min-gap 500 \
  --min-inversion 200
```

### Disabling Auto-Orientation

To see original contig orientations:
```bash
python genomeViz.py \
  --reference ref.fna \
  --assembly asm.fna \
  --gff genes.gff3 \
  --output results/ \
  --no-auto-orient
```

### Publication-Ready Plots Only

Skip interactive HTML generation:
```bash
python genomeViz.py \
  --reference ref.fna \
  --assembly asm.fna \
  --gff genes.gff3 \
  --output results/ \
  --no-interactive
```

## 📊 Example Data

The repository includes example bacterial genome data:

- **Reference**: ~1.9 Mbp circular chromosome
- **Assembly**: 4-contig draft assembly
- **Annotations**: Bakta-annotated GFF3 file

Run the example:
```bash
python genomeViz.py \
  --reference examples/input/reference.fna \
  --assembly examples/input/sample.fna \
  --gff examples/input/reference.gff3 \
  --output examples/output/
```

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## 📝 Citation

If you use GenomeViz in your research, please cite:

```
This Github Page
```

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🐛 Issues and Support

If you encounter any problems or have suggestions, please [open an issue](https://github.com/Aaron-Thiel/GenomeViz/issues).

## 📚 Dependencies

GenomeViz relies on these excellent tools:
- [minimap2/mappy](https://github.com/lh3/minimap2) - Fast sequence alignment
- [BioPython](https://biopython.org/) - Biological sequence handling
- [matplotlib](https://matplotlib.org/) - Static plotting
- [Plotly](https://plotly.com/python/) - Interactive plotting
- [NumPy](https://numpy.org/) & [Pandas](https://pandas.pydata.org/) - Data processing

## 🙏 Acknowledgments

GenomeViz was developed to make genome assembly visualization accessible to both bioinformaticians and lab biologists. Special thanks to the developers of the underlying tools that make this possible.

---

**Made with ❤️ for the genomics community**
