# GenomeViz Usage Guide

## Table of Contents
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Input Files](#input-files)
- [Output Files](#output-files)
- [Interpreting Results](#interpreting-results)
- [Advanced Options](#advanced-options)
- [Troubleshooting](#troubleshooting)
- [FAQ](#faq)

## Installation

### Requirements
- Python 3.7 or higher
- pip package manager
- mamba

### Step-by-Step Installation

1. **Clone the repository**:
```bash
git clone https://github.com/Aaron-Thiel/GenomeViz.git
cd GenomeViz
```

2. **Install dependencies**:
```bash
mamba create -n genomeviz python=3.11 -c conda-forge -c bioconda   mappy biopython numpy pandas matplotlib plotly -y
pip install kaleido
```

3. **Test installation**:
```bash
python genomeViz.py --version
```

## Quick Start

### Basic Command

```bash
python genomeViz.py \
  --reference your_reference.fasta \
  --assembly your_assembly.fasta \
  --gff your_genes.gff3 \
  --output results/
```

### Example with Provided Data

```bash
python genomeViz.py \
  --reference examples/input/reference.fna \
  --assembly examples/input/sample.fna \
  --gff examples/input/reference.gff3 \
  --output examples/output/
```

## Input Files

### Reference Genome (--reference)

**Format**: FASTA (.fasta, .fna, .fa)

The reference genome should be a complete, high-quality genome assembly. This can include:
- Chromosome(s)
- Plasmid(s)
- Both

**Example**:
```
>chromosome
ATGCGTACGTAGCTGATCG...
>plasmid_pXYZ
GCTAGCTAGCTAGCTAG...
```

### Assembly (--assembly)

**Format**: FASTA (.fasta, .fna, .fa)

Your draft assembly to compare against the reference. Can contain:
- Multiple contigs
- Scaffolds
- Complete or partial assemblies

**Example**:
```
>contig_1
ATGCGTACGTAGCTGATCG...
>contig_2
GCTAGCTAGCTAGCTAG...
```

### Gene Annotations (--gff)

**Format**: GFF3 (.gff, .gff3)

Gene annotations for the reference genome. Must match reference sequence IDs.

**Supported tools**:
- Bakta (preferred)
- Prokka
- PGAP
- Any GFF3-compliant annotator

**Example**:
```
##gff-version 3
chromosome  Bakta  gene  1  1500  .  +  .  ID=gene001;Name=dnaA
chromosome  Bakta  CDS   1  1500  .  +  0  ID=cds001;Parent=gene001
```

**Important**: Sequence IDs in GFF must match FASTA headers in reference!

## Output Files

GenomeViz creates an organized output structure with separate directories for each reference sequence:

```
output/
├── sequence_1/                           # Directory for first reference sequence
│   ├── sequence_1_circular.png          # Static circular plot
│   ├── sequence_1_interactive_circular.html  # Interactive circular plot
│   ├── sequence_1_linear.png            # Static linear plot
│   ├── sequence_1_interactive_linear.html    # Interactive linear plot
│   ├── sequence_1_gene_stats.csv        # Per-gene quality statistics
│   └── gene_alignments/                 # Detailed gene alignments (for clicking)
│       ├── gene_001.html
│       ├── gene_002.html
│       └── ...
├── sequence_2/                           # If multiple sequences (e.g., plasmids)
│   └── (similar structure)
├── contig_mapping.json                   # All contig-to-reference mappings
└── summary_report.txt                    # Overall analysis summary
```

### Static Circular Plots (*.png)

High-resolution circular plots suitable for publications.

**Rings**:
1. **Outer**: Gene quality (green → yellow → orange → red)
2. **Middle**: Alignment status (green/orange/red/gray)
3. **Inner**: Contig mapping (unique colors per contig)

**Resolution**: 300 DPI
**Size**: ~1-2 MB per plot

### Interactive Circular Plots (*_interactive_circular.html)

HTML files with interactive circular visualization:
- **Zoom**: Scroll or drag to zoom
- **Pan**: Click and drag to move around
- **Hover**: See detailed information
- **Click genes**: Click any gene to view detailed alignment
- **Help button**: Click "?" for usage instructions
- **Legend**: Click to show/hide traces
- **Export**: Download as PNG

**Features**:
- Gene information (name, position, quality)
- Contig details (name, coverage, identity)
- Alignment status (type, position, size)
- Clickable genes open detailed alignment views

**Usage**: Open in any web browser

### Interactive Linear Plots (*_interactive_linear.html)

Fully interactive linear genome viewer with multi-level zooming.

**Five Comprehensive Tracks**:
1. **Gene Quality**: Color-coded gene quality scores (0-100)
2. **Contig Mapping**: Assembly contigs aligned to reference
3. **Coverage Depth**: Read/contig coverage across genome
4. **Alignment Identity**: Sequence identity percentage
5. **Misassemblies**: Inversions, gaps, and overlaps

**Multi-Level Zoom Functionality**:
- **Genome Level** (default view):
  - Overview of entire genome
  - All tracks visible simultaneously
  - Navigate using range slider at bottom

- **Gene Level** (medium zoom):
  - Individual genes become visible
  - Detailed gene annotations appear
  - Gene boundaries and quality scores clear

- **Nucleotide Level** (maximum zoom):
  - Sequence-level resolution
  - Individual nucleotide differences visible
  - SNPs, insertions, and deletions highlighted

**Interactive Controls**:
- **Scroll** to zoom in/out smoothly
- **Click and drag** to pan along genome
- **Range slider** at bottom for quick navigation
- **Hover** over any element for detailed information
- **Click genes** to view detailed alignment
- **Help button** (?) for usage instructions
- **Click legend** items to toggle tracks on/off
- **Camera icon** to export as high-resolution PNG

**Performance**: Automatically optimizes for large genomes with intelligent downsampling

**Usage**: Open in any web browser (Chrome, Firefox, Safari, Edge)

**Gene Clicking**: Click any gene in Track 1 to open a detailed alignment view showing nucleotide-level differences, gaps, and mismatches

### Static Linear Plots (*.png)

Multi-track linear visualizations showing:
1. Gene-level quality
2. Coverage depth
3. Sequence identity
4. Alignment strand
5. Contig mapping

**Best for**: Detailed analysis of specific regions

### Gene Statistics (*.csv)

Per-gene quality metrics in CSV format.

**Columns**:
- `name`: Gene name/ID
- `seqid`: Sequence ID (chromosome/plasmid)
- `start`: Start position (bp)
- `end`: End position (bp)
- `strand`: Strand (+/-)
- `length`: Gene length (bp)
- `coverage_pct`: Coverage percentage (0-100%)
- `avg_identity`: Average identity (0-100%)
- `quality_score`: Overall quality (0-100)
- `status`: Gene status (complete/incomplete/missing/divergent)

**Usage**: Import into Excel, R, Python for further analysis

### Contig Mapping (contig_mapping.json)

Detailed alignment information in JSON format.

**Structure**:
```json
{
  "chromosome": {
    "contig_1": [
      {
        "ref_start": 0,
        "ref_end": 50000,
        "query_start": 0,
        "query_end": 50000,
        "strand": "+",
        "identity": 99.5,
        "coverage": 50000
      }
    ]
  }
}
```

### Summary Report (summary_report.txt)

Human-readable text summary with:
- Assembly and reference names
- Per-sequence statistics
- Gene quality breakdown
- Misassembly counts
- Contig mapping summary

## Interpreting Results

### Gene Quality Ring (Ring 1)

**Colors**:
- 🟢 **Green (Excellent)**: Quality ≥95
  - Gene is well-assembled with high coverage and identity
  - No issues detected
  
- 🟡 **Yellow (Good)**: Quality 85-95
  - Gene is mostly complete
  - Minor coverage or identity issues
  
- 🟠 **Orange (Fair)**: Quality 70-85
  - Significant quality issues
  - May have gaps or lower identity
  
- 🔴 **Red (Poor)**: Quality <70
  - Major assembly problems
  - Likely incomplete or missing

**What to do**:
- **Green/Yellow**: No action needed
- **Orange**: Investigate further
- **Red**: Consider reassembly or targeted sequencing

### Alignment Status Ring (Ring 2)

**Colors**:
- 🟢 **Green (Complete)**: 1x coverage, correct orientation
  - Perfect alignment
  - No issues
  
- 🟠 **Orange (Duplicated)**: >1x coverage
  - Multiple contigs cover same region
  - Could be:
    - True biological duplication
    - Assembly artifact (collapse failure)
    - Repeat regions
  
- 🔴 **Red (Inverted)**: Reverse strand alignment
  - Contig aligned in wrong orientation
  - Tool auto-corrects this (unless --no-auto-orient)
  - Check if correction was appropriate
  
- ⚪ **Gray (Missing)**: No coverage
  - Gap in assembly
  - Could be:
    - Missing sequence data
    - Assembly failure
    - Difficult-to-assemble region (repeats, GC-rich)

**What to do**:
- **Green**: Perfect!
- **Orange**: Check if duplication is expected
- **Red**: Verify auto-orientation correction
- **Gray**: Target for additional sequencing or reassembly

### Contig Mapping Ring (Ring 3)

Shows which assembly contigs align to each region.

**Interpretation**:
- **Single color throughout**: Well-assembled, few contigs
- **Many colors**: Fragmented assembly
- **Color changes**: Contig boundaries
- **Gaps between colors**: Unaligned regions

**What to do**:
- Many small contigs → Consider long-read sequencing
- Specific problematic regions → Targeted reassembly
- Overall good → Assembly is successful

## Advanced Options

### Adjusting Alignment Sensitivity

```bash
# More sensitive (for very similar sequences >99%)
python genomeViz.py ... --preset asm5

# Default (recommended for most cases)
python genomeViz.py ... --preset asm10

# Faster (for divergent sequences or quick checks)
python genomeViz.py ... --preset asm20
```

### Adjusting Detection Thresholds

For small genomes or plasmids:
```bash
python genomeViz.py ... --min-gap 500 --min-inversion 200
```

For large genomes:
```bash
python genomeViz.py ... --min-gap 2000 --min-inversion 1000
```

### Skipping Specific Outputs

```bash
# Skip static circular plots
python genomeViz.py ... --no-circular

# Skip interactive circular plots
python genomeViz.py ... --no-interactive

# Skip static linear plots
python genomeViz.py ... --no-linear

# Skip interactive linear plots
python genomeViz.py ... --no-interactive-linear

# Skip gene alignment generation (disables gene clicking feature)
python genomeViz.py ... --no-gene-alignments

# Only generate statistics (skip all visualizations)
python genomeViz.py ... --no-circular --no-interactive --no-linear --no-interactive-linear
```

**Note**: Using `--no-gene-alignments` will disable the gene clicking feature in interactive plots but can save time for large genomes with many genes.

### Using Interactive Linear Plots

The interactive linear plot is generated by default and provides powerful exploration capabilities:

```bash
# Default: generates all plots including interactive linear
python genomeViz.py \
  --reference ref.fna \
  --assembly asm.fna \
  --gff genes.gff3 \
  --output results/

# Open the interactive linear plot
# results/{seqid}_interactive_linear.html
```

**Best practices**:
- Start at genome level to identify problematic regions
- Zoom in on specific genes for detailed inspection
- Use range slider for rapid genome navigation
- Export high-resolution views for presentations
- Toggle tracks to focus on specific data types

### Disabling Auto-Orientation

To see original contig orientations:
```bash
python genomeViz.py ... --no-auto-orient
```

## Troubleshooting

### Problem: "No alignments found"

**Cause**: Assembly doesn't align to reference

**Solutions**:
1. Check that assembly and reference are related
2. Try different `--preset` option
3. Verify file formats are correct
4. Check sequence IDs don't have special characters

### Problem: "File not found"

**Cause**: Invalid file path

**Solutions**:
1. Check file paths are correct
2. Use absolute paths instead of relative
3. Verify files exist: `ls -l your_file.fasta`

### Problem: Plots look empty

**Cause**: GFF sequence IDs don't match reference

**Solutions**:
1. Check GFF sequence IDs: `grep "^>" reference.fasta`
2. Check GFF lines: `grep -v "^#" genes.gff | cut -f1 | sort -u`
3. Make sure IDs match exactly

### Problem: Too many contigs in legend

**Limitation**: Legend only shows ≤15 contigs

**Solutions**:
1. Check `contig_mapping.json` for full list
2. Consider if assembly is too fragmented
3. This is normal for draft assemblies

### Problem: Script is slow

**Cause**: Large genomes take time

**Solutions**:
1. Use `--preset asm20` for faster alignment
2. Skip plots: `--no-circular --no-interactive --no-linear`
3. Run on compute cluster if available
4. Be patient - quality takes time!

### Problem: Out of memory

**Cause**: Large genomes (>10 Mbp) or many contigs

**Solutions**:
1. Increase system RAM
2. Process sequences separately
3. Use `--no-interactive` to save memory
4. Close other applications

## FAQ

### Q: What genome sizes are supported?

**A**: GenomeViz works best with bacterial genomes (1-10 Mbp). It can handle:
- Small plasmids (few kb)
- Large chromosomes (up to ~15 Mbp)
- Multiple sequences simultaneously

### Q: Can I use this for eukaryotic genomes?

**A**: GenomeViz is optimized for circular bacterial genomes. For eukaryotes:
- Linear chromosomes work but visualization is suboptimal
- Large genomes (>20 Mbp) may be slow
- Many chromosomes create many output files
- Consider chromosome-by-chromosome analysis

### Q: How do I cite GenomeViz?

**A**: Citation information coming soon. For now, reference the GitHub repository.

### Q: Can I modify the colors?

**A**: Yes! The code uses standard color codes. Color definitions are organized in separate visualizer files:
- **Gene quality colors**:
  - Interactive: `src/interactive_circular_visualizer.py` and `src/interactive_linear_visualizer.py`
  - Static: `src/circular_visualizer.py` and `src/linear_visualizer.py`
- **Alignment status colors**:
  - Interactive: `src/interactive_circular_visualizer.py`
  - Static: `src/circular_visualizer.py`
- Search for `COLOR_MAP`, `gene_colors`, or `quality_colors` in the respective files
- Edit the color values and save

### Q: What if my GFF is from a different tool?

**A**: GenomeViz works with standard GFF3. It looks for:
- `gene` or `CDS` feature types
- `Name` or `ID` attributes
- Most tools (Bakta, Prokka, PGAP) work fine

### Q: Can I run this on a cluster?

**A**: Yes! GenomeViz is a command-line tool and works great on clusters:
```bash
#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

python genomeViz.py \
  --reference ref.fasta \
  --assembly asm.fasta \
  --gff genes.gff3 \
  --output results/
```

### Q: How accurate is the quality score?

**A**: Quality scores combine:
- 30% coverage (how much of gene is covered)
- 70% identity (how similar the sequence is)

This provides a good overall assessment but should be validated for critical genes.

### Q: What does "orientation-corrected" mean?

**A**: If >50% of a contig aligns in reverse, GenomeViz automatically reverse-complements it for optimal alignment. Original assembly is preserved; correction is temporary.

### Q: Can I compare multiple assemblies?

**A**: Not directly. Run GenomeViz separately for each assembly:
```bash
for asm in assembly1.fa assembly2.fa assembly3.fa; do
  python genomeViz.py \
    --reference ref.fa \
    --assembly $asm \
    --gff genes.gff3 \
    --output results_$(basename $asm .fa)/
done
```

Then compare the outputs manually.

### Q: Does GenomeViz support origin of replication alignment?

**A**: Yes! GenomeViz automatically detects the origin of replication (oriC) from your GFF annotations and rotates the reference sequence to start at oriC. This ensures:
- Consistent visualization starting points across different genomes
- Biological relevance (starting at the replication origin)
- Better comparability between assemblies

**Requirements**:
- Your GFF file must contain an oriC feature annotation
- Common annotation tools (Bakta, Prokka) often include oriC in their outputs
- If oriC is not found, the tool proceeds normally without rotation

**How it works**:
- Automatically searches for "origin_of_replication" or "oriC" in GFF annotations
- Rotates reference sequence to start at detected position
- Updates all gene coordinates accordingly
- No user action required - it's completely automatic!

## Getting Help

- **Issues**: https://github.com/Aaron-Thiel/GenomeViz/issues
- **Discussions**: https://github.com/Aaron-Thiel/GenomeViz/discussions
- **Email**: [aaron.chris.thiel@gmail.com]

---

Happy genome visualizing! 🧬
