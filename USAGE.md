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

### Static Circular Plots (*.png)

High-resolution circular plots suitable for publications.

**Rings**:
1. **Outer**: Gene quality (green → yellow → orange → red)
2. **Middle**: Alignment status (green/orange/red/gray)
3. **Inner**: Contig mapping (unique colors per contig)

**Resolution**: 300 DPI
**Size**: ~1-2 MB per plot

### Interactive Plots (*.html)

HTML files with interactive features:
- **Zoom**: Scroll or drag to zoom
- **Pan**: Click and drag to move around
- **Hover**: See detailed information
- **Legend**: Click to show/hide traces
- **Export**: Download as PNG

**Features**:
- Gene information (name, position, quality)
- Contig details (name, coverage, identity)
- Alignment status (type, position, size)

**Usage**: Open in any web browser

### Linear Plots (*.png)

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

# Skip interactive HTML plots
python genomeViz.py ... --no-interactive

# Skip linear plots
python genomeViz.py ... --no-linear

# Only generate statistics
python genomeViz.py ... --no-circular --no-interactive --no-linear
```

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

**A**: Yes! The script uses standard color codes. Search for color definitions in the code:
- Gene quality: Interactive Plotly definitions at `genomeViz.py` lines ~627 and static matplotlib definitions at lines ~941
- Alignment status: Interactive Plotly definitions at `genomeViz.py` line ~689 and static matplotlib definitions at line ~979
- Edit and save

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

## Getting Help

- **Issues**: https://github.com/Aaron-Thiel/GenomeViz/issues
- **Discussions**: https://github.com/Aaron-Thiel/GenomeViz/discussions
- **Email**: [aaron.chris.thiel@gmail.com]

---

Happy genome visualizing! 🧬
