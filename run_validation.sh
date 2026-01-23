#!/bin/bash

# GenomeViz batch runner for validation samples
# Runs each subfolder in input_validation as an individual sample
# Supports caching, mode selection, and visualization options
#
# CACHING BEHAVIOR:
#   - GenomeViz automatically caches analysis results using MD5 hashes of input files
#   - If input files haven't changed, subsequent runs skip expensive alignment/analysis steps
#   - Visualization options can be toggled without reprocessing alignments:
#     * Run 1: --no-interactive (fast, skips interactive plots)
#     * Run 2: Re-run without --no-interactive (uses cached analysis, regenerates plots)
#   - Use --force to completely reprocess all samples and override cache
#
# EXAMPLE WORKFLOWS:
#   # Fast initial run: analysis + static plots only
#   ./run_validation.sh -m all --no-interactive --no-gene-alignments
#
#   # Later: regenerate with interactive plots (uses cache for analysis)
#   ./run_validation.sh -m all
#
#   # Reprocess everything
#   ./run_validation.sh -m all --force

# Default values
INPUT_DIR="examples/input_validation"
OUTPUT_DIR="examples/output_validation_fixed"
MODE="all"  # contig, scaffold, comparison, or all
FORCE_REPROCESS=false
SKIP_INTERACTIVE=false
SKIP_STATIC=false
SKIP_COMPARISON=false
SKIP_GENE_ALIGNMENTS=true
PRESET="asm10"
VERBOSE=false

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            INPUT_DIR="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -m|--mode)
            MODE="$2"
            if [[ ! "$MODE" =~ ^(contig|scaffold|comparison|all)$ ]]; then
                echo "ERROR: Invalid mode '$MODE'. Choose: contig, scaffold, comparison, or all"
                exit 1
            fi
            shift 2
            ;;
        --preset)
            PRESET="$2"
            if [[ ! "$PRESET" =~ ^(asm5|asm10|asm20)$ ]]; then
                echo "ERROR: Invalid preset '$PRESET'. Choose: asm5, asm10, or asm20"
                exit 1
            fi
            shift 2
            ;;
        -f|--force)
            FORCE_REPROCESS=true
            shift
            ;;
        --no-interactive)
            SKIP_INTERACTIVE=true
            shift
            ;;
        --no-static)
            SKIP_STATIC=true
            shift
            ;;
        --no-comparison)
            SKIP_COMPARISON=true
            shift
            ;;
        --no-gene-alignments)
            SKIP_GENE_ALIGNMENTS=true
            shift
            ;;
        -v|--verbose)
            VERBOSE=true
            shift
            ;;
        -h|--help)
            echo "GenomeViz Batch Validation Runner"
            echo ""
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Input/Output Options:"
            echo "  -i, --input INPUT_DIR              Input directory with sample subfolders"
            echo "                                     (default: examples/input_validation)"
            echo "  -o, --output OUTPUT_DIR            Output directory for results"
            echo "                                     (default: examples/output_validation_fixed)"
            echo ""
            echo "Analysis Mode:"
            echo "  -m, --mode MODE                    Analysis mode to run:"
            echo "                                       contig     - Contig vs reference"
            echo "                                       scaffold   - Scaffold vs reference"
            echo "                                       comparison - Scaffold vs contig"
            echo "                                       all        - All three modes (default)"
            echo ""
            echo "Alignment Parameters:"
            echo "  --preset PRESET                    Minimap2 preset: asm5, asm10, asm20 (default: asm10)"
            echo ""
            echo "Visualization Options (skip expensive operations):"
            echo "  --no-interactive                   Skip interactive HTML plots"
            echo "  --no-static                        Skip static PNG plots"
            echo "  --no-comparison                    Skip comparison mode visualizations"
            echo "  --no-gene-alignments               Skip gene alignment files (default: enabled)"
            echo ""
            echo "Caching & Control:"
            echo "  -f, --force                        Force reprocessing, ignore cache"
            echo "  -v, --verbose                      Print detailed progress information"
            echo "  -h, --help                         Show this help message"
            echo ""
            echo "Examples:"
            echo "  # Run all modes with all visualizations"
            echo "  $0 -m all"
            echo ""
            echo "  # Run all modes, skip expensive interactive plots"
            echo "  $0 -m all --no-interactive"
            echo ""
            echo "  # Run only scaffold mode"
            echo "  $0 -m scaffold"
            echo ""
            echo "  # Re-run with interactive plots enabled (uses cached results)"
            echo "  $0 -m all"
            echo ""
            echo "  # Force reprocess all samples"
            echo "  $0 --force"
            echo ""
            echo "Caching Behavior:"
            echo "  - First run with --no-interactive: skips interactive plots, saves time"
            echo "  - Second run without --no-interactive: uses cached analysis, regenerates plots"
            echo "  - Use --force to recompute everything regardless of cache"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use -h or --help for usage information"
            exit 1
            ;;
    esac
done

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Display mode information
case $MODE in
    contig)
        echo "Mode: CONTIG (Contig vs Reference)"
        ;;
    scaffold)
        echo "Mode: SCAFFOLD (Scaffold vs Reference)"
        ;;
    comparison)
        echo "Mode: COMPARISON (Scaffold vs Contig)"
        ;;
    all)
        echo "Mode: ALL (Contig + Scaffold + Comparison)"
        ;;
esac

# Display visualization skip flags
if [[ "$SKIP_INTERACTIVE" == true ]]; then
    echo "Skipping: Interactive HTML plots"
fi
if [[ "$SKIP_STATIC" == true ]]; then
    echo "Skipping: Static PNG plots"
fi
if [[ "$SKIP_COMPARISON" == true ]]; then
    echo "Skipping: Comparison visualizations"
fi
if [[ "$SKIP_GENE_ALIGNMENTS" == true ]]; then
    echo "Skipping: Gene alignment files"
fi

# Display force flag if requested
if [[ "$FORCE_REPROCESS" == true ]]; then
    echo "Force: Reprocessing all samples"
fi

echo ""

# Count samples
SAMPLE_COUNT=$(find "$INPUT_DIR" -maxdepth 1 -type d ! -path "$INPUT_DIR" | wc -l)
echo "Found $SAMPLE_COUNT samples to process"
echo ""

# Track statistics
COMPLETED=0
FAILED=0
START_TIME=$(date +%s)

# Loop through each subfolder in input_validation
for SAMPLE_INPUT in "$INPUT_DIR"/*/; do
    SAMPLE=$(basename "$SAMPLE_INPUT")
    SAMPLE_OUTPUT="$OUTPUT_DIR/$SAMPLE"

    echo "========================================"
    echo "[$((COMPLETED + FAILED + 1))/$SAMPLE_COUNT] Processing: $SAMPLE"
    echo "========================================"

    # Build the command for this sample using --input for auto-detection
    CMD="python genomeViz.py \
        --input \"$SAMPLE_INPUT\" \
        --output \"$SAMPLE_OUTPUT\" \
        --mode \"$MODE\" \
        --preset \"$PRESET\""

    # Add visualization skip flags
    if [[ "$SKIP_INTERACTIVE" == true ]]; then
        CMD="$CMD --no-interactive"
    fi
    if [[ "$SKIP_STATIC" == true ]]; then
        CMD="$CMD --no-static"
    fi
    if [[ "$SKIP_COMPARISON" == true ]]; then
        CMD="$CMD --no-comparison"
    fi
    if [[ "$SKIP_GENE_ALIGNMENTS" == true ]]; then
        CMD="$CMD --no-gene-alignments"
    fi

    # Add force flag if requested
    if [[ "$FORCE_REPROCESS" == true ]]; then
        CMD="$CMD --force"
    fi

    if [[ "$VERBOSE" == true ]]; then
        echo "Command: $CMD"
    fi

    # Execute the command
    eval "$CMD"

    if [ $? -ne 0 ]; then
        echo "❌ ERROR: Failed to process $SAMPLE"
        ((FAILED++))
    else
        echo "✓ SUCCESS: $SAMPLE completed"
        ((COMPLETED++))
    fi

    echo ""
done

# Calculate elapsed time
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
MINUTES=$((ELAPSED / 60))
SECONDS=$((ELAPSED % 60))

echo "========================================"
echo "Batch Processing Complete!"
echo "========================================"
echo "Completed: $COMPLETED/$SAMPLE_COUNT"
if [[ $FAILED -gt 0 ]]; then
    echo "Failed:    $FAILED/$SAMPLE_COUNT"
fi
printf "Time:      %02d:%02d\n" $MINUTES $SECONDS
echo "Output:    $OUTPUT_DIR"
echo "========================================"

if [[ $FAILED -gt 0 ]]; then
    exit 1
else
    exit 0
fi
