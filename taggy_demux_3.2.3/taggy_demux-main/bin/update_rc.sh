#!/bin/bash

# PacBio Isoseq: Apply Cell Calling to Corrected BAM
# Applies cell calling results (bcstats) to BAM file with corrected barcodes
# Updates rc tags based on cell/non-cell classifications from Isoseq
# 
# Author: Bioinformatics Pipeline
# Version: 2.0

set -e

# Default values
INPUT_BAM=""
BCSTATS_TSV=""
OUTPUT_BAM=""
THREADS=4
KEEP_TEMP=false
VERBOSE=false

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to display help
show_help() {
    cat << EOF
${BLUE}PacBio Isoseq: Apply Cell Calling to AT Corrected BAM${NC}

${YELLOW}DESCRIPTION:${NC}
    Applies cell calling results (bcstats) to BAM files with AT corrected barcodes.
    Updates rc tags based on cell/non-cell classifications from Isoseq cell calling.

${YELLOW}USAGE:${NC}
    $0 -i <input.bam> -b <bcstats.tsv> -o <output.bam> [OPTIONS]

${YELLOW}REQUIRED ARGUMENTS:${NC}
    -i, --input     AT corrected BAM file with CB tags (from AT barcode correction)
    -b, --bcstats   bcstats TSV file from Isoseq cell calling
    -o, --output    Output BAM file with applied cell calling (rc tags updated)

${YELLOW}OPTIONAL ARGUMENTS:${NC}
    -t, --threads   Number of threads for samtools operations (default: 4)
    -k, --keep-temp Keep temporary files for debugging (default: false)
    -v, --verbose   Enable verbose output (default: false)
    -h, --help      Show this help message

${YELLOW}EXAMPLES:${NC}
    # Apply cell calling to AT corrected BAM
    $0 -i AT_corrected.bam -b bcstats_report.tsv -o final_with_cells.bam

    # With custom threads and verbose output
    $0 -i AT_corrected.bam -b bcstats.tsv -o final.bam -t 8 -v

    # Keep temporary files for debugging
    $0 -i AT_corrected.bam -b bcstats.tsv -o final.bam -k

${YELLOW}INPUT FILE REQUIREMENTS:${NC}
    - AT Corrected BAM: Must contain CB tags (corrected barcodes from AT correction)
    - bcstats TSV: Must be output from 'isoseq tag' with columns:
      Column 1: Barcode sequence
      Column 5: Classification (cell/non-cell)

${YELLOW}OUTPUT:${NC}
    - Updated BAM file with rc tags (rc:i:1 for cells, rc:i:0 for non-cells)
    - BAM index file (.bai)
    - Summary statistics

${YELLOW}TYPICAL WORKFLOW:${NC}
    1. Run AT barcode correction: isoseq correct --method AT
    2. Run cell calling: isoseq tag --design T-12U-16B  
    3. Apply cell calling to AT corrected BAM: $0 -i AT_corrected.bam -b bcstats.tsv -o final.bam
    4. Sort by barcode: samtools sort -t CB final.bam -o final.sorted.bam
    5. Deduplicate: isoseq groupdedup final.sorted.bam dedup.bam

EOF
}

# Function to log messages
log_message() {
    local level=$1
    shift
    local message="$@"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    
    case $level in
        "INFO")
            echo -e "${GREEN}[INFO]${NC} ${timestamp} - $message" >&2
            ;;
        "WARN")
            echo -e "${YELLOW}[WARN]${NC} ${timestamp} - $message" >&2
            ;;
        "ERROR")
            echo -e "${RED}[ERROR]${NC} ${timestamp} - $message" >&2
            ;;
        "DEBUG")
            if [[ "$VERBOSE" == "true" ]]; then
                echo -e "${BLUE}[DEBUG]${NC} ${timestamp} - $message" >&2
            fi
            ;;
    esac
}

# Function to check if command exists
check_command() {
    if ! command -v "$1" &> /dev/null; then
        log_message "ERROR" "Required command '$1' not found. Please install it."
        exit 1
    fi
}

# Function to validate input files
validate_inputs() {
    log_message "INFO" "Validating input files..."
    
    # Check AT corrected BAM
    if [[ ! -f "$INPUT_BAM" ]]; then
        log_message "ERROR" "AT corrected BAM file not found: $INPUT_BAM"
        exit 1
    fi
    
    # Check bcstats TSV
    if [[ ! -f "$BCSTATS_TSV" ]]; then
        log_message "ERROR" "bcstats TSV file not found: $BCSTATS_TSV"
        exit 1
    fi
    
    # Check if output directory exists
    OUTPUT_DIR=$(dirname "$OUTPUT_BAM")
    if [[ ! -d "$OUTPUT_DIR" ]]; then
        log_message "WARN" "Output directory doesn't exist, creating: $OUTPUT_DIR"
        mkdir -p "$OUTPUT_DIR"
    fi
    
    # Validate bcstats format
    local header_line=$(head -n 1 "$BCSTATS_TSV")
    local num_cols=$(echo "$header_line" | tr '\t' '\n' | wc -l)
    if [[ $num_cols -lt 5 ]]; then
        log_message "ERROR" "bcstats TSV must have at least 5 columns. Found: $num_cols"
        exit 1
    fi
    
    log_message "INFO" "Input validation passed"
}

# Function to extract barcode classifications
extract_barcodes() {
    log_message "INFO" "Extracting cell classifications from bcstats..."
    
    # Create temporary files
    local temp_dir=$(mktemp -d)
    CELLS_FILE="$temp_dir/cells.txt"
    NONCELLS_FILE="$temp_dir/noncells.txt"
    
    # Extract cell and non-cell barcodes (skip header)
    awk 'NR>1 && $5=="cell" {print $1}' "$BCSTATS_TSV" > "$CELLS_FILE"
    awk 'NR>1 && $5=="non-cell" {print $1}' "$BCSTATS_TSV" > "$NONCELLS_FILE"
    
    local cell_count=$(wc -l < "$CELLS_FILE")
    local noncell_count=$(wc -l < "$NONCELLS_FILE")
    
    log_message "INFO" "Found $cell_count cell barcodes"
    log_message "INFO" "Found $noncell_count non-cell barcodes"
    
    if [[ $cell_count -eq 0 && $noncell_count -eq 0 ]]; then
        log_message "ERROR" "No barcodes found in bcstats file. Check file format."
        exit 1
    fi
    
    echo "$temp_dir"
}

# Function to process BAM file
process_bam() {
    local temp_dir=$1
    
    log_message "INFO" "Processing AT corrected BAM file: $INPUT_BAM"
    log_message "DEBUG" "Using $THREADS threads for samtools operations"
    
    # Convert BAM to SAM, update rc tags, convert back to BAM
    samtools view -@ "$THREADS" -h "$INPUT_BAM" | python3 -c "
import sys
import os

# Load barcode classifications
cells_file = os.path.join('$temp_dir', 'cells.txt')
noncells_file = os.path.join('$temp_dir', 'noncells.txt')

print('Loading barcode classifications...', file=sys.stderr)
cell_barcodes = set()
noncell_barcodes = set()

# Load cell barcodes
try:
    with open(cells_file, 'r') as f:
        cell_barcodes = set(line.strip() for line in f if line.strip())
except FileNotFoundError:
    print(f'Warning: {cells_file} not found', file=sys.stderr)

# Load non-cell barcodes  
try:
    with open(noncells_file, 'r') as f:
        noncell_barcodes = set(line.strip() for line in f if line.strip())
except FileNotFoundError:
    print(f'Warning: {noncells_file} not found', file=sys.stderr)

print(f'Loaded {len(cell_barcodes)} cell barcodes and {len(noncell_barcodes)} non-cell barcodes', file=sys.stderr)

# Process SAM input
reads_processed = 0
rc_1_count = 0
rc_0_count = 0
cb_missing = 0
cb_unclassified = 0

for line in sys.stdin:
    line = line.rstrip()
    
    # Pass through headers unchanged
    if line.startswith('@'):
        print(line)
        continue
    
    # Process read lines
    fields = line.split('\t')
    
    # Find CB tag (corrected barcode)
    cb_value = None
    for field in fields[11:]:  # Tags start at field 12 (index 11)
        if field.startswith('CB:Z:'):
            cb_value = field[5:]  # Remove 'CB:Z:' prefix
            break
    
    # Determine new rc value
    new_rc = None
    if cb_value is None:
        new_rc = 'rc:i:0'  # No CB tag = non-cell
        cb_missing += 1
        rc_0_count += 1
    elif cb_value in cell_barcodes:
        new_rc = 'rc:i:1'
        rc_1_count += 1
    elif cb_value in noncell_barcodes:
        new_rc = 'rc:i:0'
        rc_0_count += 1
    else:
        new_rc = 'rc:i:0'  # Unclassified = non-cell
        cb_unclassified += 1
        rc_0_count += 1
    
    # Update existing rc tag or add new one
    rc_found = False
    for i, field in enumerate(fields):
        if field.startswith('rc:i:'):
            fields[i] = new_rc
            rc_found = True
            break
    
    # Add rc tag if not found
    if not rc_found:
        fields.append(new_rc)
    
    # Output updated line
    print('\t'.join(fields))
    
    reads_processed += 1
    if reads_processed % 1000000 == 0:
        print(f'Processed {reads_processed:,} reads...', file=sys.stderr)

print(f'Processing complete. {reads_processed:,} reads processed.', file=sys.stderr)
print(f'RC tags assigned: {rc_1_count:,} cells, {rc_0_count:,} non-cells', file=sys.stderr)
print(f'Missing CB tags: {cb_missing:,}, Unclassified CB: {cb_unclassified:,}', file=sys.stderr)
" | samtools view -@ "$THREADS" -bS - > "$OUTPUT_BAM"
    
    log_message "INFO" "BAM processing completed"
}

# Function to index and verify output
finalize_output() {
    log_message "INFO" "Indexing output BAM..."
    samtools index -@ "$THREADS" "$OUTPUT_BAM"
    
    log_message "INFO" "Verifying results..."
    
    # Count rc tags in output
    local rc_counts=$(samtools view "$OUTPUT_BAM" | grep -o 'rc:i:[01]' | sort | uniq -c || true)
    
    echo ""
    log_message "INFO" "RC tag distribution in output BAM:"
    echo "$rc_counts"
    
    # Get expected counts from bcstats
    local expected_cells=$(awk 'NR>1 && $5=="cell" {count++} END {print count+0}' "$BCSTATS_TSV")
    local expected_noncells=$(awk 'NR>1 && $5=="non-cell" {count++} END {print count+0}' "$BCSTATS_TSV")
    
    echo ""
    log_message "INFO" "Expected from bcstats:"
    echo "  $expected_cells rc:i:1 (cells)"
    echo "  $expected_noncells rc:i:0 (non-cells)"
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            INPUT_BAM="$2"
            shift 2
            ;;
        -b|--bcstats)
            BCSTATS_TSV="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_BAM="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -k|--keep-temp)
            KEEP_TEMP=true
            shift
            ;;
        -v|--verbose)
            VERBOSE=true
            shift
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            log_message "ERROR" "Unknown option: $1"
            show_help
            exit 1
            ;;
    esac
done

# Check required arguments
if [[ -z "$INPUT_BAM" || -z "$BCSTATS_TSV" || -z "$OUTPUT_BAM" ]]; then
    log_message "ERROR" "Missing required arguments"
    show_help
    exit 1
fi

# Main execution
echo -e "${BLUE}=== Apply Cell Calling to AT Corrected BAM ===${NC}"
log_message "INFO" "Input BAM: $INPUT_BAM"
log_message "INFO" "bcstats TSV: $BCSTATS_TSV"
log_message "INFO" "Output BAM: $OUTPUT_BAM"
log_message "INFO" "Threads: $THREADS"
echo ""

# Check required tools
log_message "INFO" "Checking required tools..."
check_command "samtools"
check_command "python3"
check_command "awk"

# Validate inputs
validate_inputs

# Extract barcodes
temp_dir=$(extract_barcodes)

# Process BAM file
process_bam "$temp_dir"

# Finalize output
finalize_output

# Cleanup
if [[ "$KEEP_TEMP" == "false" ]]; then
    log_message "DEBUG" "Cleaning up temporary files..."
    rm -rf "$temp_dir"
else
    log_message "INFO" "Temporary files kept in: $temp_dir"
fi

echo ""
echo -e "${GREEN}=== SUCCESS! ===${NC}"
log_message "INFO" "AT Corrected BAM: $INPUT_BAM"
log_message "INFO" "Output BAM: $OUTPUT_BAM (with applied cell calling)"
log_message "INFO" "Index:  $OUTPUT_BAM.bai"
echo ""
echo -e "${YELLOW}Next steps:${NC}"
echo "  samtools sort -t CB $OUTPUT_BAM -o ${OUTPUT_BAM%.bam}.sorted_by_cb.bam"
echo "  isoseq groupdedup ${OUTPUT_BAM%.bam}.sorted_by_cb.bam ${OUTPUT_BAM%.bam}.dedup.bam"
echo ""

