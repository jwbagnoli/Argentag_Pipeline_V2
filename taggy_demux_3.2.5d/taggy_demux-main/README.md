# ArgenTag single-cell read demultiplexing pipeline

Demultiplexing refers to the process of identifying the barcode(s) that each sequencing read from a given single-cell sequencing experiment is tagged with. Reads with common barcodes are assigned to the same cell. The customer-facing version of the ArgenTag pipeline is described here. An alternative version, used internally by the ArgenTag team to process internally- and externally-generated data, is described on a [separate page](doc/in-house.md).

In either case, the main output of the pipeline is a set of demultiplexed, trimmed reads. These can be fed to a downstream analysis pipeline, including the [Iso-Seq pipeline](#pacbio-data-and-iso-seq-downstream-analysis-pipeline) (recommended for PacBio reads) and the [FLAMES-based downstream analysis pipeline](#ont-data-and-flames-based-downstream-analysis-pipeline) (recommended for ONT reads). Downstream analysis is covered here only briefly, but users are encouraged to see the provided [examples](#Examples-and-downstream-analysis) and the documentation for their tool of choice for further details.

![Customer-facing demultiplexing pipeline](doc/img/customer-facing.png)

For the customer-facing pipeline, the entire pipeline (except for an optional [chimera splitting step](#chimera-splitting)) is consolidated into a single binary to make it more user friendly. The input is a file of basecalled reads in sam or fastq format, while the output is a set of demultiplexed, trimmed reads in one of the [supported output formats](#output-formats).

## Usage

    Usage: taggy_demux [OPTION...] input-file
    taggy_demux -- a demultiplexer for ArgenTag reads
    
      -2, --mate-file=FILE       Read mate reads from FILE.
      -c, --umi-pre=INT          Use INT as UMI pre context. [2]
      -C, --umi-post=INT         Use INT as UMI post context. [2]
      -d, --darwin               Output entity list for plotting with the darwin
                                 analysis tool.
      -D, --max-edit-d=INT/FLOAT Maximum edit distance to consider for linker
                                 alignment (-1 means no limit). If float and < 1,
                                 intepreted as a relative maximum edit dist. If
                                 float and > 1, rounded down. [3]
      -f, --in-fmt=STRING        Format for input read file. Valid values are
                                 "fastq", "sam". [fastq]
      -F, --out-fmt=STRING       Format for output read file. Valid values are
                                 "flames", "fastq", "scnano", "sam". [flames]
      -h, --keep-header          Preserve SAM header in output. Only meaningful if
                                 --in-fmt=sam and --out-fmt=sam. [FALSE]
      -k, --keep-failed          Keep failed reads and set nb tag to -1. Only
                                 meaningful if --out-fmt=sam. [FALSE]
      -m, --min-len=INT          Minimum length to output. [1]
      -M, --max-len=INT          Maximum length to output. -1 means no limit. [-1]
      -o, --output-dir=DIR       Output directory. [Current directory]
      -O, --orient=STRING        Orientation in which to output reads (one of
                                 "sense", "anti", "preserve" or "invert". [sense]
      -p, --trim-poly=STRING     Trim polyA/T from output sequences. Valid values
                                 are "none", "strict", "normal", "lenient". [none]
      -P, --preserve             Preserve tags in original sam record. Only
                                 meaningful if --in-fmt=sam and --out-fmt=sam.
                                 [FALSE]
      -r, --require-TSO          Require TSO to call a read as valid. [FALSE]
      -R, --max-r-bases=INT      Maximum number of bases from read to align (-1
                                 means no limit) [-1]
      -s, --split-chims          Split chimeric reads. [FALSE]
      -S, --sample=INT           Process only the first INT reads (-1 means no
                                 limit) [-1]
      -t, --trim-TSO             Trim TSO (if found) from output sequences. [FALSE]
                                
      -T, --num-threads=INT      Use INT parallel threads [1]
      -u, --umi-start=INT        Use INT as UMI start coordinate. [24]
      -U, --umi-end=INT          Use INT as UMI end coordinate. [39]
      -w, --whitelist=FILE       Use custom barcode whitelist file.
      -x, --presets=STRING       Presets for sequencing technology (hifi, ont,
                                 illu)and bead design (v1, v2). Can be combined via
                                 "+" (e.g. --presets=ont+v1).
      -?, --help                 Give this help list
          --usage                Give a short usage message
    
    Mandatory or optional arguments to long options are also mandatory or optional
    for any corresponding short options.

## Presets

For convenience, common combinations of flags are grouped into presets, which can be passed via the `-x` flag (or its long form, `--presets`). The following presets are available:

Bead design:

    -x v1 / --presets=v1: equivalent to -u 24 -U 39 -c 2 -C 2
    -x v2 / --presets=v2: equivalent to -u 24 -U 41 -c 4 -C 2

Sequencing platform:

    -x hifi / --presets=hifi: equivalent to --in-fmt=sam --out-fmt=sam --max-edit-d=0.1 --orient sense --trim-poly normal --trim-TSO --keep-header --preserve 
    -x ont / --presets=ont: equivalent to --in-fmt=fastq --out-fmt=flames --max-edit-d=0.15 --orient sense --trim-poly lenient --trim-TSO --split-chims
    -x illu / --presets=illu: equivalent to --in-fmt=paired --out-fmt=paired --max-edit-d=0.05 --orient sense --trim-poly strict --trim-TSO

Multiple presets can be combined via the `+` character (e.g. `--presets=ont+v2`), and individual settings can be overriden by flags appearing *after* the preset (e.g. `--presets=v1+hifi --out-fmt=fastq` will use general settings for PacBio Hi-Fi reads obtained with v1 beads, except the output will be in fastq format rather than sam (which is the default for this preset).

## Output formats

### PacBio-style sam format (`--out-fmt=sam`)
This format follows the [SAM format specification](http://samtools.github.io/hts-specs/SAMv1.pdf) maintained by the SAM/BAM Format Specification Working Group, making use of the optional tags to encode additional information relevant to barcode demultiplexing and single-cell analysis. Consistency with the [PacBio BAM format specification](https://pacbiofileformats.readthedocs.io/en/13.0/BAM.html) is maintained. In particular, the following tags are added:

| Tag	| Data type	| Description				|
| ----- | -------------	| -----------				|
| XA    | Z             | Order of tags names. Set to `"XM-CB"`. |
| gp    | i             | Specifies whether or not the barcode for the given read passes. Set to `1` for passing reads. |
| CB	| Z		| Corrected cell barcode. |
| CR	| Z		| Raw (uncorrected) cell barcode. Set equal to CB. |
| XC	| Z		| Raw cell barcode. Set equal to XC. |
| rc	| i		| Predicted real cell. Set to `1` for *all* reads because `taggy_demux` does *not* perform real cell identification (also known as background RNA filtering or elbow plot analysis). See the section on [Updating of the `rc` tag](#Updating-of-the-rc-tag) for details on how to achieve this.|
| nc	| i		| Number of candidate barcodes. Set to `1`. |
| nb	| i		| Edit distance from the barcode for the read to the barcode to which it was reassigned. Set to `0`. |
| XM    | Z             | Raw (after tag) or corrected (after correct) UMI. |

The above does not preclude the presence of other tags previously added by third-party tools, which will be preserved if the `--preserve` flag is set (only meaningful if `--in-fmt=sam` and `--out-fmt=sam`).

An example SAM entry (line) is shown below:

    molecule/0      4       *       0       255     *       *       0       0       GGCAYTCATG[...]CGATGGCTAG *       CB:Z:AACCAAGGAGGTAGAT   XA:Z:XM-CB      XM:Z:CGCGACTGTTCT       ic:i:1  im:Z:m84112_240530_215351_s2/139986042/ccs/40_4082      is:i:1  it:Z:CGCGACTGTTCTAACCAAGGAGGTAGAT       rc:i:1  RG:Z:e4927d21   zm:i:0

<!--

### PacBio-style bam format (`--out-fmt=bam`)
This is the binary version of the above [PacBio-style sam format](#pacbio-style-sam-format---out-fmtsam), and should be equivalent to using `--out-fmt=sam` followed by sam-to-bam conversion with a third-party tool.

-->

### FLAMES-style fastq format (`--out-fmt=flames`, default)
This is like the standard fastq format, except that read headers follow the following structure:

    @XXXX-YYYY-ZZZZ_UUUUUUUUUUUU#READID

* `XXXX`, `YYYY` and `ZZZZ` are the 3 barcodes which identify a specific cell.
* `UUUUUUUUUUUU` is the 12-nt UMI.
* `READID` is the original read ID.

An example could be
 
    @0076-0048-0089_ATACCGGCTACA#VH00444:319:AAFV5MHM5:1:1101:18421:23605

which would correspond to sequencing read VH00444:319:AAFV5MHM5:1:1101:18421:23605, which has been tagged with the barcode triplet (0076, 0048, 0089) and the UMI "ATACCGGCTACA".

### Fastq format with mapped barcode (`--out-fmt=fastq`)
This uses a standard fastq format, except that read headers follow the following structure:

    @BBBBBBBBBBBBBBBB_UUUUUUUUUUUU#READID

* `BBBBBBBBBBBBBBBB` is a unique 16-nt nucleotide pseudobarcode which identifies a specific cell. This results from a mapping from ArgenTag's barcode triplets to a single 16-nt sequence, which is artificial and will not appear verbatim in the original basecalled sequence.
* `UUUUUUUUUUUU` is the 12-nt UMI.
* `READID` is the original read ID.

An example could be
 
    @ACGTGCAGCAGACGGT_ATACCGGCTACA#VH00444:319:AAFV5MHM5:1:1101:18421:23605

which would correspond to sequencing read VH00444:319:AAFV5MHM5:1:1101:18421:23605, which has been tagged with the cell barcode "ACGTGCAGCAGACGGT" and the UMI "ATACCGGCTACA".

### scNanoGPS-style fastq format (`--out-fmt=scnano`)
This is like the standard fastq format, except that read headers follow the following structure:

    @READID_UUUUUUUUUUUU

* `READID` is the original read ID.
* `UUUUUUUUUUUU` is the 12-nt UMI.

An example could be
 
    @VH00444:319:AAFV5MHM5:1:1101:18421:23605_ATACCGGCTACA

which would correspond to sequencing read VH00444:319:AAFV5MHM5:1:1101:18421:23605, which has been tagged with the UMI "ATACCGGCTACA".

For this format, the barcode is not included in the content of the fastq file, but is instead provided in the file name (one file per barcode combination/cell).

# Chimera splitting

Chimera splitting refers to an optional step run before demultiplexing, whereby common chimeric reads are split into two or more subreads, which are then processed normally. This was previously done with a standalone "split.sh" script, which has since been merged to the main program as an optional flag `-s` (or its long form `--split-chims`). For split reads, each resulting subread keeps the original read ID, followed by a dash and the subread number (e.g. if read ID `@VH00444:319:AAFV5MHM5:1:1101:18421:23605` is split into two, this will result in two subreads `@VH00444:319:AAFV5MHM5:1:1101:18421-23605-1` and `@VH00444:319:AAFV5MHM5:1:1101:18421:23605-2`).

# Examples and downstream analysis

## PacBio data and Iso-Seq downstream analysis pipeline

For data generated on the PacBio platform, we recommend using the SAM input format (`--in-fmt=sam`), as well as the PacBio-compatible SAM format (`--out-fmt=sam`). Conversion from SAM to BAM and viceversa can be handled via samtools or similar (ideally in place, via process substitution), as detailed below. The input file is typically a segmented read file obtained from Skera (`segmented.bam` below), but can also be a regular CCS file.

### Example commands
    # Convert input bam file to sam
	samtools view -h segmented.bam > segmented.sam
    # Run demux binary with 32 threads from sam input
    NUM_THREADS=32
	mkdir -p "$OUT_DIR"
    bin/taggy_demux -T "$NUM_THREADS" -o "$OUT_DIR" --presets=hifi segmented.sam
    # The above command is equivalent to
    # bin/taggy_demux -T "$NUM_THREADS" -o "$OUT_DIR" --orient sense --in-fmt=sam \
    # --out-fmt=sam --max-edit-d=0.1 --orient sense --trim-poly normal --trim-TSO \
    # --keep-header --preserve segmented.sam
    #
    # Convert demultiplexed sam output back to bam
    samtools view -bS "$OUT_DIR"/demux.sam > "$OUT_DIR"/demux.bam
	
### Example commands (with in place bam-to-sam conversion)
    # Convert S-read bam file to sam in place via process substitution, and run demux binary with 32 threads
    NUM_THREADS=32
    bin/taggy_demux -T "$NUM_THREADS" -o "$OUT_DIR" --presets=hifi <(samtools view -h segmented.bam)
    # Convert demultiplexed sam output back to bam
    samtools view -bS "$OUT_DIR"/demux.sam > "$OUT_DIR"/demux.bam

### Updating of the `rc` tag

The demultiplexed read file generated above (`demux.sam/bam`) has the `rc` tag set to `1` for *all* reads because `taggy_demux` does *not* perform real cell identification (also known as background RNA filtering or elbow plot analysis). Users can leverage `isoseq bcstats` to perform cell calling based on their desired criterion and follow [this auxiliary guide](doc/update_rc.md) to update the `rc` tag as appropriate using the provided [convenience script](bin/update_rc.sh).

![RC tag update workflow](doc/img/rc-update-workflow.png)

### Iso-Seq downstream analysis pipeline

The output from the previous commands is compatible with the [Iso-Seq pipeline, starting from the deduplication step (Step 6)](https://isoseq.how/umi/cli-workflow.html#step-6---deduplication).

## ONT data and FLAMES-based downstream analysis pipeline

For data generated on the ONT platform, we recommend running the optional [chimera splitting step](#chimera-splitting) and outputting in the FLAMES-compatible fastq format (`--out-fmt=flames`) for direct compatibility with our [FLAMES-based downstream analysis pipeline](#FLAMES-based-downstream-analysis-pipeline).

### Example commands for a single fastq file

    # Run with 32 threads and optional chimera splitting step (see "Chimera splitting" above)
	NUM_THREADS=32
    mkdir -p "$OUT_DIR"
    bin/taggy_demux -T "$NUM_THREADS" -o "$OUT_DIR" --presets=ont "$INPUT_FASTQ_FILE"
    # The above command is equivalent to
    # bin/taggy_demux -T "$NUM_THREADS" -o "$OUT_DIR" --in-fmt=fastq \
    # --out-fmt=flames --max-edit-d=0.15 --orient sense \
    # --trim-poly lenient --trim-TSO --split-chims "$INPUT_FASTQ_FILE"

### FLAMES-based downstream analysis pipeline

The output from the previous commands can then be analyzed with our [FLAMES-based downstream analysis pipeline](https://github.com/argentagsw/at_flames).

![FLAMES-based downstream analysis pipeline overview](doc/img/FLAMES-based.png)
