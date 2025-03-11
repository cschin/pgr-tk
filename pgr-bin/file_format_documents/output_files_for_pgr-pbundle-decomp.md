# Output Files for pgr-pbundle-decomp

The `pgr-pbundle-decomp` tool generates several output files when decomposing principal bundles from sequence data. Each file serves a specific purpose in the analysis pipeline. Below is a comprehensive description of each output file format:

## 1. Principal Bundle BED File (`[prefix].bed`)

This file contains the principal bundle decomposition results in BED format, which identifies regions in the genome that share similar sequence patterns. For more detailed information, see [Principal Bundle BED File Format](principal_bundle_bed_file.md).

### Format:
```
<contig_name> <start> <end> <bundle_id>:<bundle_size>:<direction>:<start_pos>:<end_pos>:<repeat_flag>
```

### Fields:
- `contig_name`: Name of the contig or chromosome
- `start`: 0-based start position of the bundle on the contig (inclusive)
- `end`: 0-based end position of the bundle on the contig (exclusive)
- `bundle_id`: Unique identifier for the bundle
- `bundle_size`: Number of sequences contained in the bundle
- `direction`: Direction of the bundle (0 for forward, 1 for reverse)
- `start_pos`: Start position within the principal bundle coordinate system
- `end_pos`: End position within the principal bundle coordinate system
- `repeat_flag`: Classification flag - 'R' for repeat regions, 'U' for unique regions

The first line of the file contains a comment with the command used to run the tool.

## 2. Contig Summary File (`[prefix].ctg.summary.tsv`)

This tab-separated file provides detailed statistical information about bundle distribution across each contig, useful for quantitative analysis. This file is often used in conjunction with the principal bundle BED file for comprehensive analysis.

### Header:
```
#ctg length repeat_bundle_count repeat_bundle_sum repeat_bundle_percentage repeat_bundle_mean repeat_bundle_min repeat_bundle_max non_repeat_bundle_count non_repeat_bundle_sum non_repeat_bundle_percentage non_repeat_bundle_mean non_repeat_bundle_min non_repeat_bundle_max total_bundle_count total_bundle_coverage_percentage
```


### Fields:
- `ctg`: Contig name
- `length`: Total length of the contig in base pairs
- `repeat_bundle_count`: Number of repeat bundles identified
- `repeat_bundle_sum`: Total base pairs covered by repeat bundles
- `repeat_bundle_percentage`: Percentage of contig covered by repeat bundles
- `repeat_bundle_mean`: Mean length of repeat bundles
- `repeat_bundle_min`: Minimum length of repeat bundles
- `repeat_bundle_max`: Maximum length of repeat bundles
- `non_repeat_bundle_count`: Number of non-repeat (unique) bundles
- `non_repeat_bundle_sum`: Total base pairs covered by non-repeat bundles
- `non_repeat_bundle_percentage`: Percentage of contig covered by non-repeat bundles
- `non_repeat_bundle_mean`: Mean length of non-repeat bundles
- `non_repeat_bundle_min`: Minimum length of non-repeat bundles
- `non_repeat_bundle_max`: Maximum length of non-repeat bundles
- `total_bundle_count`: Total number of bundles (repeat + non-repeat)
- `total_bundle_coverage_percentage`: Percentage of contig covered by all bundles

### Example:
```
#ctg length repeat_bundle_count repeat_bundle_sum repeat_bundle_percentage repeat_bundle_mean repeat_bundle_min repeat_bundle_max non_repeat_bundle_count non_repeat_bundle_sum non_repeat_bundle_percentage non_repeat_bundle_mean non_repeat_bundle_min non_repeat_bundle_max total_bundle_count total_bundle_coverage_percentage
chr1 248956422 2156 42568912 17.1 19744 500 125680 5842 192458700 77.3 32943 200 258462 7998 94.4
```

## 3. MAP Graph GFA File (.mapg.gfa)

This file contains the Minimizer Anchor Profile (MAP) graph in GFA (Graphical Fragment Assembly) format. For detailed information about the GFA format used in PGR-TK, see [GFA Output Format in PGR-TK](gfa_format.md).

## 4. MAP Graph Index File (.mapg.idx)

This file contains the index for the MAP graph, enabling efficient querying and traversal. The index format is described in detail in the [GFA Output Format in PGR-TK](gfa_format.md) document.

## 5. Principal MAP Graph GFA File (.pmapg.gfa)

This file contains the principal MAP graph in GFA format, which is a simplified version focusing on the principal bundles. See [GFA Output Format in PGR-TK](gfa_format.md) for format details.

## 6. Principal Bundle Data File (`[prefix].pdb`)

This binary file contains the complete principal bundle data, essential for downstream analysis tools in the PGR-TK suite. It includes:

- SHIMMER parameters used for minimizer generation (k-mer size, window size, reduction factor)
- Bundle information with coordinate mappings between reference and bundle spaces
- Fragment boundary information
- Bundle connectivity data
- Sequence mapping metadata

## Related Documentation

- [Principal Bundle BED File Format](principal_bundle_bed_file.md) - Detailed explanation of the principal bundle BED file format
- [Input File Formats for pgr-pbundle-bed2svg](input_file_formats_for_pgr-pbundle-bed.md) - Documentation for the input files required by pgr-pbundle-bed2svg
- [GFA Output Format in PGR-TK](gfa_format.md) - Information about the GFA format used in PGR-TK
- [Contig SV BED Format](ctgsv.bed.md) - Documentation for the contig SV BED format, another important file format in the PGR-TK ecosystem

