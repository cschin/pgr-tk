# Output Files for pgr-pbundle-decomp

The `pgr-pbundle-decomp` tool generates several output files. Here's a description of their formats:

## 1. BED File (.bed)

This file contains the principal bundle decomposition results in BED format.

### Format:
```
<contig_name> <start> <end> <bundle_id>:<bundle_size>::<start_pos>:<end_pos>:<repeat_flag>
```

### Fields:
- `contig_name`: Name of the contig
- `start`: Start position of the bundle on the contig
- `end`: End position of the bundle on the contig
- `bundle_id`: Unique identifier for the bundle
- `bundle_size`: Size of the bundle
- `direction`: Direction of the bundle (0 or 1)
- `start_pos`: Start position within the bundle
- `end_pos`: End position within the bundle
- `repeat_flag`: 'R' for repeat, 'U' for unique

The first line of the file contains the command used to run the tool.

## 2. Contig Summary File (.ctg.summary.tsv)

This tab-separated file provides summary statistics for each contig.

### Header:
```
#ctg length repeat_bundle_count repeat_bundle_sum repeat_bundle_percentage repeat_bundle_mean repeat_bundle_min repeat_bundle_max non_repeat_bundle_count non_repeat_bundle_sum non_repeat_bundle_percentage non_repeat_bundle_mean non_repeat_bundle_min non_repeat_bundle_max total_bundle_count total_bundle_coverage_percentage
```


### Fields:
- `ctg`: Contig name
- `length`: Length of the contig
- `repeat_bundle_count`: Number of repeat bundles
- `repeat_bundle_sum`: Total length of repeat bundles
- `repeat_bundle_percentage`: Percentage of contig covered by repeat bundles
- `repeat_bundle_mean`: Mean length of repeat bundles
- `repeat_bundle_min`: Minimum length of repeat bundles
- `repeat_bundle_max`: Maximum length of repeat bundles
- `non_repeat_bundle_count`: Number of non-repeat bundles
- `non_repeat_bundle_sum`: Total length of non-repeat bundles
- `non_repeat_bundle_percentage`: Percentage of contig covered by non-repeat bundles
- `non_repeat_bundle_mean`: Mean length of non-repeat bundles
- `non_repeat_bundle_min`: Minimum length of non-repeat bundles
- `non_repeat_bundle_max`: Maximum length of non-repeat bundles
- `total_bundle_count`: Total number of bundles
- `total_bundle_coverage_percentage`: Percentage of contig covered by all bundles

## 3. MAP Graph GFA File (.mapg.gfa)

This file contains the MAP graph in GFA format.

## 4. MAP Graph Index File (.mapg.idx)

This file contains the index for the MAP graph.

## 5. Principal MAP Graph GFA File (.pmapg.gfa)

This file contains the principal MAP graph in GFA format.

## 6. Principal Bundle Data File (.pdb)

This binary file contains the principal bundle data, including SHIMMER parameters and bundle information.

