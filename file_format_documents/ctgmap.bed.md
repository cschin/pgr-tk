# *.ctgmap.bed File Format

The *.ctgmap.bed file is a tab-separated values (TSV) file that describes the alignment of query contigs to a reference genome. It follows a modified BED format with additional fields to provide detailed information about the alignments.

## File Structure

Each line in the file represents a single alignment between a query contig and the reference genome. The fields are separated by tabs and are ordered as follows:

1. Reference sequence name
2. Start position on the reference (0-based)
3. End position on the reference (exclusive)
4. Additional information (colon-separated)

## Fields

1. **Reference sequence name**: The name of the reference sequence (chromosome or scaffold).

2. **Start position**: The start position of the alignment on the reference sequence (0-based).

3. **End position**: The end position of the alignment on the reference sequence (exclusive).

4. **Additional information**: A colon-separated string containing the following fields:
   a. Query sequence name
   b. Start position on the query sequence
   c. End position on the query sequence
   d. Query contig length
   e. Orientation (0 for forward, 1 for reverse)
   f. Contig orientation (0 for forward, 1 for reverse)
   g. Target duplication flag (0 for unique, 1 for duplicated)
   h. Target overlap flag (0 for non-overlapping, 1 for overlapping)
   i. Query duplication flag (0 for unique, 1 for duplicated)
   j. Query overlap flag (0 for non-overlapping, 1 for overlapping)

## Example

```
chr1    1000    2000    contig1:500:1500:3000:0:0:0:0:0:0
```

This example line can be interpreted as follows:
- The alignment is on reference sequence "chr1" from position 1000 to 2000.
- The query contig name is "contig1".
- The alignment covers positions 500 to 1500 on the query contig.
- The total length of the query contig is 3000 base pairs.
- The alignment is in the forward orientation (0) for both the reference and the query.
- The alignment is unique and non-overlapping on both the reference and the query (all flags are 0).

## Usage

This file format is useful for visualizing and analyzing the alignment of query contigs to a reference genome. It can be used to identify structural variations, assess the quality of genome assemblies, and compare different genome versions or assemblies.