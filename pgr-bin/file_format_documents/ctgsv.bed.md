# ctgsv.bed File Format

The "ctgsv.bed" file is a tab-separated values (TSV) file that contains information about query contig alignments and potential structural variations. Each line in the file represents a region of interest in the query sequence.

## File Structure

Each line in the file consists of four columns:

1. Query Name
2. Start Position
3. End Position
4. Annotation

## Column Descriptions

1. **Query Name**: The name or identifier of the query contig.

2. **Start Position**: The starting position of the region in the query contig (0-based).

3. **End Position**: The ending position of the region in the query contig (exclusive).

4. **Annotation**: A string containing information about the region, with fields separated by colons. The annotation format is as follows:

   `<Type>:<Target>:<Target_Start>-<Target_End>:<Orientation>:<Ctg_Orientation>:<Additional_Info>`

   - `<Type>`: Indicates the type of region:
     - `QG`: Query Gap
     - `QD`: Query Duplicate
     - `QO`: Query Overlap

   - `<Target>`: The name of the target sequence this region aligns to.

   - `<Target_Start>-<Target_End>`: The start and end positions in the target sequence.

   - `<Orientation>`: The orientation of the alignment (0 for forward, 1 for reverse).

   - `<Ctg_Orientation>`: The overall orientation of the contig (0 for forward, 1 for reverse).

   - `<Additional_Info>`: Any additional information (may vary depending on the type).

## Example

```
contig1    0    1000    QG:BGN>chr1:500-1500:0:0:32000
contig1    1000 2000    QD:chr1>chr2:1000-2000:1:0:32000
contig1    2000 3000    QO:chr2>chr3:2000-3000:0:0:32000
contig1    3000 32000   QG:chr3>END
```

In this example:
- The first line shows a query gap at the beginning of contig1, aligning to chr1.
- The second line indicates a duplicated region in contig1, aligning to chr2 in reverse orientation.
- The third line shows an overlapping region in contig1, aligning to chr3.
- The last line represents the end of the contig, with a gap from the last alignment to the end of the sequence.

This file format provides a comprehensive view of how query contigs align to the target sequences, highlighting potential structural variations, duplications, and gaps in the query assembly.