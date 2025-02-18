# Input File Formats for pgr-pbundle-bed2svg

## 1. Principal Bundle BED File (Required)

This is the main input file containing bundle information.

- Format: Tab-separated values
- Each line represents a bundle segment
- Fields:
  1. Contig name
  2. Start position
  3. End position
  4. Bundle information (colon-separated):
     - Bundle ID
     - Bundle vertex count (not used)
     - Bundle direction (0 or 1)
     - Bundle vertex start (not used)
     - Bundle vertex end (not used)

Example: 
```
ctg1 1000 2000 1:5:0:1000:2000 ctg1 3000 4000 2:3:1:3000:4000
```

## 2. Annotation File (Optional)

Provides additional annotation text for each contig.

- Format: Tab-separated values
- Fields:
  1. Contig name
  2. Annotation text

Example: 
```
ctg1  Annotation for contig 1 
ctg2  Annotation for contig 2
```

## 3. Annotation Region BED File (Optional)

Defines regions for annotation tracks.

- Format: Tab-separated values
- Fields:
  1. Contig name
  2. Start position
  3. End position
  4. Title
  5. Color

Example:
```
ctg1 5000 6000 Region_A #FF0000 
ctg1 7000 8000 Region_B #00FF00
```


## 4. Offset File (Optional)

Provides offset values for each contig.

- Format: Tab-separated values
- Fields:
  1. Contig name
  2. Offset value (integer)

Example:
```
ctg1 1000 ctg2 -500
```


## 5. Dendrogram File (Optional)

Describes the hierarchical clustering of contigs.

- Format: Tab-separated values
- Three types of lines:
  1. Leaf nodes (L):
     - L    [node_id]    [contig_name]
  2. Internal nodes (I):
     - I    [node_id]    [child_node0]    [child_node1]    [node_size]    [node_height]
  3. Node positions (P):
     - P    [node_id]    [node_position]    [node_height]    [node_size]

Example:
```
L 1 ctg1 L 2 ctg2 I 3 1 2 2 0.5 P 1 0.0 0.0 1 P 2 1.0 0.0 1 P 3 0.5 0.5 2
```
