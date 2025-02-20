Principle Bundle bed File Format (Generated by `pgr-pbundle-decomp`)

The *.bed file is a tab-separated text file that contains information about principal bundle decomposition. Each line represents a bundle region in a contig.

File Header:
The file begins with a comment line containing the command used to generate the file:
# cmd: <command_string>

Data Columns:
Each subsequent line contains the following tab-separated fields:

1. Contig Name (string)
   The name of the contig or sequence.

2. Start Position (integer)
   The start position of the bundle region on the contig (0-based).

3. End Position (integer)
   The end position of the bundle region on the contig.

4. Bundle Information (string)
   A colon-separated string containing the following information:
   a. Bundle ID (integer)
   b. Bundle Size (integer)
   c. Direction (0 or 1)
      0 indicates forward direction, 1 indicates reverse direction.
   d. Start Position in Bundle (integer)
   e. End Position in Bundle (integer)
   f. Repeat Status (R or U)
      R indicates a repeat region, U indicates a unique (non-repeat) region.

Example:
contig1  1000  2000  42:5000:0:100:200:U

This line represents a bundle region on contig1, starting at position 1000 and ending at position 2000. The bundle ID is 42, with a size of 5000. It's in the forward direction (0), starts at position 100 and ends at position 200 within the bundle, and is a unique (non-repeat) region.