# PGR-tk: A PanGenomic Research Tool Kit

[![test_and_build](https://github.com/cschin/pgr-tk/actions/workflows/test_and_build.yml/badge.svg)](https://github.com/cschin/pgr-tk/actions/workflows/test_and_build.yml)

This repository is a project to provide Python and Rust libraries to facilitate pangenomics analysis. Several algorithms and data structures used for the Peregrine Genome Assembler are useful for Pangenomics analysis as well. This repo takes those algorithms and data structures, combining other handy 3rd party tools to expose them as a library in Python (with Rust code for those computing parts that need performance.)

## What is PGR-tk?

Research Preprint:

[Multiscale Analysis of Pangenome Enables Improved Representation of Genomic Diversity For Repetitive And Clinically Relevant Genes](https://www.biorxiv.org/content/10.1101/2022.08.05.502980v2)

PGR-TK provides pangenome assembly management, query and Minimizer Anchored Pangenome (MAP) Graph Generation.

![Pangenome Data Management and Minimizer Anchored Pangenome Graph Generation](/images/PGR_TK_Sketch_MAPG_construction.png)

With the MAP graph, we can use the "principal bundle decomposition" to study complicated structural variants and genome rearrangements in the human population.

![AMY1A Example](/images/AMY1A_example.png)


## Documentation, Usage and Examples

### Command Line Tools

PGR-TK provides a single `pgr` binary with grouped subcommands.  Run `pgr --help`
for the top-level map, or `pgr <group> --help` for group-level help.

**Sequence archive management**
- `agc-rs`: create and manage AGC sequence archives (`.agcrs`)
  - `agc-rs create <archive.agcrs> --sample <name> <sequences.fa>` — create a new archive
  - `agc-rs append <archive.agcrs> --sample <name> <sequences.fa>` — add a sample to an existing archive

**`pgr index` — build shimmer indices**
- `pgr index mdb`: build a shimmer index from an AGC archive
  ```
  pgr index mdb --agcrs-input <archive.agcrs> [--prefix <prefix>]
  ```
  Produces `<prefix>.mdbi`, `<prefix>.mdbv`, and `<prefix>.midx`.
- `pgr index shmmr-count`: count shimmer occurrences across three sequence sets

**`pgr align` — genome alignment**
- `pgr align alnmap`: align long contigs to a reference; produces `.alndb`, `.alnmap`, `.ctgmap.*`, `.svcnd.*`
- `pgr align map-coord`: map query coordinates to target via an alnmap/alndb file
- `pgr align liftover-gtf`: lift GTF transcript annotations from a reference to haplotype contigs

**`pgr query` — search and fetch**
- `pgr query seqs`: query a PGR-TK pangenome database; outputs hit summary and FASTA files
  ```
  pgr query seqs --pgr-db-prefix <prefix> --query-fastx-path <query.fa> --output-prefix <out> [options]
  ```
- `pgr query fetch`: list or fetch sequences from a PGR-TK database by region
- `pgr query cov` / `pgr query cov2`: compare shimmer-pair coverage between two sequence sets

**`pgr bundle` — principal bundle decomposition**
- `pgr bundle decomp`: generate principal bundle decomposition via MAP Graph from a FASTA file
- `pgr bundle svg`: generate an interactive SVG from a principal bundle BED file
- `pgr bundle sort`: generate a contig sort order from the bundle decomposition
- `pgr bundle dist`: compute pairwise alignment scores from a bundle BED file
- `pgr bundle offset` / `pgr bundle shmmr-dist` / `pgr bundle aln`: further bundle analysis tools

**`pgr variant` — variant calling and annotation**
- `pgr variant diploid-vcf`: merge two haplotype alndb files into a phased diploid VCF
- `pgr variant annotate-vcf`: annotate a VCF with gene names from a GTF file
- `pgr variant annotate-bed`: annotate BED regions with gene annotation features
- `pgr variant sv-analysis`: analyse SV candidates with principal bundle decomposition
- `pgr variant merge-sv`: merge SV candidate BED records from multiple haplotypes

**`pgr plot` — visualisation**
- `pgr plot chr-aln`: generate chromosome alignment SVG plots from an alnmap JSON file

For each subcommand, `pgr <group> <subcommand> --help` provides detailed usage information.

The API documentation is at https://genedx.github.io/pgr-tk/

A collection of Jupyter Notebooks is at https://github.com/genedx/pgr-tk-notebooks/

### Typical workflow

```bash
# 1. Build an AGC archive incrementally, one sample at a time
agc-rs create pangenome.agcrs --sample GRCh38 GRCh38.fa
agc-rs append pangenome.agcrs --sample HG002_mat HG002_mat.fa
agc-rs append pangenome.agcrs --sample HG002_pat HG002_pat.fa

# 2. Build the shimmer index (prefix defaults to "pangenome")
pgr index mdb --agcrs-input pangenome.agcrs

# 3. Query a region of interest
pgr query seqs --pgr-db-prefix pangenome --query-fastx-path query.fa \
    --output-prefix output --max-count 128 --min-anchor-count 10
```


## Built Binaries

Check https://github.com/genedx/pgr-tk/releases


## Build

See `docker/Dockerfile.build_env-20.04` for a build environment under Ubuntu 20.04.
With the proper build environment, just run `bash build.sh` to build all.

For example, on a Mac OS with Docker installed, you can clone the repository and build a Linux binary
within an Ubuntu 20.04 Linux distribution as follows:

1. Build the Docker image for a build environment:

```bash
git clone git@github.com:cschin/pgr-tk.git
cd pgr-tk/docker
ln -s Dockerfile.build_env-20.04 Dockerfile
docker build -t pgr-tk-build .
```

2. In the root directory of the repo `pgr-tk`:

```bash
docker run -it --rm -v $PWD:/wd/pgr-tk pgr-tk-build /bin/bash
```

3. Build `pgr-tk` inside the Docker container:

```bash
cd /wd/pgr-tk
bash build.sh
```

The built Python wheels will be in `target/wheels` and can be installed for Ubuntu 20.04 Python 3.8.


### Build Singularity image

If you have built pgr-tk in a Docker container, you can use the following steps to build a Singularity image.

**Step 1: Commit Docker container to image**

```bash
docker commit <container_id> <image_name>:<version>
```

**Step 2: Push Docker image to Docker Hub**

```bash
docker login # if not already logged in
docker push <image_name>:<version>
```

**Step 3: Build Singularity image**

```bash
singularity build ./pgr-tk.<version>.sif docker://<docker_repo>/<image_name>:<version>
```

**Step 4: Execute**

```bash
singularity exec --fakeroot -B <host_path>:/<container_path> ./pgr-tk.<version>.sif \
    pgr index mdb --agcrs-input pangenome.agcrs
```

Replace `<host_path>` with the actual path you wish to bind to the container.

The `--fakeroot` option allows you to build and run images as a "fake" root user.
