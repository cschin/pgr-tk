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

PGR-TK provides the following tools:

**Sequence archive management**
- `agc-rs`: create and manage AGC sequence archives (`.agcrs`)
  - `agc-rs create <archive.agcrs> --sample <name> <sequences.fa>` — create a new archive
  - `agc-rs append <archive.agcrs> --sample <name> <sequences.fa>` — add a sample to an existing archive

**Shimmer index creation** (sparse hierarchical minimizer)
- `pgr-mdb`: build a pgr shimmer index from an AGC archive
  ```
  pgr-mdb --agcrs-input <archive.agcrs> [--prefix <prefix>]
  ```
  If `--prefix` is omitted the stem of the archive path is used (e.g. `foo` for `foo.agcrs`).
  Produces `<prefix>.mdbi`, `<prefix>.mdbv`, and `<prefix>.midx`.

**Query**
- `pgr-query`: query a PGR-TK pangenome sequence database; outputs a hit summary and FASTA files from the target sequences
  ```
  pgr-query <db_prefix> <query.fa> <output_prefix> [options]
  ```

**MAP-graph and principal bundle decomposition**
- `pgr-pbundle-decomp`: generate principal bundle decomposition through MAP Graph from a FASTA file
- `pgr-pbundle-bed2svg`: generate SVG from a principal bundle bed file

**Auxiliary tools**
- `pgr-pbundle-bed2sorted`: generate an annotation file with a sorting order from the principal bundle decomposition
- `pgr-pbundle-bed2dist`: generate alignment scores between sequences using bundle decomposition from a principal bundle bed file

For each command, `command --help` provides detailed usage information.

The API documentation is at https://genedx.github.io/pgr-tk/

A collection of Jupyter Notebooks is at https://github.com/genedx/pgr-tk-notebooks/

### Typical workflow

```bash
# 1. Build an AGC archive incrementally, one sample at a time
agc-rs create pangenome.agcrs --sample GRCh38 GRCh38.fa
agc-rs append pangenome.agcrs --sample HG002_mat HG002_mat.fa
agc-rs append pangenome.agcrs --sample HG002_pat HG002_pat.fa

# 2. Build the shimmer index (prefix defaults to "pangenome")
pgr-mdb --agcrs-input pangenome.agcrs

# 3. Query a region of interest
pgr-query pangenome query.fa output --max-count 128 --min-anchor-count 10
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
    pgr-mdb --agcrs-input pangenome.agcrs
```

Replace `<host_path>` with the actual path you wish to bind to the container.

The `--fakeroot` option allows you to build and run images as a "fake" root user.
