# KSSD3: K-mer Space Sampling and Decomposition

KSSD3 is an alignment-free toolkit for fast comparison of genomic sequences. It uses probabilistic sampling of the k-mer space to build compact sketches and compute distances, containment metrics, ANI estimation, and metagenomic composition.

## Building from Source

The repository ships with a Makefile and a precompiled copy of XGBoost. To build the `kssd3` binary run:

```bash
make
```

The build requires GCC with OpenMP support and `zlib`. When the build completes, the executable is placed in `bin/kssd3`.

To make the bundled XGBoost library discoverable, append it to `LD_LIBRARY_PATH`:

```bash
make install_env
source ~/.bashrc
```

Optionally install the binary into `/usr/local/bin` with:

```bash
sudo make install
```

## Basic Usage

`kssd3` exposes several subcommands:

- `shuffle` – generate dimension-reduction shuffle files.
- `sketch` – create sketches from FASTA/FASTQ sequences.
- `dist` – compute distances between sketches.
- `set` – union/intersection/subtraction operations on sketches.
- `reverse` – recover k-mers from sketches.
- `composite` – estimate metagenomic composition.
- `matrix` – compute distance matrices.
- `ani` – estimate average nucleotide identity (ANI).

Display help for any subcommand with `--help` or `--usage`.

### Example: Sketching and ANI Calculation

```bash
# create shuffle file if needed
./bin/kssd3 shuffle --usedefault -o default

# sketch reference and query genomes
./bin/kssd3 sketch -o ref_sketches refs/*.fasta
./bin/kssd3 sketch -o qry_sketches queries/*.fasta

# build inverted index for reference sketches
./bin/kssd3 sketch -i ref_sketches

# estimate ANI
./bin/kssd3 ani -r ref_sketches -q qry_sketches -m 0 > ani.tsv
```

## Further Documentation

Command line options for each subcommand are provided by the built-in help system. For license details see [`LICENSE.txt`](LICENSE.txt). The `ani_models` directory contains example XGBoost models for ANI prediction.
