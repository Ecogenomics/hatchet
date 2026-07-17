# hatchet

Tools used to split the GTDB-Tk reference tree into smaller sub-trees.

`hatchet` produces the "divide-and-conquer" reference packages used by
[GTDB-Tk](https://github.com/Ecogenomics/GTDBTk) v2. Instead of placing query
genomes into a single, very large reference tree (which requires a large amount
of RAM), GTDB-Tk v2 first places genomes into a small **backbone** tree and then
into an appropriate **class-level sub-tree**. `hatchet` is the tool that builds
those trees, alignments, RED files, and pplacer reference packages from a GTDB
release.

> **Note:** `hatchet` is made available for transparency and reuse, but it is
> primarily intended for **internal use by the GTDB team**. It expects files and
> external tools from the GTDB build environment and is not a general-purpose
> phylogenetics utility.

## How it works

The full workflow (`hatchet_wf`) runs four steps:

1. **`pick`** – Select one high-quality representative genome per taxon at a
   given rank (default: family) to build the **backbone** tree, MSA and pplacer
   package. Representatives are chosen by a quality score (completeness,
   contamination, MAG/SAG status, contig count, ambiguous bases, MSA gaps).
2. **`red`** – Regenerate [RED](https://gtdb.ecogenomic.org/) (Relative
   Evolutionary Divergence) values for the pruned backbone tree from the
   original reference RED file.
3. **`split`** – Split the reference tree into class-level sub-trees. Each
   sub-tree is capped at roughly the size of the largest class (with a 10 %
   tolerance); one genome per absent phylum is added as an outgroup, and each
   sub-tree gets its own MSA, tree and pplacer package.
4. **`red_low`** – Regenerate RED values for each of the sub-trees.

The results are then packaged under `<out_dir>/to_copy/` ready to be dropped
into a GTDB-Tk reference data release.

## Requirements

**Python:** 3.6+

Python packages:

- [`dendropy`](https://dendropy.org/)

External command-line tools (must be on `$PATH`):

- [`genometreetk`](https://github.com/Ecogenomics/GenomeTreeTk) – tree stripping and re-rooting
- [`phylorank`](https://github.com/Ecogenomics/PhyloRank) – tree decoration
- [`taxit`](https://github.com/fhcrc/taxtastic) (taxtastic) – pplacer reference package creation

`biolib_lite` and `phylorank_lite` are already inside the package, so no
separate install is required for them.

## Installation

```bash
git clone https://github.com/Ecogenomics/hatchet.git
cd hatchet
pip install .
```

This installs the `hatchet` command. For development use `pip install -e .`.

## Usage

```bash
hatchet -h            # top-level help / list of commands
hatchet <command> -h  # help for a specific command
```

### Commands

| Command      | Purpose                                                                 |
|--------------|-------------------------------------------------------------------------|
| `pick`       | Pick one representative genome per taxon at a rank (builds the backbone) |
| `red`        | Regenerate RED values for the pruned backbone tree                       |
| `split`      | Split the reference tree into class-level sub-trees                      |
| `red_low`    | Regenerate RED values for each split sub-tree                            |
| `merge_logs` | Replace the last tree in a FastTree log with the original tree (pplacer) |
| `unroot`     | Unroot a tree                                                           |
| `hatchet_wf` | Run the full workflow: `pick` → `red` → `split` → `red_low`             |

### Full workflow example

```bash
hatchet hatchet_wf \
    --domain bac \
    --ref_tree      gtdb_r220_bac120.tree \
    --metadata      gtdb_r220_metadata.tsv \
    --msa           gtdb_r220_bac120_msa.faa \
    --tax           gtdb_r220_bac120_taxonomy.tsv \
    --red_file      gtdb_r220_bac120_red_values.tsv \
    --original_log  gtdb_r220_bac120_fasttree.log \
    --out_dir       hatchet_out
```

### Inputs

- **`--ref_tree`** – GTDB reference tree (Newick).
- **`--metadata`** – GTDB metadata table (TSV) with the standard GTDB columns
  (`gtdb_representative`, `gtdb_taxonomy`, `checkm_completeness`,
  `checkm_contamination`, `ncbi_genome_category`, `contig_count`,
  `ambiguous_bases`).
- **`--msa`** – Reference multiple sequence alignment (FASTA).
- **`--tax`** – Taxonomy file mapping each genome to its GTDB lineage (TSV).
- **`--red_file`** – Original reference RED values.
- **`--original_log`** – FastTree log used to build the pplacer packages.
- **`--domain`** – `bac` or `arc`.

### Outputs

Running `hatchet_wf` writes, under `--out_dir`:

- `backbone/` – backbone genome selection, pruned tree, and pplacer package.
- `class_level/` – one MSA, tree, RED file and pplacer package per sub-tree,
  plus `tree_mapping.tsv` (which taxon maps to which sub-tree).
- `to_copy/` – the subset of files packaged for a GTDB-Tk data release.

## Copyright

Copyright © 2021 Pierre-Alain Chaumeil. See [LICENSE](LICENSE) for details.
Released under the GNU General Public License v3.