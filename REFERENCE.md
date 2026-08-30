# arcasHLA references

arcasHLA references are built from a user-supplied, read-only IMGT/HLA release.
The source tree is only needed while building; runtime commands use the built
reference directory.

## Obtain IMGT/HLA

Clone a release:

```sh
git clone --depth 1 --branch 3.55.0 \
  https://github.com/ANHIG/IMGTHLA.git /opt/imgt-hla/3.55.0
```

Alternatively, unpack an IPD-IMGT/HLA release from the
[EBI download site](https://www.ebi.ac.uk/ipd/imgt/hla/download/). The source
must contain `hla.dat`, `wmda/hla_nom_g.txt`, and `wmda/hla_nom_p.txt`. Git LFS
must materialize `hla.dat` for GitHub releases 3.35.0 and newer.

## Build

```sh
arcasHLA reference build \
  --imgt /opt/imgt-hla/3.55.0 \
  --outdir /opt/arcashla-ref/3.55.0
```

Use `--version 3.55.0` when a release archive has no detectable version
metadata. Existing built output is protected unless `--force` is passed.
`--skip-customize` omits the large GRCh38 FASTAs when `customize` will not be
used.

The output includes complete and partial Kallisto indexes, conversion and
customization tables, static runtime configuration, and `manifest.json`.
`hla_transcripts.json` and the GRCh38 FASTAs are static arcasHLA assets rather
than IMGT-derived data; their provenance is recorded in the manifest.

## Minified references for testing

`reference build` can produce a much smaller reference for quick testing:

```sh
arcasHLA reference build \
  --imgt /opt/imgt-hla/3.55.0 \
  --outdir /opt/arcashla-ref/minified \
  --genes A,E --max-alleles-per-gene 3 \
  --skip-partial --skip-customize
```

`--genes` restricts the build to the listed genes, `--max-alleles-per-gene` caps
how many alleles of each are kept, and `--skip-partial` omits the partial allele
index, which is the largest artefact of a full build. Combining all three with
`--skip-customize` builds in seconds instead of minutes and produces a reference
of a few megabytes.

Alleles are selected in sorted order, so the same source and options always
produce the same reference. A minified reference does not cover enough of
IMGT/HLA to give correct genotypes; it gives reproducible ones. The selection is
recorded under `selection` in `manifest.json`, and `genotype --partial` refuses
to run against a reference built with `--skip-partial`.

## Select a reference

Set the environment variable for all commands:

```sh
export ARCASHLA_REF_DIR=/opt/arcashla-ref/3.55.0
arcasHLA genotype sample.1.fq.gz sample.2.fq.gz
```

Or select it per command:

```sh
arcasHLA genotype --ref /opt/arcashla-ref/3.55.0 \
  sample.1.fq.gz sample.2.fq.gz
```

The command-line path takes precedence over `ARCASHLA_REF_DIR`. A populated
repository `dat/` directory remains a deprecated fallback. `quant` retains its
historical `--ref` option for a customized sample index, so use `--reference`
there to select the built arcasHLA reference.

## Update

Place the new IMGT/HLA release in a new source directory, build into a new
output directory, validate it, then switch `ARCASHLA_REF_DIR`. Builds never
checkout, modify, or otherwise write to the IMGT/HLA source.

The legacy `reference --update`, `--version`, `--commit`, and `--rebuild`
options have been removed. Use `arcasHLA reference build --help`.
