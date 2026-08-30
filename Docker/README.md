# Container for arcasHLA

Installs up-to-date versions of arcasHLA and all dependencies via
[pixi](https://pixi.sh), using the pinned `pixi.toml`/`pixi.lock`.

### Build

In order to use this arcasHLA container, install Docker and build from the
repository root. The local checkout is copied into the image:

```sh
docker build -t <image_name> -f Docker/Dockerfile .
```

To bake in an IMGT/HLA release and built reference:

```sh
docker build -t <image_name> \
  --build-arg IMGT_HLA_VERSION=3.55.0 \
  --build-arg ARCASHLA_REF_DIR=/opt/arcashla-ref .
```

To bake in an IMGT/HLA release and built reference:

```sh
docker build -t <image_name> \
  --build-arg IMGT_HLA_VERSION=3.55.0 \
  --build-arg ARCASHLA_REF_DIR=/opt/arcashla-ref .
```

### Run

Interactively ("image_name" is as above):

```
docker run -it --entrypoint bash -v <path/to/files>:<docker/path/to/files> <image_name>
```

Noninteractively ("image_name" is as above), e.g. 'arcasHLA extract':

```
docker run \
	--rm \
	-v <path/to/files>:<docker/path/to/files> \
	-v <path/to/reference>:/opt/arcashla-ref \
	-e ARCASHLA_REF_DIR=/opt/arcashla-ref \
	<image_name> \
	arcasHLA extract --o docker/path/to/files/out_dir \
	[other options] docker/path/to/files/sample.bam & 
```
