# ezTracks
`ezTracks` (*easy tracks*) is a pyGenomeTracks’s wrapper for plotting a single GTF annotation followed by grouped bed files. `ezTracks` preprocess the input tracks allowing to change and render the plot faster.

---
## Setup
The program it's just the file `eztracks.py`. I suggest to setup a `conda` environment for installing `ezTracks` dependencies:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda create -n eztracks --yes python=3
conda activate eztracks
conda install -c anaconda freetype
conda install -c bioconda bedtools
pip install git+https://github.com/deeptools/HiCMatrix.git
pip install git+https://github.com/deeptools/pyGenomeTracks.git
```

*Note*: if you are having problems activating the environment, probably you can replace the line 
```bash
conda activate eztracks
````
with
```bash
source ~/anaconda3/etc/profile.d/conda.sh && conda activate eztracks
```

---
# Micro tutorial
## Introduction
`ezTracks` only needs a configuration file to work. To generate a plot, you need to do the following steps:
1. Activate the `eztracks` environment (or having all dependencies available in the PATH)
2. Run `python eztracks.py check path/to/config_file`
3. Run `python eztracks.py prepare path/to/config_file`
4. Run `python eztracks.py prepare path/to/config_file`
5. Modify `output_path/config.ini` for cosmetic changes and rerun the `pyGenomeTracks` command.
## 1. Test data
We are going to work with the test data provided in this repository. Enter this folder, and examine the configuration file:
```bash
> cat test/test_config.ini
[default]
chr = chrX
start = 73816524
end = 73859248
output_path = test_output
input_gtf = test/test_tracks/test.gtf

[tracks]
CRS = test/test_tracks/CRS.bed.gz
phastCons = test/test_tracks/phastcons
repeats = test/test_tracks/repeats

[plot]
output_img = my_tracks.pdf
gtf_title = Gencode_v34
gtf_height = 10
padding = -10000
width = 40
dpi = 100
```
