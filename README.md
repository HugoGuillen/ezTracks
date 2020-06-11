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
`ezTracks` only needs a configuration file to work. To generate a plot, you need to do the following steps:
1. Activate the `eztracks` environment (or having all dependencies available in the PATH)
2. Run `python eztracks.py check path/to/config_file`
3. Run `python eztracks.py prepare path/to/config_file`
4. Run `python eztracks.py prepare path/to/config_file`
5. Modify `output_path/config.ini` for cosmetic changes and rerun the `pyGenomeTracks` command.