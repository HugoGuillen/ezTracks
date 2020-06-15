![logo](logo.png)
# ezTracks
`ezTracks` (*easy tracks*) is a [pyGenomeTracks](https://github.com/deeptools/pyGenomeTracks) wrapper for plotting a single GTF annotation followed by grouped bed files. `ezTracks` preprocesses the input tracks allowing to change and render the plot faster.

![output plot](test_output/my_tracks.png)

---
## Setup
The program it's just the file `eztracks.py` and works on Linux and Windows (Ubuntu WSL1) (I have not checked on Mac but it should work too). 

I suggest to setup a `conda` environment for installing `ezTracks` dependencies:

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
`ezTracks` only needs a configuration file to work. The only options implemented are the same as the sample file `test/test_config.ini`. The coordinates must be in bed format. To generate the tracks shown above you just need to enter the `test` folder and run: 

```bash
conda activate eztracks
python eztracks.py check test/test_config.ini
python eztracks.py prepare test/test_config.ini
python eztracks.py draw test/test_config.ini
```

The complete tutorial is located [here](tutorial.md).
