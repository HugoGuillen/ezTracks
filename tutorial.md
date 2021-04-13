# Tutorial
## Introduction
`ezTracks` only needs a configuration file to work. The only options implemented are the same as the sample file `test/test_config.ini`. The coordinates must be in bed format.

To process to generate a plot is:
1. Activate the `eztracks` environment (or having all dependencies available in the PATH)
2. Run `python eztracks.py check path/to/config_file`
3. Run `python eztracks.py prepare path/to/config_file`
4. Run `python eztracks.py draw path/to/config_file`
5. Modify `output_path/config.ini` for cosmetic changes and rerun the `pyGenomeTracks` command.

## 1. Check the integrity of the test data
We are going to work with the test data provided in this repository (they are cropped to the region of interest of the configuration file, so you may want to try with your own data after finishing this tutorial). In your terminal enter this folder, and examine the configuration file.
```bash
# TO UPDATE
> cat test_config/config.ini
```

Now, we are going to activate the environment and verify the config file and the input data:
```bash
# TO UPDATE
> conda activate eztracks
> python eztracks.py check test_config/test_config.ini
```

## 2. Preprocess track files
If there were no errors in the last step, now we can preprocess the input tracks. Basically we will get only the portion we are interested to plot into the folder `output_path/prep`:


```bash
> python eztracks.py prepare test/test_config.ini
# TO UPDATE

```

## 3. Plot the tracks


```bash
# TO UPDATE
> python eztracks.py draw test/test_config.ini

```

And now the plot is ready in `test_output/my_tracks.pdf`! As the terminal output indicates, you can modify `test_output/config.ini` and rerun the `pyGenomeTracks` command to redraw the plot. The possible parameters for this new configuration file are located at https://pygenometracks.readthedocs.io/en/latest/content/possible-parameters.html

![output plot](test_output/test_region/my_tracks.png)