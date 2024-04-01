# **segment_intensity**

## Installation
Before you start, ensure [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) is installed on your system.

Create a new environment and pip install [napari](https://napari.org/stable/tutorials/fundamentals/installation.html#installation):
```
conda create -y -n segment_intensity-env -c conda-forge python=3.11
conda activate segment_intensity-env
python -m pip install "napari[all]"
```

Pip install EpiTools and additional dependencies:
```
python -m pip install epitools pillow
```

## Usage
Split a multi-frame TIF file into individual TIF images:

Arguments:
* -i, --input : Path to the input multi-frame TIF file
* -o, --output : Directory where the split TIF images will be saved

```
python split_tif.py -i <input/path/to/file.tif> -o <output/path/to/directory>
```

Segment and quantify cell edge intensities from IHC data (single).

Arguments:
* -i, --input : Path to input IHC TIF file.
* -s, --spot : Gaussian sigma for controlling spot detection, default 6.
* -o, --outline : Gaussian sigma for controlling outline detection, default 3.
* -t, --threshold : Threshold for background cutoff, default 0.01.

```
python segment_intensity_single.py -i <input/path/to/file.tif>
```

Segment and quantify cell edge intensities from IHC data (batch).

Arguments:
* -i, --input : Directory containing input IHC TIF files.
* -s, --spot : Gaussian sigma for controlling spot detection, default 6.
* -o, --outline : Gaussian sigma for controlling outline detection, default 3.
* -t, --threshold : Threshold for background cutoff, default 0.01.

```
python segment_intensity_batch.py -i <input/path/to/directory>
```
