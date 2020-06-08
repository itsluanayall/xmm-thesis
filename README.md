# Introduction
Personal backup of code written for my Master thesis at ESA. The code is written for XMM-Newton data.

## Data structure 
I didn't load the data onto git because they are a lot and heavy, but make sure the data for a specific target (e.g. Markarian421) follows
this structure:
```
Markarian421
├── 0099280101
│   ├── odf
│   └── rgs
├── 0099280201
│   ├── odf
│   └── rgs
├── 0099280301
│   ├── odf
│   └── rgs
├── 0099280401
│   ├── odf
│   └── rgs
├── 0099280501
│   ├── odf
│   └── rgs
├── 0099280601
│   ├── odf
│   └── rgs
├── etc.
```
If you have a lot of data and your structure is not in this form, then take a look at the `download_multiple_data.ipynb` file. 
There you will find a function that creates this structure automatically. 

## Usage
In order to use this package you must have HEASOFT and SAS installed on your computer (go to https://www.cosmos.esa.int/web/xmm-newton/sas-installation).
If you already have these softwares installed, I recommend making a conda environment. So, finally, the starting commands are:
```
heainit
sasinit
source activate env_name
```
where env_name is the name of the conda environment.

## Logic
The script was written following the object-oriented program paradigm, which is very useful when you want to perform operations on the same type of data. The objects that I defined are the Exposure object and the Observation object, whose source code can be found in
```
observation.py
```
The methods and attributes of each object will be explained thouroghly in an upcoming documentation. 
The main file which creates the different instances of the Observation object and calling each method is:
```
RGS_lc_mrk421_automatized.py
```
This script runs the pipeline for all the observations, generating lightcurves, dividing and analysing the spectra, and saving the parameters in a final FITS file.

These two scripts rely on some functions defined in a separated file called
```
tools.py
```

All the scripts were written with the idea to leave the user very little to edit. For this purpose I made a 
```
config.json
```
file to be modified by the user. The variables the user can change are the path where to find the target data, the target redshift, the fractional variability timescale, and so on.
