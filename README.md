# Introduction
Personal backup of code written for my Master thesis at ESA. The code is written for XMM-Newton data.

##Data structure 
I didn't load the data onto git because they are a lot and heavy, but make sure the data for a specific target (e.g. Markarian421) follows
this structure:

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

If you have a lot of data and your structure is not in this form, then take a look at the `download_multiple_data.ipynb` file. 
There you will find a function that creates this structure automatically. 

##Usage
In order to use this package you must have HEASOFT and SAS installed on your computer (go to https://www.cosmos.esa.int/web/xmm-newton/sas-installation).
If you already have these softwares installed, I recommend making a conda environment. So, finally, the starting commands are:
```
heainit
sasinit
source activate env_name

```
where env_name is the name of the conda environment.
