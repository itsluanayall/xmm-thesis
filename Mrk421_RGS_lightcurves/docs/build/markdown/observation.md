# observation module


### class observation.Exposure(evenli, srcli='srcli', specli='specli', bkgli='bkgli', respli='respli')
Bases: `object`

Class for an Exposure of an Observation.


#### synchronous_times(exp2)
Compares the exposure’s start and stop times to those of another exposure. If the exposures overlap,
the method returns the start time and the stop times of the overlapping interval.


* **Parameters**

    **exp2** (*class Exposure*) – exposure you want to compare self with.



* **Raises**

    **Exception** – when the two exposures are not synchronous



* **Returns**

    two floats, representing the start time and stop time for the pair of exposures



### class observation.Observation(obsid, target_dir)
Bases: `object`

Observation class for a specific target.
In order to avoid unnecessary complexity, the method names of the class recall
the SAS commands used to analyse observations.


#### bkg_lightcurve()
Generates Background lightcurve associated to each exposure. Useful to check if there are flares.
Follows tutorial on SAS thread of RGS background.


#### check_flaring_particle_bkgr()
Checks if the background lightcurve (CCD9 of RGS) affected by solar flares ecc, has significant peaks. 
In order to do this, the method calculates the mean of the background lightcurve, its standard deviation, 
and if there are datapoints that are >3 sigma from the average, it counts that value as a significant flare. 
The method collects all these flares in an array and then calls the SAS command tabgtigen and cuts on
RATE<maxr where maxr is the minimum element of the flare array.


#### cifbuild()
Generates the Calibration Index File (CIF) for the observation
in the directory of the observation, and sets the ‘SAS_CCF’ variable pointing to it.


#### create_pairs_exposures()
Defines lists of pairs of the RGS exposures of the observation. Some observations (0510610101, 0510610201, 0136540701) 
pairs are defined by hand because these observations were interrupted and so their exposures do not 
follow the standard conventions.


#### divide_spectrum()
Splits the spectrum into pieces of 1000 seconds each. To do this, I use a while loop with 2 indices:
- j represents the end_time of each piece, so after an iteration you just add 1000 seconds to it;
- i represents the start_time of each piece, so after an iteration you just set it to j.
At the end of the processing, we will have cut into pieces of 1000s each exposure. For instance, if there are 50 pieces for each exposure, the output will be 50 spectra for RGS1 and 50 spectra for RGS2. The nomenclature of the output spectra is the following:

sourcespec{instrume}_gti{k}.fits for the source spectrum
bkgspec{instrume}_gti{k}.fits for the background spectrum

where instrume can be RGS1 or RGS2, and k=1,2,3…50 is the piece we are considering.


#### property emdir()

#### property epdir()

#### epic_lightcurve(mjdref)
Makes the lightcurve plots with matplotlib. The first plot will consist of 4 panels containing:
the soft lightcurve (0.2 - 2 keV), the hard lightcurve (2 - 10 keV) and the respective background lightcurves.
The second plot consists of 3 panels: the soft and hard lightcurves and the hardness ratio calculated as
HR := (H-S)/(H+S) where H and S are the datapoints for the hard and soft lightcurves.


#### epiclccorr(pileup=False)
Extracts a source and background raw lightcurve for EPIC-pn and runs epiclccorr to correct the raw lightcurve.


#### epproc()
Runs the epproc SAS command to process and reduce EPIC-PN data.


#### filter_epic(pileup=False)
Filter an EPIC PN event list for periods of high background flaring activity.


#### fracvartest(screen=True, netlightcurve=True, instrument='rgs')
Reads the FITS file containing the RGS source and background timeseries produced by rgslccorr. 
It then calculates excess variance, normalized excess variance and fractional variability of the lightcurve,
storing all these values into a dictionary that will then be reported in the final .csv file.


* **Parameters**

    
    * **screen** (*boolean*) – if True, prints on terminal the results of the variability quantities, defaults to True


    * **netlightcurve** (*boolean*) – if True, uses the net lightcurve (i.e. background substracted) to calculate the variability quantities, defaults to True



#### lightcurve(mjdref)
Makes the lightcurve plot and saves it in the rgs directory of the current observation
The user here has two options: you can either save the lightcurve using the dsplot command
with XmGrace as plotting package in ‘.ps’ format, or you can plot it using python’s
module matplotlib and save it in the ‘.png’ format. 
You can choose the desired option setting the USE_GRACE boolean in the config.json file.


* **Parameters**

    
    * **mjdref** (*float*) – MJD corresponding to the beginning of the XMM-Newton mission. It is needed because the times are all in MET (Mission Elapsed Time)


    * **use_grace** (*boolean*) – if set to True, plots the lightcurve using the lotting interface Grace, if set to False it uses matplotlib



#### property obsdir()

#### property obsid()

#### property odfdir()

#### odfingest()
Generates the SUM.ASC file in the observation directory,
that summarizes all the observational info
and sets ‘SAS_ODF’ variable pointing to it.


#### pn_spectrum(pileup=False)
Follows the steps from the following SAS thread: 
[https://www.cosmos.esa.int/web/xmm-newton/sas-thread-pn-spectrum-timing](https://www.cosmos.esa.int/web/xmm-newton/sas-thread-pn-spectrum-timing)


#### pn_xspec(target_REDSHIFT)
XSPEC analysis of the EPIC-pn spectrum of the observation. The steps are all logged into a file called XSPECLogFile_{self.obsid}_spectrum.txt.
The fit is performed on two different models: logparabola and powerlaw. The plot of the spectra and residuals is done
using matplotlib. The plotting function is written in tools.py.
The flux and luminosity are stored, given the target_REDSHIFT as argument.
The fitted parameters and the spectrum counts are all stored into an astropy Table that is then saved as a FITS file.


* **Parameters**

    **target_REDSHIFT** (*float*) – redshift of target



#### property rgsdir()

#### rgslccorr()
Runs the rgslccorr SAS command for each pair (RGS1+RGS2) of exposures present in the observation.
The products are the lightcurves .ds files.
For the observations 0510610101, 0510610201 and 013654701 the eventlists are written manually because
their exposures do not copme in pairs.


#### rgsproc()
Runs the rgsproc SAS command to process and reduce RGS data.


#### property target_dir()

#### vaughan_panel(N, M, timescale=70, timebinsize=25)
Generates a variability plot, along the lines of those in Vaughan et al.


* **Parameters**

    
    * **N** (*int*) – number of bins to average on (from x to <x>)


    * **M** – number of bin to average on after having averaged on N bins


    * **timescale** (*int*) – how long the duration of the lightcurve must be, in kiloseconds


    * **timebinsize** (*int*) – bin size of lightcurve to insert as input parameter to rgslccorr, in seconds



#### xspec_divided_spectra(target_REDSHIFT)
If the divide_spectrum method has already been called, this method allows to perform analysis on all the pieces 
into which we have divided the spectrum using a for loop. The results are stored in the same astropy Table
of the average spectrum.


* **Parameters**

    **target_REDSHIFT** (*float*) – redshift of target



#### xspec_divided_spectra_average(target_REDSHIFT)
XSPEC analysis of the total, average spectrum (RGS1+RGS2) of the observation. The steps are all logged into a file called XSPECLogFile_average_spectrum.txt.
The fit is performed on two different models: logparabola and powerlaw. The plot of the spectra and residuals is done
using matplotlib. The plotting function is written in tools.py
The flux and luminosity are stored, given the target_REDSHIFT as argument.
The fitted parameters and the spectrum counts are all stored into an astropy Table that is then saved as a FITS file.


* **Parameters**

    **target_REDSHIFT** (*float*) – redshift of target
