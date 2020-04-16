'''
ekstest reads the FITS file containing the EPIC source and background time se-
ries’ produced by epiclccorr and when required its associated Good Timing Interval
(GTI) file produced by tabgtigen. It eliminates all bins that contain null values as
well as those that are not in the GTIs (if the file is provided). It then performs a va-
riety of variability tests, including the Kolmogorov-Smirnov probability of constancy
test, chi-squared probability of constancy test, flare tests and on the light curve to
determine whether the source is variable and writes the various statistics and proba-
bilities that the source is not variable and the number of good bins used to determine
these values into the header of a new version of the file or to the screen. No tests are
carried out if there are insufficient good bins.
'''
from astropy.io import fits
from astropy.table import Table
import pandas as pd 
import numpy as np

def excess_variance(rates, errrates, normalized=True):
    mean = np.mean(rates)
    variance = np.var(rates, ddof=1)
    N = len(rates)
    mse = np.mean(np.square(errrates))
    xs = variance - mse
    nxs = xs/(np.square(mean))
    f_var = np.sqrt(nxs)

    err_nxs = np.sqrt( np.square(np.sqrt(2/N)*mse/np.square(mean)) + np.square(np.sqrt(mse/N) *2*f_var/mean) )
    if normalized:
        return nxs, err_nxs
    else:
        return xs

def fractional_variability(rates, errrates):
    nxs, err_nxs = excess_variance(rates, errrates, True)
    
    f_var = np.sqrt(nxs)
    err_fvar = 1/(2*f_var) * err_nxs

    return f_var, err_fvar

def ekstest(timeseries, gtifile=None, screen=True):
    
    with fits.open(timeseries) as hdul:
        #Get dataset(s) and table.
        astro_table = Table(hdul[1].data)
        dataset = astro_table.to_pandas()
        print(dataset)
        #Check important keyword consistency (notably TSTART and TIMEDEL).
        
        #Call a warning or error if necessary.

        #Recover all light curves included in table :
        #Net source rates and errors and background rates and errors are recorded in arrays of dimension Nlightcurve * Nbins.
        
        #For each light curve :
        #Delete gaps in data (when the IEEE NaN constant is found) 
        dataset_cleaned = dataset.dropna()
        rates = np.array(dataset_cleaned['RATE'])
        errrates = np.array(dataset_cleaned['ERROR'])

        #and outside of GTIs.

        #Perform variability tests on rebinned counts and background or the net source counts and the cumulative time distribution:
        #Calculate mean count rate and variance.
        avgrate = np.mean(rates)
        variance = np.var(rates, ddof=1)
        

        #Test the null hypothesis (the source is not variable) with the :
        #- Kolmogorov-Smirnov test

        #- chi-squared test
        
        #- fractional variability amplitude test
        xs = excess_variance(rates, errrates, normalized=False)
        nxs, err_nxs = excess_variance(rates, errrates, normalized=True)
        
        f_var, err_fvar = fractional_variability(rates, errrates)
        
        #- flare test
        
        #- variation test
        
        #End for
        #Write variability test results into header and/or to screen.
        if screen:
            print('numpy mean', avgrate)
            print('numpy variace', variance)
            print('xs', xs)        
            print('nxs', nxs)
            print(f'Fractional Variability: {f_var}+-{err_fvar}')
            

       


if __name__ == "__main__":
    ekstest('MOS1_lccorr.lc')
