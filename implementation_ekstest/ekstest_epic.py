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

def ekstest(timeseries, gtifile=None):
    with fits.open(timeseries) as hdul:
        print(hdul[0].header)

if __name__ == "__main__":
    ekstest('MOS1_lccorr.lc')
    print('ciao')