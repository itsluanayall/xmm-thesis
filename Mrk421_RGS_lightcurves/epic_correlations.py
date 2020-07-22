import numpy as np
import os
from config import CONFIG
from astropy.io import fits
from tools import *
from scipy import signal

def cross_correlation(rates1, rates2, times1, times2, x, str1, str2):
    
    cc_simple = np.correlate(rates1, rates2, mode='same')
    cc_simple = cc_simple / (len(rates1) * rates1.std() * rates2.std())
    fig, axs = plt.subplots(2, 1, figsize=(10,8), gridspec_kw={'hspace':0.3})
    #plt.figure(figsize=(8,8))
    print(len(x), len(rates1))
    axs[0].plot(times1, rates1, color='b', marker='o', markersize=2, linestyle='', label='soft lightcurve')
    axs[0].plot(times2, rates2, color='red', marker='o', markersize=2, linestyle='', label='hard lightcurve')
    axs[1].plot(x, cc_simple, marker='o', markersize=4)
    plt.grid(alpha=0.5)
    axs[0].set_xlabel('Time[s]')
    axs[0].set_ylabel('Mean-subtracted rates [ct/s]')
    axs[0].legend()
    axs[1].set_xlabel('Lag of soft lc relative to hard lc [s]')
    axs[1].set_ylabel('Cross-correlation')
    text = "Cross-correlation between " + str1 + "and" + str2
    plt.title(text, fontsize=15)
    maxlag = x[np.argmax(cc_simple)]
    #maxlag = cc_simple.argmax() - (len(rates1) - 1)
    print("max correlation is at lag %d" % maxlag)
    
    return cc_simple

if __name__ == "__main__":
    epic_observations = ['0670920301', '0670920401', '0670920501', '0302180101', '0150498701', '0502030101']
    #epic_observations = ['0670920501']
    version = CONFIG['version']
    sas_dir = CONFIG['sas_dir']
    ccf_dir = CONFIG['ccf_dir']
    mjdref = CONFIG['MJDREF']
    target_dir = CONFIG['target_dir']
    target_REDSHIFT = CONFIG['target_REDSHIFT']
    use_grace = CONFIG['use_grace']
    timescale_fvar = CONFIG['timescale_fvar']
    setupSAS(sas_dir=sas_dir, ccf_dir=ccf_dir)

    for obs in epic_observations:
        os.chdir(os.path.join(target_dir, obs, 'epic', 'pn')) #Go to epic/pn directory, where the rates are
        
        #Get soft and hard lightcurves rates
        with fits.open('PN_soft.lc') as hdul:
            rates_soft = hdul['RATE'].data['RATE']
            erates_soft = hdul['RATE'].data['ERROR']
            times_soft = hdul['RATE'].data['TIME']
            #Drop NaN values by making a numpy mask
            mask_nan = np.invert(np.isnan(rates_soft)) 
            rates_soft = rates_soft[mask_nan]
            erates_soft = erates_soft[mask_nan]
            times_soft = times_soft[mask_nan]
            rms_soft = np.mean(rates_soft)

        with fits.open('PN_hard.lc') as hdul:
            rates_hard = hdul['RATE'].data['RATE']
            erates_hard = hdul['RATE'].data['ERROR']
            times_hard = hdul['RATE'].data['TIME']
            
            #Drop NaN values by making a numpy mask
            mask_nan = np.invert(np.isnan(rates_hard)) 
            rates_hard = rates_hard[mask_nan]
            erates_hard = erates_hard[mask_nan]
            times_hard = times_hard[mask_nan]
            rms_hard = np.mean(rates_hard)

        duration = times_soft[-1] - times_soft[0]
        print('duration', duration)
        npts = len(rates_soft)
        new_times = np.linspace(0, duration, len(rates_soft))
        lags = np.arange(-npts + 1, npts)
        lags = np.arange(-(npts/2), (npts/2))
        #Plot the cross correlation function (note that we are subtracting the RMS. See https://stackoverflow.com/questions/49742593/numpy-correlate-is-not-providing-an-offset for an explanation)
        cc = cross_correlation(rates_soft-rms_soft, rates_hard-rms_hard, times_soft, times_hard, lags, 'soft range (0.2 - 2 keV) ', ' hard range (2 - 10 keV)') 
        plt.savefig(os.path.join(target_dir, "Products", "EPIC_Lightcurves", f"{obs}_crosscorrelation.png"))