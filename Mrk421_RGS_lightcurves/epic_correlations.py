import numpy as np
import os
from config import CONFIG
from astropy.io import fits
from tools import *
import glob
import sys
import logging
from scipy import stats
from scipy.optimize import curve_fit
MJDREF =  50814.0

logging.basicConfig(level=logging.INFO)

def skew_norm_pdf(x,e=0,w=1,a=0):
    # adapated from:
    # http://stackoverflow.com/questions/5884768/skew-normal-distribution-in-scipy
    t = (x-e) / w
    return 2.0 * w * stats.norm.pdf(t) * stats.norm.cdf(a*t)

def cross_correlation(rates1, rates2, times1, times2, x, str1, str2, obs):
    
    cc_simple = np.correlate(rates1-np.mean(rates1), rates2-np.mean(rates2), mode='same')
    cc_simple = cc_simple / (len(rates1) * rates1.std() * rates2.std()) #normalization
    fig, axs = plt.subplots(2, 1, figsize=(10,8), gridspec_kw={'hspace':0.3})
    axs[0].plot(MJDREF + (times1/86400.0), rates1/max(rates1), color='b', marker='o', markersize=3, linestyle='', label=str1)
    axs[0].plot(MJDREF + (times2/86400.0), rates2/max(rates2), color='red', marker='^', markersize=4, linestyle='', label=str2)
    axs[1].plot(x, cc_simple, marker='o', markersize=1, linestyle='')
    axs[0].ticklabel_format(useOffset=False)
    axs[0].tick_params(axis='both', which='major', labelsize=11)
    axs[1].grid(alpha=0.5)
    axs[0].grid(alpha=0.5)
    axs[0].set_xlabel('Time [MJD]', fontsize=14)
    axs[0].set_ylabel('Normalized Rate [ct/s]', fontsize=14)
    axs[0].legend(fontsize=13)
    axs[1].set_xlabel('Lag', fontsize=14)
    axs[1].set_ylabel(f'CCF', fontsize=14)
    axs[1].tick_params(axis='both', which='major', labelsize=11)

    #text = "Cross-correlation between " + str1 + "and" + str2
    #plt.title(text, fontsize=15)
    maxlag = x[np.argmax(cc_simple)]
    #axs[1].plot(maxlag, max(cc_simple), marker='x', markersize=6)
    #maxlag = cc_simple.argmax() - (len(rates1) - 1)
    print('CCF max:', max(cc_simple) )
    print(f"max correlation is at lag {maxlag} for observation {obs}")
    '''
    data = []
    if obs in ['0670920501' ]:
        x = x[100:-100]
        cc_simple = cc_simple[100:-100]
    
    elif obs in ['0670920401', '0670920301']:
        x = x[150:-50]
        cc_simple = cc_simple[150:-50]
    elif obs=='0502030101':
        x = x[400:-400]
        cc_simple = cc_simple[400:-400]       
    else:
        x = x[250:-250]
        cc_simple = cc_simple[250:-250]
    
    for i in range(0,len(x)):
        #print(x[100:-100][i])
        #print(int(cc_simple[100:-100][i]*100))
        data.append([x[i]]*int(cc_simple[i]*1000))
    
    data = [item for sublist in data for item in sublist]
    
    dist = stats.skewnorm(*stats.skewnorm.fit(data))
    x = np.arange(-100,100)
    y = dist.pdf(x)

    # the pdf is normalized, so we need to scale it to match the histogram
    y = y/y.max()
    y = y*cc_simple.max()

    plt.plot(x,y,'r',linewidth=2)
    print('max skewnorm:', x[np.argmax(y)])
    
    popt,pcov = curve_fit(skew_norm_pdf, x, cc_simple, p0=[0,1,1])

    #chi2 gauss
    #chisq_gauss =(((cc_simple-skew_norm_pdf(x, popt[0], popt[1], popt[2]))/bin_top)**2).sum()
    #ndof_gauss = len(counts) - 3
    #print('Chisquare/ndof gauss = %f/%d' % (chisq_gauss, ndof_gauss))
    plt.plot(x,skew_norm_pdf(x, popt[0], popt[1], popt[2]), 'g',linewidth=2)
    '''
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
        print(f'Processing observation {obs}.')
        os.chdir(os.path.join(target_dir, obs))
        sum_odf_dir = glob.glob('*SUM.SAS')[0]
        os.environ["SAS_ODF"] = os.path.join(target_dir, obs, sum_odf_dir)
        os.environ["SAS_CCF"] = os.path.join(target_dir, obs, 'ccf.cif')

        os.chdir(os.path.join(target_dir, obs, 'epic', 'pn')) #Go to epic/pn directory, where the rates are
        try:
            with open("ds9.reg") as f:
                region = f.read()

        except FileNotFoundError as e:
            logging.ERROR("Please extract manually source and background coordinates into a .reg file using ds9 software.")
            sys.exit(0)
            
        region = region.split('\n') #divide text file into lines

        #Source coordinates
        coordinates = region[3][4:-1].split(',')
        xcenter = float(coordinates[0])
        xmax = xcenter + float(coordinates[2])/2
        xmin = xcenter - float(coordinates[2])/2

        #Background coordinates
        coordinates_bkg = region[4][4:-1].split(',')
        xcenter_bkg = float(coordinates_bkg[0])
        xmax_bkg = xcenter_bkg + float(coordinates_bkg[2])/2
        xmin_bkg = xcenter_bkg - float(coordinates_bkg[2])/2

        #Run epiclccorr on timebin of 25s
        timebinsize = 25 #s
        if not glob.glob(os.path.join(target_dir, obs, 'epic', 'pn', "PN_hard_25.lc")):
            #Soft LC
            logging.info('Extracting Soft LC...')
            evselect_source_cmmd = f"evselect table=PNclean.fits energycolumn=PI expression='#XMMEA_EP && (PATTERN<=4) && (RAWX in [{xmin}:{xmax}]) &&! (RAWX in [{xcenter-1}:{xcenter+1}]) && (PI in [600:2000])' withrateset=yes rateset='PN_soft_raw_25.lc' timebinsize=25 maketimecolumn=yes makeratecolumn=yes"
            evselect_source_status = run_command(evselect_source_cmmd)

            evselect_bkg_cmmd = f"evselect table=PNclean.fits energycolumn=PI expression='#XMMEA_EP && (PATTERN<=4) && (RAWX>={xmin_bkg}) && (RAWX<={xmax_bkg}) && (PI in [600:2000])' withrateset=yes rateset='PN_bkg_soft_raw_25.lc' timebinsize=25 maketimecolumn=yes makeratecolumn=yes"
            evselect_bkg_status = run_command(evselect_bkg_cmmd)

            logging.info(f'Running epiclccorr for soft lightcurve of EPIC-PN...')
            epiclccorr_soft_cmmd = f"epiclccorr srctslist=PN_soft_raw_25.lc eventlist=PNclean.fits outset=PN_soft_25.lc bkgtslist=PN_bkg_soft_raw_25.lc withbkgset=yes applyabsolutecorrections=yes"
            epiclccorr_soft_status = run_command(epiclccorr_soft_cmmd)

            #Hard LC
            logging.info('Extracting Hard LC...')            
            evselect_source_cmmd = f"evselect table=PNclean.fits energycolumn=PI expression='#XMMEA_EP && (PATTERN<=4) && (RAWX in [{xmin}:{xmax}]) &&! (RAWX in [{xcenter-1}:{xcenter+1}]) && (PI in [2000:10000])' withrateset=yes rateset='PN_hard_raw_25.lc' timebinsize=25 maketimecolumn=yes makeratecolumn=yes"
            evselect_source_status = run_command(evselect_source_cmmd)

            evselect_bkg_cmmd = f"evselect table=PNclean.fits energycolumn=PI expression='#XMMEA_EP && (PATTERN<=4) && (RAWX>={xmin_bkg}) && (RAWX<={xmax_bkg}) && (PI in [2000:10000])' withrateset=yes rateset='PN_bkg_hard_raw_25.lc' timebinsize=25 maketimecolumn=yes makeratecolumn=yes"
            evselect_bkg_status = run_command(evselect_bkg_cmmd)

            logging.info(f'Running epiclccorr for hard lightcurve of EPIC-PN...')
            epiclccorr_hard_cmmd = f"epiclccorr srctslist=PN_hard_raw_25.lc eventlist=PNclean.fits outset=PN_hard_25.lc bkgtslist=PN_bkg_hard_raw_25.lc withbkgset=yes applyabsolutecorrections=yes"
            epiclccorr_hard_status = run_command(epiclccorr_hard_cmmd)


        #Get soft and hard lightcurves rates
        with fits.open('PN_soft_10.lc') as hdul:
            rates_soft = hdul['RATE'].data['RATE']
            erates_soft = hdul['RATE'].data['ERROR']
            times_soft = hdul['RATE'].data['TIME']
            #Drop NaN values by making a numpy mask
            mask_nan = np.invert(np.isnan(rates_soft)) 
            rates_soft = rates_soft[mask_nan]
            erates_soft = erates_soft[mask_nan]
            times_soft = times_soft[mask_nan]
            rms_soft = np.mean(rates_soft)

        with fits.open('PN_hard_10.lc') as hdul:
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
        npts = len(rates_soft)
        new_times = np.linspace(0, duration, len(rates_soft))
        lags = np.arange(-npts + 1, npts)
        lags = np.arange(-(npts/2), (npts/2))

        #Plot the cross correlation function (note that we are subtracting the RMS. See https://stackoverflow.com/questions/49742593/numpy-correlate-is-not-providing-an-offset for an explanation)
        cc = cross_correlation(rates_hard, rates_soft, times_hard, times_soft, lags, 'soft range (0.6 - 2 keV) ', ' hard range (2 - 10 keV)', obs=obs) 
        plt.savefig(os.path.join(target_dir, "Products", "EPIC_Lightcurves", f"{obs}_crosscorrelation.png"))
    
