"""
Spectra Division Doc. To complete.
"""
import os
from astropy.io import fits, ascii
from tools import *
from config import CONFIG
import logging
import xspec
import glob
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np

logging.basicConfig(level=logging.INFO)

print(__doc__)
if __name__ == "__main__":

    version = CONFIG['version']
    sas_dir = CONFIG['sas_dir']
    ccf_dir = CONFIG['ccf_dir']
    target_dir = CONFIG['target_dir']
    target_redshift = CONFIG['target_REDSHIFT']
    logging.info(f'MRK421 Analysis - version:{version}')
    setupSAS(sas_dir=sas_dir, ccf_dir=ccf_dir)

    #Remove the .tar files from which we extracted the data
    for directory in os.listdir(target_dir):
        if directory.endswith('.tar.gz'):
            os.remove(os.path.join(target_dir, directory))
    
    #Create Products directory
    if not os.path.isdir(f'{target_dir}/Products'):
        os.makedirs(f'{target_dir}/Products')
    if not os.path.isdir(f'{target_dir}/Products/RGS_Spectra'):
        os.makedirs(f'{target_dir}/Products/RGS_Spectra')

    sample_observations = ['0099280301', '0411080101', '0502030101', '0658800701', '0791781201', '0411082701']


    for observation in sample_observations:


        os.environ["SAS_ODF"] = f"{target_dir}/{observation}/*SUM.SAS"
        os.environ["SAS_CCF"] = f"{target_dir}/{observation}/ccf.cif"
        if not os.path.isdir(f'{target_dir}/Products/RGS_Spectra/{observation}'):
            os.makedirs(f'{target_dir}/Products/RGS_Spectra/{observation}')
        if not os.path.isdir(f'{target_dir}/{observation}/rgs/divided_spectra'):
            os.makedirs(f'{target_dir}/{observation}/rgs/divided_spectra')

        os.chdir(f"{target_dir}/{observation}/rgs")
        

        evenli = glob.glob('*EVENLI0000.FIT')
        srcli = glob.glob('*SRCLI_0000.FIT')
        respli = glob.glob('*RSPMAT1003.FIT')
        pairs_events = sort_rgs_list(evenli, 'instr')
        pairs_srcli = sort_rgs_list(srcli, "instr")
        pairs_respli = sort_rgs_list(respli, "instr")

        print(pairs_events)
        print(pairs_srcli)
        print(pairs_respli)

        if not pairs_events:
            continue

        for n in range(0, len(evenli)):
            
            with fits.open(pairs_events[0][n]) as hdul:
                instrume = hdul['EVENTS'].header['INSTRUME']
                start_time = hdul['EVENTS'].header['TSTART']
                stop_time = hdul['EVENTS'].header['TSTOP']

            n_intervals = int(np.ceil((stop_time-start_time)/1000))
            
            if not glob.glob(f'divided_spectra/sourcespecRGS1*') or not glob.glob(f'divided_spectra/sourcespecRGS2*'):
                logging.info(f"Processing {pairs_events[0][n]}...")
                step = 1000 #seconds
                i = start_time
                j = i + step
                k = 1

                #Divide evenlist in parts of 1000s
                while (i<j) and (j<=stop_time):

                    tabgtigen_cmd = f"tabgtigen table={pairs_events[0][n]} expression='TIME in [{i}:{j}]' gtiset=divided_spectra/gti{instrume}_file{k}.fits"
                    status_tabgtigen = run_command(tabgtigen_cmd)

                    rgsfilter_cmd = f"rgsfilter  mergedset={pairs_events[0][n]} evlist=divided_spectra/event_gti{instrume}file{k}.fits auxgtitables=divided_spectra/gti{instrume}_file{k}.fits"
                    status_rgsfilter = run_command(rgsfilter_cmd)

                    rgsspectrum_cmd = f"rgsspectrum evlist=divided_spectra/event_gti{instrume}file{k}.fits srclist={pairs_srcli[0][n]} withspectrum=yes bkgcorrect=no spectrumset=divided_spectra/sourcespec{instrume}_gti{k}.fits withbkgset=yes bkgset=divided_spectra/bgspec{instrume}_gti{k}.fits order=1 rebin=1 edgechannels=2 spectrumbinning=lambda withfracexp=no badquality=1"
                    status_rgsspectrum = run_command(rgsspectrum_cmd)

                    # Change index
                    i = j
                    j = i + step
                    if j>stop_time:
                        j = stop_time
                    k += 1
            
        ## Spectra analysis with XSPEC
        if not glob.glob(f'{target_dir}/Products/RGS_Spectra/{observation}/*.gif'):
            
            logging.info(f"Starting spectral analysis with XSPEC for observation {observation}.")
            xspec.Xset.chatter = 10
            xspec.Xset.logChatter = 20
            logFile = xspec.Xset.openLog("divided_spectra/XSPECLogFile.txt") #Create and open a log file for XSPEC
            logFile = xspec.Xset.log
            spectra_table = Table(names=('ObsId', 'Instrument', 'Piece', 'Model',
                                        'PhoIndex', 'PhoIndex_sigma', 'Alpha', 'Alpha_sigma', 'Beta', 'Beta_sigma',
                                        'Flux[erg/cm2/s]', 'Flux_error_range',
                                        'Luminosity[e+44erg/s]', 'Lumin_error_range',
                                        'Fit_statistic'),
                                        dtype=('object', 'object','object', 'object', 'object', 'object', 'object', 'object', 'object', 'object', 'object', 'object', 'object', 'object', 'object'))
            os.chdir(f"{target_dir}/{observation}/rgs/divided_spectra")
            model_list = ['const*tbabs*zpowerlw', 'const*tbabs*zlogpar']

        
            for i in range(1, n_intervals+1):

                #Load RGS1 + RGS2 data
                xspec.AllData(f"1:1 sourcespecRGS1_gti{i}.fits 2:2 sourcespecRGS2_gti{i}.fits")

                spectrum1 = xspec.AllData(1)
                spectrum1.background = f"bgspecRGS1_gti{i}.fits"
                spectrum1.response = f"../{pairs_respli[0][0]}"

                spectrum2 = xspec.AllData(2)
                spectrum2.background = f"bgspecRGS2_gti{i}.fits"
                spectrum2.response = f"../{pairs_respli[0][1]}"
                
                xspec.AllData.ignore("bad")
                xspec.AllData.ignore('**-0.331 2.001-**')
                
                for model in model_list:
                    m1 = xspec.Model(model)
                    m2 = xspec.AllModels(2) #Retrieve the model object assigned to data group 2

                    #Freeze constant of RGS1 to 1 to account for calibration
                    m1.constant.factor = 1
                    m1.constant.factor.frozen = True

                    #Freeze redshift for both RGS1 and RGS2
                    if m1.expression=='constant*TBabs*zpowerlw':
                        m1.zpowerlw.Redshift = target_redshift
                        m1.zpowerlw.Redshift.frozen = True
                    
                        m2.zpowerlw.Redshift = target_redshift
                        m2.zpowerlw.Redshift.frozen = True
                        
                    if m1.expression=='constant*TBabs*zlogpar':
                        m1.zlogpar.Redshift = target_redshift
                        m1.zlogpar.Redshift.frozen = True
                    
                        m2.zlogpar.Redshift = target_redshift
                        m2.zlogpar.Redshift.frozen = True
                
                    #Perform Fit
                    xspec.Fit.statMethod = "cstat"     
                    xspec.Fit.nIterations = 100
                    xspec.Fit.criticalDelta = 1e-1
                    xspec.Fit.query = 'yes'
                    xspec.Fit.perform() 

                    #Plotting
                    xspec.Plot.splashPage = False   # XSPEC version and build data information will not be printed to the screen when the first plot window is initially opened
                    xspec.Plot.device = f'{observation}_{model}_{i}.gif/gif'
                    xspec.Plot.xAxis = 'keV'
                    xspec.Plot.setRebin(minSig=3, maxBins=4096) 
                    xspec.Plot.addCommand(f'label top {observation} {model} Part {i}')
                    xspec.Plot('ldata res model') 
                    os.rename(f"{target_dir}/{observation}/rgs/divided_spectra/{observation}_{model}_{i}.gif" , f"{target_dir}/Products/RGS_Spectra/{observation}/{observation}_{model}_{i}.gif")
                    
                    #use matplotlib 
                    chans1 = xspec.Plot.x(1)
                    chans1_err = xspec.Plot.xErr(1)
                    rates1 = xspec.Plot.y(1)
                    rates1_err = xspec.Plot.yErr(1)
                    chans2 = xspec.Plot.x(2)
                    chans2_err = xspec.Plot.xErr(2)
                    rates2 = xspec.Plot.y(2)
                    rates2_err = xspec.Plot.yErr(2)
                    
                    folded1 = xspec.Plot.model(1)
                    folded2 = xspec.Plot.model(2)
                    fig = plt.figure(figsize=(15,10))

                    ax1 = plt.subplot(211)
                    plt.errorbar(chans1, rates1, yerr=rates1_err, xerr=chans1_err, label='RGS1')
                    plt.errorbar(chans2, rates2, yerr=rates2_err, xerr=chans2_err, label='RGS2')
                    plt.yscale('log')
                    plt.xscale('log')
                    plt.ylim(1e-7, 150)
                    plt.title(f"{observation} {model} Part {i}", fontsize=30)
                    plt.xlabel('Energy [keV]', fontsize=17)
                    plt.ylabel('Normalized counts [s-1 keV-1]', fontsize=17)
                    ax1.legend(loc='lower right', fontsize='x-large')

                    ax2 = plt.subplot(212, sharex=ax1)
                    plt.title('Residuals', fontsize=30)
                    rates1_array = np.array(rates1)
                    rates2_array = np.array(rates2)
                    folded1_array = np.array(folded1)
                    folded2_array = np.array(folded2)
                    res1 = rates1_array - folded1_array
                    res2 = rates2_array - folded2_array
                    plt.errorbar(chans1, res1, yerr=rates1_err, linestyle='', label='RGS1')
                    plt.errorbar(chans2, res2, yerr=rates2_err, linestyle='', label='RGS2')
                    plt.hlines(0, plt.xlim()[0], plt.xlim()[1], color='m')
                    plt.xlabel('Energy [keV]', fontsize=17)
                    plt.ylabel('Normalized counts [s-1 keV-1]', fontsize=17)
                    ax2.legend(loc='lower right', fontsize='x-large')
                    plt.tight_layout(pad=4.0)
                    plt.savefig(f"{target_dir}/{observation}/rgs/divided_spectra/{observation}_{model}_{i}.png")
                    plt.close()
                    
                    #Calculate Flux and Luminosity and store their values 
                    xspec.AllModels.calcFlux('0, , err')
                    xspec.AllModels.calcLumin(f'0, , {target_redshift}, err') 
                    flux = spectrum1.flux[0] #erg/cm2/s
                    lumin = spectrum1.lumin[0] #e+44 erg/s
                    flux_range = f"{spectrum1.flux[1]} - {spectrum1.flux[2]}"
                    lumin_range = f"{spectrum1.lumin[1]} - {spectrum1.lumin[2]}"

                    #Store parameter results of fit
                    if m1.expression=='constant*TBabs*zpowerlw':
                        phoindex = m1(3).values[0]
                        phoindex_sigma = m1(3).sigma
                        alpha, beta, alpha_sigma, beta_sigma = ('', '', '', '')
                    if m1.expression=='constant*TBabs*zlogpar':
                        alpha = m1(3).values[0]
                        beta = m1(4).values[0]
                        alpha_sigma = m1(3).sigma
                        beta_sigma = m1(4).sigma
                        phoindex, phoindex_sigma = ('', '')
                    
                    fit_statistic = xspec.Fit.statistic

                    #Save output table
                    spectra_table.add_row((observation,"RGS1+RGS2", i, m1.expression, phoindex, phoindex_sigma, alpha, alpha_sigma, beta, beta_sigma, flux, flux_range, lumin, lumin_range, fit_statistic))
                    ascii.write(table=spectra_table, output=f'spectra{observation}_table.csv', format='csv', overwrite=True)
                    ascii.write(table=spectra_table, output=f'{target_dir}/Products/RGS_Spectra/{observation}/{observation}_table.csv', format='csv', overwrite=True)

            # Close XSPEC's currently opened log file.
            xspec.Xset.closeLog()
            
