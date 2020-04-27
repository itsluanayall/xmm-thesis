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
    #setupSAS(sas_dir=sas_dir, ccf_dir=ccf_dir)

    observation = '0099280301'
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
    logFile = xspec.Xset.openLog("divided_spectra/XSPECLogFile.txt") #Create and open a log file for XSPEC
    logFile = xspec.Xset.log
    spectra_table = Table(names=('ObsId', 'Instrument', 'Piece', 'Model', 'Parameters', 'Flux[erg/cm2/s]', 'Luminosity[e+44erg/s]'), dtype=('object', 'object', 'object', 'object', 'object', 'object', 'object'))
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
            
            #Calculate Flux and Luminosity and store their values 
            xspec.AllModels.calcFlux('0, , err')
            xspec.AllModels.calcLumin(f'0, , {target_redshift}, err') 
            flux = spectrum1.flux[0] #erg/cm2/s
            lumin = spectrum1.lumin[0] #e+44 erg/s

            #Store parameter results of fit
            if m1.expression=='constant*TBabs*zpowerlw':
                phoindex = m1(3).values[0]
                phoindex_sigma = m1(3).sigma
                parameter_dict = {"PhoIndex": phoindex, "PhoIndex_sigma": phoindex_sigma}
                print(parameter_dict)
            if m1.expression=='constant*TBabs*zlogpar':
                alpha = m1(3).values[0]
                beta = m1(4).values[0]
                alpha_sigma = m1(3).sigma
                beta_sigma = m1(4).sigma
                parameter_dict = {"Alpha": alpha, "Alpha_sigma": alpha_sigma, "Beta": beta, "Beta_sigma": beta_sigma}
                print(parameter_dict)

            spectra_table.add_row((observation,"RGS1+RGS2", i, m1.expression, str(parameter_dict), flux, lumin))
            ascii.write(table=spectra_table, output=f'spectra{observation}_table.csv', format='csv', overwrite=True)
            
    # Close XSPEC's currently opened log file.
    xspec.Xset.closeLog()
    
