
"""
Simple python script based on module PySpec with the scope to fit different models on the same dataset of exclusively RGS (1 and 2) data.
The output will be printed in .ps format in the working directory.
To choose between RGS1 and RGS2, run the python command with appended 'R1' or 'R2' respectively.
The results are printed in a log file in the working directory.
"""


import xspec
import matplotlib.pyplot as plt
import argparse 
import sys 


#Parser options
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('RGS', type=str, help='Type R1 or R2.')
args = parser.parse_args()

#load data 
try:
    if args.RGS == 'R1':
        logFile = xspec.Xset.openLog("RGS1LogFile.txt") #Create and open a log file for XSPEC output
        logFile = xspec.Xset.log  # Get the Python file object for the currently opened log
        s1 = xspec.Spectrum('P0658801601R1U002SRSPEC1001.FIT') 
        s1.background = 'P0658801601R1U002BGSPEC1001.FIT'
        s1.response = 'P0658801601R1U002RSPMAT1001.FIT'
        instr = 'RGS1'
    elif args.RGS == 'R2':
        logFile = xspec.Xset.openLog("RGS2LogFile.txt") 
        logFile = xspec.Xset.log  
        s1 = xspec.Spectrum('P0658801601R2U002SRSPEC1001.FIT') 
        s1.background = 'P0658801601R2U002BGSPEC1001.FIT'
        s1.response = 'P0658801601R2U002RSPMAT1001.FIT'
        instr = 'RGS2'
    else:
        raise ValueError
except ValueError:
    print('Please enter a valid instrument.')
    sys.exit(1)

#ignore bad data
s1.ignore('bad')
s1.ignore('1.9-**')  #change if necessary

#Define the models to fit on the data
model_list = ['powerlaw', 'logpar', 'phabs*logpar', 'phabs*powerlaw', 'tbabs*eplogpar', 'tbabs*powerlaw', 'wabs*zlogpar', 'wabs*logpar', 'phabs*cutoffpl', 'tbabs*srcut' ] #change if necessary

for model in model_list:

    #Define the model
    m = xspec.Model(model)
    if 'zlogpar' in m.componentNames:
         m.zlogpar.Redshift = 0.03 #otherwise it is freezed to 0.

    #Fit the model
    xspec.Fit.nIterations = 100
    xspec.Fit.criticalDelta = 1e-1
    xspec.Fit.perform()   
    
    #Calculate Flux and Luminosity
    xspec.AllModels.calcFlux('0, 1.9, err')
    xspec.AllModels.calcLumin(', , 0.03, err')  #note: 0.3 is the redshift for Mrk421

    #Plotting: set the device, set the axes and title and plot data
    xspec.Plot.device = f'{instr}_{model}.ps/cps'
    xspec.Plot.xAxis = 'keV'
    xspec.Plot.setRebin(minSig=3, maxBins=4096) 
    xspec.Plot.addCommand(f'label top {instr} {model}')
    xspec.Plot('ldata res')    

# Close XSPEC's currently opened log file.
xspec.Xset.closeLog()

    

    
