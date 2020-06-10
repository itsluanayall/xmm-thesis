import os
import matplotlib.pyplot as plt
from config import CONFIG
from astropy.table import Table
import seaborn as sns
import pandas as pd
import numpy as np

if __name__ == "__main__":
    
    target_dir = CONFIG['target_dir'] 
    products_dir = os.path.join(target_dir, "Products")
    if not os.path.isdir(os.path.join(target_dir, "Products", "Plots_spectra")):
        os.makedirs(os.path.join(target_dir, "Products", "Plots_spectra"))

    #Read nH data 
    hdul_spec = Table.read(os.path.join(products_dir, "RGS_Spectra", "spectra_table.fits"), hdu=1)
    data_spec = hdul_spec.to_pandas()
    data_spec = data_spec[data_spec['tbinid'] == 0]   #Use only average spectra

    data_spec_zlogp = data_spec[data_spec['model']=='constant*TBabs*zlogp']
    data_spec_zpowe = data_spec[data_spec['model']=='constant*TBabs*zpowe']
    data_distrib = pd.DataFrame({'nH_logpar': data_spec_zlogp['nH'].values, 'nH_powerlaw': data_spec_zpowe['nH'].values })
    nH_top_logpar = data_spec_zlogp['nH_up'].values - data_spec_zlogp['nH'].values
    nH_bot_logpar = data_spec_zlogp['nH'].values - data_spec_zlogp['nH_low'].values
    nH_top_powerlaw = data_spec_zpowe['nH_up'].values - data_spec_zpowe['nH'].values
    nH_bot_powerlaw = data_spec_zpowe['nH'].values - data_spec_zpowe['nH_low'].values


    g = (sns.jointplot("nH_logpar", "nH_powerlaw", data=data_distrib, space=0, color="skyblue", marker='.').plot_joint(sns.kdeplot,zorder=0, n_levels=6))
    plt.errorbar(data_distrib.nH_logpar, data_distrib.nH_powerlaw, linestyle='', ecolor='gray', alpha=0.4, xerr=(nH_bot_logpar, nH_top_logpar), yerr=(nH_bot_powerlaw, nH_top_powerlaw))
    
    #Plot maximum (for now doesn't work)
    '''
    x = g.lines[0].get_xdata() # Get the x data of the distribution
    y = g.lines[0].get_ydata() # Get the y data of the distribution
    maxid = np.argmax(y) # The id of the peak (maximum of y data)
    plt.plot(x[maxid],y[maxid], 'bo', ms=10)
    '''

    #Save and show
    plt.savefig(os.path.join(products_dir, "Plots_spectra", "nH_density.png"))
    plt.show()