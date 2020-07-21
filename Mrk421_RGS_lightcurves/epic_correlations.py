import numpy as np
import os

def cross_correlation(self, rates1, rates2, x, str1, str2):
    
    cc_simple = np.correlate(df_burst1, df_burst2, mode='same')
    plotcc = plt.plot(x, cc_simple, marker='o', markersize=4)
    plt.grid(alpha=0.5)
    plt.xlabel('Seconds since beginning of shifting')
    plt.ylabel('Cross-correlation')
    text = "Cross-correlation between " + str1 + "and" + str2
    plt.title(text, fontsize=15)
    
    return cc_simple