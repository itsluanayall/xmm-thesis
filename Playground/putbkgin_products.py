import os

pwd = os.getcwd()
os.chdir(pwd+'/Markarian421')

for directory in os.listdir():
    if directory != "Products":
        os.chdir(pwd+'/Markarian421/'+ directory+ '/rgs')
        for subdir in os.listdir():
            if subdir.endswith('_check_rates.fit.ps'):
                os.rename(pwd+'/Markarian421/'+ directory+ '/rgs/'+ subdir, pwd+'/Markarian421/'+ 'Products/'+'Backgrnd_LC/'+subdirq)
