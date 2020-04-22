import os
pwd = os.getcwd()
os.chdir(pwd+'/Markarian421')
for directory in os.listdir():
	os.chdir(pwd+'/Markarian421/'+ directory+ '/rgs')
	for subdir in os.listdir():
		#if subdir.endswith('RGS_rates.ds.png') or subdir.endswith('RGS_rates.ds'):
		os.remove(pwd+'/Markarian421/'+ directory+ '/rgs/'+ subdir)

