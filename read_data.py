# Modified to read physics data and plot for a full run or period (just change the path corresponding to the run or period) 
# Orginal code by Yannick to plot a single .lh5 file
import os, glob
import sys, h5py
import numpy as np
import math
import time
import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import csv
import pygama
#from pygama import lh5
from lgdo import LH5Store   #'from pygama.lgdo' changed to 'from lgdo'
from scipy.optimize import curve_fit
from scipy.stats import crystalball
from scipy.special import erfc
from scipy import integrate
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import json
from scipy.optimize import fsolve
from datetime import datetime
date = datetime.now().strftime("%Y-%m-%d")
from matplotlib import rc
plt.rcParams.update({
    'text.latex.preamble': r"\usepackage{amsmath}",
    "text.usetex": False,
    "figure.figsize": (10.4, 8.4),
    "font.size": 24, #16
    "font.family": 'serif'}) #,
    #"font.serif": ["Computer Modern Roman"]})


filter_fct = 'trapEmax_ctc_cal'  #'cuspEmax_ctc_cal' #cuspEmax_ctc_cal
#run = '003' #'027'
#tier = 'dsp' #'hit' 'dsp'
EE = np.array([])

#ch = 'ch1104000'
#det_name = ['ch1104000', 'ch1104001', 'ch1104002', 'ch1104003', 'ch1104004', 'ch1104005', 'ch1105600', 'ch1105602', 'ch1105603']
f = open("detectors.csv")
det_name = []
dir_name = []

for row in csv.reader(f):
    if row[0] == 'dets':
        continue
    det_name.append(str(row[0]))

x = 0

for row in csv.reader(f):
    dir_name.append(str(row[0]))
#path = '/global/cfs/projectdirs/legend/data/lngs/l200/public/prodenv/prod-blind/tmp/auto/generated/tier/hit/phy/p04/r002/'

#path = '/global/homes/a/arreme/training_data/l200/phy/p04/r003/r003_dsp/l200-p04-r003-phy-20230501T205951Z-tier_dsp.lh5'     #'/disk/gfs_atp/legend/data/l60/phy/r' + run + '/r' + run + '_' + tier + '/'

path = '/global/cfs/projectdirs/legend/data/lngs/l200/public/prodenv/prod-blind/tmp/auto/generated/tier/hit/phy/p04/r000/'

#for filename in glob.glob(os.path.join(path, 'l200-p04-r002-phy-*-tier_hit.lh5')): 
    #with open(os.path.join(os.getcwd(), filename), 'r') as f:    


for subdir, dirs, files in os.walk(path):
    for file in files:
        with open(os.path.join(subdir, file)) as f:  
#with open(path, 'r') as f:   
        #print(f)
            print(os.path.join(subdir, file))
            store = LH5Store() #lh5.Store()

            #get detector list. assumes all files have same channel map!
            with h5py.File(f.name, 'r') as ff:
                det_list = list(ff.keys())
                #print(ff.keys()) #print(ff.keys())
            det_list = list(filter(lambda x: x.startswith('c'), det_list)) #g
            #print(det_list)
            for det in det_list:
                if det not in det_name:
                    continue
                print(det)

                #pull out E array
                energy_name = filter_fct #cusp, trap, zac + tier + 
                #print(filter_fct)
                e_lh5, n_E = store.read_object(det + '/hit/' + energy_name, f.name)  #store.read_object(energy_name, f.name) line49
                try:
                    EE_0 = e_lh5.nda
                except(AttributeError):
                    continue
                print(len(EE_0))
                if n_E != len(EE_0):
                    print(f'got n_E ({n_E}) != len(EE_0) ({len(EE_0)})')
                    continue
                EE = np.append(EE, EE_0)


print(len(EE))
xmin = 360
xmax = 5500
#EE = EE[np.where(EE < xmax)]
#EE = EE[np.where(xmin < EE)]
print(len(EE))


binwidth = 10
plt.figure(0, figsize = (18.4, 6.4)) #14.4
(counts, bins, patches) = plt.hist(EE, bins = np.arange(int(xmin) - binwidth / 2., int(xmax) + 1 + binwidth / 2., binwidth), color = "midnightblue", histtype='step', rwidth = 0.9, label = 'Histogram data'); #bins centered around integer energy value
plt.yscale("log")
#plt.xscale("log")
plt.text(410, 160, r'$^{39}$Ar', fontsize = 28)
plt.text(900, 100, r'$2\nu\beta\beta$', fontsize = 28)
plt.text(1350, 115, r'$^{40}$K', fontsize = 28)
plt.text(1550, 260, r'$^{42}$K', fontsize = 28)
plt.text(2590, 13, r'$^{208}$Tl', fontsize = 28)
#plt.yticks([10**0, 10**1])
plt.ylabel("Events / (10 keV)") #yr
plt.xlim(xmin, xmax)
#plt.xlim(int(en - window / 2) - 2 * binwidth, int(en + window / 2) + 2 * binwidth)
#plt.xticks([100, 200, 300, 400, 500, 600, 700, 800, 900, 1000])
plt.xlabel("Energy [keV]")
plt.axvline(x = 2038.6, color = "darkred", linestyle = '--', alpha = 1, label = '$Q_{\\beta\\beta}$')
plt.legend(frameon = False)
plt.savefig("/global/homes/a/arreme/test_python_files/L200_hist_p04_r000_full.png", bbox_inches = 'tight')
plt.savefig("/global/homes/a/arreme/test_python_files/L200_hist_p04_r000_full.pdf", bbox_inches = 'tight')
plt.show()
plt.close()
