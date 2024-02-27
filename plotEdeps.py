import os
import sys
import glob
import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.optimize import curve_fit
from scipy import stats
#from scipy.signal import find_peaks
from matplotlib import rc
#plt.rcParams.update({
 #   'text.latex.preamble': r"\usepackage{amsmath}",
  #  "text.usetex": True, 
   # "figure.figsize": (10.4, 8.4),  #10.4,8.4
    #"font.size": 28,#30, #16
    #"font.family": 'serif'}) #,
    #"font.serif": ["Computer Modern Roman"]})
import pandas as pd
import ROOT


print(1)
#sys.exit()#takes ages like that
df = pd.read_csv('/global/homes/a/arreme/legend-swdev-scripts/python/testpos_orig.dat', header = 1, engine = 'python') #pd.read_csv('/disk/gfs_atp/yannmu/LGND/n-sims/deps/Ge77/LGND_1T_CDR_LNGS_Ge77_0eV.dat', header = 1, engine = 'python') #pd.read_csv('/disk/gfs_atp/yannmu/LGND/n-sims/deps/Ge77/14String_Ge77_0eV_all.dat', header = 1, engine = 'python') #pd.read_csv('/disk/gfs_atp/yannmu/LGND/n-sims/deps/calsims_16_positions_combined/14String_CalSource_16Positions_10stat.dat', header = 1, engine = 'python') #/disk/gfs_atp/yannmu/LGND/n-sims/deps/14String_n_in-det_20_MeV.dat for Ga76 /disk/gfs_atp/yannmu/LGND/n-sims/deps/Ar41/14String_Ar41_0eV.dat for Ar41 /disk/gfs_atp/yannmu/LGND/n-sims/deps/Ge77/14String_Ge77_0eV_all.dat for Ge77      Ar41/LGND_1T_CDR_LNGS_Ar41_0eV
df_save_str = 'Th228_L000_cal' #'DEP' Ar41_L1000  Ge77_L1000' #'Ge77_cuts' #DEP, Ar41, Ge77, Ge75_M2, Ga76 , ...Ge77_smear (with extend from DEP M1 analysis) 10stat
df['keV'] = df['edep']*1000
#df1 = pd.read_csv('/disk/gfs_atp/yannmu/LGND/n-sims/deps/calsims_16_positions_combined/14String_CalSource_16Positions_M1_10stat.dat', header = 0, engine = 'python')
#df1['keV'] = df1['edep']*1000
print(1)
print(df.head())
print(df.tail())
df = df.query(" vol != 'Detector' ")
df = df[df['vol'].str.contains("ActiveDet")]
df['keVround'] = df['keV'].apply(lambda x: round(x*10.)/10)
print(df.head())
print(df.tail())
df1 = df.groupby('ien').filter(lambda g: len(g) == 1) #only multiplicity cut for individual file, not for merged files!
print(1)
df2 = df.query('extent <= 1') #2.9812188602699186 1.0069959265369008
print(1)
df3 = df1.query('extent <= 1') #2.9812188602699186 1.0069959265369008
print(1)
#dfM2 = df.groupby('ien').filter(lambda g: len(g) == 2) 



#linedf = df.groupby('keV').count()
#linedf = linedf.loc[linedf.ien > linedf.ien.quantile(.75)]
#elines = linedf.index.tolist()
#elines=  []


ROI = 2050 #2050 265 #center of ROI
ROI_width = 100 #100 3.10442 5 #half of ROI window


#Ge75: M2 265 +- 5 [3.10442] keV: 4374 [4147] out of 10**6

#sample = 1
#df = df.sample(frac = sample)


#smear by GERDAII+ IC resolution (from cal GSTR)
#1 sigma at DEP: 1.09929, FWHM at DEP: 2.58863, 0.5 fwhm 1.294315, 1.5 fwhm 3.882945, 2 fwhm 5.17726
a = 0.280
b = 5.83e-4
def get_sigma(E):
    sigma = np.sqrt(a + b*E)

    #if E 511 or SEP... add broadening
    if abs(E-2103.50544) < 1e-5:
        sigma += 0.9
    if abs(E-510.99891) < 1e-5:
        sigma = np.sqrt(get_sigma(2103.50544)**2 - get_sigma(1592.50653)**2) #because sigma(DEP)**2 = sigma(SEP)**2 + sigma(broadened 511)**2
    return sigma

ft = ROOT.TF1("tail", "([1] /(2.*[2])) *exp ( x/ [2] + [0]*[0]/ 2. / ([2]*[2]) ) * erfc ( x / sqrt(2.) / [0]  + [0] /sqrt(2.) / [2] )", -60, 30)
ft.SetParameter(0, get_sigma(2614.)) #sigma at 2.6 MeV
ft.SetParameter(1, 1) #tail amplitude, just has to be non-zero
ft.SetParameter(2, 2.66) #average beta of ICs at 2.6 MeV
def get_smear(E):
    if np.random.random() < 1.01: #0.97 include or neglect tail
        return np.random.normal(loc = 0, scale = get_sigma(E))
    return ft.GetRandom()

df['keVsmear'] = df.keV.apply(lambda E: E + get_smear(E)) #REMOVE ZERO IF SMEARING NEEDED!!!
df1['keVsmear'] = df1.keV.apply(lambda E: E + get_smear(E))
print(1)
#energy threshold
df = df.query('200 < keVsmear < 5000') #5000 keVsmear
df1 = df1.query('200 < keVsmear < 5000') 
#df = df.query('keVsmear > 2575 and keVsmear<2640')
print(1)

df1['keVsmear'] = df1.keV.apply(lambda E: E + get_smear(E))
df1 = df1.query('200 < keVsmear < 5000')
df2['keVsmear'] = df2.keV.apply(lambda E: E + get_smear(E))
df2 = df2.query('200 < keVsmear < 5000')
df3['keVsmear'] = df3.keV.apply(lambda E: E + get_smear(E))
df3 = df3.query('200 < keVsmear < 5000')
#dfM2['keVsmear'] = dfM2.keV.apply(lambda E: E + get_smear(E))
#dfM2 = dfM2.query('200 < keVsmear < 5000')

#df before cuts, df1 after AC, df2 after PSD (DEP), df3 after AC & PSD (DEP)
roicount = df.loc[abs(df['keVsmear'] - ROI) <= ROI_width].groupby('vol').count().ien.sum()  #2050, 100
print(roicount)
roicount1 = df1.loc[abs(df1['keVsmear'] - ROI) <= ROI_width].groupby('vol').count().ien.sum()  
print(roicount1)
roicount2 = df2.loc[abs(df2['keVsmear'] - ROI) <= ROI_width].groupby('vol').count().ien.sum()
print(roicount2)
roicount3 = df3.loc[abs(df3['keVsmear'] - ROI) <= ROI_width].groupby('vol').count().ien.sum()
print(roicount3)
#no smear: Ge77: 434626, 346822, 13210, 8969, Ar41: 37, 34, 0, 0 --> https://aip.scitation.org/doi/pdf/10.1063/1.533433 upper limit 90% CL for b = 0, n = 0 is 2.62 (alternative constructions 2.3, 2.44)
#with smearing applied: Ge77: 434686, 346800, 3303, 2169, Ar41: 37, 33, 0, 0, (0)
#roicountM2 = dfM2.loc[abs(dfM2['keVsmear'] - ROI_width) <= ROI_width].groupby('vol').count().ien.sum()  
#print(roicountM2)


"""
#determine extent cut #plus minus x * sigma (1.09929), within 1 fwhm (2.58863): 1.294315
fig, ax3 = plt.subplots()
dfall = df.query('abs(keVsmear - 1592.5065339700002) > 5.17726').extent
dfall.hist(bins = 1500, ax = ax3, density = True, alpha = .99, label = 'other', grid = None, color = 'darkorange')#'#ff7f0e')  #<20
dfdep = df1.query('abs(keVsmear - 1592.5065339700002) <= 5.17726').extent
dfdep.hist(bins = 1500, ax = ax3, density = True, alpha = .99, label = 'DEP, M1, +-2FWHM', grid = None, color = 'darkblue') 
ax3.set_xlabel('Extent [mm]')
ax3.set_ylabel('Normalised frequency')
ax3.set_xlim(0, 15)
ax3.set_ylim(0, 1)
#ax3.axvline(x = 2.9812188602699186, linestyle = '--', color = 'darkred', label = '$r_{90}$') #2.9812188602699186
ax3.legend(frameon = False, loc = 'upper right')
fig.savefig('Plots/DEPextent_4fwhm_smear.png') #4fwhm_smear
fig.savefig('Plots/DEPextent_4fwhm_smear.pdf')     #DEPextent_4fwhm_smear
plt.show()
sys.exit()

fig, ax4 = plt.subplots()
#dfall = df.query('abs(keVsmear - 1592.5065339700002) > 5.17726').extent
dfall.hist(bins = 1500, ax = ax4, density = True, alpha = .99, label = 'other', grid = None, color = 'darkorange') #'#ff7f0e')
#dfdep = df.query('abs(keVsmear - 1592.5065339700002) <= 5.17726').extent
dfdep.hist(bins = 1500, ax = ax4, density = True, alpha = .99, label = 'DEP, M1, +-2FWHM', grid = None, color = 'darkblue')
ax4.set_xlabel('Extent [mm]')
ax4.set_ylabel('Normalised frequency')
ax4.set_xlim(0, 4)
#ax4.axvline(x = 2.9812188602699186, linestyle = '--', color = 'darkred', label = '$r_{90}$')
ax4.legend(frameon = False, loc = 'upper right')
#fig.savefig('Plots/DEPextentzoom_4fwhm_smear.png') #298_all
#fig.savefig('Plots/DEPextentzoom_4fwhm_smear.pdf')
plt.show()


#df['SSE'] = df.extent < 1

dep_mm = df1.query('abs(keVsmear - 1592.5065339700002) <= 5.17726').extent #0.01 1 * 1.09929
#other_mm = df.query('abs(keVsmear - 1592.5065339700002) > 5.17726').extent

print(str(dep_mm.quantile(.25)) + ' , ' + str(dep_mm.quantile(.5)) + ' , ' + str(dep_mm.quantile(.75)) + ' , ' + str(dep_mm.quantile(.68)) + ' , ' + str(dep_mm.quantile(.9)) + ' , ' + str(dep_mm.quantile(.95)) + ' , ' + str(dep_mm.quantile(.99)))
#df 0.01: 0.459319051248606 , 0.5888165789977128 , 0.7901888326968735 , 0.7160006376794539 , 2.9812188602699186 , 12.239739121153496 , 41.91429918431555
#df1 4fwhm smear: 0.6262186355421019 , 1.1583837567576525 , 21.731852314034434 , 14.261793731064142 , 42.7709332697818 , 54.250375247321145 , 71.54097645399506

#print(str(other_mm.quantile(.25)) + ' , ' + str(other_mm.quantile(.5)) + ' , ' + str(other_mm.quantile(.75)) + ' , ' + str(other_mm.quantile(.68)) + ' , ' + str(other_mm.quantile(.9)) + ' , ' + str(other_mm.quantile(.95)) + ' , ' + str(other_mm.quantile(.99)))
#df 0.01: 0.03843471810301332 , 1.9957440109389983 , 14.41532847385521 , 9.701077274242955 , 32.50967692733586 , 45.16050370650382 , 66.13286268451036
#df1 4fwhm smear: 1.1067270732070933 , 9.834650869946437 , 24.173092492068125 , 19.127733316281088 , 41.75836340958905 , 52.697503511442434 , 70.10146588167649
#df 4fwhm smear: 0.9244443370757032 , 9.49550548764303 , 24.14177135131984 , 18.995778236289333 , 42.018965880805816 , 53.10811474450616 , 70.76672797100419
print(1)

"""



bins = np.arange(*list(df.keVsmear.agg(['min', 'max'])), step = 5) #step = 1, 3, 5
nbins = len(bins)
fig, ax = plt.subplots()#figsize = (16.4, 8.4)) #figsize = (12.4, 6.4)) #(5, 3)) 6, 4
#plt.hist(df.query('vol == "DetectorArray_001_StringAssembly_004_DetUnit_005_ActiveDet_001"').keVsmear, bins = bins, stacked = True, log = True, histtype = 'step', color = 'red')
#ns, bins, _ = plt.hist(dfM2.keVsmear, bins = bins, stacked = True, log = True, histtype = 'step', color = "darkblue", label = "M2") #'#2f2fff')
ns, bins, _ = plt.hist(df.keVsmear, bins = bins, stacked = True, log = True, histtype = 'step', color = "darkblue", label = "before cuts") #'#2f2fff')   #before cuts
ns, bins, _ = plt.hist(df1.keVsmear, bins = bins, stacked = True, log = True, histtype = 'step', color = "darkgreen", label = "after AC") #'#2f2fff')
ns, bins, _ = plt.hist(df2.keVsmear, bins = bins, stacked = True, log = True, histtype = 'step', color = "darkorange", label = "after PSD") #'#2f2fff')
ns, bins, _ = plt.hist(df3.keVsmear, bins = bins, stacked = True, log = True, histtype = 'step', color = "darkred", label = "after AC+PSD ") #'#2f2fff')  after cuts
ax.axvspan(1950, 2150, alpha = 0.4, color = "lightblue")
ax.set_xlabel('Event energy [keV]')
ax.set_ylabel('Counts / 5 keV') #/3 keV
ax.legend(frameon = False, prop={'size': 12}, loc = 'upper right')#, prop={'size': 20})#, loc = 'lower center', prop={'size': 28})#'lower center')   'upper right'
#ax.set_xlim(0, 3000)
ax.set_yscale('log')
#ax.plot([2614], [4*10**5], marker = 'v', c= 'k')
#ax.annotate('FEP 2.6 MeV', [2700, 3.5*10**5], horizontalalignment = 'left') #10**4.8
#ax.plot([1592.50653], [2*10**4], marker = 'v', c= 'k')
#ax.annotate('DEP 1.6 MeV', [1592, 3*10**4], horizontalalignment = 'center')
#ax.plot([2103.5], [9*10**4], marker = 'v', c= 'k')
#ax.annotate('SEP 2.1 MeV', [2050, 1.8*10**5], horizontalalignment = 'center')
fig.savefig('../Plots/simspec_s44_' + df_save_str + '_smear.png', bbox_inches = 'tight')
fig.savefig('../Plots/simspec_s44_' + df_save_str + '_smear.pdf', bbox_inches = 'tight')
#plt.show()
plt.close()

#sys.exit()

bins = np.arange(1500, 2500, step = 5) #2590, 2650 1500 2500 240 290 #step=   0.1    .3    1572, 1612.5    0.5
nbins = len(bins)
#left, bottom, width, height = [0.75, 0.6, 0.2, 0.2]
fig, ax2 = plt.subplots()#figsize = (16.4, 8.4)) #figsize = (12.4, 6.4)) #(3, 3))
#df = df.query('keVsmear >= 1500 and keVsmear <= 2500') #2590 2650   1500, 2500     1572, 1612
#df1 = df1.query('keVsmear >= 1572 and keVsmear <= 1612') #2590 2650
#df2 = df2.query('keVsmear >= 1500 and keVsmear <= 2500') #2590 2650
#df3 = df3.query('keVsmear >= 1500 and keVsmear <= 2500') #2590 2650
#dfM2 = dfM2.query('keVsmear >= 1500 and keVsmear <= 2500') #2590 2650
#ax2.hist(df.query('vol == "DetectorArray_001_StringAssembly_004_DetUnit_005_ActiveDet_001"').keVsmear, bins = bins, stacked = True, log = True, histtype = 'step', color = 'red')
#ns, bins, _ = ax2.hist(dfM2.keVsmear, bins = bins, stacked = True, log = False, histtype = 'step', color = 'darkblue', label = "M2") #'#2f2fff')
ns, bins, _ = ax2.hist(df.keVsmear, bins = bins, stacked = True, log = False, histtype = 'step', color = 'darkblue', label = "before cuts") #'#2f2fff')
ns, bins, _ = ax2.hist(df1.keVsmear, bins = bins, stacked = True, log = False, histtype = 'step', color = 'darkgreen', label = "after AC") #'#2f2fff')
ns, bins, _ = ax2.hist(df2.keVsmear, bins = bins, stacked = True, log = False, histtype = 'step', color = 'darkorange', label = "after PSD") #'#2f2fff')
ns, bins, _ = ax2.hist(df3.keVsmear, bins = bins, stacked = True, log = False, histtype = 'step', color = 'darkred', label = "after AC+PSD") #'#2f2fff')   after    cuts
#ns, bins, _ = ax2.hist(df1.keVsmear, bins = bins, stacked = True, log = False, histtype = 'step', color = 'darkblue', label = "after AC")
ax2.axvspan(1950, 2150, alpha = 0.4, color = "lightblue")
ax2.set_xlabel('Event energy [keV]')
ax2.set_ylabel('Counts / 5 keV') #/0.3    0.5 keV
#ax2.set_yticks([0, 1, 2, 3, 4])
#ax2.set_ylim(0.1, 10**5) #0.5, 0.5 * 10**4
#ax2.set_yscale("log")
#ax2.set_xlim(1572, 1612)
#ax2.axvline(x = 1592.50653, linestyle = '--', color = 'darkred', label = 'DEP')
#ax2.axvline(x = 1592.50653 - 5.17726, linestyle = '--', color = 'darkorange', label = '$\pm 2$ FWHM')
#ax2.axvline(x = 1592.50653 + 5.17726, linestyle = '--', color = 'darkorange')
ax2.legend(frameon = False, prop={'size': 12}, loc = 'upper right')#, prop={'size': 20}) #
fig.savefig('../Plots/simspeczoom_s44_' + df_save_str + '_smear.png', bbox_inches = 'tight')
fig.savefig('../Plots/simspeczoom_s44_' + df_save_str + '_smear.pdf', bbox_inches = 'tight')
#plt.show()
plt.close()
print("Finished")




##'survival' fraction
#plotvals = [(x/y, (x/y)*np.sqrt(1/x+1/y)) for x, y in zip(*ns)]
#suv, suverr = zip(*plotvals)
#linx = [list(abs(bins - line)).index(min(abs(bins - line))) for line in elines]
#fig, ax = plt.subplots()
##for eline in elines:
##    ax.axvline(eline, c = 'k', alpha = .5, linestyle = '--')
#ax.errorbar(bins[0:-1], suv, yerr = suverr, marker = 'o', linestyle ='')
#ax.errorbar(bins[linx], np.array(suv)[linx], yerr = np.array(suverr)[linx], marker = 'o', linestyle ='', c = 'r')
#fig.savefig('suv.pdf')
#
#for vol, groupdf in df.groupby('vol'):
#  print(vol, len(groupdf))
#  fig, ax = plt.subplots()
#  groupdf.hist('keVsmear', bins = nbins, ax = ax, log = True)
#  ax.set_title(vol)
#  ax.set_xlabel('Event energy [keV]')
#  fig.savefig('plot/'+vol+'.pdf')
