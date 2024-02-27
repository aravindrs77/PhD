import sys
import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import itertools
import numpy as np
import pandas as pd
from matplotlib import rc
#plt.rcParams.update({
 #   'text.latex.preamble': r"\usepackage{amsmath}",
  #  "text.usetex": True, 
   # "figure.figsize": (10.4, 8.4), #10.4,8.4   8.4, 6.4
   # "font.size": 22, #28#16
    #"font.family": 'serif'}) #,
    #"font.serif": ["Computer Modern Roman"]})
from mpl_toolkits.mplot3d import Axes3D

print("test1")

#L1000 CDR Geometry
cal_x = [96.7013, 96.4898, -96.9126, -96.2279] #mm
cal_y = [121.136, -121.304, -120.967, 121.472]
cal_z = [402, 524, 624, 724]

# Simulated Calibration source
x_cal = 96.7013
y_cal = 121.136 
z_cal = 624

size = 3500 #1000 #290 #2000 #3000
nbins = 500 #0   150   #125   #4000,2000     150

detrad = 40 #38#37 #detector radius

#SOURCEDF
r = 150 #cal source Radii location
#angs = [90, 330, 210] 
angs = [90, 0, 270, 180] 
i = 1
sources = []
for ang in angs:
    x = r * np.sin(ang * 2*np.pi/360) ## WRONG WAY AROUND BECAUSE OF HOW RALPH IMPLEMENTED IT
    y = r * np.cos(ang * 2*np.pi/360)
    sources.append([i, x, y])
    #if i == 6:
    #    circ = plt.Circle((x, y), 3, color = 'red')
    #    ax.add_artist(circ)
    #ax.text(x, y, i, ha = 'center', va = 'center', fontweight = 'bold', fontsize = 16, color = 'green')
    i+= 1
sourcedf = pd.DataFrame(sources, columns = ['sourceid', 'xc', 'yc']).set_index('sourceid')

df = pd.concat([pd.read_csv(depfile, header = 1) for depfile in ['/global/homes/a/arreme/legend-swdev-scripts/python/testpos_orig.dat']]) #Ge77_0V   Ar41_0eV    #Ga76: '/disk/gfs_atp/yannmu/LGND/n-sims/deps/14String_n_in-det_20_MeV.dat'  #'/disk/gfs_atp/yannmu/LGND/n-sims/deps/n_all/14String_n_CalSource_all.dat'   /disk/gfs_atp/yannmu/LGND/n-sims/deps/Ar41/14String_Ar41_0eV.dat   /disk/gfs_atp/yannmu/LGND/n-sims/deps/Ge77/14String_Ge77_0eV_all.dat   /disk/gfs_atp/yannmu/LGND/n-sims/deps/n_all/LGND_1T_CDR_LNGS_n_CalSource_all_2.dat
df_save_str = 'L200_Th228_cal'  #   'Ge77_L1000' #'n_all_2'# Ar41_L1000 Ge77_L1000'   #'ncapture_bins'    'Ar41'    'Ge77'     #only pos_dist with electrons, 250 keV, all steps written Ge75      
#df = pd.concat([pd.read_csv(depfile, header = 1) for depfile in ['14Stringrealistic_CalSource_42_Position_-0360mm-morestats.dat']])
df['keV'] = df['edep']*1000
#df['x0'] = abs(df['x0'])
elines = list(set([round(x*1e5)/1e5 for x in df.groupby('keV').count().sort_values('edep').query('edep > 50').sort_index().index.tolist()])) 

#df = df.sample(frac = .01)

#set DEP
#df['evtype'] = df.keV.apply(lambda x: 'DEP' if abs(x-1592.50653)<1e-5 else 'SEP' if abs(x-2103.50544)<1e-5 else 'FEP' if abs(x - 2614.504353969453) < .00000000001 else 'other')
#df = df.groupby('ien').filter(lambda g: len(g) == 1) #NO MULTIPLICITY CUT

#smear by resolution (from cal GSTR)
a = 7.07e-1
b = 4.29e-4
def get_sigma(E):
    sigma = np.sqrt(a + b*E)

    #if E 511 or SEP... add broadening
    if abs(E-2103.50544) < 1e-5:
        sigma += 0.9
    if abs(E-510.99891) < 1e-5:
        sigma = np.sqrt(get_sigma(2103.50544)**2 - get_sigma(1592.50653)**2) #because sigma(DEP)**2 = sigma(SEP)**2 + sigma(broadened 511)**2

    return sigma
df['keVsmear'] = df.keV.apply(lambda E: E + np.random.normal(loc = 0, scale = get_sigma(E)))

#energy threshold
#df = df.query('keVsmear > 400 and keVsmear<5000')

#set SSE
#df['SSE'] = df.extent < 1

print(df)
print(df.head())
print(df.tail())
#df = df.query(" vol != 'Detector' ")    #!=    Detector means LAr!!!
#df = df[df['vol'].str.contains("ActiveDet")]
print(df.head())
print(df.tail())
print(df)


fig, ax = plt.subplots()  #figsize = (5,5))
ax.set_aspect('equal')
stringdata = []
rs = [235] #stringRadii
angs = [360./14] #det angles
i = 1
#ax.plot(df['x0'], df['y0'], linestyle = '', marker = '.', alpha = .2)
h = ax.hist2d(df['x0'], df['y0'], bins = 500, cmin = 1)#, norm = matplotlib.colors.LogNorm())#vmin = 10, vmax = 1000)) #cmap = 'Reds' bins = 500    int(nbins / 2 + 200)
fig.colorbar(h[3], ax=ax)
for r, ang in zip(rs, angs): #detectors
    #circ = plt.Circle((0, 0), r, fill = None, ec = 'r')
    #ax.add_artist(circ)
    for a in np.arange(0, 360, ang):
        x = r * np.cos(a * 2*np.pi/360)
        y = r * np.sin(a * 2*np.pi/360)
        #ax.text(x, y, i, ha = 'center', va = 'center', fontsize = 32, color = 'white') #'black'    darkorange
        circ = plt.Circle((x, y), detrad, fill = None, ec = 'k', linewidth = .1)
        #ax.add_artist(circ)
        stringdata.append([i, x, y])
        i+= 1
for i in [1, 2, 3, 4]:
  c = 'k'
  if i == 2: c = 'r'
  #circ = plt.Circle(tuple(sourcedf.loc[i]), 3, color = c)
  #ax.text(tuple(sourcedf.loc[i])[0]+15, tuple(sourcedf.loc[i])[1]+15, i, ha = 'center', va = 'center', fontsize = 10, color = c)
  #ax.add_artist(circ)
for cx, cy in zip(cal_x, cal_y):
    circ = plt.Circle((cx, cy), radius = 20, color = 'darkorange')
    ax.add_artist(circ)

circ = plt.Circle((x_cal, y_cal), radius = 20, color = 'red')
ax.add_artist(circ)
#ax.set_xlim([-235-detrad-20, 235+detrad+20])
#ax.set_ylim([-235-detrad-20, 235+detrad+20])
ax.set_xlabel('x [mm]')
ax.set_ylabel('y [mm]')
ax.set_ylim([-size,size])
ax.set_xlim([-size,size])
#plt.text(x = -550, y = -2500, s = 'x [mm]', fontsize = 28)     #-500, -2500
#plt.text(x = -1350, y = 100, s = 'y [mm]', fontsize = 28, rotation = 90)   #-2800, 200
plt.text(x = size * 1.4, y = 1000, s = 'No. of events', rotation = 90) #28    455, 100      #3000, 500
plt.savefig('../Plots/xy_s44_' + df_save_str + '.png')#, transparent = False, bbox_inches = 'tight')
plt.savefig('../Plots/xy_s44_' + df_save_str + '.pdf')#, transparent = False, bbox_inches = 'tight')
#plt.show()
plt.close()


#detector strings

#stringdf = pd.DataFrame(stringdata, columns = ['stringid', 'xs', 'ys']).set_index('stringid')
#df['String'] = df.vol.str.extract('StringAssembly_(.+)_Det').astype(int)
#df['localR'] = df.apply(lambda r: np.sqrt((r['x0'] - stringdf.loc[r['String']]['xs'])**2 + (r['y0'] - stringdf.loc[r['String']]['ys'])**2), axis = 1)


#fig.savefig('DEPgeo.pdf')

#dead layer

"""
#study distribution
desdf = df.groupby('evtype').localR.agg(['count', 'mean', 'std'])
desdf['meanerr'] = desdf.apply(lambda x: x['std']/np.sqrt(x['count']), axis = 1)
desdf['diff'] = desdf[['mean']]-desdf.loc['other']['mean']
desdf['differr'] = np.sqrt( desdf[['meanerr']]**2 + desdf.loc['other']['meanerr'] ) 
print(desdf[['mean', 'meanerr', 'diff', 'differr']])

nbins = 10
dl = 1
df['active'] = df.apply(lambda x: x.localR < detrad-dl, axis = 1)
sumdf = pd.DataFrame()
sumdf['total'] = df.groupby('evtype').ien.count()
sumdf['active'] = df.groupby('evtype').active.sum()
sumdf['frac'] = sumdf.apply(lambda x: 1-(x.active/x.total), axis = 1)
sumdf['fracerr'] = sumdf.apply(lambda x: x.frac*np.sqrt(1./x.total + 1./(x.total-x.active)), axis = 1)
print(sumdf[['frac', 'fracerr']]* 100)
"""

"""
#rdist figure
fig, ax = plt.subplots()
#plt.hist([v.tolist() for a, v in df.groupby('evtype').localR][::-1], bins = nbins, label = ['other', 'SEP', 'DEP', 'FEP'], normed = True)
plt.hist([{a: v.tolist() for a, v in df.groupby('evtype').localR}[l] for l in ['other', 'DEP', 'SEP', 'FEP']][::-1], bins = nbins, label = ['other', 'DEP', 'SEP', 'FEP'][::-1], density = True)
ax.legend()
ax.set_xlabel('Radius [mm]')
ax.set_ylabel('a.u.')
#ax.set_title('all dets')
fig.savefig('rdist.pdf')
"""

sourcepos =  sourcedf.loc[2, ['xc', 'yc']].tolist()
#sourcepos.append(400+ -360) #-360
sourcepos.append(-360) #-360
xc, yc, zc = sourcepos
wo = 3.2
ho = 7.95*2
wi = 2
hi = 4
#x/z projection
fig, ax = plt.subplots()
tempdf = df
#tempdf['x0'] = abs(tempdf['x0'])
#ax.plot(tempdf['x0'], tempdf['z0'], marker = '.', linestyle = '', alpha = .2)
h = ax.hist2d(tempdf['y0'], tempdf['z0'], bins = nbins, cmin = 1)#, norm=matplotlib.colors.LogNorm())#vmin = 10, vmax = 1000)) #cmap = 'Reds', bins = 500   nbin * 2
#tempdf.plot('x0', 'z0', ax = ax, marker = '.', linestyle = '')
#ax.add_patch(patches.Rectangle((yc-wo/2, zc-ho/2+4.15), wo, ho, color = 'red', zorder = 9))
#ax.add_patch(patches.Rectangle((yc-wi/2, zc-hi/2), wi, hi, color = 'yellow', zorder = 10))
for cy in cal_y:
    for cz in cal_z:
        circ = plt.Circle((cy, cz), radius = 10, color = 'darkorange')
        ax.add_artist(circ)
ax.set_xlabel('y [mm]')
ax.set_ylabel('z [mm]')
for cz in cal_z:
    circ_cal = plt.Circle((y_cal, cz), radius = 10, color = 'red')
    ax.add_artist(circ_cal)
#ax.set_ylim([-120, 420])
#ax.set_title('hits with y > 100')
#ax.set_ylim([-600, 600])
#ax.set_yticks(np.arange(-700, 101, step = 100))
#ax.set_xlim([-1000,1000])
#ax.set_ylim([-1000,1000])
ax.set_xlim([-size,size])
fig.colorbar(h[3], ax=ax)
plt.text(x = size * 1.4, y = 2500, s = 'No. of events', rotation = 90)  #445    1500
plt.savefig('../Plots/yz_s44_' + df_save_str + '.png', transparent = False, bbox_inches = 'tight')
plt.savefig('../Plots/yz_s44_' + df_save_str + '.pdf', transparent = False, bbox_inches = 'tight')
#plt.show()
plt.close()


sourcepos =  sourcedf.loc[2, ['xc', 'yc']].tolist()
#sourcepos.append(400+ -360) #-360
sourcepos.append(-360) #-360
xc, yc, zc = sourcepos
wo = 3.2
ho = 7.95*2
wi = 2
hi = 4
#x/z projection
fig, ax = plt.subplots()
tempdf = df
#tempdf['y0'] = abs(tempdf['y0'])
#ax.plot(tempdf['x0'], tempdf['z0'], marker = '.', linestyle = '', alpha = .2)
h = ax.hist2d(tempdf['x0'], tempdf['z0'], bins = nbins, cmin = 1)#, norm=matplotlib.colors.LogNorm())#vmin = 10, vmax = 1000)) #cmap = 'Reds', bins = 500
#tempdf.plot('x0', 'z0', ax = ax, marker = '.', linestyle = '')
#ax.add_patch(patches.Rectangle((yc-wo/2, zc-ho/2+4.15), wo, ho, color = 'red', zorder = 9))
#ax.add_patch(patches.Rectangle((yc-wi/2, zc-hi/2), wi, hi, color = 'yellow', zorder = 10))
for cx in cal_x:
    for cz in cal_z:
        circ = plt.Circle((cx, cz), radius = 10, color = 'darkorange')
        ax.add_artist(circ)
ax.set_xlabel('x [mm]')
ax.set_ylabel('z [mm]')
for cz in cal_z:
    circ_cal = plt.Circle((x_cal, cz), radius = 10, color = 'red')
    ax.add_artist(circ_cal)
#plt.text(x = -400, y = -850, s = 'x [mm]', fontsize = 28)
#plt.text(x = -2820, y = -300, s = 'z [mm]', fontsize = 28, rotation = 90)
#ax.set_ylim([-120, 420])
#ax.set_title('hits with y > 100')
#ax.set_ylim([-750, 150])
#ax.set_yticks(np.arange(-700, 101, step = 100))
#ax.set_ylim([-600,600])
ax.set_xlim([-size,size])
fig.colorbar(h[3], ax=ax)
plt.text(x = size * 1.4, y = 2500, s = 'No. of events', rotation = 90, fontsize = 28) #445
plt.savefig('../Plots/xz_s44_' + df_save_str + '.png')#, transparent = False, bbox_inches = 'tight')
plt.savefig('../Plots/xz_s44_' + df_save_str + '.pdf')#, transparent = False, bbox_inches = 'tight')
#plt.show()
plt.close()

#sys.exit()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
tempdf = df
ax.scatter(tempdf['x0'], tempdf['y0'], tempdf['z0'], marker = '.')
ax.grid(False)
ax.set_xticks([-4000, 4000]) #[-750, 0, 750])
ax.set_yticks([-5000, 5000]) #[-750, 0, 750])
ax.set_zticks([-4000, 4000]) #[-500, 0, 500])
ax.set_xlabel('x [mm]')
ax.set_ylabel('y [mm]')
ax.set_zlabel('z [mm]')
#ax.set_xlim([-235-detrad-20, 235+detrad+20])
#ax.set_ylim([-235-detrad-20, 235+detrad+20])
#ax.set_zlim([-750, 150])
fig.savefig('../Plots/xyz_s44_' + df_save_str + '.png', transparent = False)
fig.savefig('../Plots/xyz_s44_' + df_save_str + '.pdf', transparent = False)
#plt.show()
plt.close()
