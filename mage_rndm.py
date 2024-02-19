import numpy as np
import random
import csv
import pandas as pd
import sys
import os

# Writing header into output file
with open('/global/homes/a/arreme/python/testpos_orig.dat', 'w') as f: 
    f.write("nev 10000\n")
    f.write("%s,%s,%s,%s,%s,%s,%s\n" % ("ien", "vol", "edep", "extent", "x0", "y0", "z0"))

# Input File
path = '/global/homes/a/arreme/MaGe/deps/sim/'
file = 'LGND_200_Baseline_Calibration_wSourceNumber_44_Height_'

# Position of radioactive sources
pos = ['-96', '24', '124', '224']

# Activity of the sources corresponding to above positions
act = [1986.5, 1042.7, 119.5, 2651.7]
maxact = max(act)

en = []
vol = []
edep = []
extent = []
x0 = []
y0 = []
z0 = []
print("Started")
for i in range(len(pos)):
    with open(path + file + pos[i] + 'mm.dat', 'r') as f:
        reader = csv.reader(f)
        next(reader) #skip header row 1 in input file
        next(reader) #skip header row 2 in input file
        for row in reader:
            en.append(row[0])
            vol.append(row[1])
            edep.append(row[2])
            extent.append(row[3])
            x0.append(row[4])
            y0.append(row[5])
            z0.append(row[6])
        print(len(en)) #original length of ien

        #filter the unique event ids
        ien = []
        [ien.append(i) for i in en if i not in ien]
        #print(ien, len(ien))

        #random selection of events
        rndm = int(act[i]*len(ien)/maxact)
        #print(rndm)
        ien = random.sample(ien, rndm)  #selection is in a random order
        ien.sort()  #sort back to ascending 
        print(ien, len(ien))   
        c = 0
        with open('/global/homes/a/arreme/python/testpos_orig.dat', 'a') as f:
            for i in range(len(en)):
                if en[i] in ien:
                    c = c + 1
                    f.write("%s,%s,%s,%s,%s,%s,%s\n" % (en[i], vol[i], edep[i], extent[i], x0[i], y0[i], z0[i]))

        en.clear()
        vol.clear()
        edep.clear()
        extent.clear()
        x0.clear()
        y0.clear()
        z0.clear()
        ien.clear()
        c = 0

print("Finished")   