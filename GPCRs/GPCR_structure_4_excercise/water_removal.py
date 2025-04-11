from scipy.stats import scoreatpercentile
import numpy as np
import sys

# Works with gro files only
ifile = sys.argv[1]
ofile = sys.argv[2]

f = open(ifile, 'r')

lines = f.readlines()
atomlines = lines[2:-1]  # The lines that contain atom informations

# Figure out the positions of cholesterol oxygens
chol_positions = []
for i in atomlines:
    words = i.split()
    if words[0][-3:] == 'CHL' and words[1] == 'O1':
        chol_positions.append(float(words[-1]))

# Calculate the 5th and 95th percentile of the cholesterol O positions
l_bound = scoreatpercentile(chol_positions, 5)
u_bound = scoreatpercentile(chol_positions, 95)

# Figure out the positions of water molecules
water_positions = {}
for i in atomlines:
    words = i.split()
    if words[0][-3:] == 'SOL':
        name = words[0][:-3]
        if name not in water_positions:
            water_positions[name] = [float(words[-1])]
        else:
            water_positions[name].append(float(words[-1]))

# Identify all the water molecules between the boundaries for deletion
to_del = []
for k,v in water_positions.items():
    if np.any(v > l_bound) and np.any(v < u_bound):
        to_del.append(k)

# Actually delete them, adjust the number of atoms too (gro needs it to work)
n_atoms = int(lines[1])
reduced_atomlines = []
for i in atomlines:
    words = i.split()
    if words[0][-3:] == 'SOL' and words[0][:-3] in to_del:
        n_atoms -= 1
    else:
        reduced_atomlines.append(i)

# Report on water molecules to make adjustments of topology easy
print ('Initial number of water molecules:', len(water_positions))
print ('Number of water molecules to be deleted:', len(to_del))
print ('Final number of water molecules:', len(water_positions) - len(to_del)) 

# Write the adjusted lines into a new file
f2 = open(ofile, 'w')
f2.write(lines[0])
f2.write(str(n_atoms) + '\n')
f2.writelines(reduced_atomlines)
f2.write(lines[-1])
f2.close()




