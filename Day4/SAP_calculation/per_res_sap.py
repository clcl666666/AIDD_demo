import sys, os
import pymol
from pymol import cmd
import gzip, bz2
import re
import string
import colorsys
import numpy as np

pdbname = cmd.get_object_list()[0]

def read_pdb( fname ):
    if os.path.exists( fname + ".pdb.gz" ):
        return gzip.open(fname + ".pdb.gz" )
    elif os.path.exists( fname + ".pdb" ):
        return open( fname + ".pdb" )
    elif os.path.exist( fname + ".pdb.bz2" ):
        return open( fname + ".pdb.bz2" )
    else:
        return False

f = read_pdb( pdbname )
if not f:
    print("No pdb files found, exit!")
    exit(0)

seqs = []
saps = []

for line in f:
    if not "per_res_sap" in line:
        continue

    sp = line.strip().split()
    seqpos = int(sp[0].split("_")[-1])
    sap = float(sp[1])

    seqs.append(seqpos)
    saps.append(sap)


min_sap = 0
max_sap = 3

def sap_to_color(sap):
    rgb = colorsys.hsv_to_rgb(0, np.interp(sap, [min_sap, max_sap], [0, 1]), 1)

    r = rgb[0] * 255
    g = rgb[1] * 255
    b = rgb[2] * 255

    color = "0x%02x%02x%02x"%(int(r), int(g), int(b))

    return color


for seqpos in range(1, max(seqs)+1):
    if ( seqpos not in seqs):
        continue
    idx = seqs.index(seqpos)
    sap = saps[idx]

    color = sap_to_color(sap)

    # colorCPK("resi %i"%seqpos, color
    cmd.color(color, "resi %i"%seqpos,)

    if ( sap > 0.5 ):
        cmd.label("resi %i and name CA"%seqpos, "\"%.1f\""%sap)
