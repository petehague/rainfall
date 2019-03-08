#!/usr/bin/env python3

from astropy.io import fits
import sys
import os

filename = sys.argv[1]
inputfile = fits.open(filename)
outfilename = os.path.splitext(filename)[0] + ".txt"

print(filename, outfilename)

naxis1 = inputfile[0].header['NAXIS1']

outfile = open(outfilename, "w")
for ra in range(naxis1):
    for dec in range(naxis1):
        datum = inputfile[0].data[0,0,dec,ra]
        outfile.write("{}{}".format(datum, " " if dec<naxis1-1 else ""))
    outfile.write("\n")
outfile.close()
