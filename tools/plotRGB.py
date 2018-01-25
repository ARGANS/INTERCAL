#!/usr/bin/env python
# -*- coding: utf-8 -*-

from PIL import Image
import sys
from osgeo import gdal
import numpy as np
import re
import os

gdal.UseExceptions()


def usage():
    print 'Usage: plotRGB.py FILE where FILE are three valid files, may use a glob'


def get_data(filename):
    print filename
    ds = gdal.Open(filename)
    band = ds.GetRasterBand(1)
    data = band.ReadAsArray()
    return data


def parse_filename(filename):
    filename = filename.split('/')[-1]
    bandext = '_B[0-9]'
    regex = re.compile(bandext)
    fileid = re.split(regex, filename)[0]
    return fileid


def check_filename(args):
    fileid = parse_filename(args[1])
    if fileid != parse_filename(args[2]) or fileid != parse_filename(args[3]):
        print 'ERROR: Files not from same product'
        exit(2)
    return fileid


def main():
    if len(sys.argv) != 4:
        print 'ERROR: Incorrect arguments'
        usage()
    else:
        outfile = check_filename(sys.argv)+'.jpg'
        path = os.getcwd()

        RGB = []
        RGB.append(get_data(os.path.join(path,str(sys.argv[1]))))
        RGB.append(get_data(os.path.join(path,str(sys.argv[2]))))
        RGB.append(get_data(os.path.join(path,str(sys.argv[3]))))

        whitepix = 0.15

        rgbArray = np.zeros((len(RGB[0]), len(RGB[0][1]), 3), 'uint8')
        rgbArray[..., 0] = np.maximum(0, np.minimum(RGB[2] * 255 / (whitepix), 255))  # red
        rgbArray[..., 1] = np.maximum(0, np.minimum(RGB[1] * 255 / (whitepix), 255))  # green
        rgbArray[..., 2] = np.maximum(0, np.minimum(RGB[0] * 255 / (whitepix), 255))  # blue

        img = Image.fromarray(rgbArray)
        img.save(outfile, quality=92)


if __name__ == "__main__":
    main()
