#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
from __future__ import print_function

import sys
import argparse
import logging
import os.path
from osgeo import gdal
from osgeo import osr
import numpy as np
import fiona


gdal.UseExceptions()

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

formatter = logging.Formatter('%(levelname)s:%(asctime)s - %(message)s', '%d-%m-%Y %H:%M:%S')

fh = logging.FileHandler('intercal.log', 'w')
fh.setFormatter(formatter)
logger.addHandler(fh)

ch = logging.StreamHandler(sys.stdout)
ch.setLevel(logging.INFO)
ch.setFormatter(formatter)
logger.addHandler(ch)


def create_parser():
    parser = argparse.ArgumentParser(prog='intercal',
                                     description='Relative radiometric normalization routine',
                                     epilog='''Performs a histogram match with the reference image, by default.
                                     with type=pif will use selected areas to fit''')
    parser.add_argument('-r', '--refdir',
                        default=os.path.join(os.path.dirname(__file__), 'reference/'),
                        help="Path to ref_file images, defaults to ./reference",
                        metavar='DIR', type=lambda x: is_valid_dir(parser, x))
    parser.add_argument('-i', '--indir',
                        default=os.path.join(os.path.dirname(__file__), 'input/'),
                        help="Path to images to be calibrated, defaults to ./input",
                        metavar='DIR', type=lambda x: is_valid_dir(parser, x))
    parser.add_argument('-o', '--outdir',
                        default=os.path.join(os.path.dirname(__file__), 'output/'),
                        help="Path to output calibrated images, defaults to ./output",
                        metavar='DIR', type=lambda x: is_valid_dir(parser, x))
    parser.add_argument('-p', '--pif',
                        default=os.path.join(os.path.dirname(__file__), 'PIF/PIF.shp'),
                        help="Path to shapefile defining PIFs, defaults to ./PIF/PIF.shp",
                        metavar='FILE', type=lambda x: is_valid_file(parser, x))
    parser.add_argument('--debug', default=False, action='store_true',
                        help="Provide debugging information")
    parser.add_argument('-t', '--type', default='match', type=str,
                        choices=['match', 'pif'],
                        help="Type of calibration - '[match]','pif'")
    return parser


def is_valid_dir(parser, path):
    # helper function to check if directory exists
    if not os.path.isdir(path):
        parser.error('The directory {} does not exist!'.format(path))
    else:
        return path


def is_valid_file(parser, filename):
    # helper function to check if file exists
    if not os.path.isfile(filename):
        parser.error('The file {} does not exist!'.format(filename))
    else:
        return filename


def parse_filename(filename):
    # parses filename returns band if *_Bn.tif
    band = filename.split('_')[-1].split('.')[0]
    return band


def open_geotiff(path):
    # returns file descriptor of open GeoTiff, or skips invalid files
    try:
        ds = gdal.Open(path, gdal.GA_ReadOnly)
    except RuntimeError:
        logging.debug("I/O exception! Skipping %s", path)
        return None
    if ds.GetDriver().LongName != 'GeoTIFF':
        logging.debug("Incompatible format! Skipping %s", path)
        return None
    if ds.RasterCount != 1:
        logging.debug("Multiple bands! Skipping %s", path)
    return ds


def write_geotiff(im, meta, outfile):
    # write out geotiff using metadata provided as argument
    outfile = outfile
    logger.info("Writing output")
    driver = gdal.GetDriverByName('GTiff')
    if driver is None:
        raise ValueError("Cannot find GeoTiff driver")
    ds = driver.Create(
        outfile,
        meta['bbox']['cols'],
        meta['bbox']['rows'],
        1,
        gdal.GDT_Float32,
        options=[
            'TILED=YES',
            'BIGTIFF=IF_SAFER',
            'BLOCKXSIZE=' + str(meta['blocksize'][0]),
            'BLOCKYSIZE=' + str(meta['blocksize'][1]),
            'COMPRESS=DEFLATE'
        ]
    )
    band = ds.GetRasterBand(1)
    band.WriteArray(im)
    band.SetNoDataValue(meta['ndv'])
    band.ComputeStatistics(1)
    ds.SetGeoTransform(meta['geotransform'])
    ds.SetProjection(meta['projection'])
    band.FlushCache()
    band = None
    ds = None
    return band, ds


def get_bbox_coords(ds):
    # returns bounding box coordinates and resolution
    ulx, xres, xrot, uly, yrot, yres = ds.GetGeoTransform()
    cols = ds.RasterXSize
    rows = ds.RasterYSize
    lrx = ulx + cols * xres + rows * xrot
    lry = uly + rows * yres + cols * yrot
    tlg = brg = None  # corners in lat/lon
    if osr.SpatialReference(ds.GetProjection()).IsProjected:
        src = osr.SpatialReference()
        src.ImportFromWkt(ds.GetProjection())
        tgt = osr.SpatialReference()
        tgt.ImportFromEPSG(4326)
        tlg = coord_transform(src, tgt, ulx, uly)[0:2]
        brg = coord_transform(src, tgt, lrx, lry)[0:2]
    loc = {
        'ulx':  ulx,  'uly':  uly,
        'lrx':  lrx,  'lry':  lry,
        'xres': xres, 'yres': yres,
        'xrot': xrot, 'yrot': yrot,
        'cols': cols, 'rows': rows,
        'tlg':  tlg,  'brg':  brg
    }
    return loc


def get_statistics(ds):
    # return geotiff statistics as dictionary
    band = ds.GetRasterBand(1)
    bmin, bmax, mean, std = band.ComputeStatistics(True)
    stats = {'min': bmin,
             'max': bmax,
             'mean': mean,
             'std': std}
    return stats


def get_metadata(ds):
    # create a dictionary from geotiff metadata values
    bbox = get_bbox_coords(ds)
    blocksize = ds.GetRasterBand(1).GetBlockSize()
    geotransform = ds.GetGeoTransform()
    projection = ds.GetProjection()
    ndv = ds.GetRasterBand(1).GetNoDataValue()
    meta = {'bbox': bbox,
            'blocksize': blocksize,
            'ndv': ndv,
            'geotransform': geotransform,
            'projection': projection
            }
    return meta


def match_histogram(ref_ds, ds):
    """
    Given an input image this routine adjusts the mean and standard deviation
    to those of a reference image. Handles masked arrays.
    """
    stats = get_statistics(ref_ds)
    mean_ref = stats['mean']
    std_ref = stats['std']
    stats = get_statistics(ds)
    mean_input = stats['mean']
    std_input = stats['std']
    scale = std_ref / std_input
    offset = mean_ref - mean_input * std_ref / std_input
    band = ds.GetRasterBand(1)
    ndv = band.GetNoDataValue()
    arr = band.ReadAsArray()
    arr = np.ma.masked_equal(arr, ndv, copy=False)
    arr = arr * scale + offset
    return arr


def do_polyfit(x, y, degree):
    # perform a linear regression fit of two arrays
    results = {}
    coeffs = np.polyfit(x, y, degree)
    results['polynomial'] = coeffs.tolist()
    p = np.poly1d(coeffs)
    yhat = p(x)
    ybar = np.sum(y)/len(y)
    ssreg = np.sum((yhat-ybar)**2)
    sstot = np.sum((y - ybar)**2)
    results['determination'] = ssreg / sstot
    return results


def coord_transform(src, tgt, x, y):
    # transform points from source to target projection
    transform = osr.CoordinateTransformation(src, tgt)
    return transform.TransformPoint(x, y)


def get_point_coords(ds, opts):
    # extract list of coordinates from shape file
    proj_ref = osr.SpatialReference()
    proj_ref.ImportFromWkt(ds.GetProjection())
    with fiona.open(opts.pif, 'r') as shp:
        epsg = int(shp.crs['init'].split(':')[1])
        proj_src = osr.SpatialReference()
        proj_src.ImportFromEPSG(epsg)
        points = []
        for feature in shp:
            lon = feature['geometry']['coordinates'][0]
            lat = feature['geometry']['coordinates'][1]
            coords = coord_transform(proj_src, proj_ref, lon, lat)
            points.append(coords)
    return points


def convert_proj2pixel(coords, meta):
    # calculate the pixel index corresponding to affine coordinates
    x = coords[0]
    y = coords[1]
    xorigin = meta['bbox']['ulx']
    yorigin = meta['bbox']['uly']
    pixelwidth = meta['bbox']['xres']
    pixelheight = meta['bbox']['yres']
    xoffset = int((x - xorigin) // pixelwidth)
    yoffset = int((y - yorigin) // pixelheight)
    return xoffset, yoffset


def get_point_data(ds, points):
    # given a set of affine points return image data at points
    meta = get_metadata(ds)
    band = ds.GetRasterBand(1)
    im = band.ReadAsArray()
    data = []
    for point in points:
        x, y = convert_proj2pixel(point, meta)
        data.append(np.mean(im[y-1:y+2, x-1:x+2]))
    return data


def apply_correction(ds, fit, outfile):
    # apply scaling and offset to raster and write out
    meta = get_metadata(ds)
    band = ds.GetRasterBand(1)
    im = band.ReadAsArray()
    im = (im - fit['polynomial'][1]) / fit['polynomial'][0]
    write_geotiff(im, meta, outfile)
    return


# ALGORITHM DRIVERS #

def do_histogram_match(opts):
    # driver for the histogram matching algorithm
    logger.info("Performing histogram matching")
    for reffile in os.listdir(opts.refdir):
        ref_ds = open_geotiff(os.path.join(opts.refdir, reffile))
        if ref_ds is not None:
            band = parse_filename(reffile)
            logger.info("Processing band %s", band)
            for infile in os.listdir(opts.indir):
                if parse_filename(infile) == band:
                    in_ds = open_geotiff(os.path.join(opts.indir, infile))
                    if in_ds is not None:
                        logger.info("Correcting %s", infile)
                        outfile = os.path.join(opts.outdir, 'C_'+infile)
                        im = match_histogram(ref_ds, in_ds)
                        meta = get_metadata(in_ds)
                        write_geotiff(im, meta, outfile)
    return


def do_interp_pifs(opts):
    # driver for the point PIF interpolation algorithm
    logger.info("Performing PIF interpolation")
    for reffile in os.listdir(opts.refdir):
        ref_ds = open_geotiff(os.path.join(opts.refdir, reffile))
        if ref_ds is not None:
            logger.info("Reference %s", reffile)
            band = parse_filename(reffile)
            points = get_point_coords(ref_ds, opts)
            data_ref = get_point_data(ref_ds, points)
            for infile in os.listdir(opts.indir):
                if parse_filename(infile) == band:
                    in_ds = open_geotiff(os.path.join(opts.indir, infile))
                    if in_ds is not None:
                        logger.info("Correcting %s", infile)
                        points = get_point_coords(in_ds, opts)
                        data_input = get_point_data(in_ds, points)
                        fit = do_polyfit(data_ref, data_input, 1)
                        outfile = os.path.join(opts.outdir, 'C_'+infile)
                        apply_correction(in_ds, fit, outfile)
    return


def main():
    parser = create_parser()
    opts = parser.parse_args()
    if opts.debug:
        logger.setLevel(logging.DEBUG)
    if opts.type == 'match':
        do_histogram_match(opts)
    elif opts.type == 'pif':
        do_interp_pifs(opts)
    else:
        # TODO: Add Histogram stretch
        # TODO: Add Polygon PIFs
        # TODO: Add Auto PIF selection
        parser.error('The type method is not implemented')
    return


if __name__ == "__main__":
    print("Exiting")
    exit()
    main()

    for handler in logger.handlers:
        handler.close()
        logger.removeFilter(handler)
