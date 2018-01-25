# test_intercal.py

from intercal import *

import os
import pytest
import numpy as np
from osgeo import gdal
from osgeo import osr


@pytest.fixture()
def cli():
    parser = create_parser()
    return parser


# Test command line parameter parsing -----------------------------------------

def test_cli_default_refdir(cli):
    opts = cli.parse_args([])
    assert os.path.exists(opts.refdir)


def test_cli_default_indir(cli):
    opts = cli.parse_args([])
    assert os.path.exists(opts.indir)


def test_cli_default_outdir(cli):
    opts = cli.parse_args([])
    assert os.path.exists(opts.outdir)


def test_cli_default_shapefile(cli):
    opts = cli.parse_args([])
    assert os.path.isfile(opts.pif)


def test_cli_default_debug(cli):
    opts = cli.parse_args([])
    assert opts.debug is False


def test_cli_default_type(cli):
    opts = cli.parse_args([])
    assert opts.type == 'match'


def test_cli_refdir_notexist(cli):
    with pytest.raises(SystemExit):
        cli.parse_args(['-r', 'notexist'])


def test_cli_indir_notexist(cli):
    with pytest.raises(SystemExit):
        cli.parse_args(['-i', 'notexist'])


def test_cli_outdir_notexist(cli):
    with pytest.raises(SystemExit):
        cli.parse_args(['-o', 'notexist'])


def test_cli_shapefile_notexist(cli):
    with pytest.raises(SystemExit):
        cli.parse_args(['-p', 'notexist'])


def test_cli_type_notexist(cli):
    with pytest.raises(SystemExit):
        cli.parse_args(['-t', 'notexist'])


def test_cli_debug_is_set(cli):
    opts = cli.parse_args(['--debug'])
    assert opts.debug is True

# Test access to GeoTIFF files ------------------------------------------------


def test_geotiff_not_correct_format():
    testfile = os.path.join(os.path.dirname(__file__), 'files/small_tif_notags.png')
    assert open_geotiff(testfile) is None


def test_geotiff_multiband_skip():
    testfile = os.path.join(os.path.dirname(__file__), 'files/small_tif_2band.png')
    assert open_geotiff(testfile) is None


def test_geotiff_is_valid_format():
    testfile = os.path.join(os.path.dirname(__file__), 'files/small_tif_tags.tif')
    assert open_geotiff(testfile) is not None


def test_parse_filename():
    assert parse_filename('C_32TLP_20160823_B1.TIF') == 'B1'
    assert parse_filename('T51JXM_20161229T014652_B02.jp2') == 'B02'


def test_bbox_coordinates_latlon():
    testfile = os.path.join(os.path.dirname(__file__), 'files/small_tif_tags.tif')
    coords = get_bbox_coords(open_geotiff(testfile))
    assert pytest.approx(coords['ulx']) == 150.9100000
    assert pytest.approx(coords['uly']) == -34.1700000
    assert pytest.approx(coords['lrx']) == 150.9491667
    assert pytest.approx(coords['lry']) == -34.2300000


def test_bbox_coordinates_utm():
    testfile = os.path.join(os.path.dirname(__file__), 'files/small_tif_proj.tif')
    coords = get_bbox_coords(open_geotiff(testfile))
    assert pytest.approx(coords['ulx']) == -(117+38.0/60.0+30.24/3600.0)
    assert pytest.approx(coords['uly']) == 33+56.0/60.0+37.8/3600.0
    assert pytest.approx(coords['lrx']) == -(117+18.0/60.0+31.15/3600.0)
    assert pytest.approx(coords['lry']) == 33+39.0/60.0+54.26/3600.0


def test_geotiff_write():
    array = np.array(((0.1, 0.2, 0.3, 0.4),
                      (0.2, 0.3, 0.4, 0.5),
                      (0.3, 0.4, 0.5, 0.6),
                      (0.4, 0.5, 0.6, 0.7),
                      (0.5, 0.6, 0.7, 0.8)))

    lat = np.array(((10.0, 10.0, 10.0, 10.0),
                    (9.5, 9.5, 9.5, 9.5),
                    (9.0, 9.0, 9.0, 9.0),
                    (8.5, 8.5, 8.5, 8.5),
                    (8.0, 8.0, 8.0, 8.0)))

    lon = np.array(((20.0, 20.5, 21.0, 21.5),
                    (20.0, 20.5, 21.0, 21.5),
                    (20.0, 20.5, 21.0, 21.5),
                    (20.0, 20.5, 21.0, 21.5),
                    (20.0, 20.5, 21.0, 21.5)))

    xmin, ymin, xmax, ymax = [lon.min(), lat.min(), lon.max(), lat.max()]
    nrows, ncols = np.shape(array)
    xres = (xmax - xmin) / float(ncols)
    yres = (ymax - ymin) / float(nrows)
    geotransform = (xmin, xres, 0, ymax, 0, -yres)

    output_raster = gdal.GetDriverByName('GTiff').Create('/tmp/test.tif', ncols, nrows, 1,
                                                         gdal.GDT_Float32)
    output_raster.SetGeoTransform(geotransform)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    output_raster.SetProjection(srs.ExportToWkt())
    output_raster.GetRasterBand(1).WriteArray(array)
    output_raster = None

    ds = open_geotiff('/tmp/test.tif')
    stats = get_statistics(ds)
    assert pytest.approx(stats['min']) == 0.1
    assert pytest.approx(stats['max']) == 0.8
    assert pytest.approx(stats['mean']) == 0.45
    assert pytest.approx(stats['std']) == 0.18027756


# Test functions --------------------------------------------------------------


def test_polyfit():
    centigrade = range(0, 100, 2)
    fahrenheit = (9/5.0 * np.array(centigrade)) + 32
    result = do_polyfit(centigrade, fahrenheit, 1)
    assert pytest.approx(result['polynomial'][0]) == 1.8
    assert pytest.approx(result['polynomial'][1]) == 32
    assert pytest.approx(result['determination']) == 1


def test_match_histogram():
    reference = os.path.join(os.path.dirname(__file__), 'files/reference_B1.tif')
    target = os.path.join(os.path.dirname(__file__), 'files/target_B1.tif')
    ref_ds = open_geotiff(reference)
    tgt_ds = open_geotiff(target)
    array1 = ref_ds.GetRasterBand(1).ReadAsArray()
    array2 = tgt_ds.GetRasterBand(1).ReadAsArray()
    correct = match_histogram(ref_ds, tgt_ds)
    assert np.mean(array2) != np.mean(array1)
    assert np.std(array1) != np.std(array2)
    assert np.mean(correct) == np.mean(array1)
    assert pytest.approx(np.std(correct)) == np.std(array1)


def test_nodata_reference():
    reference = os.path.join(os.path.dirname(__file__), 'files/nodataval_B1.tif')
    target = os.path.join(os.path.dirname(__file__), 'files/target_B1.tif')
    ref_ds = open_geotiff(reference)
    tgt_ds = open_geotiff(target)
    array1 = ref_ds.GetRasterBand(1).ReadAsArray()
    array1 = np.ma.masked_equal(array1, -99)
    array2 = tgt_ds.GetRasterBand(1).ReadAsArray()
    array2 = np.ma.masked_equal(array2, -99)
    correct = match_histogram(ref_ds, tgt_ds)
    assert np.mean(array1) != np.mean(array2)
    assert np.std(array1) != np.std(array2)
    assert pytest.approx(np.mean(correct)) == np.mean(array1)
    assert pytest.approx(np.std(correct)) == np.std(array1)


def test_nodata_target():
    reference = os.path.join(os.path.dirname(__file__), 'files/reference_B1.tif')
    target = os.path.join(os.path.dirname(__file__), 'files/nodataval_B1.tif')
    ref_ds = open_geotiff(reference)
    tgt_ds = open_geotiff(target)
    array1 = ref_ds.GetRasterBand(1).ReadAsArray()
    array1 = np.ma.masked_equal(array1, -99)
    array2 = tgt_ds.GetRasterBand(1).ReadAsArray()
    array2 = np.ma.masked_equal(array2, -99)
    correct = match_histogram(ref_ds, tgt_ds)
    assert np.mean(array1) != np.mean(array2)
    assert np.std(array1) != np.std(array2)
    assert pytest.approx(np.mean(correct)) == np.mean(array1)
    assert pytest.approx(np.std(correct)) == np.std(array1)
