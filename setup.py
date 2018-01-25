from distutils.core import setup

setup(
    name='InterCal',
    version='0.1.0',
    author='Dr S. M. Emsley (ARGANS)',
    author_email='semsley@argans.co.uk',
    license='LICENSE',
    description='Optical sensor inter-calibration tool',
    long_description=open('README.md').read(),
    requires=['numpy', 'fiona', 'gdal']
)
