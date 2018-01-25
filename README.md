# INTERCAL - Optical Sensor Inter-calibration Tool
This tool is developed in the scope of the CoReSyF project (co-resyf.eu),
which is a European Commission Horizon 2020 funded project to provide a
Coastal Waters Research Synergy Framework by creating an online platform
to support research applications that use Earth Observation (EO) data for
coastal water modelling.

Image inter-calibration is required for images due to differences in 
viewing and illumination characteristics, and due to differences in the
spectral response function of different sensors. Despite sensors being
calibrated pre-launch, when in use they degrade over time i.e. drift.
Recalibration against an on-board reference standard suffers similar
problems i.e. the reference unit degrades over time. Consequently, many
on-the-fly calibration methods rely on analysis of over PICS, targets
of known and invariant radiance i.e. sensor verification is required
for absolute sensor calibration. This is a challenging task.

This tool does not attempt absolute calibration but attempts to rescale
sensor reflectances (top-of-atmosphere) against a selected reference
image.

## Use
Use `python intercal.py -h` to view the command line options.

## Notes
The current version requires Sentinel-2 images as single band GeoTIFF
files that are of the same area of interest (AOI) and co-registered.

## Collaborators
Romain Serra (ACRI)
