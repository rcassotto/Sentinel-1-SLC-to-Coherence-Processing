# Sentinel-1-SLC-to-Coherence-Processing

### Python script / wrapper to process Sentinel-1 Coherence images from SLC zipped files
This script is written as a python wrapper around ESA's (European Space Agency) open-sourced Sentinel Application Platform (SNAP). It relies heavily on  the graph process tool (gpt) to pre-process Sentinel-1 Single Look Complex (SLC) files into coherence images. There are also external calls to gdal and a shell script to unzip embedded map overlay kml files that are used to check for sufficient overlap between coherence image pairs.

It is recommended to execute from within a dedicated python environment, such as conda. 

### _Python dependencies_
- geopandas
- shapely
- pandas
- fiona
<br/>


### _Non-Python dependencies_
- ESA's SNAP toolkit (https://step.esa.int/main/download/snap-download/)
- gdal
<br/>

### _Initial Setup Requirements_
The 
