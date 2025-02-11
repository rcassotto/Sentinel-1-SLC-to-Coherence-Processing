# Sentinel-1-SLC-to-Coherence-Processing

### Python script / wrapper to process Sentinel-1 Coherence images from SLC zipped files
This script is written as a python wrapper around ESA's (European Space Agency) open-sourced Sentinel Application Platform (SNAP). It relies heavily on SNAP's graph process tool (gpt) to pre-process Sentinel-1 Single Look Complex (SLC) files into coherence images. There are also external calls to gdal and a shell script to unzip embedded map overlay kml files that are used to check for sufficient overlap between coherence image pairs.

### _Initial Setup Requirements_
1) Install ESA's SNAP program, if not already installed (see below).
2) Install gdal, if not already installed (see below).
3) Configure a python environment and install necessary dependencies (see below).
4) Create a dedicated directory for precise and restituted orbits. 
5) Open a python editor, and perform a search and replace for hard-coded paths for the following variables: 
    - **aux_poe**: precise orbit files.
    - **res_orb**: restituted orbit files.
    - **gdal_ext**: full path for the location of gdal_transate.
    - **base_snap_exe**: full pathway for ESA SNAP's gpt binary.
    - **workflow_dir**: Workflow directory where placed the _extract_map_overlays.sh_ shell script. 
<br/>

### _Non-Python dependencies_
- ESA's SNAP toolkit (https://step.esa.int/main/download/snap-download/)
- gdal(https://gdal.org/en/stable/download.html)
<br/>

### _Python dependencies_
It is **_strongly_** recommended to create a dedicated python environment, install the necessary dependencies, and execute the scripts from within that environment.  

- geopandas
- shapely
- pandas
- fiona
<br/>
 
### _Executing Script_
To run the script, perform the following:    

1) Activate the dedicated python environment

2) Use your favorite text editor to modify the input file (e.g. _FIREDpy_process_coh_input_asc.txt_) for your desired inputs. 
    - **_slc_file_loc_**: location of SLC zip files.
    - **_roi_polygon_**: for future use. Leave as blank quotes. 
    - **_roi_path_**: for future use. Leave as blank quotes. 
    - **_output_resolution_m_**: desired output resolution in meters. **NOTE higher resolutions will result in longer processing times.**
    - **_sys_index_var_**: starting "Python" index for SLC to Coh processing.  Set to 0 to start from the beginning. Set to a higher integer to continue processing a stack of images, if or when the processing is interrupted. 

3) Execute the command: **_python3_** **_Process_SLC2Coh_wSNAP_v1.0.py_** **_FIREDpy_process_coh_input_asc.txt_**

4) The script will take considerable time to process (hours-days) depending on the output resolution and number of files. It will process several folders and files. The final coherence images will appear in subdirectories witin the "Processed_Data" folder.  The subdirectories for each image will be named with the reference and secondary images as yyyymmddTHHMMSS_yyyymmddTHHMMSS (e.g. 20200806T141458_20200818T141458). 
       
