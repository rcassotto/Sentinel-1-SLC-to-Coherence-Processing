import sys
import os

sys.path.append(".")

#from Main.helper_functions import general_functions
from SNAP_preProcessing_functions.helper_functions import general_functions
from pathlib import Path
import pathlib
import zipfile
import shutil
import datetime
from datetime import timedelta
from datetime import datetime
import glob

import shapely
import geopandas as gpd 
import numpy as np
import pandas as pd 

## Add driver support
from fiona.drvsupport import supported_drivers
supported_drivers['LIBKML'] = 'rw'

overlapping_threshold_prct=40

 ## calculate overlap between a list containing image pairs - RKC
def calculate_S1_scene_overlap(infilenames_and_paths):
     kml_polygons = read_kml_coords(infilenames_and_paths) 
     print(kml_polygons)
     poly_coords_0 = kml_polygons[0].get_coordinates(ignore_index=True)
     poly_coords_1 = kml_polygons[1].get_coordinates(ignore_index=True)
    
     ## Calcualte S1 scene areas
     deltaX_0_deg = np.max(np.diff(poly_coords_0.x))
     deltaY_0_deg = np.max(np.diff(poly_coords_0.y))
     deltaX_1_deg = np.max(np.diff(poly_coords_1.x))
     deltaY_1_deg = np.max(np.diff(poly_coords_1.y))
     poly_0_meanLat = np.mean(poly_coords_0.y)
     poly_1_meanLat = np.mean(poly_coords_1.y)
     poly_0_area = (deltaX_0_deg * np.cos(np.deg2rad(poly_0_meanLat))) * (deltaY_0_deg) *np.square(111)
     poly_1_area = (deltaX_1_deg * np.cos(np.deg2rad(poly_1_meanLat))) * (deltaY_1_deg) *np.square(111)
     print("Scene 0 area (km^2): ", poly_0_area)   
     print("Scene 1 area (km^2): ", poly_1_area)
     
     ## Calculate the shared geometry between polygons via shapely's intersection
     overlapping_region = shapely.intersection(kml_polygons[0], kml_polygons[1])
     for nn in range(0, len(overlapping_region)):
         print("Overlapping region: ", overlapping_region[nn]); print("")

     ## Use geopandas to convert overlap to a geopandas polygon, and calculate overlapping areas
     overlapping_polygon = gpd.GeoDataFrame(index=[0], crs='epsg:4326', geometry=overlapping_region)       
     test4NoOverlap=(overlapping_region.is_empty).tolist()
     
     ## check for empty overlapping region
     if test4NoOverlap[0] == True: 
         area_check_prct = np.zeros([1,])
         print("The scenes do not overlap. Moving on........")
     else:
         ## Calculate area of overlapping region
         poly_coords = overlapping_polygon.get_coordinates(ignore_index=True)
         mean_latitude = np.mean(poly_coords.y)
         overlap_area = list(overlapping_region.area) 
         overlap_area_km2 = np.float64(overlap_area)*np.square(111)*np.cos(np.deg2rad(mean_latitude)) 
         print("Mean Latitude:", mean_latitude)    
         print("Overlapping area: (deg)", overlap_area[0])
         print("Overlapping area: (km^2)", overlap_area_km2[0])
    
         ## Determine if overlapping region is more than n-% of original scenes
         area_check_prct = overlap_area_km2 / np.min([poly_0_area, poly_1_area])*100
     return area_check_prct


 ## function to read Sentinel-1 scene end points from *_map_overlay.kml files - RKC
def read_kml_coords(infilenames_and_paths):
     kml_polygons = []  # initialize a list of kml_polygons
     for n in range(0,len(infilenames_and_paths)):    
         print("Reading in: ", infilenames_and_paths[n])
         geo_df = gpd.read_file(infilenames_and_paths[n],driver='KML')  # READ KML file to a geopandas dataframe       
         df= pd.DataFrame(geo_df)         # Create Pandas Dataframe from GeoPandas 
         geodf = gpd.GeoDataFrame(df, geometry='geometry')
         poly1 = geodf.geometry; 
         kml_polygons.append(poly1)
     return kml_polygons
         




class set_up_for_insar_run():
#    def __init__(self, json_file_path, base_path, max_delta_t, sys_arg): # Original 
    def __init__(self, json_file_path, base_path, max_delta_t, sys_arg, data_path): # RKC mod to pass data_path
        self.sys_arg = sys_arg
        self.json_file_path = json_file_path
        self.base_path = base_path
        self.data_path = data_path
        self.max_delta_t = int(max_delta_t)
        self.slc_pairs_dict = None
        self.kml_frame_coords_kml_path = None
        

    def make_pair_list(self):
        print('Data path in make_pair_list: ', self.data_path)
        slc_list = glob.glob(self.data_path + "/*.zip")  # this is a list
        slc_list__ = list(glob.glob(self.data_path + "/*.zip")) # this is a list. Why redundant?
        slc_date_list = []
        max_delta_t_ob = timedelta(self.max_delta_t)
        for SLC in slc_list:
            slc_date = general_functions.slc_path_2_date(SLC)
            slc_date_list.append(slc_date)
        slc_list_ = [x for y, x in sorted(zip(slc_date_list, slc_list__))]
        slc_date_list_sorted = sorted(slc_date_list)
        slc_len = len(slc_date_list_sorted)
        j, k = 0, 0
        slc_pair_dict = {}
        for slc_date in slc_date_list_sorted:
            for i in range(j, slc_len):
                slc_pair_iset_dict = {}
                if i + 2 <= slc_len:
                    #t1, t2 = (datetime.strptime(slc_date, '%Y%m%d'),
                     #         datetime.strptime(slc_date_list_sorted[i + 1], '%Y%m%d'))  # original
                    t1, t2 = (datetime.strptime(slc_date, '%Y%m%dT%H%M%S'),
                              datetime.strptime(slc_date_list_sorted[i + 1], '%Y%m%dT%H%M%S'))  # RKC mod for time
                    t_delta = t2 - t1

#### code to check over overlap needs to occur here or within the next 'if' statement, before the generaton of Pairs_##.json                  
                    
                   #if t_delta <= max_delta_t_ob:  # original
                    if t_delta <= max_delta_t_ob and t_delta > timedelta(days=5):    # RKC edit to only include interferograms >11 day (to avoid adjacent scenes)
                        pair_name = 'Pair' + '_' + str(k)
                        slc_list_temp = [str(slc_list_[j]), str(slc_list_[i + 1])]
                        slc_pair_dict[pair_name] = slc_list_temp
                        slc_pair_iset_dict[pair_name] = slc_list_temp
                        
                        # print("HI THERE")
                        # print("length of slc_list_temp:", len(slc_list_temp))
                        print("SLC list temp:" , slc_list_temp)
                        
                        ### Use slc_list_temp, parse S1 filename, rename with *overlay.kml file, then run the overlap check. 
                        ### If overlap exceeds 40%, continue with code; if not, skip to next iteration of for loop (continue)
                        ### Create new list with overlay files. Pass this through the calculate overlap functions (minus geopackage writing steps)
                        
                        ## Create a list of *overlay.kml files
                        overlay_list_temp = []
                        for nOverlay in range(0, len(slc_list_temp)):
                            zipfile_path, zipfile_name = os.path.split(slc_list_temp[nOverlay])
                            overlay_fname = zipfile_name.split('_')[5] + "_map_overlay.kml"
                            overlay_fname_and_path = os.path.join(zipfile_path, overlay_fname)
                            overlay_list_temp.append(overlay_fname_and_path)
                            print("overlay base: ", overlay_fname_and_path)

                        
                        ## Check for sufficient overlap
                        area_check_prct = calculate_S1_scene_overlap(overlay_list_temp)
                        print("Area Overlap (%): ", area_check_prct)
                        
                        if (area_check_prct > overlapping_threshold_prct):
                            print("The two scenes overlap by at least ", overlapping_threshold_prct, "%. Continuing with interferogram pairing.")
                            
                            ## Create json files for each pair
                            pairs_dir = Path(self.data_path + '/Pairs')
                            if not pairs_dir.exists():
                                pairs_dir.mkdir()
                            slc_pair_iset_name = self.data_path + ("/Pairs/" + pair_name + ".json") # Creates Pair_##.json within /Pairs folder
                            general_functions.write_json_from_dict(slc_pair_iset_dict, slc_pair_iset_name)
                            k = k + 1
                            
                        else:
                            print("The scenes only overlap by", int(area_check_prct[0]),"%, which is less than the " , overlapping_threshold_prct,"%. Skipping this pairing.") 
                            continue


#                         ## Create json files for each pair
#                         pairs_dir = Path(self.data_path + '/Pairs')
# #                        print(type(pairs_dir))
# #                        print("Pairs dir: ", pairs_dir, type(pairs_dir))
#                         if not pairs_dir.exists():
#                             pairs_dir.mkdir()
#                         slc_pair_iset_name = self.data_path + ("/Pairs/" + pair_name + ".json") # Creates Pair_##.json within /Pairs folder
#                         general_functions.write_json_from_dict(slc_pair_iset_dict, slc_pair_iset_name)
#                         k = k + 1
            j = j + 1
            
        ## Create Pairs_list.json to iterate over
        print("Set_up data path: ", self.data_path)
        json_save_name = Path(self.data_path + "/Pairs_List.json")  # Ryan
        general_functions.write_json_from_dict(slc_pair_dict, json_save_name)
        self.slc_pairs_dict = general_functions.open_json_file(json_save_name)
        

    def create_pair_directories(self):
        jj = 0
        for slc_list in self.slc_pairs_dict.values():
            slc_1 = general_functions.slc_path_2_date(slc_list[0])
            slc_2 = general_functions.slc_path_2_date(slc_list[1])
            path_processed_data = Path(self.data_path + '/Processed_Data')
            path_cache_directory = Path(self.data_path + '/Cache')
            
            pair_directory_name = slc_1 + '_' + slc_2
            pair_cache_directory_name = pair_directory_name + '_cache'
            pair_cache_directory_path = Path(self.data_path + '/Cache/' + pair_cache_directory_name)
            pair_cache_directory_path_java = Path(self.data_path + '/Cache/' + pair_cache_directory_name + '/java_tmp')
            self.kml_frame_coords_kml_path = Path(self.data_path + '/kml_frame')  # RKC

            
            pair_directory_path = Path(str(path_processed_data) + '/' + pair_directory_name)# RKC
            if not path_processed_data.exists():
                path_processed_data.mkdir()
            if not path_cache_directory.exists():
                path_cache_directory.mkdir()
            if not pair_cache_directory_path.exists():
                pair_cache_directory_path.mkdir()
            if not pair_cache_directory_path_java.exists():
                pair_cache_directory_path_java.mkdir()
            if not self.kml_frame_coords_kml_path.exists():
                self.kml_frame_coords_kml_path.mkdir()
            if not pair_directory_path.exists():
                pair_directory_path.mkdir()
            if jj == 0 and self.sys_arg == 0:
                pair_cache_directory_path_aux = Path(self.data_path + '/Cache/' + pair_cache_directory_name + '/.esa_snap')
                
                if not pair_cache_directory_path_aux.exists():
                    pair_cache_directory_path_aux.mkdir()
                jj = jj + 1

    def extract_meta_data(self):
        meta_data_directory = Path(self.data_path + '/Meta_Data')
        if not meta_data_directory.exists():
            meta_data_directory.mkdir()
#        s1_data_dir_names = self.data_path.glob('S1*.zip')
        s1_data_dir_names = glob.glob(self.data_path + '/S1*.zip')
        ii = 0
        for zip_safe_dir in s1_data_dir_names:
            with zipfile.ZipFile(zip_safe_dir, "r") as zf:
                if ii == 0:
                    preview_path = [s for s in zf.namelist() if 'preview' and '.kml' in s][0]
                    parent_direct_kml = self.kml_frame_coords_kml_path / Path(preview_path).parent.parent
                    zf.extract(preview_path, path=self.kml_frame_coords_kml_path)
                    path_to_extracted = self.kml_frame_coords_kml_path / preview_path
                    shutil.copy2(path_to_extracted, self.kml_frame_coords_kml_path)
                    shutil.rmtree(parent_direct_kml)
                ii = ii + 1
                list_annotation_xml = ([s for s in zf.namelist() if "annotation" and "xml" in s
                                        and ("annotation" and "/s1a-" in s or
                                             "annotation" and "/s1b-" in s)])
                paths_list = []
                parent_direct_meta = meta_data_directory / Path(list_annotation_xml[0]).parent.parent
                for annotation_xml_path in list_annotation_xml:
                    zf.extract(annotation_xml_path, path=meta_data_directory)
                    path1 = meta_data_directory / annotation_xml_path
                    paths_list.append(path1)
                for item in paths_list:
                    item_path = Path(item)
                    shutil.copy2(item_path, meta_data_directory)
                shutil.rmtree(parent_direct_meta)
            zf.close()

    def run_functions(self):
        self.make_pair_list()
        self.create_pair_directories()
        self.extract_meta_data()
