# 讀取檔案用
import os
import glob
from osgeo import gdal
from pathlib import Path

# 建立資料夾是否存在的判斷函數
def dir_maker(dir_path):
    dir_path_list = Path(dir_path).parts
    if Path(dir_path).is_dir() == False:
        for i in range(len(dir_path_list)):
            now_dir = './' + '/'.join([dir_path_list[j] for j in range(i + 1)])
            if Path(now_dir).is_dir() == False:
                os.mkdir(now_dir)
                if i == len(dir_path_list) - 1:
                    print(f"新增此路徑: {now_dir}")
    else:
        pass

# 期作別判斷
def Which_Period(date):
    year = int(date[:4])
    month = int(date[4:6])
    if year == 2024:
        if month < 7:
            return '113-1'
        elif month > 6:
            return '113-2'
        
def merge_tif_files(indicator, input_folder, output_folder, date):
    # 獲取指定指標的所有 TIFF 文件
    tif_files = [os.path.join(input_folder, f) for f in os.listdir(input_folder) if f.startswith(indicator) and f.endswith('.tif')]
    
    # 獲取第一個 TIFF 文件的投影和地理變換信息
    src_ds = gdal.Open(tif_files[0])
    projection = src_ds.GetProjection()
    geotransform = src_ds.GetGeoTransform()
    src_ds = None
    
    # 使用 gdal.Warp 合併 TIFF 文件
    output_file = os.path.join(output_folder, f'{date}_{indicator}_merged.tif')
    gdal.Warp(output_file, tif_files, format='GTiff', srcSRS=projection, dstSRS=projection, xRes=geotransform[1], yRes=-geotransform[5])