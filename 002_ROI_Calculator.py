'''
這個程式用來切出ROI的範圍，並且計算ROI區域範圍內的植生指標
'''

from osgeo import gdal, osr, ogr
import shapefile

from skimage import io, filters
from skimage.color import rgb2hsv, rgb2lab
import skimage
import numpy as np 

import glob
from pathlib import Path 

import csv 
import os
import threading


from tqdm import tqdm, trange
import time


FieldName = "Field"
fieldnames = ['site', 'date', 'Period', 'Field', 'Id', 'GC', 'NDVI', 'NDVI_sum', 'NDVI_max', 'NDVI_min', 'NDVI_range', 'NDVI_std', 'NDRE', 'NDRE_sum', 'NDRE_max', 'NDRE_min', 'NDRE_range', 'NDRE_std']

st_time = time.time()

#期作別判斷
def Which_Period(date):
    year = int(date[:4])
    month = int(date[4:6])
    if year == 2022:
        if month < 7:
            return '111-1'
        elif month > 6:
            return '111-2'
    elif year == 2023:
        if month < 7:
            return '112-1'
        elif month > 6:
            return '112-2'

#建立資料夾是否存在的判斷函數
def dir_maker(dir_path):
    dir_path_list = Path(dir_path).parts
    if Path(dir_path).is_dir()==False:
        for i in range(len(dir_path_list)):
            now_dir = "./" + "/".join([dir_path_list[j] for j in range(i + 1)])
            print(i, now_dir)
            if Path(now_dir).is_dir()==False:
                os.mkdir(now_dir)
    else:
        pass

csv_file_path = f'./04_Result/Data/ROI/'
csv_file_name = 'ROI_Result.csv'
dir_maker(csv_file_path)

with open(f'{csv_file_path}{csv_file_name}','w',encoding='utf-8',newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()

    #讀取6band圖片檔案
    for image_path in glob.glob('./01_Tif/6band/BG/*.tif'):
        #print出圖片路徑
        print(image_path)

        #建立資料容器
        data = {}

        #讀取場域跟日期
        site = Path(image_path).parts[-2]
        date = Path(image_path).parts[-1][:10].replace('-','')
        Period = Which_Period(date)
        print(site, date, Period)

        if Period != '112-1':
            continue

        #初始資料寫入
        data['site'] = site
        data['date'] = date
        data['Period'] = Period

        #讀取6band的原圖資訊
        dataset = gdal.Open(image_path,gdal.GA_ReadOnly)
        img  = dataset.ReadAsArray()/ 65535
        img[img==1.0]=0

        ds = dataset.GetGeoTransform()
        
        B = img[0,:,:]
        G = img[1,:,:]
        R = img[2,:,:]
        RE = img[3,:,:]
        NIR = img[4,:,:]

        RGB_img = np.dstack((R,G,B))

        HSV = rgb2hsv(RGB_img)
        Lab = rgb2lab(RGB_img)

        H = HSV[:,:,0]
        a = Lab[:,:,1]

        #計算NDVI跟NDRE
        NDVI = (NIR - R) / (NIR + R + (10**(-10)))
        NDRE = (NIR - RE) / (NIR + RE + (10**(-10)))
        GBGAP = G - B
        GRGAP = G - R
        B_MEAN = np.nanmean(B)

        #分別使用大津演算法二值化
        thresh1 = filters.threshold_otsu(NIR)
        thresh2 = filters.threshold_otsu(NDVI)
        thresh3 = filters.threshold_otsu(NDRE)
        thresh6 = filters.threshold_otsu(a)

        #建立mask
        mask1 = (NIR >= thresh1) * 1.0
        mask2 = (NDVI >= thresh2) * 1.0
        mask3 = (NDRE >= thresh3) * 1.0
        mask4 = (GBGAP >= 0) * 1.0
        mask5 = (GRGAP >= 0) * 1.0
        mask6 = (a < thresh6) * 1.0
        mask7 = (B <= B_MEAN) * 1.0
        
        if date == '20220315':
            mask = mask2 * mask6 * mask7

        else:
            mask = mask1 * mask2 * mask3 * mask4 * mask5

        #建立放大圖mask的資料夾
        mask_dir = f'./04_Result/Image/{site}/ROI/{date}/mask'
        dir_maker(mask_dir)

        #將大圖mask的tif寫入
        tif_writer = gdal.Translate(f'{mask_dir}/{date}.tif', dataset,bandList=[1],outputType=gdal.GDT_Float32) 
        tif_writer.GetRasterBand(1).WriteArray(mask)
        del tif_writer

        print('大圖切割完畢', '花費時間',round(time.time() - st_time,2),"秒")        
        #碼表重新計算時間
        st_time = time.time()

        #把ROI的shapefile讀出來
        shapefile_path = glob.glob(f'./02_Shapefile/{site}/ROI/{Period}_ROI.shp')[0]

        #讀取整個shapefile
        input_shapefile = ogr.Open(shapefile_path)
        input_layer = input_shapefile.GetLayer(0)

        #針對每個日期每個區域，建立用shapefile切下來的6band檔案
        dir_maker(f'./04_Result/Shapefile/{site}/ROI/{Period}')
        dir_maker(f'./04_Result/Image/{site}/ROI/{date}/tif')

        #在這個shapefile裡面依據每個區塊儲存小塊shapefile
        for i in range(input_layer.GetFeatureCount()):
            try:
                feature = input_layer.GetFeature(i)
                Id = feature.GetField('Id')
                Field = feature.GetField('Field')

                shp_output_path = f'./04_Result/Shapefile/{site}/ROI/{Period}'
                #建立儲存路徑
                save_shp_path = os.path.join(shp_output_path, (f"{Field}_{Id}.shp"))
                if os.path.isfile(save_shp_path) == True or Field == 'temp':
                    continue

                driver = ogr.GetDriverByName('ESRI Shapefile') ## a driver
                datasource = driver.CreateDataSource(save_shp_path)
                layer = datasource.CreateLayer('layerName',geom_type=ogr.wkbPolygon)
                layer.CreateFeature(feature)
                geometry = feature.GetGeometryRef()
            except:
                pass


        #開啟要切成小塊tif的圖檔
        input_raster = gdal.Open(image_path)
        
        #刷新碼錶
        print('小塊shapefile切割完畢', '花費時間',round(time.time() - st_time,2),"秒")
        st_time = time.time()

        #用每個小shapefile切割tif
        for shp_path in glob.glob(f'./04_Result/Shapefile/{site}/ROI/{Period}/*.shp'):
            try:
                Field = Path(shp_path).parts[-3]

                tif_file_path = f'./04_Result/Image/{site}/ROI/{date}/tif/' + Path(shp_path).stem + '.tif'
                if os.path.isfile(tif_file_path) == True:
                    continue
                #os.system('%s 02_arcpy_cuter.py   %s  %s  %s'%(arcgis_python_path, image_path, shp_path, Plant_tif_file_path))

                st_time = time.time()

                #獲得輸入shapefile的左上右下座標資料
                r = shapefile.Reader(shp_path)
                minX, minY, maxX, maxY = r.bbox

                #裁切tif檔案：6band
                ds = gdal.Warp(
                    tif_file_path,#輸出raster路徑
                    input_raster,#輸入raster檔案，用gdal.Open()開啟
                    format = 'GTiff',#格式
                    cutlineDSName = shp_path,#shapefile路徑
                    outputBounds = [minX, minY, maxX, maxY],
                    dstNodata = 0,
                )

                del ds

            except Exception as e:
                print(shp_path)
                quit()
                pass

        print('shapefile切割結束', '花費時間',round(time.time() - st_time,2),"秒") 
        st_time = time.time()

        #讀取剛切完的小塊tif列表
        tif_file_path_list = glob.glob(f'./04_Result/Image/{site}/ROI/{date}/tif/*.tif')

        for i in trange(len(tif_file_path_list)):
            tif_file_path = tif_file_path_list[i]
            tif_name = Path(tif_file_path).stem

            Field = tif_name.split("_")[0]
            Id = tif_name.split("_")[1]

            data['Field'] = Field
            data['Id'] = Id

            #讀取原圖tif檔的資訊
            dataset = gdal.Open(tif_file_path,gdal.GA_ReadOnly)
            #讀取tif檔案原圖圖片資訊
            image  = dataset.ReadAsArray()/ 65535
            image[image==1.0]=0
            ds = dataset.GetGeoTransform()

            #讀取所有band
            B = image[0,:,:]
            G = image[1,:,:]
            R = image[2,:,:]
            RE = image[3,:,:]
            NIR = image[4,:,:]


            RGB_img = np.dstack((R,G,B))

            HSV = rgb2hsv(RGB_img)
            Lab = rgb2lab(RGB_img)

            H = HSV[:,:,0]
            a = Lab[:,:,1]

            #計算NDVI跟NDRE
            NDVI = (NIR - R) / (NIR + R + (10**(-10)))
            NDRE = (NIR - RE) / (NIR + RE + (10**(-10)))
            GBGAP = G - B
            GRGAP = G - R

            #建立mask
            mask1 = (NIR >= thresh1) * 1.0
            mask2 = (NDVI >= thresh2) * 1.0
            mask3 = (NDRE >= thresh3) * 1.0
            mask4 = (GBGAP >= 0) * 1.0
            mask5 = (GRGAP >= 0) * 1.0
            mask6 = (a < thresh6) * 1.0
            mask7 = (B <= B_MEAN) * 1.0

            if date == '20220315':
                mask = mask2 * mask6 * mask7

            else:
                mask = mask1 * mask2 * mask3 * mask4 * mask5

            mask_path = f'./04_Result/Image/{site}/ROI/{date}/mask'

            #將小圖mask的tif寫入
            tif_writer = gdal.Translate(f'{mask_path}/{tif_name}.tif', dataset,bandList=[1],outputType=gdal.GDT_Float32) 
            tif_writer.GetRasterBand(1).WriteArray(mask)
            del tif_writer

            GC_area = np.count_nonzero(image[0,:,:] * mask)
            total_area = np.count_nonzero(image[0,:,:])
            GC_rate =(GC_area / total_area)

            VIs = ['NDVI', 'NDRE']

            for VI in VIs:
                img = eval(VI) * mask
                VI_path = f'./04_Result/Image/{site}/ROI/{date}/tif/{VI}'
                dir_maker(VI_path)

                tif_writer = gdal.Translate(f'{VI_path}/{tif_name}.tif', dataset,bandList=[1],outputType=gdal.GDT_Float32) 
                tif_writer.GetRasterBand(1).WriteArray(img)
                del tif_writer

            #計算小塊的NDVI與其他植生指標數值
            NDVI = NDVI * mask
            NDRE = NDRE * mask
            NDVI[np.isnan(NDVI)] = 0
            NDRE[np.isnan(NDRE)] = 0
            NDVI = NDVI[NDVI!=0]
            NDRE = NDRE[NDRE!=0]

            data['GC'] = GC_rate
            data['NDVI'] = np.nanmean(NDVI)

            if str(data['NDVI']) == 'nan':
                data['NDVI'] = 0
                data['NDVI_sum'] = 0
                data['NDVI_max'] = 0
                data['NDVI_min'] = 0
                data['NDVI_range'] = 0
                data['NDVI_std'] = 0

            else:
                data['NDVI_sum'] = np.nansum(NDVI)
                data['NDVI_max'] = np.nanmax(NDVI)
                data['NDVI_min'] = np.nanmin(NDVI)
                data['NDVI_range'] = np.nanmax(NDVI) - np.nanmin(NDVI)
                data['NDVI_std'] = np.nanstd(NDVI)


            data['NDRE'] = np.nanmean(NDRE)
            if str(data['NDRE']) == 'nan':
                data['NDRE'] = 0
                data['NDRE_sum'] = 0
                data['NDRE_max'] = 0
                data['NDRE_min'] = 0
                data['NDRE_range'] = 0
                data['NDRE_std'] = 0

            else:
                data['NDRE_sum'] = np.nansum(NDRE)
                data['NDRE_max'] = np.nanmax(NDRE)
                data['NDRE_min'] = np.nanmin(NDRE)
                data['NDRE_range'] = np.nanmax(NDRE) - np.nanmin(NDRE)
                data['NDRE_std'] = np.nanstd(NDRE)

            writer.writerow(data)

print('分析結束', '花費時間',round(time.time() - st_time,2),"秒")
