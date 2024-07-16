"""
GDAL安裝的說明文章
https://cloud.tencent.com/developer/article/1621201

執行方法：
python 01_VIs_Calculator.py

1.建立要放NDVI、NDRE結果的csv，寫入欄位名稱
2.建立要寫入csv的資料容器，是由{key:value}所組成
3.以for迴圈形式將要處理的大圖都讀取出來(01_TIF)
4.找出大圖當中的Otsu閾值，找出植體的位置
5.把單株的shapefile炸裂，變成多個單元
6.依照炸裂後的shapefile單元，將大圖的tif檔案進行切割
7.將切割後的tif檔案，進行NDVI、NDRE與其統計值的運算

"""

#讀取檔案用
import os
import glob
from pathlib import Path

#分析影像用
from osgeo import gdal, ogr
from skimage import io, filters # 用於Otsu大津演算法
from skimage.color import rgb2hsv, rgb2lab
import shapefile
import numpy as np  # 矩陣運算

#圖片顯示用
import matplotlib.pyplot as plt

#輸出csv用
import csv

#顯示處理進度與計時用
from tqdm import tqdm, trange
import threading
import time
from datetime import datetime

from xx_VIs_transfomator import *
from T01_NoDataKill import NoData_kill

#建立植生指標名稱list，用來建立圖層儲存迴圈
VIs = ['NDVI', 'NDRE','NDI','RGI','ExG']
# VIs = ["gray","HSV-H","HSV-S","HSV-V","Lab-L","Lab-a","Lab-b","YCbCr-Y","YCbCr-Cb","YCbCr-Cr","YDbDr-Y","YDbDr-Db","YDbDr-Dr","YIQ-Y","YIQ-I","YIQ-Q","YPbPr-Y","YPbPr-Pb","YPbPr-Pr","YUV-Y","YUV-U","YUV-V","HSL-H","HSL-S","HSL-L","GRVI","NDI","RGI","ExG","ExR","ExGR","GDVI","GNDVI","GWDRVI","CIg","MSR_G","GSAVI","MSAVI-2","GOSAVI","GRDVI","NGI","NREI-1","NNIR","MDD-1","MNDI-1","MEVI-2","MNDRE","NRI-2","MCARI1","MCARI2","NDVI","RVI","DVI","RDVI","WDRVI","SAVI","OSAVI","MSAVI-1","TNDVI","MSR","VIopt","NLI","MNLI","NDVI_x_RVI","SAVI_x_SR","MTCI","DATT","NNIRI","NREI-2","NRI-1","MDD-2","MRESR","MNDI-2","MEVI-1","MNDRE2","RETVI","MCARI3","MCARI4","MTCARI","MRETVI","MCCCI","MCARI1/OSAVI","MCARI2/OSAVI","MTCARI/OSAVI","MCARI1/MRETVI","MTCARI/MRETVI","NDRE","RERVI","REDVI","RERDVI","REWDRVI","RESAVI","REOSAVI","MRESAVI","REVIopt","CIre","MSR_RE","RENDVI","RESR","MREDVI"]

#輸出csv欄位名稱
fieldnames = ['Id', 'site', 'date', 'B', 'G', 'R', 'RE', 'NIR', 'GC'] + [f'{VI}' for VI in VIs] + [f'{VI}_sum' for VI in VIs] + [f'{VI}_max' for VI in VIs] + [f'{VI}_min' for VI in VIs] + [f'{VI}_range' for VI in VIs] + [f'{VI}_std' for VI in VIs]

#給定開始時間
st_time = time.time()

#建立資料夾是否存在的判斷函數
def dir_maker(dir_path):
    dir_path_list = Path(dir_path).parts
    if Path(dir_path).is_dir()==False:
        for i in range(len(dir_path_list)):
            now_dir = "./" + "/".join([dir_path_list[j] for j in range(i + 1)])
            if Path(now_dir).is_dir()==False:
                os.mkdir(now_dir)
    else:
        pass

#判斷是否存在data與images資料夾
dir_maker(f'./04_Result/Data')
dir_maker(f'./04_Result/Image')
now = str(datetime.now())[:18].replace(':','').replace('/','').replace('-','').replace(' ','')
with open(f'./04_Result/Data/{now}_ROI_Result.csv','w',encoding='utf-8-sig',newline='') as f:
    writer = csv.DictWriter(f,fieldnames=fieldnames)
    #寫入欄位名稱，key值
    writer.writeheader()


#開啟csv檔案(新增資料格式)
with open(f'./04_Result/Data/{now}_ROI_Result.csv','a',encoding='utf-8-sig',newline='') as f:
    writer = csv.DictWriter(f,fieldnames=fieldnames)

    #把需要分析的原圖都先load出來，標定site跟data
    for image_path in glob.glob(f'./01_Tif/6band/BG/*band*.tif'):
        print(f'目標路徑：{image_path}')
        #建立要寫入csv的資料容器
        data = {}

        #TIF檔名稱
        tif_name = Path(image_path).stem
        #調查場域
        site = Path(image_path).parts[-2]
        #調查日期
        date = tif_name[:10].replace('-','')

        data['site'] = site
        data['date'] = date

        if int(date[:4]) != 2023:
            continue

        if date == '20231101':
            continue

        #辨識期作別
        if int(date[:4]) == 2022:
            if int(date[4:6]) < 7:
                Period = '111-1'
            else:
                Period = '111-2'
        elif int(date[:4]) == 2023:
            if int(date[4:6]) < 7:
                Period = '112-1'
            else:
                Period = '112-2'



        #讀取tif原圖為dataset
        dataset = gdal.Open(image_path,gdal.GA_ReadOnly)

        #讀取dataset當中的座標資訊
        ds = dataset.GetGeoTransform()

        #讀取tif原圖的影像資料，並且正規化
        img  = dataset.ReadAsArray() / 65535

        #將數值為1(即原本為65535之值)設定為0
        img[img==1.0]=0

        #讀取影像資料中的各波段，以進行植生指標運算
        B = img[0,:,:]
        G = img[1,:,:]
        R = img[2,:,:]
        RE = img[3,:,:]
        NIR = img[4,:,:]

        #計算NDVI跟NDRE               
        NDVI = (NIR - R) / (NIR + R + (10**(-10)))
        NDRE = (NIR - RE) / (NIR + RE + (10**(-10)))
        NDI = (G - R) / (G + R + (10)**(-10))
        RGI =  R / (G + (10)**(-10))
        ExG = (2 * G - B - R) / (R + B + G + (10)**(-10))
        GBGAP = G - B
        NIRGRATIO = (NIR) / (G + (10)**(-10))

        #分別使用大津演算法二值化
        thresh1 = filters.threshold_otsu(NIR)
        thresh2 = filters.threshold_otsu(NDVI)
        thresh3 = filters.threshold_otsu(NDRE)
        thresh4 = filters.threshold_otsu(NDI)
        thresh5 = filters.threshold_otsu(RGI)
        thresh6 = filters.threshold_otsu(ExG)
        thresh7 = filters.threshold_otsu(GBGAP)
        thresh8 = filters.threshold_otsu(NIRGRATIO)

        #建立mask
        mask1 = (NIR >= thresh1) * 1.0
        mask2 = (NDVI >= thresh2) * 1.0
        mask3 = (NDRE >= thresh3) * 1.0
        mask4 = (NDI >= 0) * 1.0
        mask5 = (RGI >= thresh5) * 1.0
        mask6 = (GBGAP >= 0) * 1.0
        mask7 = (NIRGRATIO > 0) * (NIRGRATIO <= 2.8) * 1.0
        mask8 = (R != 0) * (G !=0) * (B != 0) * 1.0

        if date == '20220310':
            mask2 = (NDVI >= 0.4) * 1.0
            mask = mask2 * mask6 * mask8
        else:
            mask = mask1 * mask2 * mask3 * mask6 * mask8

        dir_maker(f'./04_Result/Image/{site}/ROI/{date}/')

        new_dataset = gdal.Translate(f'./04_Result/Image/{site}/ROI/{date}/mask.tif', dataset,bandList=[1],outputType=gdal.GDT_Float32) 
        new_dataset.WriteArray(mask)
        new_dataset = None

        day = date[4:8]

        #把植株定位後的Shapefile讀出來，每個shapefile的polygon當作一個分析單元，以每個場域為單元
        shapefile_path = glob.glob(f'./02_Shapefile/{site}/ROI/{Period}*.shp')[0]

        #讀取整個shapefile
        input_shapefile = ogr.Open(shapefile_path)
        input_layer = input_shapefile.GetLayer(0)

        #建立調查區域與日期的檔案
        dir_maker(f'./04_Result/Image/{site}/{date}')

        print(f'{site} {date}切割小塊的shapefile')
        #依照shapefile裡面的數量執行迴圈內容
        print(input_layer.GetFeatureCount())

        for s in trange(input_layer.GetFeatureCount()):

            #指定要輸出各別單株shapefile的檔案路徑位置
            shp_output_path = f'./04_Result/Shapefile/{site}/ROI/'

            #組合儲存路徑與名稱

            feature = input_layer.GetFeature(s)
            Id = feature.GetField('Field')
            save_shp_path = os.path.join(shp_output_path,(f"{Id}.shp"))

            if os.path.isfile(save_shp_path):
                continue

            #針對每個日期每個區域，建立用shapefile切下來的檔案資料夾
            dir_maker(f'./04_Result/Shapefile/{site}/ROI/')
            dir_maker(f'./04_Result/Image/{site}/ROI/{date}/tif')

            driver = ogr.GetDriverByName('ESRI Shapefile') ## a driver
            datasource = driver.CreateDataSource(save_shp_path)
            print(save_shp_path)
            layer = datasource.CreateLayer('layerName',geom_type=ogr.wkbPolygon)
            layer.CreateFeature(feature)

        print(f'{site} {date}小塊的shapefile切割完畢', '花費時間',round(time.time() - st_time,2),"秒")        

        #計時&重新計時
        print(f'{site} {date}切割小塊的tif檔案')
        st_time = time.time()

        #讀取該場域該日期的所有shapefile檔案
        shapefile_path_list = glob.glob(f'./04_Result/Shapefile/{site}/ROI/*.shp')

        #建立小塊tif檔案切割進度條
        for i in trange(len(shapefile_path_list)):
            #指定的shapefile路徑
            shp_path = shapefile_path_list[i]

            if Path(shp_path).stem == 'temp':
                continue

            #小塊tif檔案的輸出路徑
            tif_file_output_path = f'./04_Result/Image/{site}/ROI/{date}/tif/' + Path(shp_path).stem + '.tif'

            #獲得輸入shapefile的左上右下座標資料
            r = shapefile.Reader(shp_path)
            minX, minY, maxX, maxY = r.bbox

            #裁切tif檔案：6band
            # ds = gdal.Warp(
            #     tif_file_output_path,#輸出raster路徑
            #     image_path,#輸入raster檔案，用gdal.Open()開啟
            #     format = 'GTiff',#格式
            #     cutlineDSName = shp_path,#shapefile路徑
            #     cropToCutline=True,
            #     dstNodata = 0,
            # )

            ds = gdal.Warp(
                tif_file_output_path,#輸出raster路徑
                image_path,#輸入raster檔案，用gdal.Open()開啟
                format = 'GTiff',#格式
                cutlineDSName = shp_path,#shapefile路徑
                outputBounds = [minX, minY, maxX, maxY],
                dstNodata = 0,
            )

            NoData_kill(tif_file_output_path, tif_file_output_path)


        #計時&重新計時
        print(f'{site} {date}小塊tif切割結束', '花費時間',round(time.time() - st_time,2),"秒")
        st_time = time.time()

        #建立該場域該日期的所有樣區的tif檔案路徑list
        tif_file_path_list = glob.glob(f'./04_Result/Image/{site}/ROI/{date}/tif/*.tif')

        #建立分析植生指標的進度條
        for i in trange(len(tif_file_path_list)):
            #指定要分析的tif檔案
            tif_file_path = tif_file_path_list[i]

            #該tif檔案名稱
            tif_name = Path(tif_file_path).stem

            data['Id'] = tif_name
            #用檔案名稱把資訊讀出來
            Id = tif_name.split("_")[0]

            #讀取原圖tif檔的資訊
            dataset = gdal.Open(tif_file_path,gdal.GA_ReadOnly)
            ds = dataset.GetGeoTransform()

            #讀取tif檔案原圖圖片資訊
            image  = dataset.ReadAsArray()/ 65535
            image[image==1.0]=0

            #讀取所有band
            B = image[0,:,:]
            G = image[1,:,:]
            R = image[2,:,:]
            RE = image[3,:,:]
            NIR = image[4,:,:]

            #計算NDVI跟NDRE                
            NDVI = (NIR - R) / (NIR + R + (10**(-10)))
            NDRE = (NIR - RE) / (NIR + RE + (10**(-10)))
            NDI = (G - R) / (G + R + (10)**(-10))
            RGI =  R / (G + (10)**(-10))
            ExG = (2 * G - B - R) / (R + B + G + (10)**(-10))
            GBGAP = G - B
            NIRGRATIO = (NIR) / (G + (10)**(-10))

            #建立mask
            mask1 = (NIR >= thresh1) * 1.0
            mask2 = (NDVI >= thresh2) * 1.0
            mask3 = (NDRE >= thresh3) * 1.0
            mask4 = (NDI >= 0) * 1.0
            mask5 = (RGI >= thresh5) * 1.0
            mask6 = (GBGAP >= 0) * 1.0
            mask7 = (NIRGRATIO > 0) * (NIRGRATIO <= 2.8) * 1.0
            mask8 = (R != 0) * (G !=0) * (B != 0) * 1.0

            #建立植體遮罩
            if date == '20220310':
                mask2 = (NDVI >= 0.4) * 1.0
                mask = mask2 * mask6 * mask8
            else:
                mask = mask1 * mask2 * mask3 * mask6 * mask8

            dir_maker(f'./04_Result/Image/{site}/ROI/{date}/mask')
            # print(tif_name)
            new_dataset = gdal.Translate(f'./04_Result/Image/{site}/ROI/{date}/mask/{tif_name}.tif', dataset,bandList=[1],outputType=gdal.GDT_Float32) 
            new_dataset.WriteArray(mask)
            new_dataset = None

            #計算冠層覆蓋率，運算邏輯：計算原圖(image)的隨便一個波段所有非零像素植體的數量，再除以原圖該波段的所有非零像素植體數量
            GC_area = np.count_nonzero(image[0,:,:] * mask)
            total_area = np.count_nonzero(image[0,:,:])
            GC_rate =(GC_area / total_area)

            B = B * mask
            G = G * mask
            R = R * mask
            RE = RE * mask
            NIR = NIR * mask

            data['B'] = np.nanmean(B)
            data['G'] = np.nanmean(G)
            data['R'] = np.nanmean(R)
            data['RE'] = np.nanmean(RE)
            data['NIR'] = np.nanmean(NIR)

            #將GC_rate作為冠層覆蓋濾放入data
            data['GC'] = GC_rate

            # VIs = VIs_tranfer(B, G, R, RE, NIR)
            dir_maker(f'./04_Result/Image/{site}/ROI/{date}/VI/')

            #計算NDVI數列的非nan平均值，如果還是nan值，則平均值與其他統計量一樣設定為0
            for VI in VIs:
                VI_img = eval(VI)
                VI_img_nozero = VI_img[VI_img!=0.0]
                if VI == "RGI":
                    VI_img_nozero = VI_img_nozero[VI_img_nozero<=1.0]

                data[f'{VI}'] = np.nanmean(VI_img_nozero)
                if str(data[f'{VI}']) == 'nan':
                    data[f'{VI}'] = 0
                    data[f'{VI}_sum'] = 0
                    data[f'{VI}_max'] = 0
                    data[f'{VI}_min'] = 0
                    data[f'{VI}_range'] = 0
                    data[f'{VI}_std'] = 0

                else:
                    data[f'{VI}_sum'] = np.nansum(VI_img_nozero)
                    data[f'{VI}_max'] = np.nanmax(VI_img_nozero)
                    data[f'{VI}_min'] = np.nanmin(VI_img_nozero)
                    data[f'{VI}_range'] = np.nanmax(VI_img_nozero) - np.nanmin(VI_img_nozero)
                    data[f'{VI}_std'] = np.nanstd(VI_img_nozero)

                #確認輸出路徑資料夾是否存在
                try:
                    #設定要儲存的tif檔案路徑，以原圖的dataset設定，並放入一個band，輸出格式為float32
                    new_dataset = gdal.Translate(f'./04_Result/Image/{site}/ROI/{date}/VI/{VI}_{tif_name}.tif', dataset,bandList=[1],outputType=gdal.GDT_Float32) 
                    new_dataset.WriteArray(VI_img)
                    new_dataset = None
                except:
                    pass

            #將一小塊tif檔案的分析結果寫入csv檔案中
            writer.writerow(data)

    #分析完畢，計時
    print(f'{site} {date}分析結束', '花費時間',round(time.time() - st_time,2),"秒")