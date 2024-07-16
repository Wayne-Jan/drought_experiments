import cv2
from osgeo import gdal
import scipy.interpolate
import numpy as np
 
def read_img(filename):
    dataset=gdal.Open(filename)
 
    im_width = dataset.RasterXSize
    im_height = dataset.RasterYSize
 
    im_geotrans = dataset.GetGeoTransform()
    im_proj = dataset.GetProjection()
    im_data = dataset.ReadAsArray(0,0,im_width,im_height)
 
    del dataset 
    return im_proj,im_geotrans,im_width, im_height,im_data
 
 
def write_img(filename, im_proj, im_geotrans, im_data):
    if 'int8' in im_data.dtype.name:
        datatype = gdal.GDT_Byte
    elif 'int16' in im_data.dtype.name:
        datatype = gdal.GDT_UInt16
    else:
        datatype = gdal.GDT_Float32
 
    if len(im_data.shape) == 3:
        im_bands, im_height, im_width = im_data.shape
    else:
        im_bands, (im_height, im_width) = 1,im_data.shape 
 
    driver = gdal.GetDriverByName("GTiff")
    dataset = driver.Create(filename, im_width, im_height, im_bands, datatype)
 
    dataset.SetGeoTransform(im_geotrans)
    dataset.SetProjection(im_proj)
 
    if im_bands == 1:
        dataset.GetRasterBand(1).WriteArray(im_data)
    else:
        for i in range(im_bands):
            dataset.GetRasterBand(i+1).WriteArray(im_data[i])
 
    del dataset
 
def NoData_kill(in_path, out_path):
    # data = cv2.imread("cq_test.tif")
    print(in_path)
    im_proj,im_geotrans,im_width, im_height,im_data = read_img(in_path)
    mask = np.isnan(im_data)
    c, w, h = mask.shape
    mask_list = []
    for i in range(c):
        if mask[i].__contains__(True):
            mask_list.append(mask[i])
 
    for m in mask_list:
        m = m + 0
        m = np.uint8(m)
        inpainted_img = cv2.inpaint(im_data, m, inpaintRadius=3, flags=cv2.INPAINT_TELEA)
        im_data = inpainted_img
    # cv2.imwrite('./fixed2.tif', data)
    write_img(out_path, im_proj, im_geotrans, im_data)
 
#too slow
def NoData_kill2(in_path, out_path):
    im_proj,im_geotrans,im_width, im_height, data = read_img(in_path)
    data = data.transpose(2,1,0)
 
    # a boolean array of (width, height) which False where there are missing values and True where there are valid (non-missing) values
    mask = ~( (data[:,:,0] == 255) & (data[:,:,1] == 255) & (data[:,:,2] == 255) )
 
    # array of (number of points, 2) containing the x,y coordinates of the valid values only
    xx, yy = np.meshgrid(np.arange(data.shape[1]), np.arange(data.shape[0]))
    xym = np.vstack( (np.ravel(xx[mask]), np.ravel(yy[mask])) ).T
 
    # the valid values in the first, second, third color channel,  as 1D arrays (in the same order as their coordinates in xym)
    data0 = np.ravel( data[:,:,0][mask] )
    data1 = np.ravel( data[:,:,1][mask] )
    data2 = np.ravel( data[:,:,2][mask] )
 
    # three separate interpolators for the separate color channels
    interp0 = scipy.interpolate.NearestNDInterpolator( xym, data0 )
    interp1 = scipy.interpolate.NearestNDInterpolator( xym, data1 )
    interp2 = scipy.interpolate.NearestNDInterpolator( xym, data2 )
 
    # interpolate the whole image, one color channel at a time    
    result0 = interp0(np.ravel(xx), np.ravel(yy)).reshape( xx.shape )
    result1 = interp1(np.ravel(xx), np.ravel(yy)).reshape( xx.shape )
    result2 = interp2(np.ravel(xx), np.ravel(yy)).reshape( xx.shape )
 
    # combine them into an output image
    result = np.dstack( (result0, result1, result2) )
    result = result.transpose(2,1,0)
    write_img(out_path, im_proj, im_geotrans, result)
 
if __name__ == "__main__":
    NoData_kill('./003_Field/ARI/20200929/tif/D1_N4_CP_TNG71.tif', './test.tif')