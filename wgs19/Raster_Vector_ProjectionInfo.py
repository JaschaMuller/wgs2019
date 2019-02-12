import os, sys, re, math
from osgeo import gdal, ogr, osr
import numpy as np
import pandas as pd

gdal.UseExceptions()
ogr.UseExceptions()
osr.UseExceptions()

def wkt2epsg(wkt, forceProj4 = False, cstypeStat = False):
    '''
    Function Def: wkt2epsg(wkt, forceProj4 = False, cstypeStat = False)

    Adapted from source:
    source: https://gis.stackexchange.com/questions/20298/is-it-possible-to-get-the-epsg-value-from-an-osr-spatialreference-class-using-th
    
    wkt         ~  Porjection information in the wkt format [REQUIRED]
    forceProj4  ~  This is if there is a problem and you want to force a Proj4 text as a return value [OPTIONAL]
    cstypeStat  ~  If this is specified as True, the function will only return if the projection is a Geographic or a projected coordinate system. [OPTIONAL]

    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
    
    returns either:    ~ list containing epsg code, if found, and a confidence ranking (1 to 6, with 1 being the highes confidence) .e.g ['EPSG:32734', 1]
                       ~ None, if the EPSG code could not be found.
                       ~ cstype is True, the function will only return 'GEOGCS' or 'PROJCS'                 
    '''
                       
    rel_dir = os.path.dirname(os.path.realpath(__file__))
    rel_dir = rel_dir.replace('\\', '/')
    rel_dir = rel_dir + '/'
    epsg = rel_dir + 'epsg'
    esri = rel_dir + 'esri'
    code= None
    
    p_in = osr.SpatialReference()
    s = p_in.ImportFromWkt(wkt)

    if s == 5:  # invalid WKT
        return None
    if p_in.IsLocal() == 1:             # this is a local definition
        return p_in.ExportToWkt()
    if p_in.IsGeographic() == 1:        # This is a geographic srs
        cstype = 'GEOGCS'
    else:
        cstype = 'PROJCS'               # This is a projected srs

    if cstypeStat == True:
        return cstype
    
    an = p_in.GetAuthorityName(cstype)
    ac = p_in.GetAuthorityCode(cstype)

    if an is not None and ac is not None:   # return EPSG code
        return ['%s:%s' % (p_in.GetAuthorityName(cstype), p_in.GetAuthorityCode(cstype)),1]
    else: # try brute force appproach by grokking proj epsg definition file
        p_out = p_in.ExportToProj4()
        
        if p_out:
            if forceProj4 is True:
                return [p_out,6]

            f = open(epsg)
            for line in f:       
                if line.find(p_out) != -1:
                    m = re.search('<(\\d+)>', line)
                    if m:
                        code = m.group(1)
                        break                     
            if code:
                return ['EPSG:%s' % code,2]
            else:               # did not find the spatial reference in the epsg file, now systematically searching the esri file which has a slightly different format
                p_out_ellps = None
                p_out_datum = None
                   
                if p_out.find('+ellps=') != -1:
                    elp_res = re.search(r'\+ellps=(\S+)', p_out)
                    p_out_ellps = elp_res.group()
                    p_out_ellps = p_out_ellps.strip()

                if p_out.find('+datum=') != -1:
                    dt_res = re.search(r'\+datum=(\S+)', p_out)
                    p_out_datum = dt_res.group()
                    p_out_datum = p_out_datum.strip()
                    
                if p_out_ellps == None and p_out_datum == None:
                    if p_out.find('+no_defs') != -1:
                        nodef_1 = re.search(r'(.+)(?=\+no_def)', p_out)
                        new_p_out = nodef_1.group(1)
                        new_p_out = new_p_out.strip()
                    else:
                        nodef_1 = re.search(r'(.+)(?=no_defs)', p_out)
                        new_p_out = nodef_1.group(1)
                        new_p_out = new_p_out.strip()
                        
                elif p_out_ellps == None and p_out_datum != None:
                    rslt = re.search(r'(.+)(?=\+datum=)', p_out)
                    new_p_out = rslt.group(1)
                    new_p_out = new_p_out.strip()
                else:
                    rslt = re.search(r'(.+)(?=\+ellps=)', p_out)
                    new_p_out = rslt.group(1)
                    new_p_out = new_p_out.strip()
                
                f.close()
                f = open(esri)
                for line in f:

                    if  line[0] != '#' and (line[0:8] == '<102022>' or line[0:6] == '<4326>' or line[0:7] == '<32734>'): 

                        line_ellps = None
                        line_datum = None
                        spatial_code = 'None'   
                        if line.find('+ellps=') != -1:
                            elp_res = re.search(r'\+ellps=(\S+)', line)
                            line_ellps = elp_res.group()

                        if line.find('+datum=') != -1:
                            dt_res = re.search(r'\+datum=(\S+)', line)
                            line_datum = dt_res.group()

                        if line_ellps == None and line_datum == None:
                            if line.find('+no_defs') != -1:
                                nodef_1 = re.search(r'<\d+>(.+)(?=\+no_def)', line)
                                new_line = nodef_1.group(1)
                                new_line = new_line.strip()
                            else:
                                nodef_1 = re.search(r'<\d+>(.+)(?=no_defs)', line)
                                new_line = nodef_1.group(1)
                                new_line = new_line.strip()
                                
                        elif line_ellps == None and line_datum != None:
                            rslt = re.search(r'<\d+>(.+)(?=\+datum=)', line)
                            new_line = rslt.group(1)
                            new_line = new_line.strip()
                        else:
                            rslt = re.search(r'<\d+>(.+)(?=\+ellps=)', line)
                            new_line = rslt.group(1)
                            new_line = new_line.strip()

                        if len(new_line) == 0:      # This is not 100% water tight becuase the know EPSG code are usually catched in either the first of second phases, however more logical checks can be placed here
                            new_line = 'NoThInG'
                        elif new_line == '+proj=longlat':
                            if line_ellps != p_out_ellps and line_datum != p_out_datum:
                                new_line = 'NoThInG'
                                                    
                        if new_line == new_p_out:
                            if (p_out_datum == line_datum and p_out_datum != None and line_datum != None) and \
                            (p_out_ellps == line_ellps and p_out_ellps != None and line_ellps != None):
                                m = re.search('<(\\d+)>', line)
                                if m:
                                    code = m.group(1)
                                    spatial_code = 'D&E'
                                    break
                            elif (p_out_datum == line_datum and p_out_datum != None and line_datum != None) or \
                            (p_out_ellps == line_ellps and p_out_ellps != None and line_ellps != None):
                                m = re.search('<(\\d+)>', line)
                                if m:
                                    code = m.group(1)
                                    spatial_code = 'D|E'
                                    break
                            else:
                                m = re.search('<(\\d+)>', line)
                                if m:
                                    code = m.group(1)
                                    spatial_code = '!D!E'
                                    break
                            
                if code:
                    if spatial_code == 'D&E':                       
                        return ['EPSG:%s' %code, 3]
                    elif spatial_code == 'D|E':
                        return ['EPSG:%s' %code, 4]
                    elif spatial_code == '!D!E':
                        return ['EPSG:%s' %code, 5]
                else:
                    return None
        else:
            return None
                
def wkt2cstype(wkt, forceProj4 = False, cstypeStat = False):
    
    p_in = osr.SpatialReference()
    s = p_in.ImportFromWkt(wkt)

    if s == 5:  # invalid WKT
        return 'None'
    if p_in.IsLocal() == 1:             # this is a local definition
        return p_in.ExportToWkt()
    if p_in.IsGeographic() == 1:        # This is a geographic srs
        cstype = 'GEOGCS'
    else:
        cstype = 'PROJCS'               # This is a projected srs

    return cstype



def get_raster_wkt(input_raster):
    '''
    Function Def: get_raster_wkt(input_raster)
    
    input_raster ~ the path and name of the input raster image (e.g. .tif) [REQUIRED]

    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    returns     ~ The function returns the spatial reference object of input_raster
                  in an well known text (WKT) format
    '''
    input_layer = gdal.Open(input_raster)
    prj_ras = input_layer.GetProjection()
    srs_ras = osr.SpatialReference(wkt=prj_ras)
    return srs_ras.ExportToWkt()

def get_raster_epsg(inraster):
    '''
    Function Def: get_raster_epsg(inraster)

    inraster     ~ path to a raster file [REQUIRED]
    - - - - - - - - - - - - - - - - - - - -

    returns      ~ list containing epsg code, if found, and a confidence ranking (1 to 6, with 1 being the highes confidence) .e.g ['EPSG:32734', 1]
                 ~ None, if the EPSG code could not be found.
    
    '''
    raster_epsg = wkt2epsg(get_raster_wkt(inraster), forceProj4 = False, cstypeStat = False)
    return raster_epsg


def project_raster(inRas, outRas, outEpsg, outPixelsize=None, NoDataValue = -9999, snap_clip=False, Snap_Raster=None, resample_alg = 'Nearest', data_type = None):   
    '''
    Function Def: project_raster(inRas, outRas, outEpsg, outPixelsize=None, NoDataValue = -9999, snap_clip=False, Snap_Raster=None, resample_alg = 'Nearest')
    
    inRas        ~ the path and name of the input raster (e.g. .tif) [REQUIRED]
    outRas       ~ the path and name of the desired output projected raster [REQUIRED]
    outEpsg      ~ the output coordinate system in EPSG formate.g. 'EPSG:4326', this will be overwritten if snap is True [REQUIRED]
    outPixelsize ~ the pixel size of outRas in meters (the default is its own pixel size), if snap is true, the pixel size will be overwritten with snap raster pixel size [OPTIONAL]
                 ~ if you specify the pixelsize to somehting other than its owns pixel size, the function can serve as a resampling fucntion
    NoDataValue  ~ The Nodata value of the output raster
    snap_clip    ~ weather or not to run the snap raster process, the default is FALSE, the pixel size and extent of the snap raster will overwrite outPixelsize [OPTIONAL]
    resample_alg ~ What algorithm to use for the resampling process, the default is Nearest Neigbour, can be, exclusively: 'Nearest', 'Bilinear', 'Cubic', 'Cubic_spline' [OPTIONAL]

    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    returns      ~ This function does not return anything, but outputs a projected raster, either snapped or not at outRas
   
    '''
    
    ePsG_list = get_raster_epsg(inRas)

    if ePsG_list == None:
        raise Exception('PROJECTION ERROR: The input coordinate system (EPSG code) could not be defined. The file most likely contains a custom projection made in a geospatial program without an known EPSG code\n'  + \
                        'Please re-project the file: ' + inRas + ' to a defined coordinate system then try again')
        
    
    epsg_from = ePsG_list[0]
    epsg_from_num = int(epsg_from[5:len(epsg_from)])
    from_coord = osr.SpatialReference()
    from_coord.ImportFromEPSG(epsg_from_num)

    if snap_clip == True:                           ##### STILL SORT OUT DATA TYPE 8BIT 32BIT 
        SnP_epsg = get_raster_epsg(Snap_Raster)
        SnP_epsg = SnP_epsg[0]
        epsg_to_num = int(outEpsg[5:len(SnP_epsg)]) 
        to_coord = osr.SpatialReference()
        to_coord.ImportFromEPSG(epsg_to_num)
        coord_trans = osr.CoordinateTransformation(from_coord, to_coord)

        inData = gdal.Open(inRas)
        gt = inData.GetGeoTransform()
        xRS = inData.RasterXSize
        yRS = inData.RasterYSize

        snapRaster = gdal.Open(Snap_Raster)
        SN_xRS = snapRaster.RasterXSize
        SN_yRS = snapRaster.RasterYSize
        gtSNAP = snapRaster.GetGeoTransform()   # Geotransformation # ulx[0], xres[1], xskew[2], uly[3], yskew[4], yres[5]         
        
        outPixelsize_x = gtSNAP[1]
        outPixelsize_y = gtSNAP[5]
        
        proj_ras = gdal.GetDriverByName( 'MEM' ).Create( "", SN_xRS, SN_yRS, 1, gdal.GDT_Float32)
        proj_ras.SetGeoTransform(gtSNAP)
        proj_ras.SetProjection(to_coord.ExportToWkt())
        proj_ras_band = proj_ras.GetRasterBand(1)       
        proj_ras_band.SetNoDataValue(NoDataValue)   

        nodata_array = np.full((SN_yRS, SN_xRS), NoDataValue)   # NoData HACK
        proj_ras_band.WriteArray(nodata_array)                  # NoData HACK
        
        if resample_alg == 'Nearest':
            res = gdal.ReprojectImage(inData, proj_ras, from_coord.ExportToWkt(), to_coord.ExportToWkt(), gdal.GRA_NearestNeighbour) # GRA_NearestNeighbour, GRA_CubicSpline, GRA_Bilinear, GRA_Cubic, GRA_Lanczos  
        elif resample_alg == 'Bilinear':
            res = gdal.ReprojectImage(inData, proj_ras, from_coord.ExportToWkt(), to_coord.ExportToWkt(), gdal.GRA_Bilinear)
        elif resample_alg == 'Cubic':
            res = gdal.ReprojectImage(inData, proj_ras, from_coord.ExportToWkt(), to_coord.ExportToWkt(), gdal.GRA_Cubic)
        elif resample_alg == 'Cubic_spline':
            res = gdal.ReprojectImage(inData, proj_ras, from_coord.ExportToWkt(), to_coord.ExportToWkt(), gdal.GRA_CubicSpline)
            
        
        gdal.GetDriverByName('GTiff').CreateCopy(outRas, proj_ras)
        inData=None
        proj_ras=None
        
    
    else:
        epsg_to_num = int(outEpsg[5:len(outEpsg)])  
        to_coord = osr.SpatialReference()
        to_coord.ImportFromEPSG(epsg_to_num)

        inData = gdal.Open(inRas)
        gt = inData.GetGeoTransform()
        xRS = inData.RasterXSize
        yRS = inData.RasterYSize

        coord_trans = osr.CoordinateTransformation(from_coord, to_coord)

        (ulx, uly, ulz) = coord_trans.TransformPoint(gt[0], gt[3])
        (lrx, lry, lrz) = coord_trans.TransformPoint(gt[0] + gt[1]*xRS, gt[3] + gt[5]*yRS)

        inCStype = wkt2epsg(get_raster_wkt(inRas), cstypeStat = True)
        outCStype = wkt2epsg(to_coord.ExportToWkt(), cstypeStat = True)
        if outPixelsize is None:    
            if inCStype == outCStype:
                outPixelsize = gt[1]
            elif inCStype == 'GEOGCS' and outCStype == 'PROJCS':
                outPixelsize = int((lrx - ulx)/xRS)
            else:
                outPixelsize = (lrx - ulx)/xRS
        else:
            if outCStype == 'GEOGCS':
                if inCStype == 'GEOGCS':

                    ulx_trns = gt[0]
                    uly_trns = gt[3]

                    inSpatial_ref = osr.SpatialReference()
                    inSpatial_ref.ImportFromEPSG(epsg_from_num)

                    utm_band = str((math.floor((ulx_trns + 180)/6)%60)+1)
                    if  len(utm_band) == 1:
                        utm_band = '0'+ utm_band
                    if uly_trns >= 0:
                        epsg_code = '326' + utm_band
                    else:
                        epsg_code = '327' + utm_band

                    epsg_code = int(epsg_code)
                    outSpatial_ref = osr.SpatialReference()
                    outSpatial_ref.ImportFromEPSG(epsg_code)

                    Geo_to_Proj = osr.CoordinateTransformation(inSpatial_ref, outSpatial_ref)        

                    (UlX, UlY, UlZ) = Geo_to_Proj.TransformPoint(gt[0], gt[3])
                    (LrX, LrY, LrZ) = Geo_to_Proj.TransformPoint(gt[0] + gt[1]*xRS, gt[3] + gt[5]*yRS)

                    num_of_pixels = int(round((LrX - UlX)/outPixelsize, 0))
                    outPixelsize = abs(((gt[0] + gt[1]*xRS) - ulx_trns)/num_of_pixels)
                    
                elif inCStype == 'PROJCS':
                    
                    num_of_pixels = ((gt[0] + gt[1]*xRS) - ulx_trns)/pixel_size 
                    outPixelsize = (lrx - ulx)/num_of_pixels
               
                
        if data_type == '8bit':
            #NoDataValue = 0                                                                                                                         #8BIT HERE
            proj_ras = gdal.GetDriverByName( 'MEM' ).Create( "", int((lrx-ulx)/outPixelsize), int((uly-lry)/outPixelsize), 1, gdal.GDT_Byte)        #8BIT HERE
        elif data_type == '32bit':
            proj_ras = gdal.GetDriverByName( 'MEM' ).Create( "", int((lrx-ulx)/outPixelsize), int((uly-lry)/outPixelsize), 1, gdal.GDT_Float32)     #32BIT HERE
        else:
            proj_ras = gdal.GetDriverByName( 'MEM' ).Create( "", int((lrx-ulx)/outPixelsize), int((uly-lry)/outPixelsize), 1, gdal.GDT_Float32)     #32BIT HERE

        new_gt = (ulx, outPixelsize, gt[2], uly, gt[4], -outPixelsize)
        proj_ras.SetGeoTransform(new_gt)
        proj_ras.SetProjection(to_coord.ExportToWkt())
        proj_ras_band = proj_ras.GetRasterBand(1)
        proj_ras_band.SetNoDataValue(NoDataValue)
        
        nodata_array = np.full((int((uly-lry)/outPixelsize), int((lrx-ulx)/outPixelsize)), NoDataValue)     # NoData HACK
        proj_ras_band.WriteArray(nodata_array)                                                              # NoData HACK

        if resample_alg == 'Nearest':
            res = gdal.ReprojectImage(inData, proj_ras, from_coord.ExportToWkt(), to_coord.ExportToWkt(), gdal.GRA_NearestNeighbour) # GRA_NearestNeighbour, GRA_CubicSpline, GRA_Bilinear, GRA_Cubic, GRA_Lanczos  
        elif resample_alg == 'Bilinear':
            res = gdal.ReprojectImage(inData, proj_ras, from_coord.ExportToWkt(), to_coord.ExportToWkt(), gdal.GRA_Bilinear) # GRA_NearestNeighbour, GRA_CubicSpline, GRA_Bilinear, GRA_Cubic, GRA_Lanczos  
        elif resample_alg == 'Cubic':
            res = gdal.ReprojectImage(inData, proj_ras, from_coord.ExportToWkt(), to_coord.ExportToWkt(), gdal.GRA_Cubic) # GRA_NearestNeighbour, GRA_CubicSpline, GRA_Bilinear, GRA_Cubic, GRA_Lanczos  
        elif resample_alg == 'Cubic_spline':
            res = gdal.ReprojectImage(inData, proj_ras, from_coord.ExportToWkt(), to_coord.ExportToWkt(), gdal.GRA_CubicSpline) # GRA_NearestNeighbour, GRA_CubicSpline, GRA_Bilinear, GRA_Cubic, GRA_Lanczos  
            

        gdal.GetDriverByName('GTiff').CreateCopy(outRas, proj_ras)
        inData = None
        proj_ras = None
        
    #if snap != False:
    #    print('\tSnapping ...')
    #    snap_raster_prj(snap, proj_ras, outRas, to_coord.ExportToWkt())
    #    proj_ras=None
    #else:
    #    gdal.GetDriverByName('GTiff').CreateCopy(outRas, proj_ras)
    #    proj_ras=None


# Alternative relative path
rel_dir = os.path.dirname(os.path.realpath(__file__))
rel_dir = rel_dir.replace('\\', '/')
rel_dir = rel_dir + '/'

print('RelDir 1: ', rel_dir1)

rel_dir = os.path.dirname(__file__) # path of the epsg is always relative to this file
rel_dir = rel_dir + '/'

print('RelDir 2', rel_dir)
