##:::::::::::::::::::::::::::::::::
## Project: Receiving Sentinel 2 images and finding  various spatial indices 
## Author: Aleksandar Tendjer
## Description:
##Different earth surfaces reflect the solar radiation differently and each raster layer represents how much incident solar radiation is reflected at a particular wavelength bandwidth.
##Bands 6 to 9 are in the Near Infrared Range (NIR)
##10 m resolution band 2, band 3, band 4 and band 8
##20 m resolution band 5, band 6, band 7, band 11 and band 12
##60 m resolution band 1, band 9 and band 10
##LINK to the image from scihub : https://scihub.copernicus.eu/dhus/odata/v1/Products('c92416ce-4477-496c-8552-d96f2fb10519')/$value
##:::::::::::::::::::::::::::::::::
#install.packages(c("sp", "rgdal", "raster", "viridis", "rasterVis"))

### Load packages
library(sp)
library(rgdal)
library(raster)
library(ggplot2)
#library that can represent images for people with low eye sight
library(viridis)
library(rasterVis)
library(sen2r)



setwd("C:\\Users\\WhyteRabyte\\Desktop\\FERIPR~2\\1SEMES~1\\ALGORI~1\\vaje\\")
path_Base_Pictures_10m<-"C:\\Users\\WhyteRabyte\\Desktop\\FERIPR~2\\1SEMES~1\\ALGORI~1\\vaje\\slika\\S2A_MS~2\\S2A_MS~1.SAF\\GRANULE\\L2A_T3~1\\IMG_DATA\\R10m\\"
path_Base_Pictures_20m<-"C:\\Users\\WhyteRabyte\\Desktop\\FERIPR~2\\1SEMES~1\\ALGORI~1\\vaje\\slika\\S2A_MS~2\\S2A_MS~1.SAF\\GRANULE\\L2A_T3~1\\IMG_DATA\\R20m\\"
path_Base_Pictures_60m<-"C:\\Users\\WhyteRabyte\\Desktop\\FERIPR~2\\1SEMES~1\\ALGORI~1\\vaje\\slika\\S2A_MS~2\\S2A_MS~1.SAF\\GRANULE\\L2A_T3~1\\IMG_DATA\\R60m\\"

###########################################################################
###########################################################################
###                                                                     ###
###                           SECTION 1:                                ###
###                           FUNCTION DEFINITION                       ###
###                                                                     ###
###########################################################################
###########################################################################
##########HELP FUNCTIONS##########

raster2matrix<-function(input_raster){
     img_dim=dim(input_raster)
     #create array from raster
       img_array=as.array(values(input_raster))
       #new matrix
         img_matrix=matrix(as.numeric(0),nrow =img_dim[1],ncol=img_dim[2])
         k=img_dim[1]
         #go from last row to the 1st row w column goin from 1 to last 
           for(i in img_dim[1]:1){
               l=0
               for(j in 1:img_dim[2])
                 {
                     img_matrix[l,k]=img_array[(img_dim[1]*(i-1))+j]
                     #}else
                       # {
                       #img_matrix[i,j]=img_array[j]
                       
                       #}
                       l=l+1
                     }
               k=k-1
             }
         return(img_matrix)
       }
#FUNCTION TO CHANGE RESOLUTION-NEAREST NEIGHBOUR 
upscaling_fun<-function(img_matrix,new_dim){
  old_img=img_matrix
  old_dim=dim(img_matrix)
  
  
  # calculate row and column ratios
  
  old_dim=dim(old_img)
  row_ratio=old_dim[1]/new_dim[1]
  col_ratio=old_dim[2]/new_dim[2]
  #initialize new matrix
  new_img=matrix(as.numeric(0), nrow = new_dim[1], ncol = new_dim[2])
  #do nearest neighboor scaling
  for( i in 1:new_dim[1])
  {
    for( j in 1:new_dim[2])
    {
      srcX=ceiling(i*(row_ratio))
      srcY=ceiling(j*(col_ratio))
      
      if(srcX==as.numeric(0)){
        
        srcX=1
        
      }
      if(srcY==as.numeric(0))
      {
        srcY=1
      }
      
      new_img[i,j]=old_img[srcX,srcY]
    }
  }
  return(new_img)
}
##########FUNCTIONS FOR SPATIAL INDEX FINDING##########
#9,8
ndmi_sentinel<-function(band9,band8){
  #new calculations
  band9=raster::brick(band9)
  band8=raster::brick(band8)
  ndmi=overlay(band9, band8, fun=function(x,y){(x-y)/(x+y)})
  output_name = "ndmi.tiff"#paste( output_name,"ndvi",'.tif',sep='')
  
  #export the image to the working directory
  raster::writeRaster(ndmi, filename = output_name)
  
  return(ndmi)
  
}
ndmi_sentinel2<-function(band9,band8){
  ndmi=(band9-band8)/(band9+band8);
  return(ndmi)
}
#Modified soil-adjusted vegetation index (MSAVI and MSAVI2)
#4,8
msavi_sentinel<-function(red,nir){
  bands<-c(red,nir)
  print(bands)
  
  #new calculations
  red=raster::brick(red)
  nir=raster::brick(nir)
  
 # MSAVI = (2*NIR+1-√((2*NIR+1)^2 - 8*(NIR-RED)))/2 
  # MSAVI2 = (B08 + 1) - 0.5 * sqrt((2 * B08 - 1) ^ 2 + 8 * B04)) other Formula
  msavi_2 = overlay(nir, red, fun=function(x,y){(x+1)-0.5*sqrt((2*x-1)^2+8*y)})
  
  #setting up the export name
  #output_name <-list.files(path=path_Base_Pictures,full.names = FALSE,recursive=TRUE, pattern="B04_10m.jp2")
  #calculating the total length of the path so that we can count from the end
  #my_length<-nchar(output_name)
  #add that its ndvi
  #output_name<-substr(output_name,my_length-30,my_length-7)
  output_name = "msavi2.tiff"#paste( output_name,"ndvi",'.tif',sep='')
  
  #export the image to the working directory
  raster::writeRaster(ndvi, filename = output_name)
  
  return(msavi_2)
}

#4,8 respectively
ndvi_sentinel<-function(red,nir){
  bands<-c(red,nir)
  print(bands)
  
  #new calculations
  red=raster::brick(red)
  nir=raster::brick(nir)
  
  
  ndvi=(nir-red)/(red+nir)
  
  #setting up the export name
  #output_name <-list.files(path=path_Base_Pictures,full.names = FALSE,recursive=TRUE, pattern="B04_10m.jp2")
  #calculating the total length of the path so that we can count from the end
  #my_length<-nchar(output_name)
  #add that its ndvi
  #output_name<-substr(output_name,my_length-30,my_length-7)
  output_name = "ndvi.tiff"#paste( output_name,"ndvi",'.tif',sep='')
  
  #export the image to the working directory
  raster::writeRaster(ndvi, filename = output_name)
  
  return(ndvi)
}
#2,4,8 respectively
evi_sentinel<-function(band2,band4,band8){
  
  #new calculations
  evi=(band8-band4)/(band8+6*band4-7.5*band2+1)
  
  
  return(evi)
  
}
#function for ndbi index from bands 11 and 8
ndbi_sentinel<-function(band8,band11){
  band8=raster::brick(band8)
  band11=raster::brick(band11)
  ndbi=(band11-band8)/(band11+band8)
  output_name= "ndbi.tif"#paste( output_name,"ndvi",'.tif',sep='')
  
  #export the image to the working directory
  raster::writeRaster(ndbi, filename = output_name)
  
  return(ndbi)
  
}
msi_sentinel<-function(band8,band11){
  band8=raster::brick(band8)
  band11=raster::brick(band11)
  msi=overlay(band11, band8, fun=function(x,y){(x/y)})
  output_name= "msi.tiff"
  
  #export the image to the working directory
  raster::writeRaster(msi, filename = output_name)
  
  return(msi)
}
#function finding ndwi index from bands 3 and 11
ndwi_sentinel<-function(green,nir){
  #new calculations
  
  green=raster::brick(green)
  nir=raster::brick(nir)
  
  #find index
  ndwi=(green-nir)/(green+nir)
  
  
  output_name <-"ndwi.tiff" #list.files(path=image_folder,full.names = FALSE,recursive=TRUE, pattern="B03_10m.jp2")
  #calculating the total length of the path so that we can count from the end
  #my_length<-nchar(output_name)
  #output_name<-substr(output_name,my_length-27,my_length-12)
  #output_name <- paste('ndwi', output_name,'.tif',sep='')
  
  #export the image to the working directory
  raster::writeRaster(ndwi, filename = output_name)
  return(ndwi)
}
gndvi<-function(band8,band3){
  band8=raster::brick(band8)
  band3=raster::brick(band3)
  
  #find index
  gndvi=overlay(band8, band3, fun=function(x,y){(x-y)/(x+y)})
  
  
  output_name <-"gndvi.tiff" #list.files(path=image_folder,full.names = FALSE,recursive=TRUE, pattern="B03_10m.jp2")
  #calculating the total length of the path so that we can count from the end
  #my_length<-nchar(output_name)
  #output_name<-substr(output_name,my_length-27,my_length-12)
  #output_name <- paste('ndwi', output_name,'.tif',sep='')
  
  #export the image to the working directory
  raster::writeRaster(gndvi, filename = output_name)
  return(gndvi)
}
###########################################################################
###########################################################################
###                                                                     ###
###                           SECTION 1:                                ###
###                           DATA INPUT                                ###
###                                                                     ###
###########################################################################
###########################################################################

  b1_pic<-"T34UDV_20190922T094031_B01_60m.jp2"
  b2_pic<-"T34UDV_20190922T094031_B02_10m.jp2"
  b3_pic<-"T34UDV_20190922T094031_B03_10m.jp2"
  b4_pic<-"T34UDV_20190922T094031_B04_10m.jp2"
  b4_pic_10<-"T34UDV_20190922T094031_B04_10m.jp2"
  b5_pic<-"T34UDV_20190922T094031_B05_20m.jp2"
  b6_pic<-"T34UDV_20190922T094031_B06_20m.jp2"
  b7_pic<-"T34UDV_20190922T094031_B07_20m.jp2"
  b8_pic<-"T34UDV_20190922T094031_B08_10m.jp2"
  b8A_pic<-"T34UDV_20190922T094031_B8A_20m.jp2"
  b9_pic<-"T34UDV_20190922T094031_B09_60m.jp2"
  #b10_pic<-"T33TWK_20181010T095031_B10.jp2"
  b11_pic<-"T34UDV_20190922T094031_B11_20m.jp2"
  b12_pic<-"T34UDV_20190922T094031_B12_20m.jp2"
  
  b1_pic=paste(path_Base_Pictures_60m,b1_pic,sep="")
  b2_pic=paste(path_Base_Pictures_10m,b2_pic,sep="")
  b3_pic=paste(path_Base_Pictures_10m,b3_pic,sep="")
  b4_pic=paste(path_Base_Pictures_10m,b4_pic,sep="")
  b4_pic_10=paste(path_Base_Pictures_10m,b4_pic_10,sep="")
  b5_pic=paste(path_Base_Pictures_20m,b5_pic,sep="")
  b6_pic=paste(path_Base_Pictures_20m,b6_pic,sep="")
  b7_pic=paste(path_Base_Pictures_20m,b7_pic,sep="")
  b8_pic=paste(path_Base_Pictures_10m,b8_pic,sep="")
  b8A_pic=paste(path_Base_Pictures_20m,b8A_pic,sep="")
  b9_pic=paste(path_Base_Pictures_60m,b9_pic,sep="")
  #b10_pic=paste(path_Base_Pictures_60m,b10_pic,sep="")
  b11_pic=paste(path_Base_Pictures_20m,b11_pic,sep="")
  b12_pic=paste(path_Base_Pictures_20m,b12_pic,sep="")
  
  #read jp2 picture and transform into a raster pic 
  b1_raster<-raster(readGDAL(b1_pic))
  
  b2_raster<-raster(readGDAL(b2_pic))
  b3_raster<-raster(readGDAL(b3_pic))
  b4_raster<-raster(readGDAL(b4_pic))
  b4_raster_10<-raster(readGDAL(b4_pic_10))
  b5_raster<-raster(readGDAL(b5_pic))
  b6_raster<-raster(readGDAL(b6_pic))
  b7_raster<-raster(readGDAL(b7_pic))
  b8_raster<-raster(readGDAL(b8_pic))
  b8A_raster<-raster(readGDAL(b8A_pic))
  b9_raster<-raster(readGDAL(b9_pic))
  #b10_raster<-raster(readGDAL(b10_pic))
  b11_raster<-raster(readGDAL(b11_pic))
  b12_raster<-raster(readGDAL(b12_pic))
  #or translate  gdal_translate("T33XVJ_20170803T125711_B01.jp2", "B01.tif")
  #and then read gdal 
  ###########################################VISUALISATION MANIPULATION####################################
  summary(raster_b1)
  plot(raster_b1)
  image(raster_b1, col= viridis_pal(option="D")(10), main="Sentinel 2-Czorsztyn area")
  
  #STACK RGB- show rgb picture
  windows()
  pic_RGB <- stack(list(b4_raster, b3_raster, b2_raster))              # creates raster stack
  plotRGB(pic_RGB, axes = TRUE, stretch = "lin", main = "Sentinel RGB colour composite")
  
  windows()
  gplot(b8_raster) +
    geom_raster(aes(x = x, y = y, fill = value)) +
    # value is the specific value (of reflectance) each pixel is associated with
    scale_fill_viridis_c() +
    coord_quickmap() +
    ggtitle("sea and mountains") +
    xlab("Longitude") +
    ylab("Latitude") +
    theme_classic() +   					    # removes defalut grey background
    theme(plot.title = element_text(hjust = 0.5),             # centres plot title
          text = element_text(size=20),		       	    # font size
          axis.text.x = element_text(angle = 90, hjust = 1))  # rotates x axis text
  
  
  t<-stack(b1_raster,b2_raster,b3_raster,b4_raster,b5_raster,b6_raster,b7_raster,b8_raster,b9_raster,b10_raster,b11_raster,b12_raster)
  windows()
  gplot(t) +
    geom_raster(aes(x = x, y = y, fill = value))+
    scale_fill_viridis_c() +
    facet_wrap(~variable) +
    coord_quickmap()+
    ggtitle("Sentinel 2 image, all bands, raster plots") +
    xlab("Longitude") +
    ylab("Latitude") +
    theme_classic() +
    theme(text = element_text(size=20),
          axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(plot.title = element_text(hjust = 0.5))
  #############################################################################################
  
  ###########################################################################
  ###########################################################################
  ###                                                                     ###
  ###                           SECTION 3:                                ###
  ###                           INDEX FiNDING                             ###
  ###                                                                     ###
  ###########################################################################
  ###########################################################################
  
  
  #Normalise difference water INDEX
  ndwi=ndwi_sentinel(b3_raster,b11_raster)
  windows()
  plot(ndwi)
  
  #Normalise difference vegetative INDEX
  # Values of 1 show vegetation, 0 bare and -1 water
  ndvi=ndvi_sentinel(b4_raster_10,b8_raster)
  #show the vegetation index index 
  windows()
  plot(ndvi)
  
  
  #since (10m res)1830*6=1098(60m res ) we need to interpolate 
  # we will do nereast neighbour interpolation
  
  #transform raster to matrix
  b4_matrix=raster2matrix(b4_raster)
  b4_matrix[1830,1830]
 # b4_matrix_10=raster2matrix(b4_raster_10)
  
  
  new_dim=dim(b4_raster_10)
  #FUNCTION TO CHANGE RESOLUTION-NEAREST NEIGHBOUR 
  # calculate row and column ratios
  transformed_b4=upscaling_fun(b4_matrix,new_dim)
  
  old_img[2][5]
  new_img[2][2]
  gc()
  new_dim=dim(b8_raster)
  b1_matrix=raster2matrix(b1_raster)
 
  b1_matrix_10m=upscaling_fun(b1_matrix,new_dim)

  array=c(1:dim(b11_raster)[1])
  
  b1_raster_10=raster(xmn=xmin(b8_raster),xmx=xmax(b8_raster),ymn=ymin(b8_raster),ymx=ymax(b8_raster),res=c(10,10),crs=crs(b1_raster))
  b8_raster
  b4_raster
  values(b1_raster_10)=c(b1_matrix_10m)

  #evi_sentinel(b8_raster,b4_raster,b1_raster_10)
  ############TRANSLATING 20m to 10m resolution  
  new_dim=dim(b8_raster)
  
  b11_matrix_2=raster2matrix(b11_raster)
  
  b11_matrix_10_2=upscaling_fun(b11_matrix_2,new_dim)
  
  dim(b11_matrix_10_2)
  
  b11_raster_10=raster(xmn=xmin(b8_raster),xmx=xmax(b8_raster),ymn=ymin(b8_raster),ymx=ymax(b8_raster),res=c(10,10),crs=crs(b11_raster))
  
  values(b11_raster_10)=c(b11_matrix_10_2)
  #delete the raster since we dont need it 
  rm(b11_raster)
  msi=msi_sentinel(b8_raster,b11_raster_10)
  
  values(b11_raster_10)
  values(b11_raster)
  windows()
  plot(msi)
  
  #NDBI

  new_dim=dim(b8_raster)
  b11_matrix=raster2matrix(b11_raster)
  rm(b11_raster)  
  gc()
  b11_matrix_10=upscaling_fun(b11_matrix,new_dim)
  rm(b11_matrix)
  gc()
  
  
  b11_raster_10=raster(xmn=xmin(b8_raster),xmx=xmax(b8_raster),ymn=ymin(b8_raster),ymx=ymax(b8_raster),res=c(10,10),crs=crs(b8_raster))
  
  values(b11_raster_10)=c(b11_matrix_10)
  
 ndbi=ndbi_sentinel(b8_raster,b11_raster_10)
  rm(ndbi)
  
 #GNDVI
  gndvi=gndvi(b8_raster,b3_raster)
  rm(b3_raster)
  rm(gndvi)
  gc()
  #NDMI
  b9_matrix=raster2matrix(b9_raster)
  new_dim=dim(b8_raster)
  b9_matrix_10=upscaling_fun(b9_matrix,new_dim)
  b8_matrix=raster2matrix(b8_raster)
  #get the ndmi by another way 
  ndmi=ndmi_sentinel2(b9_matrix_10,b8_matrix)
  ndmi_raster=raster(xmn=xmin(b8_raster),xmx=xmax(b8_raster),ymn=ymin(b8_raster),ymx=ymax(b8_raster),res=c(10,10),crs=crs(b9_raster))
  values(ndmi_raster)=c(ndmi)
  output_name= "ndmi.tif"
  
  #export the image to the working directory
  raster::writeRaster(ndmi_raster, filename = output_name)
  
  rm(b9_pic)
  rm(b9_matrix_10)
  rm(b8_matrix)
  rm(b9_raster)
  rm(b9_matrix)
  gc()
#EVI
  b2_matrix=raster2matrix(b2_raster)
  b8_matrix=raster2matrix(b8_raster)
  b4_matrix=raster2matrix(b4_raster)
  rm(b2_raster)
  rm(b8_raster)
  rm(b4_raster)
  gc()
  #evi=evi_sentinel(b2_matrix,b4_matrix,b8_matrix)
  evi=(b8_matrix-b4_matrix)/(b8_matrix+6*b4_matrix-7.5*b2_matrix+1)
  evi_raster=raster(xmn=xmin(b8_raster),xmx=xmax(b8_raster),ymn=ymin(b8_raster),ymx=ymax(b8_raster),res=c(10,10),crs=crs(b8_raster))
  values(evi_raster)=c(evi)
  output_name= "evi.tif"
  rm(b2_matrix)
  rm(b8_matrix)
  rm(b4_matrix)
  rm(b2_pic)
  rm(b8_pic)
  rm(b4_pic)
  rm(new_dim)
  rm(ndmi)
  rm(ndmi_raster)
  gc()
  #export the image to the working directory
  raster::writeRaster(evi_raster, filename = output_name, format="GTiff")
  
  #################################5 CUSTOM INDICES#################################
  #MSAVI2
  msavi_2=msavi_sentinel(b4_raster_10,b8_raster)
  windows()
  plot(msavi_2)
  values(b4_raster)
  dim(b4_raster)
  #Normalized Difference Glacier Index (NDGI)
  b4_matrix=raster2matrix(b4_raster)
  ndgi=(b3_matrix-b4_matrix)/(b3_matrix+b4_matrix)
  ndgi_raster=raster(xmn=xmin(b4_raster),xmx=xmax(b4_raster),ymn=ymin(b4_raster),ymx=ymax(b4_raster),res=c(10,10),crs=crs(b4_raster))
  values(ndgi_raster)=c(ndgi)
  output_name= "ndgi.tif"
  rm(b4_raster)
  rm(b3_matrix)
  raster::writeRaster(ndgi_raster, filename = output_name, format="GTiff")
  rm(ndgi)
  rm(ndgi_raster)
  ############Normalized Difference Snow Index (NDSI)############
  rm(b8_raster)
  b11_matrix_10=raster2matrix(b11_raster_10)
  rm(b11_raster_10)
  gc()
  b3_matrix=raster2matrix(b3_raster)
  ndsi=(b3_matrix-b11_matrix_10)/(b3_matrix+b11_matrix_10)
  ndsi_raster=raster(xmn=xmin(b3_raster),xmx=xmax(b3_raster),ymn=ymin(b3_raster),ymx=ymax(b3_raster),res=c(10,10),crs=crs(b3_raster))
  values(ndsi_raster)=c(ndsi)
  output_name= "ndsi.tif"
  rm(b11_matrix_10)
  rm(b3_raster)
  rm(ndsi)
  raster::writeRaster(ndsi_raster, filename = output_name, format="GTiff")
  
  ###########Normalized Pigment Chlorophyll Ratio Index (NPCRI)-for crops ############
  #################https://www.geo.university/pages/spectral-indices-with-multispectral-satellite-data
  b4_matrix=raster2matrix(b4_raster)
  rm(b4_raster)
  b2_matrix=raster2matrix(b2_raster)
  rm(b4_matrix)
  npcri=(b4_matrix-b2_matrix)/(b4_matrix+b2_matrix)
  npcri_raster=raster(xmn=xmin(b2_raster),xmx=xmax(b2_raster),ymn=ymin(b2_raster),ymx=ymax(b2_raster),res=c(10,10),crs=crs(b2_raster))
  values(npcri_raster)=c(npcri)
  output_name= "npcri.tif"
  raster::writeRaster(npcri_raster, filename = output_name, format="GTiff")
  rm(npcri)
  rm(npcri_raster)
  ###########Structure Insensitive Pigment Index (SIPI) f SIPI = (NIR – Blue) / (NIR – Red)############
  b4_matrix=raster2matrix(b4_raster)
  rm(b4_raster)
  b2_matrix=raster2matrix(b2_raster)
  rm(b2_raster)
  b8_matrix=raster2matrix(b8_raster)
  rm(b8_raster)
  sipi=(b8_matrix-b2_matrix)/(b8_matrix-b4_matrix)
  sipi_raster=raster(xmn=xmin(b2_raster),xmx=xmax(b2_raster),ymn=ymin(b2_raster),ymx=ymax(b2_raster),res=c(10,10),crs=crs(b2_raster))
  values(sipi_raster)=c(sipi)
  rm(sipi)
  rm(sipi_raster)
  gc()
  output_name= "sipi.tif"
  raster::writeRaster(sipi_raster, filename = output_name, format="GTiff")
  ###########OSAVI is the soil-adjusted vegetation index optimized for agricultural monitoring.############
  b4_matrix=raster2matrix(b4_raster)
  
  b8_matrix=raster2matrix(b8_raster)
  rm(b8_matrix)
  osavi=(b8_matrix-b4_matrix)/(b8_matrix+b4_matrix+0.16)
  osavi_raster=raster(xmn=xmin(b4_raster),xmx=xmax(b4_raster),ymn=ymin(b4_raster),ymx=ymax(b4_raster),res=c(10,10),crs=crs(b4_raster))
  rm(b4_raster)
  values(osavi_raster)=c(osavi)
  rm(b4_matrix)
  gc()
  rm(osavi)
  output_name= "osavi.tif"
  raster::writeRaster(osavi_raster, filename = output_name, format="GTiff")
  gc()
  rm(osavi_raster)