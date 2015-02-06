writeRasterCT=function(raster,cols,outfile){
  system(paste0("gdalbuildvrt -overwrite ",outfile," ",raster))
  vrt=scan(outfile,what="char")                        
  hd=c("<ColorInterp>Palette</ColorInterp>","<ColorTable>")
  ft="</ColorTable>"
  colR=cols
  cols=data.frame(t(col2rgb(colR)))
  ct=paste("<Entry c1=\"",cols$red,"\" c2=\"",cols$green,"\" c3=\"",cols$blue,"\" c4=\"255\"/>")
  cti=grep("ColorInterp",vrt)  # get index of current color table
  vrt2=c(vrt[1:(cti-1)],hd,ct,ft,vrt[(cti+1):length(vrt)])
  ## update missing data flag following http://lists.osgeo.org/pipermail/gdal-dev/2010-February/023541.html
#  csi=grep("<ComplexSource>",vrt2)  # get index of current color table
#  vrt2=c(vrt2[1:csi],"<NODATA>0</NODATA>",vrt2[(csi+1):length(vrt2)])
  write.table(vrt2,file=outfile,col.names=F,row.names=F,quote=F)
}
