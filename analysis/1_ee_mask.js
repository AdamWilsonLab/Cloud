// MCD09_MaskAlbedo

// specify what happens
var verbose=true       // print info about collections along the way (slows things down)
var drawmap=true       // add image to map
var exportDrive=true  // add exports to task window
var DownloadURL=false  // add exports to task window

//  Specify destination and run name
var driveFolder="EarthEngineOutput";
var run="GlobalSRTM_slopeLTE2minalbedo30"

// Get current date as string for file metadata
var currentdate = new Date();
var date= currentdate.getFullYear()+''+("0" + (currentdate.getMonth() + 1)).slice(-2)+''+currentdate.getDate();
print(date)


//  MODIS water mask
var water=ee.Image('MODIS/MOD44W/MOD44W_005_2000_02_24').select("water_mask")
var watermask=water.mask(water.eq(1))

// Import Cloud frequency data for comparison (not used in generation of mask)
var cf_01=  ee.Image('GME/images/04040405428907908306-05137678721221478767')
var cf_02=  ee.Image('GME/images/04040405428907908306-04053480542821848562')
var cf_03=  ee.Image('GME/images/04040405428907908306-15722273339382880469')
var cf_04=  ee.Image('GME/images/04040405428907908306-01763722745770023372')
var cf_05=  ee.Image('GME/images/04040405428907908306-17401395946472583516')
var cf_06=  ee.Image('GME/images/04040405428907908306-14795255083217681318')
var cf_07=  ee.Image('GME/images/04040405428907908306-15789600467113128339')
var cf_08=  ee.Image('GME/images/04040405428907908306-14214644652936744608')
var cf_09=  ee.Image('GME/images/04040405428907908306-06867970932243171500')
var cf_10=  ee.Image('GME/images/04040405428907908306-13259991075269553771')
var cf_11=  ee.Image('GME/images/04040405428907908306-07710157189606598984')
var cf_12=  ee.Image('GME/images/04040405428907908306-01845458384947302296')
// CF derived layers
var cf=ee.ImageCollection([cf_01,cf_02,cf_03,cf_04,cf_05,cf_06,cf_07,cf_08,cf_09,cf_10,cf_11,cf_12])
var intra= ee.Image('GME/images/04040405428907908306-16043517108931311115').divide(100);
var inter=ee.Image('GME/images/04040405428907908306-00562669341254644056').divide(100);
// annual average
var annual=cf.mean()

print(cf_01.getInfo())

///////////////////////////////////////////////////////////////////////////
// MCD43B3 Albedo
var alb=ee.ImageCollection('MODIS/MCD43B3').select('Albedo_WSA_vis')
//print(alb.getInfo())

var alb_01=alb.filter(ee.Filter.calendarRange(1,1,"month")).mean()
var alb_02=alb.filter(ee.Filter.calendarRange(2,2,"month")).mean()
var alb_03=alb.filter(ee.Filter.calendarRange(3,3,"month")).mean()
var alb_04=alb.filter(ee.Filter.calendarRange(4,4,"month")).mean()
var alb_05=alb.filter(ee.Filter.calendarRange(5,5,"month")).mean()
var alb_06=alb.filter(ee.Filter.calendarRange(6,6,"month")).mean()
var alb_07=alb.filter(ee.Filter.calendarRange(7,7,"month")).mean()
var alb_08=alb.filter(ee.Filter.calendarRange(8,8,"month")).mean()
var alb_09=alb.filter(ee.Filter.calendarRange(9,9,"month")).mean()
var alb_10=alb.filter(ee.Filter.calendarRange(10,10,"month")).mean()
var alb_11=alb.filter(ee.Filter.calendarRange(11,11,"month")).mean()
var alb_12=alb.filter(ee.Filter.calendarRange(12,12,"month")).mean()
var alb_mean=alb.mean()
var alb_sd=alb.reduce(ee.call("Reducer.sampleStdDev"))
var alb_min=alb.min();
var alb_max=alb.max();
var alb_range=alb_max.subtract(alb_min)

//seasonality
var seas=ee.Image("GME/images/04040405428907908306-00317830000645337584")

// Elevation
var dem=ee.Image('CGIAR/SRTM90_V4')//GME/images/04040405428907908306-08319720230328335274')

print(dem.getInfo())//
var slope=ee.Algorithms
                      .Terrain(dem)
                      .select("slope")
                      .focal_median(500,"square","meters")
                      .reproject('EPSG:4326',[0.00083333333, 0, -180, 0, -0.00083333333,60]);
                    
// load maximum slope slope_mx_GMTED2010_md
//var slope=ee.Image('GME/images/04040405428907908306-11428276147442013580').focal_max(2000,"circle","meters")

//centerMap(-90.79994, 44.21912, 2);
  var palette="08306b,0d57a1,2878b8,4997c9,72b2d7,a2cbe2,c7dcef,deebf7,f7fbff"
  
addToMap(cf_01,{min: 0, max: 100,palette:palette}, 'January',0);
addToMap(cf_02,{min: 0, max: 100,palette:palette}, 'February',0);
addToMap(cf_03,{min: 0, max: 100,palette:palette}, 'March',0);
addToMap(cf_04,{min: 0, max: 100,palette:palette}, 'April',0);
addToMap(cf_05,{min: 0, max: 100,palette:palette}, 'May',0);
addToMap(cf_06,{min: 0, max: 100,palette:palette}, 'June',0);
addToMap(cf_07,{min: 0, max: 100,palette:palette}, 'July',0);
addToMap(cf_08,{min: 0, max: 100,palette:palette}, 'August',0);
addToMap(cf_09,{min: 0, max: 100,palette:palette}, 'September',0);
addToMap(cf_10,{min: 0, max: 100,palette:palette}, 'October',0);
addToMap(cf_11,{min: 0, max: 100,palette:palette}, 'November',0);
addToMap(cf_12,{min: 0, max: 100,palette:palette}, 'December',0);


// add albedo to map
var apalette="000000,ffffff"
var arange=[0,700]
addToMap(alb_01,{min: arange[0], max: arange[1],palette:apalette}, 'Albedo_01',0);
addToMap(alb_02,{min: arange[0], max: arange[1],palette:apalette}, 'Albedo_02',0);
addToMap(alb_03,{min: arange[0], max: arange[1],palette:apalette}, 'Albedo_03',0);
addToMap(alb_04,{min: arange[0], max: arange[1],palette:apalette}, 'Albedo_04',0);
addToMap(alb_05,{min: arange[0], max: arange[1],palette:apalette}, 'Albedo_05',0);
addToMap(alb_06,{min: arange[0], max: arange[1],palette:apalette}, 'Albedo_06',0);
addToMap(alb_07,{min: arange[0], max: arange[1],palette:apalette}, 'Albedo_07',0);
addToMap(alb_08,{min: arange[0], max: arange[1],palette:apalette}, 'Albedo_08',0);
addToMap(alb_09,{min: arange[0], max: arange[1],palette:apalette}, 'Albedo_09',0);
addToMap(alb_10,{min: arange[0], max: arange[1],palette:apalette}, 'Albedo_10',0);
addToMap(alb_11,{min: arange[0], max: arange[1],palette:apalette}, 'Albedo_11',0);
addToMap(alb_12,{min: arange[0], max: arange[1],palette:apalette}, 'Albedo_12',0);
addToMap(alb_mean,{min: arange[0], max: arange[1],palette:apalette}, 'Albedo_mean',0);
addToMap(alb_sd,{min: 0, max: 100,palette:apalette}, 'Albedo_sd',0);
addToMap(alb_min,{min: arange[0], max: arange[1],palette:apalette}, 'Albedo_min',0);
addToMap(alb_max,{min: arange[0], max: arange[1],palette:apalette}, 'Albedo_max',0);
addToMap(alb_range,{min: arange[0], max: arange[1],palette:apalette}, 'Albedo_range',0);


// build mask to eliminate hotspots
//region to process water mask
var geodesic = ee.Geometry.Rectangle(-180, -60, 180, 60);
var maskbox = ee.Geometry(geodesic, null, false);

//start with image of cloud frequency to define the grid 
var mask_base=cf_06.gte(0).not()
// Identify areas with high albedo and high variability of albedo (SD) 
var mask_alb=alb_min.gte(180).and(alb_sd.gte(50)).and(slope.lte(2))//.clip(maskbox) // 150:50:0  .and(alb_sd.lte(200))
// Identify water pixels with high albedo and high variability and buffer by 2km 
var mask_water=alb_min.gte(25).and(alb_sd.gte(10)).and(water.eq(1))//.clip(maskbox)
// Combine the masks
var mask=mask_base.where(mask_alb,1)
                  .where(mask_water,2)
                  .where(mask_water.and(mask_alb),3)
                  .focal_max(3000,"circle","meters")
//                  .clip(maskbox)
                  .int8()

  var palette2="0000ff,00ff00,ff0000"
addToMap(annual,{min: 0, max: 100,palette:palette}, 'Mean Annual',0);
addToMap(annual.mask(mask.eq(0)),{min: 0, max: 100,palette:palette}, 'Mean Annual (Masked)',0);
addToMap(intra,{min: 0, max: 25,palette:palette2}, 'IntraannualVariability',0);
addToMap(inter,{min: 0, max: 40,palette:palette2}, 'Interannual Variability',0);
addToMap(watermask,{min: 0, max: 1,palette:"0000ff"}, 'Water',0);

addToMap(seas,{}, 'Seasonality',1);

//add the mask  .mask(mask.eq(0))
addToMap(mask.mask(mask.gte(1)),{min:1,max:3,palette:"ff0000,00ff00,0000ff"},"mask",0)
addToMap(mask_alb.mask(mask_alb.eq(1)),{min:.25,max:1.5,palette:"000000,ff0000"},"maskalb",0)
addToMap(mask_water.mask(mask_water.eq(1)),{min:.25,max:1.5,palette:"000000,ff0000"},"maskwater",0)
addToMap(dem,{min:0,max:500},"DEM",0)
addToMap(slope,{min:0,max:5,palette:"ff0000,000000"},"Slope",0)

//print(mask.getInfo())
//print(annual.getInfo())

if(exportDrive){
  //define export filename
  var filename='mask_'+run+'_'+date;
  if(verbose){  print('Exporting to: '+filename)}

  exportImage(mask,filename,
    {'maxPixels':1000000000,
    'driveFolder':driveFolder,
    'crs': 'EPSG:4326',
//    'crs_transform': '[0.0083333333, 0, -180, 0, -0.0083333333,90]',
//    'dimensions':'[43200,21600]'
    'crs_transform': '[0.0083333333, 0, -180, 0, -0.0083333333,60]',
    'dimensions':'[43200,14400]'
//    'crs_transform': '[0.0083333333, 0, -180, 0, -0.0083333333,90]',
//    'dimensions':'[43200,21600]'
  });
}

if(DownloadURL){

var url=mask.getDownloadURL(
    {'maxPixels':1000000000,
    'crs': 'EPSG:4326',
//    'crs_transform': '[0.0083333333, 0, 54, 0, -0.0083333333,39]',
//    'dimensions':'[840,840]'
    'crs_transform': '[0.0083333333, 0, -180, 0, -0.0083333333,90]',
    'dimensions':'[43200,21600]'
  });
  print(url)
}
