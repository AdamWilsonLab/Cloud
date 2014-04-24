// MCD09CF
///////////////////////////////////////////////////////////
// Adam.wilson@yale.edu
// Generate a cloud frequency climatology from M*D09 data
// This script does the following
// 1) Extracts daily cloud flag from MOD09GA and MYD09GA and create single collection for each month
// 2) Calculate mean cloud frequency for each month
// 3) Calculate overall mean and sd of monthly mean cloud frequencies for each month
// 4) Export the monthly mean/sd cloud frequency separately (as separate tasks) for each month
////////////////////////////////////////////////////////////

//  Specify destination and run name
var driveFolder="ee_mcd09cf";
var run="g3"

// limit overall date range  (only dates in this range will be included)
var datestart=new Date("2000-01-01")  // default time zone is UTC
var datestop=new Date("2014-03-31")

// specify which months (within the date range above) to process
var monthstart=1
var monthstop=12

//  Sensors to process
var mcols = ['MOD09GA','MYD09GA'];

// specify what happens
var verbose=false       // print info about collections along the way (slows things down)
var drawmap=true       // add image to map
var test1=false         // report on all images requested via dates above, but only add image below to map
var exportDrive=!test1  // add exports to task window

// set testing sensor and month, these only apply if test1==t above
var testsensor='MOD09GA'
var testmonth=1

// define regions
var globe = '[[-180, -89.95], [-180, 89.95], [180, 89.95], [180, -89.95]]';  
var sahara = '[[-18, 0], [-18, 30], [15, 30], [15, 0]]';  
// choose region to export (if exportDrive==true)
var region = globe


/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
// Generate some filtering variables based on input above and date

// Get current date as string for file metadata
var currentdate = new Date();
var date= currentdate.getFullYear()+''+("0" + (currentdate.getMonth() + 1)).slice(-2)+''+currentdate.getDate();
print(date)



/////////////////////////////////////////////////////////////////////////////
/// Functions

//function to select terra or aqua and subset by date
function getMCD09(modcol,datestart,datestop){
  var getdates=ee.Filter.date(datestart,datestop);
  return(ee.ImageCollection(modcol).filter(getdates));
}

//function to subset a collection by year and month
function fyearmonth(collection,year,month){
  // define filters
  var getyear=ee.Filter.calendarRange(year,year,"year");
  var getmonth=ee.Filter.calendarRange(month,month,"month");
  // extract values
  var output=collection.filter(getmonth).filter(getyear);
  return(output);
}

///////////////////////////////////////////////////////////////////////////
// Loop through months and get cloud frequency for each month in year range

var monthmean=function(collection,month){
  //For this month, build array to hold means for every year
  var mod = new Array (nYears);
  // loop over years and and put the mean monthly CF in the array
    for (var i=0; i<nYears; i ++) {
      mod[i]= fyearmonth(collection,years[i],month). // filter by year-month
      map(getcf).                              // extract cloud frequency
      mean();                                  // take the mean
      }
  if(verbose){print('Processing '+nYears+' years of data for Month '+month)}
  // build an image collection for all years for this month using the array
  return(ee.ImageCollection(mod));
}
/**
 * Returns an image containing just the specified QA bits.
 *
 * Args:
 *   image - The QA Image to get bits from.
 *   start - The first bit position, 0-based.
 *   end   - The last bit position, inclusive.
 *   name  - A name for the output image.
 */
var getQABits = function(image, start, end, newName) {
    // Compute the bits we need to extract.
    var i, pattern = 0;
    for (i = start; i <= end; i++) {
       pattern += Math.pow(2, i);
    }
    return image.select([0], [newName])
                  .bitwise_and(pattern)
                  .right_shift(start);
};



// function to extract MOD09 internal cloud flag
var getcf=function(img){ 
  var img2=getQABits(img.select(['state_1km']),10, 10, 'cloud'); // extract cloud flag from bit 10
  return(img2.
        mask(img2.gte(0)).           // mask out NAs
        gte(1).                      // Convert to logical
        multiply(ee.Image(10000)))};   // Multiply by 10000 to facilate taking mean as an integer 0-10000
  
// function to return the total number of observations in the collection
var getcount=function(img) {
  var img2=img.select(['num_observations_1km']);
  return img2.
          mask(img2.gte(0)).           // mask out NAs
          select([0],['nObs']).        // rename to nObs
          multiply(ee.Image(10000))};    // multiply to 100 to facilate taking mean as an integer

// function to return the proportion of days with at least one observation in the collection
var getprop=function(img) {
  var img2=img.select(['num_observations_1km']);
  return img2.
          mask(img2.gte(0)).           // mask out NAs
          gte(1).                      // Convert to logical to get percentage of days with n>0
          multiply(ee.Image(10000)).      // multiply to 100 to facilate taking mean as an integer
          select([0],['pObs'])};        // rename to nObs


//function to return year+1
var yearplus=function (x, y) { return yearstart +y ; }

//////////////////////////////////////////////////////////////////////////////////////////////////////
//  Clip high-latitude nonsense
//  For an unknown reason there are pixels at high latitudes with nobs>0 in winter (dark) months.
//  These values were identified by viewing the monthly CF and nObs datasets and identifying the latitude
//  at which there was a stripe of fully missing data (nObs=0) followed by an erroneous region at higher latitudes
//  This returns a polygon for any month which can be used to clip the image
var monthbox=function(month) {
  // limits for each month (12 numbers correspond to jan-december limits)
  var ymin=[-90,-90,-90,-90,-70,-70,-70,-77,-77,-90,-90,-90];
  var ymax=[ 74, 84, 90, 90, 90, 90, 90, 90, 90, 90, 77, 69];
  // draw the polygon, use month-1 because month is 1:12 and array is indexed 0:11
  var ind=month-1;
  var box = ee.Geometry.Rectangle(-180,ymin[ind],180,ymax[ind]);
return (box)};


////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////
//  Start processing the data



// loop over collections and months
  for (var c=0;c<mcols.length;c++) {
    var mcol=mcols[c];
  for (var tmonth = monthstart; tmonth <= monthstop; tmonth ++ ) {

// if test1 is true, use only test settings 
if(test1){
  var mcol=testsensor
  var tmonth=testmonth
}

if(verbose){
  print('Starting processing for:'+ mcol+' month '+tmonth);
}

// identify start and stop years
var yearstart=datestart.getUTCFullYear();
var yearstop=datestop.getUTCFullYear();

// update startdates based on terra/aqua start dates
// otherwise missing month-years will cause problems when compiling the mean
// Terra (MOD) startdate: February 24, 2000
if(mcol=='MOD09GA'&tmonth==1&yearstart==2000) yearstart=2001
/// Aqua (MYD) startdate: July 4, 2002
if(mcol=='MYD09GA'&tmonth<7&yearstart<=2002) yearstart=2003
if(mcol=='MYD09GA'&tmonth>=7&yearstart<=2002) yearstart=2002

// get array of years to process
var years=Array.apply(0, Array(yearstop-yearstart+1)).map(yearplus);
var nYears=years.length;
if(verbose){print('Processing '+years)}



// build a combined M*D09GA collection limited by start and stop dates
var MCD09all=getMCD09(mcol,datestart,datestop);

//////////////////////////////////////////////////////////////////////////////////////////////////////
//  Process cloud frequency

var MCD09m=monthmean(MCD09all,tmonth);
if(verbose){
  print('MCD09m info:');
  print(MCD09m.getInfo());
}

// take overall mean and SD across all years for this month
var MCD09_mean=MCD09m.mean();
var MCD09_sd=MCD09m.reduce(ee.call("Reducer.sampleStdDev"));

/////////////////////////////////////////////////////////////////////////////////////////////////////
///  Process MODIS observation frequency/proportion

// get number of obs and proportion of days with at least one obs
var MCD09_pObs = MCD09all.filter(ee.Filter.calendarRange(tmonth,tmonth,"month")).map(getprop).mean();
var MCD09_nObs = MCD09all.filter(ee.Filter.calendarRange(tmonth,tmonth,"month")).map(getcount).mean();

//////////////////////////////////////////////////////////////////////////////////////////////////////
//  Process monthly cloud trends
//var trend= MCD09m.formaTrend(null,2)
//if(verbose) print(trend.getInfo())

/////////////////////////////////////////////////////////////////////////////////////////////////////
// Build a single 8-bit image with all bands
// mask by nobs (only keep pixels where nobs > 0.75)
var MCD09=MCD09_mean.
          addBands(MCD09_sd).
          addBands(MCD09_nObs).
          addBands(MCD09_pObs).
          mask(MCD09_nObs.gte(25)).
          clip(monthbox(tmonth)).
          int16();

if(exportDrive){
  //define export filename
  var filename=date+'_'+run+'_'+yearstart+yearstop+'_'+mcol+'_'+tmonth;
  if(verbose){  print('Exporting to: '+filename)}

  exportImage(MCD09,filename,
    {'maxPixels':1000000000,
    'driveFolder':driveFolder,
    'crs': 'SR-ORG:6974', //4326
    'scale': '926.625433055833',
    'region': region
  });
}

if(test1) break;  // if running testing, dont complete the loop

} // close loop through months
} // close loop through collections

// Draw the map?
if(drawmap) {
//  var palette="000000,00FF00,FF0000";
  var palette="000000,FFFFFF";

  addToMap(MCD09.select([0]),{min:0,max:10000,palette:palette}, "mean");
  addToMap(MCD09.select([1]),{min:0,max:5000,palette:palette}, "sd");
  addToMap(MCD09.select([2]),{min:0,max:15000,palette:palette}, "nObs");
  addToMap(MCD09.select([3]),{min:0,max:10000,palette:palette}, "pObs");
//  addToMap(trend.select('long-trend'), {min:-10, max:10, palette:palette}, 'trend');
}
