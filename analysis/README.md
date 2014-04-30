Analysis Scripts
==============

This analysis involves large datasets (~250TB) and both cloud-based (Google Earth Engine) and local processing and thus cannot simply be 're-run'.  However, the scripts are organized sequentially and sould be useful to replicate the analysis.  

* [setup.R](setup.R) loads necessary libraries, identifies useful directories, and loads some common data.  It is typically run at the beginning of each script.
* [1_ee_MCD09CF.js](1_ee_MCD09CF.js) contains the Google Earth Engine (GEE) javascript used to process and export the 'raw' MOD09GA and MYD09GA products into monthly climatologies.  The resulting files were then exported to a personal Google Drive account and downloaded to a local server.
* [2_compile.R](2_compile.R) takes the files from GEE and mosaics/reprojects them to a standard WGS84 grid
* [3_bias.R](3_bias.R) performs the VSNR artifact removal procedure to account for MODIS' orbital effects
* [4_combine.R](4_combine.R) calculates the combined (terra+aqua) mean cloud frequency
* [5_validate.R](5_validate.R) performs the station validation
* [6_Biome.R](6_Biome.R) calculates the biome summaries
* [6_Figures.R](6_Figures.R) generates the figures used in the manuscript