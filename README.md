High-resolution Cloud Climatology for Global Land Areas
=====

## Adam M. Wilson and Walter Jetz

Code associated with the generation of  a Global Cloud Climatology derived from MODIS MOD09 cloud flags.



# Remaining Issues

* Missing data in inter/intra variability
* Albedo issues (Australia)
   *  beaches (high albedo sites with shallow water)
* Albedo  over high elevation areas
* Fine grain cloud variation 
* Validation: find region with high station density to do fine-grain validation example. 
* Timeseries
* Detailed images of interesting regions
* Comparison to coarser products in appendix
* Could you develop correlograms that characterizes the spatial autocorrelation for the main global precipitation or could related products and our cloud climatology (mean or intra-annual variation), ideally including the full (or 95%CI) global variation around the mean estimate for a given lag, and potentially separating the tropics from non-tropics. So this could include Worldclim (which also starts at 1km lag) and other cloud or prec products (CRU etc) that start at say 5 or 50km smallest lag. Sounds easy of course, but the computationally challenging (and also additionally novel) part would be to truly get all 1km pixels represented.  This could be a great way to not only characterize the world’s climate variation in an integrative way, but to bring out the particular scale and associated processes that empirical models to date have failed to address over large extents. Sound right?
* add quantiles to biome plot


## Four stories
 * Amazing detail of cloud cover
    * confirms detail in some regions with real differences
 * Getting it wrong in landscape  (interpolated precipitation layers)
    *   biodiversity research
 * global picture
    * inter vs. intra annual variability
    * hotspots of cloud extremes
    * 
    

New lens to look at the planet and monitor change
  cloud forest-environment niches could be shifting/disappearing
  
## Data distribution
Environmental layers on standardized grid

### Use-cases
1.  Regional (rectangular) subset
2.  Regional (polygons) uploaded or drawn subset with ID.  For each ID, person wants average, min, max, other summaries
3.  Points (single pixel or with buffer) with ID.
4.  






