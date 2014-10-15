
 # install_github("pingles/redshift-r")
getebird=function(con, sptaxon, nulltaxon=NULL,region){
  print(paste("Extracting data, this can take a few minutes..."))
  if(!is.null(nulltaxon)){  #return only species in list as 'non-detection'
  dbGetQuery(conn, 
     paste(   
     "WITH ebird_subset as (SELECT all_species_reported,taxonomic_order,latitude,longitude,observation_date,sampling_event_identifier,group_identifier,effort_distance_km,effort_area_ha,duration_minutes",
        "FROM ebird",
        "WHERE latitude BETWEEN ",paste(bbox(region)["y",],collapse=" AND "),
          "AND longitude BETWEEN ",paste(bbox(region)["x",],collapse=" AND "),
          "AND floor(taxonomic_order) IN (",paste(c(sptaxon,nulltaxon),collapse=","),")),",
      "presence as (SELECT DISTINCT latitude,longitude,observation_date,sampling_event_identifier,group_identifier,effort_distance_km,effort_area_ha,duration_minutes,1 AS presence ",
        "FROM ebird_subset",
        "WHERE floor(taxonomic_order) IN (",paste(sptaxon,collapse=","),")),",
      "absence as (SELECT DISTINCT latitude,longitude,observation_date,sampling_event_identifier,group_identifier,effort_distance_km,effort_area_ha,duration_minutes,0 AS presence",
        "FROM ebird_subset",
        "WHERE all_species_reported='t'",
           "AND sampling_event_identifier NOT IN (SELECT sampling_event_identifier FROM presence))",
     "SELECT latitude,longitude,observation_date,presence,effort_distance_km,effort_area_ha,duration_minutes FROM presence",
     "UNION",
     "SELECT latitude,longitude,observation_date,presence,effort_distance_km,effort_area_ha,duration_minutes FROM absence"))
}
 if(is.null(nulltaxon)){  ## return all species as 'non-detections'
  dbGetQuery(conn, 
             paste(   
               "WITH ebird_subset as (SELECT all_species_reported,taxonomic_order,latitude,longitude,observation_date,sampling_event_identifier,group_identifier,effort_distance_km,effort_area_ha,duration_minutes",
               "FROM ebird",
               "WHERE latitude BETWEEN ",paste(bbox(region)["y",],collapse=" AND "),
               "AND longitude BETWEEN ",paste(bbox(region)["x",],collapse=" AND "),"),",
               "presence as (SELECT DISTINCT latitude,longitude,observation_date,sampling_event_identifier,group_identifier,effort_distance_km,effort_area_ha,duration_minutes,1 AS presence ",
               "FROM ebird_subset",
               "WHERE floor(taxonomic_order) IN (",paste(sptaxon,collapse=","),")),",
               "absence as (SELECT DISTINCT latitude,longitude,observation_date,sampling_event_identifier,group_identifier,effort_distance_km,effort_area_ha,duration_minutes,0 AS presence",
               "FROM ebird_subset",
               "WHERE all_species_reported='t'",
               "AND sampling_event_identifier NOT IN (SELECT sampling_event_identifier FROM presence))",
               "SELECT latitude,longitude,observation_date,presence,effort_distance_km,duration_minutes,effort_area_ha FROM presence",
               "UNION",
               "SELECT latitude,longitude,observation_date,presence,effort_distance_km,duration_minutes,effort_area_ha FROM absence"))
 }
}
  
  
#presabs=  getebird(con=conn,sptaxon=sptaxon,nulltaxon=nulltaxon,region=reg)
#plot(latitude~longitude,data=presabs[presabs$presence==0,],col="black",pch=".")
#points(latitude~longitude,data=presabs[presabs$presence==1,],col="red",pch=16,cex=.5)

#dbGetQuery(conn, 
#           paste(   
#             "SELECT all_species_reported,taxonomic_order,latitude,longitude,observation_date,sampling_event_identifier,group_identifier",
#             "SELECT * ",
#             "FROM ebird LIMIT 10"))
             
