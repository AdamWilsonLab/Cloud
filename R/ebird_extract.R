
 # install_github("pingles/redshift-r")
  
require(redshift)
conn <- redshift.connect("jdbc:postgresql://mol-points.c98tkbi1cfwj.us-east-1.redshift.amazonaws.com:5439/mol?tcpKeepAlive=true",username="mol",password="Apu5apu5")

getebird=function(con, sptaxon, nulltaxon,region){
  print(paste("Extracting data, this can take a few minutes..."))
  dbGetQuery(conn, 
     paste(   
     "WITH ebird_subset as (SELECT all_species_reported,taxonomic_order,latitude,longitude,observation_date,sampling_event_identifier,group_identifier",
        "FROM ebird",
        "WHERE latitude BETWEEN ",paste(bbox(region)["y",],collapse=" AND "),
          "AND longitude BETWEEN ",paste(bbox(region)["x",],collapse=" AND "),
          "AND floor(taxonomic_order) IN (",paste(c(sptaxon,nulltaxon),collapse=","),")),",
      "presence as (SELECT DISTINCT latitude,longitude,observation_date,sampling_event_identifier,group_identifier,1 AS presence ",
        "FROM ebird_subset",
        "WHERE floor(taxonomic_order) IN (",paste(sptaxon,collapse=","),")),",
      "absence as (SELECT DISTINCT latitude,longitude,observation_date,sampling_event_identifier,group_identifier,0 AS presence",
        "FROM ebird_subset",
        "WHERE all_species_reported='t'",
           "AND sampling_event_identifier NOT IN (SELECT sampling_event_identifier FROM presence))",
     "SELECT latitude,longitude,observation_date,presence FROM presence",
     "UNION",
     "SELECT latitude,longitude,observation_date,presence FROM absence"))
}
  
  
#presabs=  getebird(con=conn,sptaxon=sptaxon,nulltaxon=nulltaxon,region=reg)
#plot(latitude~longitude,data=presabs[presabs$presence==0,],col="black",pch=".")
#points(latitude~longitude,data=presabs[presabs$presence==1,],col="red",pch=16,cex=.5)
  
