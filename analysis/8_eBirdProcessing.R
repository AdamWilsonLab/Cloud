
### using standard reference data
## which colnames are in our list
tsp=grep(paste(c(idcols,spcols),collapse="|"),hdr,value=T)
d=read.csv.sql(f, sql = paste0("select * FROM file WHERE ",
                               "COMMON.NAME name REGEXP '",paste(gsub("_"," ",spcols),collapse="|"),"' AND ",
                               "LATITUDE  BETWEEN ",bbox(reg)["y","min"]," AND ",bbox(reg)["y","max"]," AND ",
                               "LONGITUDE BETWEEN ",bbox(reg)["x","min"]," AND ",bbox(reg)["x","max"]," AND ",
                               "PRIMARY_CHECKLIST_FLAG=1"," LIMIT 5 ;"),sep="\t")



## get colnames in table
hdr=colnames(read.csv(f,nrows=1))
## which colnames are in our list
tsp=grep(paste(c(idcols,spcols),collapse="|"),hdr,value=T)
d=system.time(read.csv.sql(f, sql = paste0("select \"",paste(tsp,collapse="\",\""),"\" FROM file WHERE ",
                                           "LATITUDE  BETWEEN ",bbox(reg)["y","min"]," AND ",bbox(reg)["y","max"]," AND ",
                                           "LONGITUDE BETWEEN ",bbox(reg)["x","min"]," AND ",bbox(reg)["x","max"]," AND ",
                                           "PRIMARY_CHECKLIST_FLAG=1",";")))
colnames(d)=gsub("\"","",colnames(d))
## loop through species and clean up "X" flags
for(i in tsp){
  x=d[,i]
  x[x=="X"]="1"
  x=as.numeric(x)
  x[x>0]=1
  d[,i]=x
  #  print(i)
}
## calculate number of trials and presences
d$trials=rowSums(d[,tsp])>0
d=d%.%filter(trials>0)
d$presence=d[,paste(sp,collapse="_")]
## subset to columns of interest
d=d[,c("LATITUDE","LONGITUDE","YEAR","MONTH","presence","trials")]
print(f)
return(d)
}




### try ff
colClasses=c("factor","integer",rep("factor",18),rep("numeric",2),"Date",rep("factor",8),rep("numeric",4),"logical","numeric","character","character","character")
d2=read.csv.ffdf(file=f,sep="\t",nrows=200000,colClasses="factor")


##### exploring fread



library("RSQLite")
field.types <- list(
  GLOBAL.UNIQUE.IDENTIFIER="TEXT",
  TAXONOMIC.ORDER="TEXT",
  CATEGORY="TEXT",
  COMMON.NAME="TEXT",
  SCIENTIFIC.NAME="TEXT",
  SUBSPECIES.COMMON.NAME="TEXT",
  SUBSPECIES.SCIENTIFIC.NAME="TEXT",
  OBSERVATION.COUNT="TEXT", 
  BREEDING.BIRD.ATLAS.CODE="TEXT",
  AGE.SEX="TEXT",
  COUNTRY="TEXT",
  COUNTRY.CODE="TEXT",
  STATE="TEXT",
  STATE.CODE="TEXT",
  COUNTY="TEXT",
  COUNTY.CODE="TEXT",
  IBA.CODE="TEXT",
  BCR.CODE="TEXT",
  LOCALITY="TEXT",
  LOCALITY.ID="TEXT",
  LOCALITY.TYPE="TEXT",
  LATITUDE="REAL",
  LONGITUDE="REAL",
  OBSERVATION.DATE="TEXT",
  TIME.OBSERVATIONS.STARTED="TEXT",
  TRIP.COMMENTS="TEXT",
  SPECIES.COMMENTS="TEXT",
  OBSERVER.ID="TEXT",
  FIRST.NAME="TEXT",
  LAST.NAME="TEXT",
  SAMPLING.EVENT.IDENTIFIER="TEXT",
  PROTOCOL.TYPE="TEXT",
  PROJECT.CODE="TEXT",
  DURATION.MINUTES="TEXT",
  EFFORT.DISTANCE.KM="TEXT",
  EFFORT.AREA.HA="TEXT",
  NUMBER.OBSERVERS="TEXT",
  ALL.SPECIES.REPORTED="TEXT",
  GROUP.IDENTIFIER="TEXT",
  APPROVED="TEXT",
  REVIEWED="TEXT",
  REASON="TEXT",
  X="TEXT")    

db <- dbConnect(SQLite(), dbname="ebird.sqlite") ## will make, if not present
dbGetQuery(db,"PRAGMA encoding")

10501620

dbWriteTable(conn=db, name="ebird", value=f,
             row.names=FALSE, header=TRUE, field.types=field.types,sep="\t")
dbGetQuery(db, "CREATE INDEX IF NOT EXISTS idx_species ON ebird (SCIENTIFIC.NAME)")
dbGetQuery(db, "CREATE INDEX IF NOT EXISTS lat ON ebird (LATITUDE)")
dbGetQuery(db, "CREATE INDEX IF NOT EXISTS lon ON ebird (LONGITUDE)")

dbDisconnect(db)

#dbWriteTable(conn=db, name="ebird", value=f, row.names=FALSE, header=TRUE, field.types=field.types,sep="\t")
#Error in try({ : 
#                 RS-DBI driver: (RS_sqlite_import: /mnt/data2/projects/mol/points_Aug_2014/ebd_relMay-2014.txt line 10501620 expected 43 columns of data but found 27)
#               [1] FALSE
#               > system(paste("sed -n '10501620p' ",f," "))
#               URN:CornellLabOfOrnithology:EBIRD:OBS227847206  17672	species	Loggerhead Shrike	Lanius ludovicianus			X			United States	US	Arizona	US-AZ	Cochise	US-AZ-003	US-AZ_902	34	Carr Canyon	L128935	H	31.4421997	-110.2855988	1997-05-03			Patagonia Nature Conservacy¿¿¿/? ¿	obsr88953	Brenda	Tekin	S16533503	eBird - Casual Observation	EBIRD					1		1	0		#



field.types <- list(
  GLOBAL.UNIQUE.IDENTIFIER="character",
  TAXONOMIC.ORDER="character",
  CATEGORY="character",
  COMMON.NAME="character",
  SCIENTIFIC.NAME="character",
  SUBSPECIES.COMMON.NAME="character",
  SUBSPECIES.SCIENTIFIC.NAME="character",
  OBSERVATION.COUNT="character", 
  BREEDING.BIRD.ATLAS.CODE="character",
  AGE.SEX="character",
  COUNTRY="character",
  COUNTRY.CODE="character",
  STATE="character",
  STATE.CODE="character",
  COUNTY="character",
  COUNTY.CODE="character",
  IBA.CODE="character",
  BCR.CODE="character",
  LOCALITY="character",
  LOCALITY.ID="character",
  LOCALITY.TYPE="character",
  LATITUDE="numeric",
  LONGITUDE="numeric",
  OBSERVATION.DATE="character",
  TIME.OBSERVATIONS.STARTED="character",
  TRIP.COMMENTS="character",
  SPECIES.COMMENTS="character",
  OBSERVER.ID="character",
  FIRST.NAME="character",
  LAST.NAME="character",
  SAMPLING.EVENT.IDENTIFIER="character",
  PROTOCOL.TYPE="character",
  PROJECT.CODE="character",
  DURATION.MINUTES="character",
  EFFORT.DISTANCE.KM="character",
  EFFORT.AREA.HA="character",
  NUMBER.OBSERVERS="character",
  ALL.SPECIES.REPORTED="character",
  GROUP.IDENTIFIER="character",
  APPROVED="character",
  REVIEWED="character",
  REASON="character",
  X="character")    


sfn=function(f,species,lat,lon){
  d=data.table(scan(f,what=field.types,nlines=1000,sep="\t",skip=1),stringsAsFactors=F)
  d%.%filter(LATITUDE>=bbox(reg)["y","min"],
             LATITUDE<=bbox(reg)["y","max"],
             LONGITUDE>=bbox(reg)["x","min"],
             LONGITUDE<=bbox(reg)["x","max"],
             ALL.SPECIES.REPORTED==1)}
  
  grep(paste(sub("_"," ",spcols),collapse="|"),as.character(d2$SCIENTIFIC.NAME))
  
  
  