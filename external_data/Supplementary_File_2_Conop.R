#Note: This script uses results extracted from the excel file produced by rGrapCor to assign the border effect on each Site


# CONOP formating # Including border effect

# Set working space

setwd("D:/CONOP_Daye/CONOP_2024_01_B")

# Load hits format:
# Hist  have the columns: subject,	depth,	hits,	site
hits <- read.csv("Neogene_zonation_Oct3_2023_filtered.csv")

# All Holocene samples in 1 group # LADs = 0 ; all FADs in border

holosamp <- c("L1_Holocene","L2_Holocene","L3_Holocene","L5_Holocene",
  "L6_Holocene","L7_Holocene","L8_Holocene","Amacayacu")

holos <- which(hits$site %in% holosamp)

hits$site[holos] <- "HOLOCENE"
hits$depth[holos] <- 0.0


# SITES (.sct file)
usites <- unique(hits$site)

# Abreviated sites (optional)
abb <- paste0(ifelse(grepl("AM", usites),"A",
  ifelse(grepl("Holocene", usites),"HO",substr(usites,1,3))), 
  gsub("[^0-9.]+", "", gsub("1AS","",usites)))


sctfile <- data.frame(A = 1:length(usites), B = abb, C = 1:length(usites), 
  D = usites, E = 1)

# change to single Quotation 
sctfile$B <- paste0(" '", sctfile$B, "'")
sctfile$D <- paste0(" '", sctfile$D, "'")

# Explore and make changes (leave only 3 characters)
fix(sctfile)

# chack no repeated numbers or abbreviations

any(duplicated(sctfile$A)) # must be FALSE
any(duplicated(sctfile$B)) # must be FALSE


# Save file

write.table(sctfile ,"section.sct",row.names = FALSE, col.names = FALSE, 
  sep = "\t", quote = FALSE)


# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

# 2. TAXA .EVT file
# labels your taxa and assigns them a number and nickname. 

usub <- sort(unique(hits$subject))


# Abbreviate names function (same in rgraphcor)

AbbNames <- function(vector) {
  abbreviations <- character(length = length(vector))
  splitnames <- strsplit(vector, " ")
  for (i in 1:length(vector)){
    # One name
    if(length(splitnames[[i]]) == 1){
      vector[i] <- paste(substr(splitnames[[i]][1],1,2),
                         substr(splitnames[[i]][2],1,5), sep = "")
    }
    # Two names
    else if(length(splitnames[[i]]) == 2){
      vector[i] <- if(splitnames[[i]][2] == "^sp") {
        paste(substr(splitnames[[i]][1],1,3), splitnames[[i]][2], sep = "")
      }
      else {
        paste(substr(splitnames[[i]][1],1,1),
              substr(splitnames[[i]][2],1,5), sep = "")
      }
    }
    # Three names
    else if(length(splitnames[[i]])==3){
      vector[i] <- paste(substr(splitnames[[i]][1],1,1),
                         substr(splitnames[[i]][2],1,3),
                         substr(splitnames[[i]][3],1,5), sep = "")
      # Assuming that the unwanted word is always in the middle 
    }
  }
  return(vector)
} # abb names


abbn <- AbbNames(usub)

# cut the names into 20 characters only

snames <- strsplit(usub, split = " ")

usnames <- unlist(lapply(snames, function(x)  
  paste(substr(x[1],1,6),x[2],x[3])))

unames20 <- substr(gsub(" NA", "", usnames),1,20)
any(duplicated(unames20)) # must be FALSE

evtfile <- data.frame(A = 1:length(usub), B = abbn, C = unames20 ) 

# check names and make necessary changes 
fix(evtfile)

any(duplicated(evtfile$B)) # yes but the numbers are different, so no problem

# change to Single Quotation 
evtfile$B <- paste0(" '", evtfile$B, "'")
evtfile$C <- paste0(" '", evtfile$C, "'")

# write format evt
write.table(evtfile ,"events.evt", row.names = FALSE, col.names = FALSE, 
 sep = "\t", quote = FALSE)

# Save a key for complete subject names for future reference (optional)
sub.key <- data.frame(Subject = usub, Abbreviation = abbn, 
  Abb20char = unames20, Number = 1:length(usub) ) 
write.table(sub.key,"SubjectKey.csv", row.names = FALSE, col.names = TRUE,
  sep = ",")

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

## .DAT file (fads and lads)
# fad & lads

# All negative numbers because all sites are wellcores
hits$depthN <-  hits$depth * -1

flh <- aggregate(depthN ~ subject + site, data = hits, FUN = range)

# fad/ lad of each taxa on each site
fadladsite <- data.frame(subject = flh$subject, site = flh$site, 
  FAD = flh$depth[,1], LAD = flh$depth[,2])


library(reshape2)

# reshape dataframe to get Fads and lads by depth separated
mfl <- melt(fadladsite, id.vars = c("subject","site"),
  variable_name = "event",  value.name = "depth")

# level within sites
levs <- aggregate(rep(1,nrow(mfl)) ~ depth + site , data = mfl, FUN = max)

# sort by site, then value
levs1 <- levs[order(levs$depth, decreasing = F),]
levs <- levs1[order(levs1$site),]


for ( i in 2:nrow(levs)){
  levs[i,3] <- ifelse(levs$site[i] == levs$site[i - 1],
    levs[i,3] + levs[i -1,3], 1)
}

# set codes Site + Depth
levs$code <- paste(levs$site,levs$depth)
mfl$code <- paste(mfl$site,mfl$depth)

# fix quotation 
mfl$site <- paste0(" '", mfl$site,"'")

df1 <- data.frame(
  Nsubject = sub.key$Number[match(mfl$subject, sub.key$Subject)],
  Event = ifelse(mfl$variable == "FAD", 1,2),
  Nsite = sctfile$C[match(mfl$site, sctfile$D)],
  D = mfl$depth,
  Nlevel = levs[match(mfl$code, levs$code) , 3])

df1$SITE <- sctfile$D[match(df1$Nsite,sctfile$C )]
df1$SITE <- gsub(" '","", df1$SITE)
df1$SITE <- gsub("'","", df1$SITE)

# match sites with Borders

# Load borders assigned from rGraphCor output:
library(readxl)
mdata <- read_excel("CS_2023-12-20.xlsx", sheet = "Metadata")

df1$bfad <- mdata$BorderFAD[match(df1$SITE, mdata$SiteName)]
df1$bfad <- ifelse(is.na(df1$bfad), 0, df1$bfad)

df1$blad <- mdata$BorderLAD[match(df1$SITE, mdata$SiteName)]
df1$blad <- ifelse(is.na(df1$blad), 0, df1$blad)

# number of levels on each site
snlevs <- aggregate(levs[,3] ~ levs[,2], FUN = max)
df1$nlevs <- snlevs[match(df1$SITE, snlevs[,1]),2]

# number of level above which is considered a border
df1$blad_lev <- df1$nlevs - df1$blad

# mark borders

df1$ISBORDER <- "NO"

# FADS in Holocene samples
df1$ISBORDER <- ifelse(df1$SITE =="HOLOCENE" & 
  df1$Event == 1, "YES", df1$ISBORDER)

# border fads
df1$ISBORDER <- ifelse(df1$Nleve <= df1$bfad & df1$Event == 1,
  "YES", df1$ISBORDER )

df1$ISBORDER <- ifelse(df1$Nleve > df1$blad_lev & df1$Event == 2,
  "YES",df1$ISBORDER )

# keep only non borders events
df2 <- df1[which(df1$ISBORDER == "NO"), ]

levsite <- aggregate(df2$Nlevel ~ df2$Nsite, FUN = c)

nwlevsite <- lapply(levsite[,2], function(x) 
  data.frame(oldlevel = sort(unique(x)), newlevel = order(sort(unique(x)))))

# Set new levels (consecutive after  filtering borders
df2$NwLev <- NA

for (i in 1: length(nwlevsite)){
  ies <- which(df2$Nsite == i)
  df2$NwLev[ies] <- nwlevsite[[i]][match(df2$Nlevel[ies], nwlevsite[[i]][,1]),2]
}


# Eliminate subjects with incomplete  range ( only fad or only Lad)

incrange <- aggregate(df2$Event ~ df2$Nsubject, FUN = unique)
incs <- which(lapply(incrange[,2],length) < 2)

df3 <- df2[df2$Nsubject %in% incs == FALSE, ]

# Re enumerate Taxa 
ntaxa <- data.frame(oldN = sort(unique(df3$Nsubject)))
ntaxa$newN <- order(ntaxa$oldN)

df3$Newsubject <- ntaxa$newN[match(df3$Nsubject, ntaxa$oldN)]

# check levels again

levsite2 <- aggregate(df3$NwLev ~ df3$Nsite, FUN = c)

nwlevsite2 <- lapply(levsite2[,2], function(x) 
  data.frame(oldlevel = sort(unique(x)), newlevel = order(sort(unique(x)))))

# Set new levels (consecutives after  filtering borders
df3$NwLev2 <- NA

for (i in 1: length(nwlevsite)){
  ies <- which(df3$Nsite == i)
  df3$NwLev2[ies] <- nwlevsite2[[i]][match(df3$Nlevel[ies], nwlevsite2[[i]][,1]),2]
}


datfile <- data.frame(
  A = df3$Newsubject,
  B = df3$Event,
  C = df3$Nsite,
  D = df3$D,
  E = df3$NwLev2,
  F = df3$Event,
  G = 1,
  H = 1)

#  avoid zeros by changing for 0.001 (negative because is a wellcore)
datfile$D <- ifelse (datfile$D == 0, -0.0001, datfile$D)

# SORT BY A, B, and C, in that order. 
# group together all of the first occurrences of the first taxon in order
# by section. Next will be the last appearances for that same taxon,
# in the same order by section.

datfile <- datfile[order(datfile$C),]
datfile <- datfile[order(datfile$B),]
datfile <- datfile[order(datfile$A),]

head(datfile, 30)

# write file
write.table(datfile,"loadfile.dat", row.names = FALSE, col.names = FALSE, sep = "\t")


# ........
# for reference in CFG file 
# n sites
sort(unique(datfile$C))

# n taxa
sort(unique(datfile$A))
