4/14/2020 - Jacob LaFontaine (jlafonta@usgs.gov)

# Information
Processing of datasets to produce segment-based parameters for the USGS National Hydrologic Model with Geospatial Fabric version 1.1 (GF_v11)

- Parameters computed:
	- `seg_width`
	- `seg_depth`
	- `mann_n`

- Base datasets:
	- Existing PRMS parameters from the NHM paramdb version 1.1:
		- `seg_slope`

- Existing datasets:
	- GF_v11 stream segment shapefile
	- River network shapefiles for North America and Central America from the HydroSHEDS HydroRIVERS website (https://www.hydrosheds.org/page/hydrorivers)

# Requirements
The following external code and datasets are required to complete the workflow below.
- R (https://www.r-project.org)
- River network shapefiles for North America and Central America from the HydroSHEDS HydroRIVERS website (https://www.hydrosheds.org/page/hydrorivers)

# Workflow
Processing:
1. HydroRIVERS shapefiles of North and Central America were merged into one shapefile in ArcGIS.
2. The merged shapefile of HydroRIVERS was converted to a point coverage in ArcGIS.
3. A spatial join was done using the nearest neighbor method to transfer the HydroRIVERS attributes to the GF_v11 stream segments.
4. The new attribute table for GF_v11 stream segments was output to a table.
5. Code was written in R to read the GF_v11 attribute table and output attributes of stream width and depth to .csv files for inclusion in the PRMS parameter database.
6. Manning's n was computed using the `seg_slope` parameter from the PRMS parameter database. The equation for computing Manning's n from segment slope was obtained from Bray (1979). That equation is n = 0.1 * seg_slope ** 0.18
  $$
  n = 0.1 * seg\_slope^{0.18}
  $$
7. The parameter `mann_n` was output to a .csv file for inclusion in the PRMS parameter database.


# Scripts
R code used to compute parameters:
```R
mainDir <- "E:/NHM/"
setwd(mainDir)

dataDir <- "G:/GIS/GeospatialFabric/gf_v_1_1_f/shapefiles/"
paramDir <- "E:/NHM/paramdbs/v1_1/paramdb-master-20200413/paramdb-master/"
outDir <- "E:/NHM/paramdbs/v1_1/paramdb-master-20200413/new_params/"

seg.length <- read.table(paste(paramDir, "seg_length.csv", sep=""), sep=",", header=T,
                         stringsAsFactors=F)
seg.slope <- read.table(paste(paramDir, "seg_slope.csv", sep=""), sep=",", header=T,
                        stringsAsFactors=F)
seg.area <- read.table(paste(paramDir, "seg_cum_area.csv", sep=""), sep=",", header=T,
                       stringsAsFactors=F)
hydro.data <- read.table(paste(dataDir, "gf_v11f_segments_geom.txt", sep=""), sep=",", header=T,
                         stringsAsFactors=F)
hydro.data.sort <- hydro.data[order(hydro.data$nsegment_v),]

nsegment <- length(hydro.data.sort$nsegment_v)
seg.index <- 1:nsegment

seg.data <- matrix(-999, nrow=nsegment, ncol=11)
seg.data.df <- data.frame(seg.data)
colnames(seg.data.df) <- c("seg_id", "seg_length", "seg_slope", "seg_cum_area", "seg_width", "seg_width5",
                           "seg_width95", "seg_depth", "seg_depth5", "seg_depth95", "seg_mann_n")

seg.data.df[,1] <- hydro.data.sort$nsegment_v # segment id
seg.data.df[,2] <- seg.length$seg_length # segment length
seg.data.df[,3] <- seg.slope$seg_slope # segment slope
seg.data.df[,4] <- seg.area$seg_cum_area # segment drainage area

for (i in seg.index) {
  seg.data.df[i,5] <- max(hydro.data.sort$WIDTH[i], hydro.data.sort$a_WIDTH[i]) # segment width
  seg.data.df[i,6] <- max(hydro.data.sort$WIDTH5[i], hydro.data.sort$a_WIDTH5[i]) # segment width lower
  seg.data.df[i,7] <- max(hydro.data.sort$WIDTH95[i], hydro.data.sort$a_WIDTH95[i]) # segment width upper
  seg.data.df[i,8] <- max(hydro.data.sort$DEPTH[i], hydro.data.sort$a_DEPTH[i]) # segment depth
  seg.data.df[i,9] <- max(hydro.data.sort$DEPTH5[i], hydro.data.sort$a_DEPTH5[i]) # segment depth lower
  seg.data.df[i,10] <- max(hydro.data.sort$DEPTH95[i], hydro.data.sort$a_DEPTH95[i]) # segment depth upper
}

seg.data.df[,11] <- 0.1 * (seg.data.df$seg_slope^0.18)

write.table(seg.data.df, file=paste(paramDir, "seg_params.csv", sep=""), row.names=F, col.names=T,
            quote=F, sep=',')

# Individual parameter files for paramdb
seg.width <- data.frame(seg.data.df$seg_id, seg.data.df$seg_width)
colnames(seg.width) <- c("$id","seg_width")
write.table(seg.width, file=paste(outDir, "seg_width.csv", sep=""), row.names=F, col.names=T,
            quote=F, sep=',')

seg.depth <- data.frame(seg.data.df$seg_id, seg.data.df$seg_depth)
colnames(seg.depth) <- c("$id", "seg_depth")
write.table(seg.depth, file=paste(outDir, "seg_depth.csv", sep=""), row.names=F, col.names=T,
            quote=F, sep=',')

seg.mann.n <- data.frame(seg.data.df$seg_id, seg.data.df$seg_mann_n)
colnames(seg.mann.n) <- c("$id", "mann_n")
write.table(seg.mann.n, file=paste(outDir, "mann_n.csv", sep=""), row.names=F, col.names=T,
            quote=F, sep=',')
```
