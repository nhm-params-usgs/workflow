mainDir <- "E:/NHM/"
setwd(mainDir)

dataDir <- "G:/GIS/GeospatialFabric/gf_v_1_1_f/shapefiles/"
paramDir <- "E:/NHM/paramdbs/v1_1/paramdb-master-20200413/paramdb-master/"
outDir <- "E:/NHM/paramdbs/v1_1/paramdb-master-20200413/new_params/"

seg.length <- read.table(paste0(paramDir, "seg_length.csv"), sep=",", header=T, stringsAsFactors=F)
seg.slope <- read.table(paste0(paramDir, "seg_slope.csv"), sep=",", header=T, stringsAsFactors=F)
seg.area <- read.table(paste0(paramDir, "seg_cum_area.csv"), sep=",", header=T, stringsAsFactors=F)
hydro.data <- read.table(paste0(dataDir, "gf_v11f_segments_geom.txt"), sep=",", header=T, stringsAsFactors=F)
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

write.table(seg.data.df, file= paste0(paramDir, "seg_params.csv"), row.names=F, col.names=T,
            quote=F, sep=',')

# Individual parameter files for paramdb
seg.width <- data.frame(seg.data.df$seg_id, seg.data.df$seg_width)
colnames(seg.width) <- c("$id","seg_width")
write.table(seg.width, file= paste0(outDir, "seg_width.csv"), row.names=F, col.names=T, quote=F, sep=',')

seg.depth <- data.frame(seg.data.df$seg_id, seg.data.df$seg_depth)
colnames(seg.depth) <- c("$id", "seg_depth")
write.table(seg.depth, file= paste0(outDir, "seg_depth.csv"), row.names=F, col.names=T, quote=F, sep=',')

seg.mann.n <- data.frame(seg.data.df$seg_id, seg.data.df$seg_mann_n)
colnames(seg.mann.n) <- c("$id", "mann_n")
write.table(seg.mann.n, file= paste0(outDir, "mann_n.csv"), row.names=F, col.names=T, quote=F, sep=',')