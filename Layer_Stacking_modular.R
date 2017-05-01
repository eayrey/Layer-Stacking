###############################################################################
#                             LAYER STACKING                                  #
#                                                                             #
#                                                                             #
# Layer stacking is a novel algorithm for segmenting individual trees from a  # 
# LiDAR point cloud. Layer stacking works best on small-medium sized areas    #
# (<10ha), and with Leaf Off High resolution LiDAR (+5ppm). It is also        #        
# important to normalize the LiDAR data using a digital elevation model prior #
# to running layer stacking.                                                  #
#                                                                             #
# As with most segmentation algorithms, there are parameters which may        #
# require region and LiDAR specific fine tuning.                              #
#                                                                             #
###############################################################################

#Please specify the folder in which LiDAR files are stored. Only ".las" files 
# can be stored here
fpath <- "C:\\location\\location"

#Please specify the folder in which the processed csv files will be stored.
Output <- "C:\\location\\location"

#Number of processor cores to use for parallel processesing
#Default is the number of cores you have in your machine minus 1.
n= detectCores(all.tests = FALSE, logical = TRUE)-1

#If this is a predominantly hardwood stand (trees without leaves/needles), 
# set this switch to true.
hw=FALSE

#Consider using this switch if you're working with very high density LiDAR, or small tree crowns, or you feel too few trees were segmented
#Essentially runs a finer filter over the LiDAR and may detect trees very small or close to one another
d=FALSE

#Specify a cutoff value for filtering tree crowns, which referes to the
#number standard deviations below which an abnormally large tree crown 
#segment will be allowed. A small number may be needed in dense forests with
#small tree crowns, a larger number may be needed in forests with large trees
c=2
if (hw==TRUE){
  c=2.5
}

#Threshold (meters) for removing small trees. Generally trees shorter than this threshold
#will be removed.
t=3

####### Caution, these are 'last resort' parameters #########################

#Refers to the width of buffers places around each cluster. A larger value may be 
#needed with low density LiDAR if the points are spaced far apart from one another.
#0.5m was ideal for the high density LiDAR tested.
buf_width=0.6

#Refers to the width of the tree's core. Clusters touching this core
#will be considered part of that tree. A larger value may be 
#needed with low density LiDAR if the points are spaced far apart from one another
#0.6m was ideal for the high density LiDAR tested.
core_width=0.8

#the following packages are required
library(rLiDAR)
library(doParallel)
library(raster)
library(sp)
library(rgeos)
library(prevR)
library(fpc)

#################Buffer function############################
buffer.groups<-function(group, buf_width=.6){
  pts <- SpatialPoints(group[,1:2])
  buffered_layer=gBuffer(pts, width=buf_width)
  buffered_layer=disaggregate(buffered_layer)
  for (i in 1:length(buffered_layer)){
    poly=buffered_layer[i]
    outerRings = Filter(function(f){f@ringDir==1},poly[1]@polygons[[1]]@Polygons)
    outerBounds = SpatialPolygons(list(Polygons(outerRings,ID=1)))
    if (i==1){buffers=outerBounds}else{buffers=rbind.SpatialPolygons(buffers, outerBounds, makeUniqueIDs =TRUE)}
  }
  return(buffers)
}

##########Layer Stacking###################################

layer_stacker=function(inFile, buffer.groups,Output, p, n, t, c,d, hw, buf_width, core_width){
  #seperate into layers at 1m intervals
  inFile=data.frame(inFile)
  inFile=unique(inFile)
  inFile$layer= round(inFile$Z)
  inFile=subset(inFile,inFile$layer>0)
  layers=split(inFile, inFile$layer)

  #####################Local MAxima########################
  #Generates the number of local maxima for different LiDAR smoothing levels. 

  high_pts=inFile
  e=extent(min(high_pts[,1]), max(high_pts[,1]),min(high_pts[,2]),max(high_pts[,2]))
  #1m Raster, .5m Raster, and 2m Raster
  r1 <- raster(e, ncol=(e[2]-e[1]), nrow=(e[4]-e[3]))
  r_forth<- raster(e, ncol=(e[2]-e[1])*4, nrow=(e[4]-e[3])*4)
  x1 <- rasterize(inFile[,1:2], r1, inFile[,3], fun=max)
  smoothed1=focal(x1, w=matrix(1, nrow=3, ncol=3),fun=function(x){mean(x,na.rm=TRUE)})
  localmax1 <- focal(smoothed1,w=matrix(1,nrow=3,ncol=3), fun = max, pad=TRUE, padValue=NA)
  trueLM1<-(smoothed1==localmax1)*x1
  trueLM1[trueLM1==0]<-NA

  #################Initial clusters###########################
  cl=makeCluster(n)
  registerDoParallel(cl)
  rasters=foreach(i=1:(length(layers)), .errorhandling="pass") %dopar% {
    library(raster)
    library(sp)
    library(rgeos)
    library("fpc")
    layer=data.frame(layers[i])
    
    if (i<=3){
      if (length(layer)>1){
        db=dbscan(layer[,1:2], eps=1.5, MinPts=5)
        lowlayer=layer
        lowlayer$group<-db$cluster
        lowlayer=subset(lowlayer, group==0)
        layer<-lowlayer
      }
    }
    
    LM1=trueLM1
    LM1[LM1<i]<-NA
    points1=rasterToPoints(LM1)
    
    if (nrow(points1)==0 || nrow(layer)<=nrow(points1) ){treecut1=kmeans(layer[,1:2],1)
    }else{treecut1=tryCatch({kmeans(layer[,1:2],points1[,1:2])}, error=function(e){kmeans(layer[,1:2],length(points1[,1]), nstart=15)})}
    
    layer$group1<-treecut1$cluster
    lay_groups1=split(layer, layer$group1)
    buffered_layer1<-lapply(lay_groups1,buffer.groups)
    #for each layer, bind them
    buffered_layer1s=do.call(rbind.SpatialPolygons, c(buffered_layer1, makeUniqueIDs = TRUE))
    bbb=buffered_layer1s
    #bbb=rbind.SpatialPolygons(buffered_layer1s,buffered_layer_halfs, buffered_layer_3qs,buffered_layer_1hcs,buffered_layer_halfhcs,buffered_layer_3qhcs, makeUniqueIDs = TRUE)
    coords=sapply(bbb@polygons, function(x) coordinates(x@Polygons[[1]]))
    if (length(bbb)>1 && length(bbb)==length(coords)){
      bbb=bbb[!duplicated(coords)]}
    
    x <- rasterize(bbb, r_forth)
    layer_raster=x>=0
    #Adding additional weight to polygons near the top
    if (hw==FALSE){
      if (i/length(layers)>=.7){
        layer_raster7=x>=0
        layer_raster=layer_raster+layer_raster7
      }
      if (i/length(layers)>=.8){
        layer_raster8=x>=0
        layer_raster=layer_raster+layer_raster8
      }
      if (i/length(layers)>=.9){
        layer_raster9=x>=0
        layer_raster=layer_raster+layer_raster9
      }
    } 
    layer_raster[is.na(layer_raster[])] <- 0 
    c(layer_raster)
  }
  #stopCluster(cl)
  
  #Construct the overlap map from individual rasters
  for (r in 1:length(rasters)){
   #plot(rasters[[r]][[1]])
   if (r==1){overlap_map_forth=rasters[[r]][[1]]}else{overlap_map_forth=overlap_map_forth+rasters[[r]][[1]]}
  }

  #Local maxima on the overlap maps
  r_half<- raster(e, ncol=(e[2]-e[1])*2, nrow=(e[4]-e[3])*2)
  overlap_map_half=resample(overlap_map_forth, r_half)
  overlap_map1=resample(overlap_map_forth, r1)

  smoothed1=focal(overlap_map1, w=matrix(1, nrow=3, ncol=3),fun=function(x){mean(x,na.rm=TRUE)})
  smoothed_half=focal(overlap_map_half, w=matrix(1, nrow=3, ncol=3),fun=function(x){mean(x,na.rm=TRUE)})
  smoothed_forth=focal(overlap_map_forth, w=matrix(1, nrow=3, ncol=3),fun=function(x){mean(x,na.rm=TRUE)})

  localmax1 <- focal(smoothed1,w=matrix(1,nrow=3,ncol=3), fun = max, pad=TRUE, padValue=NA)
  localmax_half <- focal(smoothed_half,w=matrix(1,nrow=3,ncol=3), fun = max, pad=TRUE, padValue=NA)
  localmax_forth <- focal(smoothed_forth,w=matrix(1,nrow=3,ncol=3), fun = max, pad=TRUE, padValue=NA)

  trueLM1<-(smoothed1==localmax1)*overlap_map1
  trueLM_half<-(smoothed_half==localmax_half)*overlap_map_half
  trueLM_forth<-(smoothed_forth==localmax_forth)*overlap_map_forth

  trueLM1[trueLM1==0]<-NA
  trueLM_half[trueLM_half==0]<-NA
  trueLM_forth[trueLM_forth==0]<-NA

  #####Clustering for polygons
  #cl=makeCluster(n)
  #registerDoParallel(cl)
  polys=foreach(i=1:(length(layers)), .errorhandling="pass") %dopar% {
   library(raster)
   library(sp)
   library(rgeos)
   library("fpc")
  
   layer=data.frame(layers[i])
   layer=unique(layer)
  
    if (i<=3){
      if (length(layer)>1){
        db=dbscan(layer[,1:2], eps=1.5, MinPts=5)
        lowlayer=layer
       lowlayer$group<-db$cluster
        lowlayer=subset(lowlayer, group==0)
       layer<-lowlayer
     }
    }
  
    points1=rasterToPoints(trueLM1)
    points_half=rasterToPoints(trueLM_half)
    
    points_forth=rasterToPoints(trueLM_forth)
    
    if (nrow(points1)==0 || nrow(layer)<nrow(points1) ){treecut1=kmeans(layer[,1:2],1)
    }else{treecut1=tryCatch({kmeans(layer[,1:2],points1[,1:2])}, error=function(e){kmeans(layer[,1:2],(length(points1[,1])-1), nstart=5)})}
  
    if (nrow(points_half)==0 || nrow(layer)<nrow(points_half)  ){treecut_half=kmeans(layer[,1:2],1)
    }else{treecut_half=tryCatch({kmeans(layer[,1:2],points_half[,1:2])}, error=function(e){kmeans(layer[,1:2],(length(points_half[,1])-1), nstart=5)})}


    if (nrow(points_forth)==0 || nrow(layer)<nrow(points_forth)  ){treecut_forth=kmeans(layer[,1:2],1)
    }else{treecut_forth=tryCatch({kmeans(layer[,1:2],points_forth[,1:2])}, error=function(e){kmeans(layer[,1:2],(length(points_forth[,1])-1), nstart=5)})}

    layer$group_forth<-treecut_forth$cluster
    lay_groups_forth=split(layer, layer$group_forth)
    buffered_layer_forth<-lapply(lay_groups_forth,buffer.groups)
    buffered_layer_forth=do.call(rbind.SpatialPolygons, c(buffered_layer_forth, makeUniqueIDs = TRUE))
     

    layer$group1<-treecut1$cluster
    layer$group_half<-treecut_half$cluster
  
    lay_groups1=split(layer, layer$group1)
    lay_groups_half=split(layer, layer$group_half)
  
    buffered_layer1<-lapply(lay_groups1,buffer.groups)
    buffered_layer_half<-lapply(lay_groups_half,buffer.groups)
  
    buffered_layer1s=do.call(rbind.SpatialPolygons, c(buffered_layer1, makeUniqueIDs = TRUE))
    buffered_layer_halfs=do.call(rbind.SpatialPolygons, c(buffered_layer_half, makeUniqueIDs = TRUE))

    bbb=rbind.SpatialPolygons(buffered_layer1s,buffered_layer_halfs, makeUniqueIDs = TRUE)
    
    coords=sapply(bbb@polygons, function(x) coordinates(x@Polygons[[1]]))
    bbb=bbb[!duplicated(coords)]
    c(bbb)
  }
  #stopCluster(cl)

  #tree core identification
  trueLM_half[trueLM_half<=t]<-NA
  cores=rasterToPoints(trueLM_half, spatial=TRUE)
  cores@data=data.frame(extract(overlap_map_half, cores))
  names(cores@data)="overlaps"
  cores=cores[rev(order(cores$overlaps)),]

  #Create buffers around cores, dissolve them. 
  if (d==TRUE){
    for (co in 1:length(cores)){
      if(co==1){buf_cores=gBuffer(cores[co,1], width=core_width, quadsegs=2)}else{buf_cores=rbind.SpatialPolygons(buf_cores,gBuffer(cores[co,1], width=core_width, quadsegs=3),makeUniqueIDs = TRUE)}
    }
  }else{
    for (co in 1:length(cores)){
      if(co==1){buf_cores=gBuffer(cores[co,1], width=1, quadsegs=2)}else{buf_cores=rbind.SpatialPolygons(buf_cores,gBuffer(cores[co,1], width=1, quadsegs=3),makeUniqueIDs = TRUE)}
    }
  }

  aggregated=aggregate(buf_cores)
  bcs=disaggregate(aggregated)
  
  for (co in 1:length(bcs)){
    if(co==1){buf_cores=gBuffer(gCentroid(bcs[co,1]), width=core_width, quadsegs=3)}else{buf_cores=rbind.SpatialPolygons(buf_cores,gBuffer(gCentroid(bcs[co,1]), width=core_width, quadsegs=3),makeUniqueIDs = TRUE)}
  }

  #cores are scrambled from aggregation, need to identify which is which 
  #and assign an overlap value
  reorganizer=function(core, buff_core){
   if (point.in.SpatialPolygons(core[2], core[3], buff_core)){
      core[1]
    }
  }

  dfc <- as.data.frame(cores)
  #cl=makeCluster(n)
  #registerDoParallel(cl)
  overlapn2=foreach(bc=1:length(buf_cores), .errorhandling="pass") %dopar% {
   library(prevR)
   bc=buf_cores[bc,1]
   overlapn=apply(dfc, 1,reorganizer, buff_core=bc)
   max(unlist(overlapn))
  }
  #stopCluster(cl)

  overlapn2=t(data.frame(overlapn2))
  row.names(overlapn2)=sapply(slot(buf_cores, "polygons"), function(x) slot(x, "ID"))
  overlapn2=data.frame(overlapn2)
  buf_cores=SpatialPolygonsDataFrame(buf_cores, data=overlapn2)
  buf_cores=buf_cores[rev(order(buf_cores$overlapn2)),]

  #cl=makeCluster(n)
  #registerDoParallel(cl)
  polys2=foreach(l=1:(length(layers)), .errorhandling="pass") %dopar% {
   layer_ps=polys[[l]][[1]]
   intersections=lapply(seq(1:length(layer_ps)), function(x){sum(na.omit(over(buf_cores, layer_ps[x,])))})
   layer_ps=layer_ps[!intersections>1]
   layer_ps
  }
  #stopCluster(cl)

  ####################CLIP OUT TREE POLYGONS##############################
  #cl=makeCluster(n)
  #registerDoParallel(cl)
  tree_polys=foreach(c=1:length(buf_cores), .errorhandling="pass") %dopar% {
   library(sp)
   library(rgeos) 
   library(prevR)
   core=buf_cores[c,1]
   poly_attributes=slot(core,"polygons")
   coords=matrix(sapply(core@polygons, function(x) coordinates(x@Polygons[[1]])),ncol=2,byrow=F)
  
   poly_areas=NULL
   for (l in 1:(length(polys2))){
     layer_ps=polys2[[l]]
     if (typeof(layer_ps)=="S4"){
       ids=sapply(slot(layer_ps, "polygons"), function(x) slot(x, "ID"))
       polysin=lapply(ids, function(x) {point.in.SpatialPolygons(coords[,1],coords[,2], layer_ps[match(x,ids)])})
       polysin=unlist(lapply(polysin,any))==T
       layer_in=layer_ps[polysin]
       layer_in=SpatialPolygonsDataFrame(layer_in, data=data.frame(row.names=sapply(slot(layer_in, "polygons"), function(x) slot(x, "ID")),rep(l,length(layer_in))))
       colnames(layer_in@data)=c("height")
       #nnn<-1:length(layer_ps)
       #results[[l]][[2]]=layer_ps[-nnn[polysin]]
       if (l==1){in_polys=layer_in
       }else{in_polys=rbind(in_polys,layer_in,makeUniqueIDs = TRUE)}
     }
   }
  
   if (length(in_polys)>1){
     #ANOTHER somewhat ARBITRARY THRESHOLD FOR REMOVING POLYS NOT CENTERED AT THE CORE OF THE TREE
      center_pt=gCentroid(core)
      in_poly_centroids=lapply(seq(to=nrow(in_polys)), function(x){gCentroid(in_polys[x,])})
      in_poly_centroids=matrix(coordinates(in_poly_centroids), ncol=2, byrow=TRUE)
      dists=spDistsN1(in_poly_centroids, center_pt)
      if (length(dists)>0){
       sd_cutoff=sd(dists[dists>core_width])*c
       goods=dists<median(dists)+sd_cutoff
       in_polys=in_polys[goods,]
      }
    
    #ANOTHER somewhat ARBITRARY THRESHOLD FOR REMOVING ABNORMALLY LARGE POLYGONS
      poly_areas=gArea(in_polys, byid=TRUE)
      if (length(poly_areas)>0){
       sd_cutoff=sd(poly_areas[poly_areas>.8])*c
       goods=poly_areas<median(poly_areas)+sd_cutoff
       in_polys=in_polys[goods,]
      }
    #Add core to the tree
     core_poly=SpatialPolygonsDataFrame(core, data.frame(row.names=sapply(slot(core, "polygons"), function(x) slot(x, "ID")), c+1000000 ))
      colnames(core_poly@data)=c("height")
      in_polys=rbind(in_polys,core_poly,makeUniqueIDs = TRUE)
      c(in_polys)
    }else{NULL}
  }
  #stopCluster(cl)
  
  tree_polys=Filter(Negate(function(x) is.null(unlist(x))), tree_polys)
#Need to remove trees with few layers (they're fragments)
  large_trees=function(x){
    if(length(x[[1]][[1]])>=t*2){
      TRUE
    }else{FALSE}
  }
  bigts=unlist(lapply(tree_polys, large_trees))
  tree_polys=tree_polys[bigts]

  ###################CLIP LIDAR POINTS######################

  #cl=makeCluster(n)
  #registerDoParallel(cl)
  segmented=foreach(l=1:length(layers), .combine='rbind', .errorhandling="pass") %dopar% {
    library(sp)
    library(prevR)
    layer=data.frame(layers[l])
    layer$TreeNumber=NA
    colnames(layer)=c("X","Y","Z","Intensity","ReturnNumber","Layer","TreeNumber")
    for (t in 1:length(tree_polys)){
      tree=tree_polys[t][[1]][[1]]
    #tree=in_polys
      polyin=tree@data==l
      nnn<-1:length(polyin)
      num_polys=seq(1:sum(polyin))
    
      core_poly=tree[tree@data$height>50,]
      pointsinc=point.in.SpatialPolygons(layer[,1],layer[,2], core_poly)
      if (length(nnn[polyin])>0){
        polyin=tree[nnn[polyin],]
        pointsin=(lapply(num_polys, function(x) {point.in.SpatialPolygons(layer[,1],layer[,2], polyin[x,])}))
        pointsin=Reduce("|",pointsin)
        pointsin=pointsin | pointsinc
        layer[pointsin,]$TreeNumber=t
        
        if (sum(pointsin>0)){
          layer[pointsin,]$TreeNumber=t
        }
      }
    }
    layer
  }
  stopCluster(cl)

  write.csv(segmented, file=paste(Output, substr(a[p], 1, nchar(a[p])-4), "_Done", ".csv",sep=""))
}

a <- list.files(fpath)

for (p in 1:length(a)){
  inFile=readLAS(paste(fpath,a[p],sep="\\"))
  ptm <- proc.time()
  layer_stacker(inFile, buffer.groups,Output, p, n, t, c,d, hw, buf_width, core_width)
  proc.time() - ptm
  print (p/length(a)*100)
} 
