
tool_exec <- function(in_params, out_params){
  # check/load require packages---------------------------------------------
  if(!requireNamespace("arcgisbinding", quietly = TRUE))
    install.packages("arcgisbinding", quiet = TRUE, dependencies = TRUE)  
  if(!requireNamespace("raster", quietly = TRUE))
    install.packages("raster", quiet = TRUE, dependencies = TRUE)
  if(!requireNamespace("sf", quietly = TRUE))
    install.packages("sf", quiet = TRUE, dependencies = TRUE)  
  if(!requireNamespace("terra", quietly = TRUE))
    install.packages("terra", quiet = TRUE, dependencies = TRUE)
  if(!requireNamespace("ggplot2", quietly = TRUE))
    install.packages("ggplot2", quiet = TRUE, dependencies = TRUE)
  if(!requireNamespace("stars", quietly = TRUE))
    install.packages("stars", quiet = TRUE, dependencies = TRUE)
  if(!requireNamespace("FNN", quietly = TRUE))
    install.packages("FNN", quiet = TRUE, dependencies = TRUE)
  if(!requireNamespace("geosphere", quietly = TRUE))
    install.packages("geosphere", quiet = TRUE, dependencies = TRUE)
  if(!requireNamespace("dplyr", quietly = TRUE))
    install.packages("dplyr", quiet = TRUE, dependencies = TRUE)
  if(!requireNamespace("RANN", quietly = TRUE))
    install.packages("RANN", quiet = TRUE, dependencies = TRUE)
  
  require(arcgisbinding)
  # arc.check_product()
  
  require(raster)         # may need to use install.packages("raster")
  require(sf)
  require(terra)
  require(ggplot2)
  require(stars)
  require(FNN)
  require(geosphere)
  require(dplyr)
  require(RANN)
  
  # define input and output parameters---------------------------------------
  raster_route1 <- in_params[[1]]      # use the TIFF file of proper resolution and extent
  optimum <- in_params[[2]]
  
  radius <- in_params[[3]]
  origin_lat1 <- in_params[[4]]
  origin_lon1 <- in_params[[5]]
  origin_lat2 <- in_params[[6]]
  origin_lon2 <- in_params[[7]]
  
  dest_lat1 <- in_params[[8]]
  dest_lon1 <- in_params[[9]]
  dest_lat2 <- in_params[[10]]
  dest_lon2 <- in_params[[11]]
  
  origin_step_lon2_lon1 <- in_params[[12]]
  origin_step_lat2_lat1 <- in_params[[13]]
  dest_step_lon2_lon1 <- in_params[[14]]
  dest_step_lat2_lat1 <- in_params[[15]]
  
  move_max1 <- in_params[[16]]
  move_min1 <- in_params[[17]]
  dist_range1 <- in_params[[18]]
  opt_range1 <- in_params[[19]]
  
  subTracks <- in_params[[20]]
  
  outputCSV <- out_params[[1]]
  outFc <- out_params[[2]]
  
  # load data--------------------------------------------------------------
  raster_route <- arc.open(raster_route1)
  raster_colA <- arc.raster(raster_route)
  raster_colB <- terra::rast(raster_colA)
  raster_col <- raster_colB - optimum
  
  # To improve efficiency, transform from raster to data frame, and remove N/A values (over water)
  rast <- as.data.frame(raster_col, xy = TRUE)  # xy = TRUE helps add x, y coordinates into the dataframe
  rast$cell <- 1:nrow(rast) # Add column with cell number
  sum(is.na(rast$Band_1)) # check if there is any na values; it would be long, so we applied sum function
  rast2 <- rast[!is.na(rast$Band_1),] #subset without NA cell (over water)
  
  getOptimumLine <- function(id, r5, df) { #returns line number in df of optimum cell
    #id = cell id for central candidate cell
    #r5 = data frame with cell template
    #df = data frame with all candidate cells
    cells <- r5$cell + id #position cell template by adding central cell ID; adds value to this column values
    cellInd <- which(df$cell %in% cells) #reduce to no NA cells in raster (returns indices w/ TRUE)
    vals = df$Band_1[cellInd]  #create vector of cell values (based on the index from above, returns Band_1 values)
    opt <- which.min(abs(vals))  #find the minimum value (index of the first min value of Band_1 values in numeric vector), might have to use which.min(abs(val))
    return(which(df$cell == df$cell[cellInd][opt])) #return the line number (index?) of the optimum cell in large data frame (here, [][] indicates values or subset inside larger set)
  } 
  
  #----Generate cell template index-----
  # generate a data frame that specifies relative cell headings and distances for a subset of raster cells 
  # only used for choosing local optimums (does not need to be very extensive)
  med_x <- median(rast$x) # Find approximate midpoint
  med_y <- median(rast$y) # Find approximate midpoint
  mid_ID <- which.min(abs(rast$x - med_x) + abs(rast$y - med_y)) # get ID of central cell (index)
  mid_x <- rast$x[mid_ID] # get x coordinate of central cell; can check the row using rast[139199, ]
  mid_y <- rast$y[mid_ID] # get y coordinate of central cell
  radius_km <- radius # km ... how big does this need to be - only a few km; multiplied by 1000 to get the km distance; changed to 150 from 30 for Lesser Yellowlegs
  rtExt <- destPoint(p=c(mid_x, mid_y), b = 90, d=radius_km*1000) #get rightmost extent of sampled area
  ltExt <- destPoint(p=c(mid_x, mid_y), b = 270, d=radius_km*1000) #get leftmost extent of sampled area
  upExt <- destPoint(p=c(mid_x, mid_y), b = 0, d=radius_km*1000) #get upper extent of sampled area
  loExt <- destPoint(p=c(mid_x, mid_y), b = 180, d=radius_km*1000) #get upper extent of sampled area
  r2 <- rast[which(rast$x >= ltExt[1] & rast$x <= rtExt[1] & rast$y >= loExt[2] & rast$y <= upExt[2]),] #cut larger raster down to sampling area
  r2m <- matrix(c(r2$x, r2$y), ncol = 2) #make a coordinate matrix to use with distGeo (calculates the shortest distance between two points in geodesic)  
  r2$dist <- round(distGeo(p1 = c(mid_x, mid_y), p2= r2m)) #use distGeo to get all distances from mid_x and mid_y
  r3 <- r2[which(r2$dist <= radius_km*1000),] #cut data frame down to cells within range threshold (radius_km); generates a circle not the rectangle raster! 
  r3m <- matrix(c(r3$x, r3$y), ncol = 2) #make matrix for bearing function
  r3$head <- bearing(p1 = c(mid_x, mid_y), p2= r3m) #get bearings for all cells
  r3$head <- ifelse(r3$head < 0, 360+r3$head,  r3$head) #correct for 360 degree wraparound
  r3$cell <- r3$cell - mid_ID #subtract middle cell ID to create cell index.
  
  #-----Initial parameters for simulation-----
  num_birds <- 200
  maxDays <- 60
  origin_area <- data.frame(lat1 = origin_lat1, lat2 = origin_lat2,  lon1 = origin_lon1, lon2 = origin_lon2) #origins generated here
  dest_area <- data.frame(lat1 = dest_lat1, lat2 = dest_lat2,  lon1 = dest_lon1, lon2 = dest_lon1) #destinations generated here
  
  #-----Initializing bird agents-----
  migrants <- data.frame(id = 1:num_birds, 
                         originlon = sample(size = num_birds, x = seq(from = origin_area$lon2, to = origin_area$lon1, by = origin_step_lon2_lon1), replace = T),
                         originlat = sample(size = num_birds, x = seq(from = origin_area$lat2, to = origin_area$lat1, by = origin_step_lat2_lat1), replace = T),
                         destlon = sample(size = num_birds, x = seq(from = dest_area$lon2, to = dest_area$lon1, by = dest_step_lon2_lon1), replace = T),
                         destlat = sample(size = num_birds, x = seq(from = dest_area$lat2, to = dest_area$lat1, by = dest_step_lat2_lat1), replace = T),
                         total_dist = 0,
                         sigma = 50,      # do we have to increase the sigma when we changed 300 to 3000; sigma is stochastic or diffusion parameter! 
                         mot=0.8,
                         wedge_size = 40, #half size of search wedge
                         move_max = move_max1*1000, #maximum movement distance; changed to 1500 from 300
                         move_min = move_min1*1000, #minimum movement distance; changed to 150 from 30
                         dist_range = dist_range1, #how many meters of range for choosing candidate cells; changed to 15000 from 3000
                         opt_range = opt_range1)  #diameter of search area surrounding candidate cell for choosing optimum cell; changed to 15000 from 3000  
  migrants$total_dist <- round(distGeo(p1<-matrix(c(migrants$originlon,migrants$originlat), ncol=2),
                                       p1<-matrix(c(migrants$destlon,migrants$destlat), ncol=2))) 
  
  #-----generate output tempate-----
  allTrack <- data.frame(id = migrants$id,
                         day = 0,
                         lon = migrants$originlon,
                         lat = migrants$originlat,
                         cell_number = 0,
                         cell_val = 0,
                         stepDist = 0,
                         destBear = 0,
                         randBear = 0,
                         arrived = FALSE)
  
  #-----generate output tempate-----
  allTrack <- data.frame(id = migrants$id,
                         day = 0,
                         lon = migrants$originlon,
                         lat = migrants$originlat,
                         cell_number = 0,
                         cell_val = 0,
                         stepDist = 0,
                         destBear = 0,
                         randBear = 0,
                         arrived = FALSE)
  
  #-----Generate parameters used in model-----
  mr = matrix(c(rast2$x, rast2$y), ncol = 2) #matrix of all xy for rast2
  mtrack = matrix(c(migrants$originlon, migrants$originlat), ncol = 2) #matrix of initial lat and lon. 
  nearest = get.knnx(data = mr, query = mtrack, k=1, algorithm = "kd_tree") #only do before loop; k=1 for nearest 1 neighbor
  allTrack$cell_number = rast2$cell[nearest$nn.index]  #look up cell numbers
  allTrack$lon <- rast2$x[nearest$nn.index]   #set track lat and lon to nearest cell coordinates (ensures start points are on land)
  allTrack$lat <- rast2$y[nearest$nn.index]   #set track lat and lon to nearest cell coordinates
  startMat <- matrix(c(allTrack$lon, allTrack$lat), ncol=2) #matrix of new start coordinates
  destMat <- matrix(c(migrants$destlon, migrants$destlat), ncol=2) #matrix of destination coordinates - keep for multiple uses
  totalDist <- round(distGeo(p1=startMat , p2=destMat)) #recalculate total distance to account for changes to coordinates in migrants
  r4 <- r3[which(r3$dist<=opt_range1),] # changed to 15000 from 3000
  
  plot=FALSE #plot tracks during model run (slows things down)
  for (day in 1:maxDays) {  
    track <- allTrack[allTrack$day == day-1,] # reduce working dataframe to past day; begin with day 0
    cont <- which(!track$arrived)             # index of birds that are still moving (not arrived)
    track <- track[cont,]                     # subset data frame to moving birds
    keepInd <- which(migrants$id %in% track$id) #index of moving birds for subsetting migrants and parameter vectors
    print(paste0("day: ", day, "  tracks complete: ", nrow(migrants) - nrow(track))) #print progress
    dMat <- destMat[keepInd,] #subset destination matrix as tracks get completed
    if(nrow(track)==0) { #exit loop if all migrations are complete
      print("Simulation complete") 
      break #break loop
    }
    
    curLocMat <- matrix(c(track$lon, track$lat), ncol=2) #get matrix for all current location coordinates
    destDist <- round(distGeo(p1=curLocMat , p2=dMat)) #get current distances to destinations
    distRatio <- destDist/totalDist[keepInd]                        #ratio of remaining distance to starting distance
    bear <- bearing(p1 = curLocMat, p2=dMat)           #get current bearings to destinations; bearing measured against north=0 and clockwise angle
    rand <- runif(n=nrow(track), min=-1, max=1)        #generate a random number for each agent from -1 to 1
    randBear <- bear + rand * (migrants$wedge_size[keepInd]) #generate a bearing - bearing to destination plus some randomness
    randBear <- ifelse(randBear < 0, 360+randBear, randBear) #correct for 360 degree wraparound; used min=-1 above and thus, leads to negative value
    stepDist <- round((migrants$mot[keepInd] * migrants$move_max[keepInd] * distRatio) + rnorm(n = nrow(track)) * migrants$sigma[keepInd]) #generate a step distance similar to abmR
    stepDist <- pmax(stepDist, migrants$move_min[keepInd])   #keep step distance above minimum threshold for each agent
    stepPt <- destPoint(p=curLocMat, b=randBear, d=stepDist) #get coordinates of candidate location
    done <- which(stepDist > destDist) #if bird is close to destination, just use destination coordinates as candidate cell center
    stepPt[done ,1] = migrants$destlon[keepInd][done] #for birds close to destination set lon to desination coordinate
    stepPt[done ,2] = migrants$destlat[keepInd][done] #for birds close to destination set lat to desination coordinate
    stepDist[done] = destDist[done]
    
    #Find nearby cells and choose an optimum value
    dpx1 <- min(stepPt[,1]) - 1; dpx2 <- max(stepPt[,1]) + 1  #narrow x range to speed calculations
    dpy1 <- min(stepPt[,2]) - 1; dpy2 <- max(stepPt[,2]) + 1  #narrow y range to speed calculations
    subInd <- which(rast2$x>dpx1 & rast2$x<dpx2 & rast2$y>dpy1 & rast2$y<dpy2) #Index of subset of reduced rast2
    r2Sub <- rast2[subInd,] #subset of large data frame
    mrSub <- matrix(c(r2Sub$x,r2Sub$y), ncol=2)  #x y matrix of subseted data frame
    nn = RANN::nn2(data = mrSub, query = stepPt, k = 1, treetype = "kd") #get nearest cell to candidate x y location
    ptIDs <- r2Sub$cell[nn$nn.idx] #index of cell IDs for candidate points; just 1 neighbor or 1 point in this case!
    optLines <- mapply(getOptimumLine, id = ptIDs, MoreArgs = list(r5=r4, df=r2Sub)) #get the line number for the optimum cell value
    r2Opt <- r2Sub[optLines,] #optimums for all agents
    
    #update track data frame and append to allTracks
    track$day <- day       #update day 
    track$lon <- r2Opt$x   #update x
    track$lat <- r2Opt$y   #update y
    track$cell_number <- r2Opt$cell  #update cell id
    track$cell_val <- r2Opt$Band_1  #update cell value
    track$stepDist <- stepDist
    track$destBear <- bear
    track$randBear <- randBear
    track$arrived[done] = TRUE       #update agent status
    allTrack <- rbind(allTrack, track)  #append data
    row.names(allTrack) <- 1:nrow(allTrack) #control annoying row names
  }
  
  #Summarize data 
  sumTracks <- allTrack %>% group_by(id) %>%
    dplyr::summarise(dist = sum(stepDist), cellVal = mean(cell_val), days = n())
  
  sumTracks <- sumTracks[order(sumTracks$cellVal),] #sort data by cell values, smallest value (best track) first; e.g. smallest mean DEM value or temperature value! 
  subTrk <- sumTracks[1:subTracks,] #subset data to 10 best tracks; changed from 10 to 20 best tracks
  tm <- 0:max(subTrk$days)  #time steps to plot
  
  #subset to lines and points for plotting
  plt <- allTrack[allTrack$id %in% subTrk$id,]
  plt <- plt[plt$day %in% tm,]
  plt$id <- as.factor(plt$id)
  
  # save the file in the given folder as CSV
  write.csv(plt, outputCSV)
  
  # writes the feature class as a shapefile
  arc.write(outFc, plt, coords = c("lon", "lat"), shape_info = list(type = "Point", WKID = 4326), overwrite = TRUE)
  
  return(out_params)
}
