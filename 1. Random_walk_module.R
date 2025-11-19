rm(list=ls())
gc()

setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # Set the working directory to the directory in which the code is stored
td<-tempdir()
dir.create(td,showWarnings = FALSE)
#
#
###//\/\/\/\/><\/\////\/\/\/\/\/\/\///\\\\\\//\/\/\/\////////////////////////></////##-#
##                                 Psudo-RandomWalk tutorial                               ##-#  
###///\/\/\/><\/\/\////\/\/\//\/\/\/\/\/\///\\\\\\\\\..\.\\\\.\\\\\\\\\\\\\\\><\\\\\##-#
#'
#'
# 0. Load the packages----
library(terra) ; library(dplyr) ; library(sf) ; library(magick)

# 0.2 Load the functions----
functions <- "./Functions" %>% list.files(recursive = FALSE,pattern = ".R$",full.names = TRUE)
lapply(functions,function(x) source(x))

par(bg=NA)

# 1. Create a random study area for our species to move
# a. Study area
study_area<-matrix(rep(1,times=100*100),ncol=100,nrow=100,byrow = TRUE) 
study_area[upper.tri(study_area)]<-NA 
study_area <- study_area[,c(100:1),drop=FALSE]
study_area <- rast(study_area)

# b. Probability of presence
prob_rast <- rnorm(0,1,n=prod(dim(study_area)))
hist(prob_rast)

# Scale between 0-1 to use as probs
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
prob_rast <- range01(prob_rast)

prob_rast <- sort(prob_rast)
prob_rast <- matrix(prob_rast,ncol=ncol(study_area),nrow=nrow(study_area))

# Add some NA values
prob_rast <- prob_rast %>% t()%>% rast()
prob_rast[sample(1:prod(dim(prob_rast)),size=prod(dim(prob_rast))/4)] <- NA
plot(prob_rast)

# c. get a look at the data for the study
par(mfrow=c(1,3))
plot(study_area,col="grey",main="Study region")
study_area[is.na(prob_rast)]<-NA # Add the gaps to the study area

plot(study_area,col="grey",main="Study region with NAs")
plot(prob_rast,main="Species suitability or\nprobability of occurrence",cex.main=0.75)

#
# 2. Select a sample of suitable starting points----
# We need to specify a set of rows and columns----
condition <- study_area==1 & !is.na(prob_rast)
cells_ids<-which(values(condition)==TRUE)

# get a 100 points
cells_ids <- cells_ids %>% sample(100)

starting_points<-data.frame(startingID=1:length(cells_ids),
                  row=rowFromCell(study_area,cells_ids),
                    col=colFromCell(study_area,cells_ids))
start<-study_area
start[cells_ids]<-999

par(mfrow=c(1,2))
plot(condition,main="Suitable areas")
plot(start,col=c("grey","tomato"),main="Starting locations")

# 2. Run the function----
# 2.a Example only with the study area matrix----
par(mfrow=c(1,1),mar=c(1,1,1,1))

for(i in 1:4){
  ColRowDisp<-random_walk(ij=starting_points[i,c("row","col")] %>% unlist() %>% c(),#coordinates in c(row,col) format of the dispersal starging
               m=study_area, # matrix or surface on which the species is moving, if m is a raster it would be transformed into a matrix
                cells_dist = 500, # Number of cells we allow the species to move
                  return_cell = FALSE, # Are species able to return to previous cells
                    jump_cells = 2*i, # Are species allow to "jump" cells? How many cells can they skip?
                      duplicated_cells = TRUE # Are species able to occupy previously used cells?
                        )
  
  points_Disp <- xy_fromRowCol(m=study_area,row=ColRowDisp[,1],col=ColRowDisp[,2])
  col_fun<-colorRampPalette(c("purple","skyblue","gold","firebrick"))
  
  points_Disp$color<-col_fun(nrow(points_Disp))
    
  # Display the information:
  plot(study_area,col="grey80",main=paste("Jump dist=",2*i))
  
  points(x=points_Disp$xcord,y=points_Disp$ycord,type="l",pch=19,
         col="grey25",lwd=2)
  
  points(x=points_Disp$xcord,y=points_Disp$ycord,type="p",pch=19,
         col=points_Disp$color %>% adjustcolor(alpha.f = 0.5))
  
  }
 
# 2.b Example with a species suitability surface----
par(mfrow=c(1,2))
  
  for(i in 1:4){
    ColRowDisp<-random_walk(ij=starting_points[i,c("row","col")] %>% unlist() %>% c(),#coordinates in c(row,col) format of the dispersal starging
                            m=study_area, # matrix or surface on which the species is moving, if m is a raster it would be transformed into a matrix
                            prob_rast = prob_rast,
                            cells_dist = 500, # Number of cells we allow the species to move
                            return_cell = TRUE, # Are species able to return to previous cells
                            jump_cells = 2*i, # Are species allow to "jump" cells? How many cells can they skip?
                            duplicated_cells = TRUE # Are species able to occupy previously used cells?
    )
    
    points_Disp <- xy_fromRowCol(m=study_area,row=ColRowDisp[,1],col=ColRowDisp[,2])
    col_fun<-colorRampPalette(c("purple","skyblue","gold","firebrick"))
    
    points_Disp$color<-col_fun(nrow(points_Disp))
    
    # Display the information:
    plot(study_area,col="grey80",main=paste("Jump dist=",2*i))
    
    points(x=points_Disp$xcord,y=points_Disp$ycord,type="p",pch=19,
           col=points_Disp$color %>% adjustcolor(alpha.f = 0.5))
    
    points(x=points_Disp$xcord,y=points_Disp$ycord,type="l",pch=19,
           col="grey25",lwd=2)
    
    points(x=points_Disp$xcord[1],y=points_Disp$ycord[1],pch="*",cex=4,
           col="green")
    
    plot(prob_rast %>% mask(study_area))
  }

# 3. Make more realistic scenarios----
# Previous examples assume that the probability of any jum happening is uniform, all values have the same probability of happening
# But in real life some scenarios are more likely than other. To simulate this we are going to change the "type" argument of the function
# from "range" to "probs" and feed the algorithm the probability of each happening.
#
# For a species that can jump a maximun of 5 cells but present a dispersal rate closer to 2.5 cells
values <- 1:5

# Compute probabilities using normal density distribution of mean 2.5 and sd of 1
probs <- dnorm(values, mean = mean(values), sd = 1)

# Normalize to sum to 1
probs <- probs / sum(probs)

# Sample 200 values according to these probabilities
set.seed(123)  # for reproducibility
freqs <- rbind(
                xtabs(~sample(values, size = 5000, replace = TRUE, prob = probs)),
                xtabs(~sample(values, size = 5000, replace = TRUE))
                )

par(mfrow=c(1,1))
barplot(freqs,beside=T,col=c("skyblue","purple"),ylab="Frequency",xlab="Jump distance")
legend(legend=c("ProbFunctino","Uniform"),col=c("skyblue","purple"),"topright",
       pch=15,pt.cex=2,title="Sampling type")

# jump_data
jumps <- data.frame(jump=values,probs=probs)
print(jumps)

# 3.1 Compare random jumps with probability driven ones ----

par(mfrow=c(1,2))
# Random jumps
ColRowDisp_random<-random_walk(ij=starting_points[4,c("row","col")] %>% unlist() %>% c(),#coordinates in c(row,col) format of the dispersal starging
                        m=study_area, # matrix or surface on which the species is moving, if m is a raster it would be transformed into a matrix
                        prob_rast = prob_rast,
                        cells_dist = 500, # Number of cells we allow the species to move
                        return_cell = TRUE, # Are species able to return to previous cells
                        jump_cells = 5, # Are species allow to "jump" cells? How many cells can they skip?
                        duplicated_cells = TRUE # Are species able to occupy previously used cells?
)

points_Disp <- xy_fromRowCol(m=study_area,row=ColRowDisp_random[,1],col=ColRowDisp_random[,2])
col_fun<-colorRampPalette(c("purple","skyblue","gold","firebrick"))

points_Disp$color<-col_fun(nrow(points_Disp))

# Display the information:
plot(prob_rast %>% mask(study_area),main=paste("Max jump dist=",5),col=gray(0:100 / 100))
  points(x=points_Disp$xcord,y=points_Disp$ycord,type="l",pch=19,
         col="black",lwd=2)
  points(x=points_Disp$xcord,y=points_Disp$ycord,type="p",pch=19,
         col=points_Disp$color %>% adjustcolor(alpha.f = 0.9),lwd=2)
  points(x=points_Disp$xcord[1],y=points_Disp$ycord[1],pch="*",cex=4,
         col="green")

# Normal distribution jumps
  ColRowDisp_normal<-random_walk(ij=starting_points[4,c("row","col")] %>% unlist() %>% c(),#coordinates in c(row,col) format of the dispersal starging
                                 m=study_area, # matrix or surface on which the species is moving, if m is a raster it would be transformed into a matrix
                                 prob_rast = prob_rast,
                                 cells_dist = 500, # Number of cells we allow the species to move
                                 return_cell = TRUE, # Are species able to return to previous cells
                                 jump_cells = jumps, # Are species allow to "jump" cells? How many cells can they skip?
                                 type="probs",
                                 duplicated_cells = TRUE # Are species able to occupy previously used cells?
  )
  
  points_Disp <- xy_fromRowCol(m=study_area,row=ColRowDisp_normal[,1],col=ColRowDisp_normal[,2])
  col_fun<-colorRampPalette(c("purple","skyblue","gold","firebrick"))
  
  points_Disp$color<-col_fun(nrow(points_Disp))
  
  # Display the information:
  plot(prob_rast %>% mask(study_area),main=paste("Mean jump dist=",2.5),col=gray(0:100 / 100))
  points(x=points_Disp$xcord,y=points_Disp$ycord,type="l",pch=19,
         col="black",lwd=2)
  points(x=points_Disp$xcord,y=points_Disp$ycord,type="p",pch=19,
         col=points_Disp$color %>% adjustcolor(alpha.f = 0.9),lwd=2)
  points(x=points_Disp$xcord[1],y=points_Disp$ycord[1],pch="*",cex=4,
         col="green")

# 4. Trasnform dispersal into density kernells ----
# We can transform the random-walks into continuous density kernells using the d.fromXY function----
  walk_points <- st_as_sf(points_Disp, coords = c("xcord","ycord")) 
  
  par(mfrow=c(1,3))
  # Transform the points directly
  kdensity <- d.fromXY(x=walk_points)
  plot(study_area,legend=F,axis=F,col="grey")
  plot(kdensity$lyr.1,add=T,legend=F) # We are missing the rest of the study area
  
  # If we provided a reference surface we can expand the kernell
  kdensity <- d.fromXY(x=walk_points,pol.ref=st_as_sfc(st_bbox(study_area)),mask.p = NULL)                     
  plot(study_area,legend=F,axis=F,col="grey")
  plot(kdensity$lyr.1,add=T,legend=F)
  
  # with this approach we can include our whole study area, however we know there are areas that are not
  # feasible for the species. We can add a mask polygon/raster to fix this.
  kdensity <- d.fromXY(x=walk_points,pol.ref=st_as_sfc(st_bbox(study_area)),mask.p = study_area)                     
  plot(kdensity$lyr.1,legend=F)
  
# 5. Simulations and dispersals ---- 
# Thorough multiple simulations and random walks we can simulate the dispersal of a species over terrain given a set of 
# pre-calculated starting points.
# Run some examples with pre-trained SDM models:
# model <- "F:/Work_Docs/OneDrive/Projects/Disease_prevalence_niche_centrality/Results/Distance_metrics/Atherurus africanus/Dist_num.rds" %>% readRDS()
  polygons_w <- "E:/OneDrive/Projects/Disease_prevalence_niche_centrality/Results/Distance_metrics/Atherurus africanus/Dist_polygons.rds" %>% readRDS()
  points <- "E:/OneDrive/Projects/Disease_prevalence_niche_centrality/Data/sp_points/clean/Atherurus africanus.csv" %>% read.csv()
  present <- "E:/OneDrive/Projects/Disease_prevalence_niche_centrality/Results/Distance_metrics/Atherurus africanus/Distance_metrics.tif" %>% rast()
  
# We are only going to use one model
  r_present <- present[[1]]
  range_sp <- polygons_w[[1]]
  sp_points <- points %>% st_as_sf(coords = c("decimalLongitude","decimalLatitude"))

# We are going to use the presence points of the species as our starting points
  starting_points <- data.frame(id=1:nrow(sp_points),
                                row=rowFromY(r_present,y=points$decimalLatitude),
                                col=colFromX(r_present,x=points$decimalLongitude))

# 5.a Run the simulations allowing for 500 movements for each starting point ----
# The maximun dispersal distance for this species its belive to be around 2.6 km per night, 
# which correspond to less than 1 cell in our raster. So we are going to assume a dispersal distance of 1
# for each event
  dispersal_points <- list()
  
  for(i in 1:nrow(starting_points)){
    skip_to_next<-FALSE
    
    tryCatch(
    points_w<-random_walk(ij=starting_points[i,c("row","col")] %>% unlist() %>% c(),#coordinates in c(row,col) format of the dispersal starging
                          m=r_present, # matrix or surface on which the species is moving, if m is a raster it would be transformed into a matrix
                          prob_rast = r_present,
                          cells_dist = 500, # Number of cells we allow the species to move
                          return_cell = TRUE, # Are species able to return to previous cells
                          jump_cells = 5, # Are species allow to "jump" cells? How many cells can they skip?
                          duplicated_cells = TRUE # Are species able to occupy previously used cells?
                          ),error= function(e) {skip_to_next<-TRUE})
    if(skip_to_next) next
    
    dispersal_points[[i]] <- points_w  
    print(i)
    }

# 5.b Display and animate the information ----
  coordinates_points <- lapply(dispersal_points, function(x)  xy_fromRowCol(m=r_present,
                                                                            row=x[,"row"],
                                                                            col=x[,"col"]))
  
  par(mfrow=c(1,1))
  plot(r_present)
  lapply(coordinates_points,function(x) x %>% st_as_sf(coords = c("xcord","ycord")) %>% 
           st_geometry() %>% 
           plot(pch=19,add=T,col=grey(0:100/100) %>% adjustcolor(alpha.f = 0.25),cex=0.5))

# 5.c Export the information to do not repeat the simulations----
  dir.create("./Data")
  coordinates_points %>% saveRDS("./Data/Dispersal_info.rds")

# 6. Generate some animations for presentation purpose ----
# Load the data  
  coordinates_points <- readRDS("./Data/Dispersal_info.rds")

# Display the general information    
  plot(r_present, col=viridis::viridis(250))
  plot(range_sp %>% st_geometry(),border="white",add=T)

# Display all points
  lapply(coordinates_points,function(x) x %>% st_as_sf(coords = c("xcord","ycord")) %>% 
           st_geometry() %>% 
           plot(pch=19,add=T,col="white",cex=0.5))
  
  
# Create a list of polygons to engulf the points and display the movement in detail
  polygons_w <- lapply(coordinates_points,function(x) x %>% st_as_sf(coords = c("xcord","ycord")) %>% 
           st_bbox() %>% st_as_sfc() %>% st_buffer(dist=0.5) %>% st_sf())
  
  polygons_w<-do.call("rbind",polygons_w)

# b. Generate the individual plots and save the results into individual folders
dir.create("./Figures/Animation",recursive=T)
  
  for(i in 37:nrow(polygons_w)){
  exit_route<-paste("./Figures/Animation",paste0("Pol_",i),sep="/")
  dir.create(exit_route,recursive = T)
  rp<-crop(r_present,polygons_w[i,] %>% vect())
  
  # Build the plots
    for(j in 1:300){
      png(paste(exit_route,paste0(j,".png"),sep="/"),width = 10,height = 8,units="cm",res=300)
      par(bg=NA)
      plot(rp,col=viridis::viridis(250),breaks=seq(0,1,length.out=250),legend=F)
      index<-1:j
      color_p <- ifelse(j==1,"white","grey")
      points(coordinates_points[[i]][index:j,c("xcord","ycord")],col=c(rep("grey",times=j-1),"white"),pch=19)
      lines(coordinates_points[[i]][index:j,c("xcord","ycord")],col="grey",lty=3)
      dev.off()
    } 
  
   # Generate the gif with the maps
    frames <- data.frame(index=list.files(exit_route,pattern=".png") %>% gsub(pattern="([0-9]+).*$",replacement="\\1") %>% as.numeric(),
               routes=list.files(exit_route,pattern=".png",full.names = T))
    frames <- frames[order(frames$index),]

    frames_g <- image_read(frames$routes)
    image_join(frames_g) %>% image_animate(fps=20) %>% image_write(paste("./Figures",paste0("Run",i,".gif"),sep="/"))

  # Remove the maps
   unlink(exit_route,recursive = T,force=T)
   gc()
  }
  
# c. Combine all the points and rasterize them----  
rasts_w <- lapply(coordinates_points,function(x) x %>% st_as_sf(coords = c("xcord","ycord")) %>% 
            vect() %>% rasterize(r_present,fun="sum"))

rast_w <- rasts_w %>% rast()
rast_w <- sum(rast_w,na.rm=T)
plot(rast_w)


png(paste("./Figures",paste0("total_points",".png"),sep="/"),width = 10,height = 8,units="cm",res=600,bg="transparent")
  par(bg=NA)
    plot(r_present,col=viridis::viridis(250),breaks=seq(0,1,length.out=250),legend=F,alpha=0.5,mar=c(1,1,1,5))
    plot(rast_w,add=T,plg=list(title="count"),col=viridis::inferno(200),alpha=0.5)
    plot(range_sp %>% st_geometry(),add=T,border="white")
dev.off()






  