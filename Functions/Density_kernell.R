#
# 0                                                 0
#  \                                               /
#   \ /----0                            \   --0---0
#    0      \                          \     / \
#            0---0 Density Kernell 0---0  0-0   0-----
#           / \   0  0--0             \/            \  
#          0   \  |\   /                           000
#               0-0 0-0--------0                  /
#
# Transforms points into a density kernell with a smooth distance based on the size of the matrix/raster

d.fromXY<-function(x, # Points used to build the density kernell
                   y.size=c(1000,1000), # size of the density kernell matrix
                   inv=FALSE, # should the final matrix need to be inverse?
                   scale01=TRUE, # should the final matrix need to be scaled between 0_1?
                   pol.ref=NULL, # reference polygon to create the extension of the kernell raster
                   mask.p=NULL # Does the initial kernell need to be masked?
){
  # 0. load the needed libraries----
  list.of.packages<-c("dplyr","terra","MASS","sf","KernSmooth")
  
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  lapply(list.of.packages,require,character.only=TRUE)
  rm(list.of.packages,new.packages)
  
  # Run the function
  if(is.null(pol.ref)){
    coord <- x %>% st_coordinates()
    
    # Check the quantiles of our coordinates to avoid errors in kde2d
    if(bandwidth.nrd(coord[,1])==0){
      coord[,1]<-coord[,1]+rnorm(n=length(coord[,1]),mean=0.05,sd=0.0025)
    }
    
    if(bandwidth.nrd(coord[,2])==0){
      coord[,2]<-coord[,2]+rnorm(n=length(coord[,2]),mean=0.05,sd=0.0025)
    }
    
    bx <- diff(range(coord[,"X"])) %>% abs() * 0.1
    by <- diff(range(coord[,"Y"])) %>% abs() * 0.1
    
    k<-KernSmooth::bkde2D(coord,bandwidth = c(bx,by),range.x = list(range(coord[,"X"]),range(coord[,"Y"])),
                          gridsize = y.size,truncate=TRUE)
    
    k_mat <- apply(t(k$fhat), 2, rev)
    
    
  }else{
    coord <- x %>% st_coordinates()
    
    # Check the quantiles of our coordinates to avoid errors in kde2d
    if(bandwidth.nrd(coord[,1])==0){
      coord[,1]<-coord[,1]+rnorm(n=length(coord[,1]),mean=0.05,sd=0.0025)
    }
    
    if(bandwidth.nrd(coord[,2])==0){
      coord[,2]<-coord[,2]+rnorm(n=length(coord[,2]),mean=0.05,sd=0.0025)
    }
    
    bx <- diff(st_bbox(pol.ref)[c(1,3)]) %>% abs() * 0.1 # Maybe we can add an extra argument to the function here to control the size of the buffer
    by <- diff(st_bbox(pol.ref)[c(2,4)]) %>% abs() * 0.1
    
    k<-KernSmooth::bkde2D(coord,bandwidth = c(bx,by),range.x = list(st_bbox(pol.ref)[c(1,3)],
                                                                    st_bbox(pol.ref)[c(2,4)]),
                          gridsize = y.size,truncate=TRUE)
    
    k_mat <- apply(t(k$fhat), 2, rev)
    
  }
  
  if(inv==TRUE){
    k_mat<-max(k_mat)-k_mat
  }
  # Create the reference raster
  if(is.null(pol.ref)){
    r1 = rast(k_mat,extent=c(min(k$x1,na.rm=T),max(k$x1,na.rm=T),min(k$x2,na.rm=T),max(k$x2,na.rm=T)))
  }else{
    r1 = rast(k_mat,extent=st_bbox(pol.ref),crs=crs(pol.ref %>% vect()))  
  }
  
  if(!is.null(mask.p)){
    if("SpatRaster" %in% class(mask.p)){
      mask.p <- resample(mask.p,r1)
      r1 <- r1 %>% mask(mask.p)
        }else{
      r1 <- r1 %>% mask(mask.p %>% vect())
        }
    }
  
  # Scale the matrix between 0-1
  scale.f<-function(x,...){
    
    min.m <- min(x,na.rm = TRUE)
    max.m <- max(x,na.rm = TRUE)
    
    rescaled.k_mat <- (x - min.m)/(max.m - min.m)
    return(rescaled.k_mat)
  }
  
  # Process the density matrix
  if(scale01==FALSE){
    print(paste0(ifelse(inv==TRUE,"Inverse ",""),"density of observations returned!"))
    r2 <- r1
  }
  
  if(scale01==TRUE){
    print(paste0(ifelse(inv==TRUE,"Inverse ",""),"density of observations, between 0-1, returned!"))
    r2 <- r1 %>% app(fun=scale.f,na.rm=TRUE)
  }
  return(r2)
  
}
