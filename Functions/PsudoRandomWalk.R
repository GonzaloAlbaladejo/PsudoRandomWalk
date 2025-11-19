#
# 0                     0
#  \                   /
#   \ /----0          0
#    0      \        / \
#            0------0   0-0   Pseudo_random walk algorithm
#           / \   0  0--0
#          0   \  |\   /
#               0-0 0-0--------0
#

# a.2 Run the function----
random_walk<-function(ij,#coordinates in c(row,col) format of the dispersal starging
                      m, # matrix or surface on which the species is moving, if m is a raster it would be transformed into a matrix
                      prob_rast=NULL, # a probability raster that shows how likely are the cells to be used by the species
                      connect_rast=NULL, # A connectivity raster that represent how accessible are the cells to the species
                      cells_dist = 100, # Number of cells we allow the species to move
                      return_cell = FALSE, # Are species able to return to previous cells
                      jump_cells = NULL, # Are species allow to "jump" cells? How many cells can they skip?
                      type="range", # Type of jump to be used, for "range" (default) the value (single valule) is used to establish the maximun 
                                    # value of the jump and the range between 1-jump_cells is used for the sample. For type "probs"
                                    # jumps are sampled according to their associated probabilities. In this case Jump_cells needs to be 
                                    # a data.frame or matrix with a jumps (cells to jump) and probs (associated probabilities) columns
                      duplicated_cells = TRUE # Are species able to occupy previously used cells?
){
  
  #  Check the starting point
    if(is.na(m[ij[1],ij[2]])|is.null(m[ij[1],ij[2]])) stop("ERROR: The dispersal origin is not valid!")
  
  # a identify the surrounding destinations for the species
  init_coords <- t(matrix(ij)) ; colnames(init_coords)<-c("row","col") # Initial coordinates
  ij_x <- NULL
  try_p <- 0
  
  while(cells_dist>0){
    # a.1 Identify the dispersal origin----
    if(is.null(ij_x)){
      ij_x <- ij
    }
    
    # a.2 Is the species jumping cells?----
    if(!is.null(jump_cells)){
      # Some jumps can lead to areas of the map/matrix with no information NA, in these cases 
      # we are going to keep jumping until we fall into back into the study region
      if(type=="range") xp <- sample(1:jump_cells,1,replace=T) # Is the species jumping?
      
      if(type=="probs"){ 
        xp <- sample(jump_cells$jump,1,prob=jump_cells$probs,replace=TRUE)
        }
      
    }else{
      xp <- 1 # The species is not jumping, therefore is moving cell by cell
    } 
    
    # a.3 Define the potential destinations:
    # Check that the point falls within the limits of the matrix
      dom <- Check_domain(m=m,xp=xp,x=ij_x["col"],y=ij_x["row"])
    
    # b. Define the range of potential destinations with the values of the matrix----
      data_mx<-expand.grid(row=dom[["rows"]][1]:dom[["rows"]][2],col=dom[["cols"]][1]:dom[["cols"]][2])
    
    if(xp!=1){
      data_mx<-data_mx[data_mx$row == max(data_mx$row) | data_mx$row == min(data_mx$row) | 
                         data_mx$col == max(data_mx$col)| data_mx$col == min(data_mx$col),]
    }
    
    # b.1 Extract the values from column and row
    vals <- apply(data_mx,1,function(vi) m[vi[1],vi[2]]) %>% unlist()
    data_mx$vals <- unlist(vals)
    
    rm(vals)
    
    # If we get the probability raster alongside the study area ----
    if(!is.null(prob_rast)){
      vals2 <- apply(data_mx,1,function(vi) prob_rast[vi[1],vi[2]]) %>% unlist()
      data_mx$probs <- unlist(vals2) 
      
      rm(vals2)
      }
    
    # c. Select the destinations
    # c.1 Can the individual go back to the origin cell (equal to non-dispersal)
    if(return_cell==FALSE){
      # Remove the initial cell from the list of potential destinations
      data_mx <- data_mx[!c(data_mx$row==ij_x[1] & data_mx$col==ij_x[2]),]
    }
    
    data_mx <- data_mx[!is.na(data_mx$vals),] # remove observations outside of the study area (cells with NA's)
    
    # c.2 Are still destinations available for the individual?
    if(nrow(data_mx)<1){
      
        if(jump_cells==1){
          message("No suitable cells for the individual to jump")
          print(paste("Dead end at dispersal event",nrow(init_coords),"try a different set of parameters"))
          break()
        }
        try_p<-0
        message("No suitable cells for the individual to jump, trying other routes")
        
        while (try_p<20 & nrow(data_mx)<1){
          print(nrow(data_mx))
          xp <- sample(1:jump_cells,1) # Change the jump parameter
          
          # Check the data again
          dom <- Check_domain(m=m,xp=xp,x=ij_x["col"],y=ij_x["row"])
          
          # b. Define the range of potential destinations with the values of the matrix----
          data_mx<-expand.grid(row=dom[["rows"]][1]:dom[["rows"]][2],col=dom[["cols"]][1]:dom[["cols"]][2])
          
          if(xp!=1){
            data_mx<-data_mx[data_mx$row == max(data_mx$row) | data_mx$row == min(data_mx$row) | 
                               data_mx$col == max(data_mx$col)| data_mx$col == min(data_mx$col),]
          }
          
          # b.1 Extract the values from column and row
          vals <- apply(data_mx,1,function(vi) m[vi[1],vi[2]]) %>% unlist()
          data_mx$vals <- unlist(vals)
          
          rm(vals)
          
          # If we get the probability raster alongside the study area ----
          if(!is.null(prob_rast)){
            vals2 <- apply(data_mx,1,function(vi) prob_rast[vi[1],vi[2]]) %>% unlist()
            data_mx$probs <- unlist(vals2) 
            
            rm(vals2)
          }
          
          # c. Select the destinations
          # c.1 Can the individual go back to the origin cell (equal to non-dispersal)
          if(return_cell==FALSE){
            # Remove the initial cell from the list of potential destinations
            data_mx <- data_mx[!c(data_mx$row==ij_x[1] & data_mx$col==ij_x[2]),]
          }
          
          data_mx <- data_mx[!is.na(data_mx$vals),]
          
          # Check the try_p paramter
          try_p <- try_p + 1
        }
        
        if(nrow(data_mx<1)){
          print(paste("Dead end at dispersal event",nrow(init_coords),"try a different set of parameters"))
          return(init_coords)
          break()
        } else {
          init_coords <- rbind(init_coords,ij_x,deparse.level = 0)
          }
      
        }else{
        init_coords <- rbind(init_coords,ij_x,deparse.level = 0)
       }
    
    # Are used cells available for later interactions?
    if(duplicated_cells==FALSE){
      m[ij_x[1],ij_x[2]]<-NA
      }
    
    # Select a destination and return its coordinates
    # a. depending on the combination of raster we can have different probability scenarios
      if("probs" %in% names(data_mx)){
        ij_y <- data_mx[sample(1:nrow(data_mx),1,prob = data_mx[,"probs"]),c("row","col")] # Coordinates of the new destination
        
      }else{
        ij_y <- data_mx[sample(1:nrow(data_mx),1),c("row","col")] # Coordinates of the new destination
          }
      
    ij_x <- ij_y
      
    cells_dist <- cells_dist-1
  }
  
  return(init_coords)
  
}  

# Check cells
# Given a domain or jump and a set of origin coordinates, this funcion will return the indexes of rows and columns that are
# of interest for the sampling. It would also return a message if the parameters are incompatible with the size of the matrix/study area

Check_domain <- function(m, # Matrix
                       x, # X-cord origin
                       y, # Y-cord Origin
                       xp # domain (jump)
                       ){
                        
  # Check the rows
    if((y-xp)<1 | (y + xp)> nrow(m)){ # Are there any conflicts with the borders of the matrix?
                          
    # Is the dispersal distance greather than the matrix domain
    if((y-xp)<1 & (y + xp)> nrow(m)){
      rows_tx <- c(1,nrow(m)) 
      stop("Dispersal distance out of domain! Reconsider a lower dispersal distance")
        }else{
      # Define the destination rows
        if((y-xp)<1){
          rows_tx<-c(1,(y+xp))
        }else{
          rows_tx<- c(nrow(m)-xp,y)
            }
          }
            }else{
          rows_tx <- c(y-xp,y+xp) # If there are no problems destination rows are calculated straigthforward
              }
  # Columns
     if((x-xp)<1 | (x + xp)> ncol(m)){
        if((x-xp)<1 & (x + xp)> ncol(m)){
          columns_tx <- c(1,1) 
          stop("Dispersal distance out of domain! Reconsider a lower dispersal distance")
            }else{
          if((x-xp)<1){
              columns_tx<-c(1,(x+xp))
                }else{
              columns_tx<- c((x-xp),ncol(m))}
                  }
            }else{
              columns_tx <- c(x-xp,x+xp) 
                        }             
  
  # Return the index of valid columns and rows
  return(list(rows=rows_tx %>% unlist(),cols=columns_tx %>% unlist()))
  
  }
# End of the script