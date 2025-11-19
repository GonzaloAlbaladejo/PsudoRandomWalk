#
# 0              0    
#  \            /     
#   \ /----0   /      
#    0      \ /       
#            0------0 0-0 Plot dispersal results
#           / \   
#          0   \  
#               0
#
# Function to plot the results of the random walk algorithm
# a. Transform the column and row indexes into coordinates and cell indexes (for rast calculations) 
xy_fromRowCol <- function(m,row,col){
  require(terra)
  # Retrieve the cell, x and y coordinates from a set of colum and row coordinates
  p<-data.frame(cell=cellFromRowCol(m,row=row,col=col),
                xcord = xFromCol(m,col),
                ycord =yFromRow(m,row))
  
  return(p)
}

# b. Function to plot the results
#  For matrix and rast data, run the required data transformations
function(m, # Matrix or rast object over which the dispersion took place
         points, # selected cells as returned by the PsudoRandomWalk functino
         kernell=F, # Returns a density kernell of the dispersal events following the format of the "m" matrix/rast
         occurrences=F, # Returns the number of dispersion events per pixel following the format of the "m" matrix/rast
         color=NULL){
  
}







# #c. Generate the gif with the maps
#   frames <- data.frame(index=list.files(exit_route,pattern=".png") %>% gsub(pattern="([0-9]+).*$",replacement="\\1") %>% as.numeric(),
#              routes=list.files(exit_route,pattern=".png",full.names = T))
#   frames <- frames[order(frames$index),]
# 
#   frames_g <- image_read(frames)
#   image_join(frames_g) %>% image_animate(fps=20) %>% image_write(paste(".gif",sep="/"))
# # 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # ramp_col<-colorRampPalette(c("tomato3","gold","skyblue4"))
# # plot(y=init_coords[,1],x=init_coords[,2],xlim=c(1,ncol(m)),ylim=c(nrow(m),1),
# #      type="b",pch=19,col=rev(ramp_col(nrow(init_coords))))
# # 
# # # Need to retunr the cells index, the matrix coordinates, and the order in which the steps are taken
# 
# #init_coords$color<-rev(ramp_col(nrow(init_coords)))