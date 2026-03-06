![Overview of the SDM analysis](logo/LogoRWRepository.png)
Rather than a formal package, **PsudoRandomWalk** contains a function to perform random-walk simulations using a sampling of the spatial (random) or ecological (pseudo-random) space of a given species species. Althought simple, this function provides multiple arguments to modify they way the sample or new locations "steps" is internally processed by the algorithm. For a detailled description of the functionallity of this repository follow the Rmd tutorial. 
![Overview of the SDM analysis](Figures/Run46.gif)

## Dowload the repository directly from R
```{r}
# 0. Load/install the needed packages
  list.of.packages<-c("httr","tidyverse")
  
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  lapply(list.of.packages,require,character.only=TRUE)
  # conflicts_prefer(dplyr::filter)
  rm(list.of.packages,new.packages)

# 1. Connect to the AutoMaxent GitHub repository
  git_hub <- "https://api.github.com/repos/GonzaloAlbaladejo/PsudoRandomWalk/git/trees/main?recursive=1"
  WalkRepo <- GET(git_hub) # Extract the repo information
  WalkRepo

# 2. Get the route to the functions
  file_path <- data.frame(unlist(lapply(content(WalkRepo)$tree, function(x) x$path)))
  colnames(file_path) = c('Path')
  head(file_path)

# Extract routes
  file_path <- file_path %>%
    separate(Path,c('folder','filename'),'/') %>%
    filter(folder == 'Functions') %>%
    filter(str_detect(filename,'.R'))

# 3. Configure the routes, download, and export scripts
  raw_route <- "https://raw.githubusercontent.com/GonzaloAlbaladejo/PsudoRandomWalk/main" #This is the raw route to the gitHub repository
  MyRoute <- paste(getwd(),"AutoMaxent",sep="/") ; dir.create(MyRoute,recursive=T,showWarnings = F)
  
  for(i in 1:nrow(file_path)){
    write_lines(content(GET(paste(raw_route,file_path$folder[i],file_path$filename[i],sep="/"))),
                paste(MyRoute,file_path$filename[i],sep="/"))
  }

# 4. Load the functions
  functions <- MyRoute %>% list.files(recursive = FALSE,pattern = ".R$",full.names = TRUE)
  lapply(functions,function(x) source(x))
```


 
