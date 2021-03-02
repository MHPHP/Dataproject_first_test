# calculate a naive shanon index using cover
simple_shanon_index <- function(data_frame_input){
  data_frame_input$simple_shanon_index <- rowSums(-data_frame_input[, grep("Cover", names(data_frame_input))]/
                                                    (rowSums(data_frame_input[,grep("Cover", names(data_frame_input))])) *
                                                    log(data_frame_input[, grep("Cover", names(data_frame_input))]/
                                                          (rowSums(data_frame_input[,grep("Cover", names(data_frame_input))]))), na.rm = TRUE  )
  return(data_frame_input)
}


# find number for a species
find_art_nr <- function(name, latin = FALSE) {
  if (latin == TRUE) {
    artsliste <- read.csv("../Data/Raw data/artsliste.csv")
    return(artsliste[grep(name, artsliste$LatArt),])
  }
  else {
    return(artsliste[grep(name, artsliste$NavnDansk),])
  }
}


convert_to_longlat <- function(data_frame_input, UTMx, UTMy) {
  #Load Libraries
  library(sp)
  
  # Create a new data frame with only the UTM data
  df <- data.frame(UTMx, UTMy)
  # turn na to 0 for calculations and save which points are na
  df_na <- df
  df[is.na(df)] <- 0
  
  
  #Make it a SP object and specify the projection
  coordinates(df) <-  ~ UTMx + UTMy
  proj4string(df) <- CRS("+proj=utm +zone=32 +ellps=intl +units=m +no_defs +datum=WGS84")
  
  #Transform the data to long lat
  df1 <- spTransform(df, CRS("+init=epsg:4326"))
  
  # Write the long lat
  data_frame_input$latitude <-  df1@coords[,2]
  data_frame_input$longtitude <-  df1@coords[,1]
  
  # Turn na back to na
  data_frame_input$latitude[is.na(df_na[1]) ] <- NA 
  data_frame_input$longtitude[is.na(df_na[2]) ] <- NA 
  
  return(data_frame_input)
}





# Make animation
animation_map <- function(data_frame_in, longtitude, latitude, index, time, duration_time = 15) {
  library(ggmap)
  library(ggplot2)
  library(gifski)
  library(png)
  
  #Load map of Denmark
  world.map <- map_data ("world")
  map1 <- world.map %>% filter(region == "Denmark")
  
  # change input values to strings, so it is easier to use.
  index <- deparse(substitute(index))
  time <- deparse(substitute(time))
  latitude <- deparse(substitute(latitude))
  longtitude <- deparse(substitute(longtitude))
  
  # plot the map
  map_ani <-  ggplot(map1)+ 
    geom_map(map = map1, aes(long, lat, map_id = region), fill="white", colour = "black") +
    coord_map()
  
  
  # Add points
  map_ani <-  map_ani +
    geom_point(data = data_frame_in, aes(y= .data[[latitude]], x = .data[[longtitude]], colour = .data[[index]]), alpha = .9, size =5) + 
    
    scale_color_gradientn(colours = rainbow(5)) +
    transition_states(factor(data_frame_in[[time]]), transition_length = 1, state_length = 1) +
    labs(title = "Year: {next_state}")
  
  #create animation
  animate(map_ani,  nframes = 2*length(unique(data_frame_in[[time]])), duration = duration_time)
  
}






# input is cover data and frekvens data as two data frames, where the rows match for the same observation, and the names of the species match.
shanon_index_v1 <- function(cover_data, frekvens_data) {
  
  #Load functions 
  library(tidyverse)
  library(fitdistrplus)
  
  #create data frame to hold the fitted values for each species
  beta_fit <- data.frame(matrix(ncol = 3, nrow = 0))
  x <- c("species","a", "b")
  colnames(beta_fit) <- x
  
  # for Each species calculate the shape parameter for the fitted beta distribution and save them in a data frame.
  for (specie in colnames(cover_data)) {
    beta_data <- cover_data[,specie]/16
    beta_data <- beta_data[!(beta_data == 0)]
    
    beta_data <-  ifelse(beta_data == 1, 0.9999, beta_data)
    
    if (length(unique( beta_data)) > 1) {
      beta_data_fitted <- fitdist(beta_data, "beta")
      beta_fit[nrow(beta_fit) + 1,] <- c(specie, beta_data_fitted$estimate[1], beta_data_fitted$estimate[2])
      
    }
  }
  
  # for each plot create a list with all the cover data and a total of how may obersevations.
  shanon_list <- c()
  
  for (row in 1:nrow(cover_data)) {
    total_cover_obs <- sum(cover_data[row,]) #is not currently used
    # Create a list with all the Cover data greater than 0 for a given row
    all_obs <- cover_data[row,] 
    all_obs <- all_obs[all_obs > 0]/16
    
    # for a given row, find out what species is found in frekvens
    species_spotted_in_frekvens <- colnames(frekvens_data[c(frekvens_data[row,]  == 1)])
    
    #For each species spotten in frekvens but not in cover, appends its posterior cover to the cover data for that row
    for (species_spotted in species_spotted_in_frekvens ) {
      if (length(cover_data[[species_spotted]][row]) == 1 && cover_data[[species_spotted]][row] == 0) {
        
        all_obs <- append(all_obs, (as.numeric(beta_fit[beta_fit$species == species_spotted,]$a) + 1 -1)/
                            (as.numeric(beta_fit[beta_fit$species == species_spotted,]$a) +
                               as.numeric(beta_fit[beta_fit$species == species_spotted,]$b) + 17 -2))
        
        total_cover_obs <- total_cover_obs +1
      }
      
      
    }
    
    #Calculate the shanon index value and append it to the list
    total_cover <- sum(all_obs)
    shanon_value <- -sum(all_obs/total_cover * log((all_obs/total_cover)))
    
    shanon_list <- append(shanon_list,shanon_value)  
    
  }
  return(shanon_list)
}

