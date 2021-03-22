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


#Find name for a specie
find_art_name <- function(nr) {
  return(artsliste$NavnDansk[artsliste$ArtID == nr])
}

convert_to_longlat <- function(data_frame_input, UTMx, UTMy, projection = "+proj=utm +zone=32 +ellps=intl +units=m +no_defs +datum=WGS84") {
  #Load Libraries
  library(sp)
  
  # Create a new data frame with only the UTM data
  df <- data.frame(UTMx, UTMy)
  # turn na to 0 for calculations and save which points are na
  df_na <- df
  df[is.na(df)] <- 0
  
  
  #Make it a SP object and specify the projection
  coordinates(df) <-  ~ UTMx + UTMy
  proj4string(df) <- CRS(projection)
  
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




shanon_index_v2 <- function(cover_data, frekvens_data) {
  
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
    
    #remove all plots with 0 in cover or remove all plots with 0 in frekvens.
    #beta_data <- beta_data[!(beta_data == 0)]  
    beta_data <- beta_data[!(frekvens_data[[specie]] == 0)]
    
    
    if (length(unique( beta_data)) > 1) {
      beta_data_fitted <- fitdist(beta_data, "beta", method = "mme")
      beta_fit[nrow(beta_fit) + 1,] <- c(specie, beta_data_fitted$estimate[1], beta_data_fitted$estimate[2])
      
    }
  }
  
  # for each plot create a list with all the cover data and a total of how may obersevations.
  shanon_list <- c()
  
  for (row in 1:nrow(cover_data)) {
    total_cover_obs <- sum(cover_data[row,]) #is not currently used
    # Create an empty list for a given row
    all_obs <- cover_data[1,]
    all_obs <- 0
    
    # for a given row, find out what species is found in frekvens
    species_spotted_in_frekvens <- colnames(frekvens_data[c(frekvens_data[row,]  == 1)])
    
    #For each species spotten in frekvens, appends its posterior cover to the cover data for that row
    for (species_spotted in species_spotted_in_frekvens ) {
        
        alpha_post <- as.numeric((as.numeric(beta_fit[beta_fit$species == species_spotted,]$a) +
                                    as.numeric(cover_data[[species_spotted]][row]) ))
        beta_post <-  as.numeric(beta_fit[beta_fit$species == species_spotted,]$b) + 16 - as.numeric(cover_data[[species_spotted]][row])
        
        all_obs <- append(all_obs, (alpha_post)/(alpha_post+beta_post))
        
      
    }
    
    #Calculate the shanon index value and append it to the list after normalizing and removing zeroes
    total_cover <- sum(all_obs)
    all_obs <- all_obs[!(all_obs == 0)]
    shanon_value <- -sum(all_obs/total_cover * log((all_obs/total_cover)))
    
    shanon_list <- append(shanon_list,shanon_value)  
    
  }
  return(shanon_list)
}


site_aggregate <- function(data, aggregate_list, aggregate_function = mean) {
  # Calculate mean site position 
  mean_site_position <-  data %>%
    group_by(site) %>% 
    summarise_at(vars(longtitude, latitude), mean, na.rm = TRUE)
  
  # rename coordinates
  #mean_site_position <-  rename(mean_site_position, longtitude_mean_site = longtitude, latitude_mean_site = latitude)
  
  # Aggrgate the selected data
  aggregate_list <- c("site", "year", aggregate_list)
  
  site_data <- data[,aggregate_list ]
  site_data2 <- aggregate.data.frame(site_data, by = list(site_data$site, site_data$year), FUN = aggregate_function)[,aggregate_list]
  
  #Count how many plots there is in each site for each year
  count_data <- site_data %>% count(site, year)
  count_data <-  rename(count_data, observations = n )
  
  # Merge the aggregated data with the position of the site and the count data
  final_data <-left_join(site_data2 ,mean_site_position, by = "site")
  final_data <- left_join(final_data, count_data, by = c("site", "year"))
  return (final_data)
}



create_terhab_data <- function(terhabtype, cover =NULL, frekvens = NULL, abiotiske = NULL, artslist = NULL, traits = NULL, cover_taxa = NULL, frekvens_taxa = NULL  ) {
  
  # if no data is provided, Load the data
  if (is.null(cover)) {
    cover <- read.csv("../../Raw data/alledata-cover.csv")
  }
  if (is.null(frekvens)) {
    frekvens <- read.csv("../../Raw data/alledata-frekvens.csv")
  }
  if (is.null(abiotiske)) {
    abiotiske<- read.csv("../../Raw data/alledata-abiotiske.csv")
  }
  if (is.null(artslist)) {
    artslist <-  read.csv("../../Raw data/artsliste.csv")
  }
  if (is.null(traits)) {
    traits <- read.csv("../../Raw data/traits.csv")
  }
  if (is.null(cover_taxa)) {
    cover_taxa <- read.csv("../../Raw data/alledata-cover-samletaxa.csv")
  }
  if (is.null(frekvens_taxa)) {
    frekvens_taxa <- read.csv("../../Raw data/alledata-frekvens-samletaxa.csv")
  }
  
  
  #Merge All the data together
  
  frekvens_cover <- merge(frekvens,cover, by = c("plot","site","year"),suffixes = c("Frekvens","Cover"))
  
  frekvens_cover <- setNames(frekvens_cover,paste0(names(frekvens_cover),ifelse(names(frekvens_cover) %in% setdiff(names(frekvens),names(cover)),"frekvens","")))
  
  frekvens_cover <- setNames(frekvens_cover,paste0(names(frekvens_cover),ifelse(names(frekvens_cover) %in% setdiff(names(cover),names(frekvens)),"cover","")))
  
  df1 <- merge(abiotiske, frekvens_cover, by = c("plot","site","year"))
  df1 <- merge(df1, frekvens_taxa, by = c("plot","site","year"))
  df1 <- merge(df1, cover_taxa, by = c("plot","site","year"))
  
  
  #change terhabtype input to string and the correct format
  
  terhabtype <- paste0("{",toString(terhabtype),"}")
  
  # select the choosen terhabtype
  df2 <- df1[df1$terhabtype == terhabtype ,]
  
  # remove species that were not observed or variables with only missing value
  df2 <- df2[, colSums(df2 != 0) > 0]
  df2 <- df2[, colSums(df2 != "mv") > 0]
  
  return(df2)
}

