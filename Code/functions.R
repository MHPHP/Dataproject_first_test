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
    artslist <- read.csv("../../Raw data/artsliste.csv")
    return(artslist[grep(name, artslist$LatArt),])
  }
  else {
    return(artslist[grep(name, artslist$NavnDansk),])
  }
}


convert_to_longlat <- function(data_frame_input, UTMx, UTMy) {
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