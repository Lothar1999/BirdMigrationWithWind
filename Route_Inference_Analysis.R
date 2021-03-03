library(rworldmap); newmap <-getMap(resolution = "low")
library(maptools); data(wrld_simpl)
library(psych)
library(rgdal)
library(sp)


ID_list_WM<- list.files("C:/Users/Mike_Werfeli/Documents/Arbeit/VoWa/Results/Wind")
ID_list_NW<- list.files("C:/Users/Mike_Werfeli/Documents/Arbeit/VoWa/Results/NoWind")

if(any((ID_list_WM == ID_list_NW) == F) == F){
  ID_list<- ID_list_WM
  rm(ID_list_NW)
  rm(ID_list_WM)
} else{
  print("File paths do not contain the same folders")
}

read_dta<- function(ID_list){
  #' This function reads the relevant data of the different inference summary files. It reads as well the WM data as the NW data
  #' Please check the description document that comes with this skript to check for the storage structure this function requires
  #' @param ID_list list which contains the name of the last folder where the summary file is stored
  #' @return list of parameters read from the different files. general structure is: Mean track, mean track sd, mean track cv and hm, where everywhere
  #' the WM data is first and the NW data second.
 
  #mean track positions
  mean_tracks_wind<- lapply(ID_list, FUN= function(x){
    df<- read.csv(paste0("C:/Users/Mike_Werfeli/Documents/Arbeit/VoWa/Results/Wind/",x, "/",x,"_sumary_WM.csv"), sep = ',')
    lat_m<- df$Lat.50.
    lon_m<- df$Lon.50.
    
    mean_coords<- cbind(lon_m,lat_m)
  })
  
  mean_tracks_Nowind<- lapply(ID_list, FUN= function(x){
    df<- read.csv(paste0("C:/Users/Mike_Werfeli/Documents/Arbeit/VoWa/Results/NoWind/",x, "/",x,"_sumary_NW.csv"), sep = ',')
    lat_m<- df$Lat.50.
    lon_m<- df$Lon.50.
    
    mean_coords<- cbind(lon_m,lat_m)
  })
  
  #mean track positions standard deviation
  mean_tracks_sd_wind<- lapply(ID_list, FUN= function(x){
    df<- read.csv(paste0("C:/Users/Mike_Werfeli/Documents/Arbeit/VoWa/Results/Wind/",x, "/",x,"_sumary_WM.csv"), sep = ',')
    lat_sd<- df$Lat.sd
    lon_sd<- df$Lon.sd
    
    mean_coords<- cbind(lon_sd,lat_sd)
  })
  
  mean_tracks_sd_Nowind<- lapply(ID_list, FUN= function(x){
    df<- read.csv(paste0("C:/Users/Mike_Werfeli/Documents/Arbeit/VoWa/Results/NoWind/",x, "/",x,"_sumary_NW.csv"), sep = ',')
    lat_sd<- df$Lat.sd
    lon_sd<- df$Lon.sd
    
    mean_coords<- cbind(lon_sd,lat_sd)
  })
  
  #Mean track positions convidence intervall
  CV_Vals_wind<- lapply(ID_list, FUN= function(x){
    df<- read.csv(paste0("C:/Users/Mike_Werfeli/Documents/Arbeit/VoWa/Results/Wind/",x, "/",x,"_sumary_WM.csv"), sep = ',')
    latCV_lt<- df$Lat.2.5.
    latCV_ut<- df$Lat.97.5.
    lonCV_lt<- df$Lon.2.5.
    lonCV_ut<- df$Lon.97.5.
    
    CV_vals<- cbind(latCV_lt, latCV_ut,lonCV_lt,lonCV_ut)
  })
  
  CV_Vals_Nowind<- lapply(ID_list, FUN= function(x){
    df<- read.csv(paste0("C:/Users/Mike_Werfeli/Documents/Arbeit/VoWa/Results/NoWind/",x, "/",x,"_sumary_NW.csv"), sep = ',')
    latCV_lt<- df$Lat.2.5.
    latCV_ut<- df$Lat.97.5.
    lonCV_lt<- df$Lon.2.5.
    lonCV_ut<- df$Lon.97.5.
    
    CV_vals<- cbind(latCV_lt, latCV_ut,lonCV_lt,lonCV_ut)
  })
  
  #Harmonic mean values
  hm_wind<- lapply(ID_list, FUN= function(x){
    df<- read.csv(paste0("C:/Users/Mike_Werfeli/Documents/Arbeit/VoWa/Results/Wind/",x, "/",x,"_sumary_WM.csv"), sep = ',')
    hm_p<- df$hm_positions[1]
    hm_b<- df$hm_behaviour[1]
    hm_i<- df$hm_inference[1]
    hm<- c(hm_p, hm_b, hm_i)
  })
  
  hm_Nowind<- lapply(ID_list, FUN= function(x){
    df<- read.csv(paste0("C:/Users/Mike_Werfeli/Documents/Arbeit/VoWa/Results/Nowind/",x, "/",x,"_sumary_NW.csv"), sep = ',')
    hm_p<- df$hm_positions[1]
    hm_b<- df$hm_behaviour[1]
    hm_i<- df$hm_inference[1]
    hm<- c(hm_p, hm_b, hm_i)
    
  })
  
  spd_wind<- lapply(ID_list, FUN=function(x){
    df<- read.csv(paste0("C:/Users/Mike_Werfeli/Documents/Arbeit/VoWa/Results/Wind/",x, "/",x,"_sumary_WM.csv"), sep = ',')
    spd<- df$spd_z0
    })
  
  spd_nowind<- lapply(ID_list, FUN=function(x){
    df<- read.csv(paste0("C:/Users/Mike_Werfeli/Documents/Arbeit/VoWa/Results/Nowind/",x, "/",x,"_sumary_NW.csv"), sep = ',')
    spd<- df$spd_z0
  })
  
  #Bind data before it is returned
  dta<-c(mean_tracks_wind, mean_tracks_Nowind, mean_tracks_sd_wind, mean_tracks_sd_Nowind, CV_Vals_wind, CV_Vals_Nowind, hm_wind, hm_Nowind,spd_wind, spd_nowind)
  
  return(dta)
  
}

get_idx<- function(ID_list){
  #' function to get the idx needed to separate the BI_dat list elements
  #' @param  ID_list list of names found in the result folders
  #' @return idx where the last element is found belonging to one of the results extracted form the different files
  idx<- length(ID_list)
  return(idx)
}

get_path_dta<- function(idx, BI_dta){
  #' getter function to extract the bird migration routes from BI_dta
  #' @param idx parameter to indicate the number of routes to be analyzed. It also indicates the index position of the last 
  #' element to be extracted
  #' @param BI_dta the data file which contains all relevant measures read from the respective summary files
  #' @return All lat and lon coordinates of the respective migration route inferences to be analyzed (WM and NW)
  path_dta_WM<- BI_dta[1:idx]
  path_dta_NW<- BI_dta[(idx+1):(2*idx)]
  return(c(path_dta_WM,path_dta_NW))
}

get_sd_dta<- function(idx,BI_dta){
  #' getter function to extract the standard deviations of the positions found along the migration route
  #' @param idx parameter to indicate the number of routes to be analyzed. It also indicates the index position of the last 
  #' element to be extracted
  #' @param BI_dta the data file which contains all relevant measures read from the respective summary files
  #' @return All standard deviation for lat and lon respecitvely (for WM and NW)
  sd_WM<- BI_dta[((2*idx)+1):(3*idx)]
  sd_NW<- BI_dta[((3*idx)+1):(4*idx)]
  return(c(sd_WM, sd_NW))
}

get_cv_dta<- function(idx, BI_dta){
  #' getter function to extract the confidence intervals of the positions found along the migration route
  #' @param idx parameter to indicate the number of routes to be analyzed. It also indicates the index position of the last 
  #' element to be extracted
  #' @param BI_dta the data file which contains all relevant measures read from the respective summary files
  #' @return All confidence intervals for lat and lon respecitvely (for WM and NW)
  cv_WM<- BI_dta[((4*idx)+1):(5*idx)]
  cv_NW<- BI_dta[((5*idx)+1):(6*idx)]
  return(c(cv_WM,cv_NW))
}

get_hm_dta<- function(idx, BI_dta){
  #' getter function to extract the harmonic means of the positions found along the migration route, the bird behaviour between those positions
  #' and the entire route inference.
  #' @param idx parameter to indicate the number of routes to be analyzed. It also indicates the index position of the last 
  #' element to be extracted
  #' @param BI_dta the data file which contains all relevant measures read from the respective summary files
  #' @return All harmonic means for the positions, the bird behaviour and the entire inference for each route inference respecitvely (for WM and NW)
  hm_wm<- BI_dta[((6*idx)+1):(7*idx)]
  hm_nw<- BI_dta[((7*idx)+1):(8*idx)]
  return(c(hm_wm,hm_nw))
}

get_spd_dta<- function(idx, BI_dta){
  #' getter function to extract the speed values the bird did probably fly given the inferred positions and the activity data
  #' @param idx parameter to indicate the number of routes to be analyzed.
  #' @param BI_dta the data file which contains all relevant measures read from the respcetive summary files
  #' @return all speed values on each segment for each route analyzed respectively.
  spd_wm<- BI_dta[((8*idx)+1):(9*idx)]
  spd_nw<- BI_dta[((9*idx)+1):(10*idx)]
  return(c(spd_wm,spd_nw))
}

trajectory_separator<- function(ID_list){
  #' function which sepparates the migration routes into three stages. First, the fall route where the bird fly from the breeding 
  #' grounds to the non-breeding locations of residence. Second the movement at the non-breeding locations of residence and finally,
  #' the spring route where the birds migrate from the non-breeding locations of residence back to the breeding grounds.
  #' @param ID_list list of available inference results
  #' @return list containing the respective route segments for each route inference available in the given folder structure
  idx<- get_idx(ID_list)
  BI_dta<- read_dta(ID_list) 
  paths<- get_path_dta(idx = idx, BI_dta = BI_dta)
  WM<- paths[1:idx]
  NW<- paths[(idx+1):(2*idx)]
  
  #Extract coordinate pairs of non-breeding location
  coord_pairs_non_B_loc_Wind<- lapply(WM, FUN=function(x){
    #Find non-breeding locations
    #Find most southerly position and find points up to 3 degree north of it
    #Compare indexes of possible positions --> find all longitudes which match latitudes
    min_lats<- which(x[,2] < min(x[,2])+3)
    min_lons<- which(x[,1] < x[,2][which(x[,2] == min(x[,2]))]+3)
    
    matching_c_lon<- match(min_lats,min_lons)
    matching_c_lat<- match(min_lons,min_lats)
    
    matching_c_lon<-matching_c_lon[!is.na(matching_c_lon)]
    matching_c_lat<-matching_c_lat[!is.na(matching_c_lat)]
    
    pairs<-cbind(x[,2][min_lats[matching_c_lat]],x[,1][min_lons[matching_c_lon]])
    
    
    
  })
  #Extract coordinate pairs of non-breeding location
  coord_pairs_non_B_loc_NoWind<- lapply(NW, FUN=function(x){
    #Find non-breeding locations
    #Find most southerly position and find points up to 3 degree north of it
    #Compare indexes of possible positions --> find all longitudes which match latitudes
    min_lats<- which(x[,2] < min(x[,2])+3)
    min_lons<- which(x[,1] < x[,2][which(x[,2] == min(x[,2]))]+3)
    
    matching_c_lon<- match(min_lats,min_lons)
    matching_c_lat<- match(min_lons,min_lats)
    
    matching_c_lon<-matching_c_lon[!is.na(matching_c_lon)]
    matching_c_lat<-matching_c_lat[!is.na(matching_c_lat)]
    
    pairs<-cbind(x[,2][min_lats[matching_c_lat]],x[,1][min_lons[matching_c_lon]])
    
    
  })
  
  #start and end idx of positions at non - breeding locations
  idx<- c(1:length(coord_pairs_non_B_loc_Wind))
  index_nb_wind<- lapply(idx, FUN=function(x){
    start_pos<- cbind(which(WM[[x]][,2]==coord_pairs_non_B_loc_Wind[[x]][1,1]),which(WM[[x]][,1]==coord_pairs_non_B_loc_Wind[[x]][1,2]))
    end_pos<- cbind(which(WM[[x]][,2]==coord_pairs_non_B_loc_Wind[[x]][length(coord_pairs_non_B_loc_Wind[[x]][,1]),1]),which(WM[[x]][,1]==coord_pairs_non_B_loc_Wind[[x]][length(coord_pairs_non_B_loc_Wind[[x]][,1]),2]))
    saePos<- rbind(start_pos, end_pos)
    
  })
  
  idx<- c(1:length(coord_pairs_non_B_loc_NoWind))
  index_nb_NoWind<- lapply(idx, FUN=function(x){
    start_pos<- cbind(which(NW[[x]][,2]==coord_pairs_non_B_loc_NoWind[[x]][1,1]),which(NW[[x]][,1]==coord_pairs_non_B_loc_NoWind[[x]][1,2]))
    end_pos<- cbind(which(NW[[x]][,2]==coord_pairs_non_B_loc_NoWind[[x]][length(coord_pairs_non_B_loc_NoWind[[x]][,1]),1]),which(NW[[x]][,1]==coord_pairs_non_B_loc_NoWind[[x]][length(coord_pairs_non_B_loc_NoWind[[x]][,1]),2]))
    saePos<- rbind(start_pos, end_pos)
    
  })
  
  #Fall routes
  #extract all positions from start to the first index of the non breeding location
  idx<- seq(1:length(coord_pairs_non_B_loc_Wind))
  fr_idx_Wind<- lapply(idx, FUN=function(x){
    fallR<- WM[[x]][1:index_nb_wind[[x]][1],]
  })
  
  idx<- seq(1:length(coord_pairs_non_B_loc_NoWind))
  fr_idx_NoWind<- lapply(idx, FUN=function(x){
    fallR<- NW[[x]][1:index_nb_NoWind[[x]][1],]
  })
  
  #Spring routes
  idx<- seq(1:length(coord_pairs_non_B_loc_Wind))
  sr_idx_Wind<- lapply(idx, FUN=function(x){
    springR<- WM[[x]][index_nb_wind[[x]][2]:length(WM[[x]][,1]),]
  })
  
  idx<- seq(1:length(coord_pairs_non_B_loc_NoWind))
  sr_idx_NoWind<- lapply(idx, FUN=function(x){
    springR<- NW[[x]][index_nb_NoWind[[x]][2]:length(NW[[x]][,1]),]
  })
  
  wind_r<- c(fr_idx_Wind, coord_pairs_non_B_loc_Wind, sr_idx_Wind)
  
  nowind_r<- c(fr_idx_NoWind, coord_pairs_non_B_loc_NoWind, sr_idx_NoWind)
  
  r_dta<- list(wind_r= wind_r,nowind_r= nowind_r)
  
  return(r_dta)
}

routes_plotter<- function(ID_list,results_path, allRoutesFilePath){
  #' Function which plots the different migration route results on (1) a plot of the route indicating fall and spring routes with error bars,
  #' (2) fall and spring routes on separate sub plots 
  #' @param ID_list list of all data folders available
  #' @return plots saved in file structure 

  list_idx<- seq(1:length(ID_list))
  WGS84 = CRS("+init=epsg:4326")
  birdMap<- spTransform(wrld_simpl, CRSobj = WGS84)
  xText<- "Longitude [°]"
  yText<- "Latitude [°]"
  
  idx<- get_idx(ID_list)
  BI_dta<- read_dta(ID_list) 
  sd_dta<- get_sd_dta(idx, BI_dta)
  sd_WM<- sd_dta[1:idx]
  sd_NW<- sd_dta[(idx+1):(2*idx)]
  paths<- get_path_dta(idx = idx, BI_dta = BI_dta)
  WM<- paths[1:idx]
  NW<- paths[(idx+1):(2*idx)]
  rm(BI_dta)
  rm(paths)
  
  route_segm_lists<- trajectory_separator(ID_list)
  i<- c(1:idx)
  fr_idx_Wind<- lapply(i, FUN = function(x){
    route_segm_lists$wind_r[[x]]
    })
  fr_idx_NoWind<- lapply(i, FUN = function(x){
    route_segm_lists$nowind_r[[x]]
  })
  
  i<- c((idx+1) : (2*idx))
  
  coord_pairs_non_B_loc_Wind<- lapply(i, FUN = function(x){
    route_segm_lists$wind_r[[x]]
  })
  
  coord_pairs_non_B_loc_NoWind<- lapply(i, FUN = function(x){
      route_segm_lists$nowind_r[[x]]
    })
  
  i<- c((2*idx+1) : (3*idx))
  sr_idx_Wind<- lapply(i, FUN = function(x){
    route_segm_lists$wind_r[[x]]
  })
  
  sr_idx_NoWind<- lapply(i, FUN = function(x){
    route_segm_lists$nowind_r[[x]]
    })
  
  lapply(list_idx, FUN = function(x){
    result_storage_path<- paste0(results_path, ID_list[x],"_Routes.png")
    
    png(file=result_storage_path, width = 650,height = 650,units = "px", pointsize = 12)
    
    plot(birdMap, xlim=c(-20,20), ylim=c(0,50), xlab=xText, ylab=yText, main = paste0("Migration routes ", ID_list[[x]]))
    
    #Error bars
    arrows(WM[[x]][,1]+sd_WM[[x]][,1],WM[[x]][,2], WM[[x]][,1]-sd_WM[[x]][,1], WM[[x]][,2], length = 0, lwd = 2, col = "gray50")
    arrows(WM[[x]][,1],WM[[x]][,2]+sd_WM[[x]][,2], WM[[x]][,1], WM[[x]][,2]-sd_WM[[x]][,2],  length = 0,lwd = 2, col = "gray50")
    
    
    arrows(NW[[x]][,1]+sd_NW[[x]][,1],NW[[x]][,2], NW[[x]][,1]-sd_NW[[x]][,1], NW[[x]][,2], length = 0, col = "gray")
    arrows(NW[[x]][,1], NW[[x]][,2]+sd_NW[[x]][,2], NW[[x]][,1], NW[[x]][,2]-sd_NW[[x]][,2],length = 0, col = "gray")
    
    #Fall route
    points(fr_idx_Wind[[x]][,1],fr_idx_Wind[[x]][,2], pch=20, col="darkorange")
    points(fr_idx_NoWind[[x]][,1],fr_idx_NoWind[[x]][,2],pch=20, col="darkorange")
    lines(fr_idx_Wind[[x]][,1],fr_idx_Wind[[x]][,2], lwd=2, col="darkorange")
    lines(fr_idx_NoWind[[x]][,1],fr_idx_NoWind[[x]][,2], col="darkorange")
    
    #Non breeding location of residence
    points(coord_pairs_non_B_loc_Wind[[x]][,2],coord_pairs_non_B_loc_Wind[[x]][,1], pch=20)
    points(coord_pairs_non_B_loc_NoWind[[x]][,2],coord_pairs_non_B_loc_NoWind[[x]][,1], pch=20)
    lines(coord_pairs_non_B_loc_Wind[[x]][,2],coord_pairs_non_B_loc_Wind[[x]][,1], lwd=2)
    lines(coord_pairs_non_B_loc_NoWind[[x]][,2],coord_pairs_non_B_loc_NoWind[[x]][,1])
    
    #Spring route
    points(sr_idx_Wind[[x]][,1],sr_idx_Wind[[x]][,2], pch=20, col="blue")
    points(sr_idx_NoWind[[x]][,1],sr_idx_NoWind[[x]][,2], pch=20, col="blue")
    lines(sr_idx_Wind[[x]][,1],sr_idx_Wind[[x]][,2], lwd=2, col="blue")
    lines(sr_idx_NoWind[[x]][,1],sr_idx_NoWind[[x]][,2], col="blue")
    
    #Add legend
    legend(15, 15, legend = c("Fall", "Fall no wind", "Spring", "Spring no wind"),lty = 1, lwd = 2:1, col = c("darkorange","darkorange","blue", "blue"), cex=0.8, title = "Migration route segments")
    axis(side = 1, at = seq(-60,60, 5))
    axis(side = 2)
    box()
    dev.off()
    
  })
 
  
  
  png(filename = allRoutesFilePath,width = 650,height = 650,units = "px", pointsize = 12)
  par(mfrow=c(2,2))
  
  col_pal<- c("green", "blue","orange")
  
  #Plot spring routes with wind on subplot top left
  WGS84 = CRS("+init=epsg:4326")
  birdMap<- spTransform(wrld_simpl, CRSobj = WGS84)
  xText<- "Longitude [°]"
  yText<- "Latitude [°]"
  srBirdPts<- function(x){
    birdposdf<-data.frame(sr_idx_Wind[[x]][,1],sr_idx_Wind[[x]][,2])
    birdpos<-SpatialPointsDataFrame(birdposdf,data = birdposdf, proj4string = WGS84)
  }
  
  srBirdLns<- function(x){
    birdLines<- Lines(list(Line(cbind(sr_idx_Wind[[x]][,1],sr_idx_Wind[[x]][,2]))),ID=as.character(1))
    birdSpatialLines<- SpatialLines(list(birdLines), proj4string =WGS84)
  }
  
  plot(birdMap, xlim=c(-20,20), ylim=c(0,50), xlab=xText, ylab=yText, main = "Spring routes wind")
  lapply(list_idx, FUN = function(x){
    birdpos<- srBirdPts(x)
    birdLines<- srBirdLns(x)
    points(birdpos, col=col_pal[x],pch=20)
    lines(birdLines, col=col_pal[x],lwd=2)
  })
  legend(20, 15, legend = ID_list[list_idx],lty = 1, lwd = 2:1, col = col_pal, cex=0.8, title = "Spring routes")
  axis(side = 1, at = seq(-60,60, 5))
  axis(side = 2)
  box()
  
  #Plot fall routes with wind on subplot top right
  frBirdPts<- function(x){
    birdposdf<-data.frame(fr_idx_Wind[[x]][,1] ,fr_idx_Wind[[x]][,2])
    birdpos<-SpatialPointsDataFrame(birdposdf,data = birdposdf, proj4string = WGS84)
  }
  
  frBirdLns<- function(x){
    birdLines<- Lines(list(Line(cbind(fr_idx_Wind[[x]][,1] ,fr_idx_Wind[[x]][,2]))),ID=as.character(1))
    birdSpatialLines<- SpatialLines(list(birdLines), proj4string =WGS84)
  }
  
  
  plot(birdMap, xlim=c(-20,20), ylim=c(0,50), xlab=xText, ylab=yText, main = "Fall routes wind")
  lapply(list_idx, FUN = function(x){
    birdpos<- frBirdPts(x)
    birdLines<- frBirdLns(x)
    points(birdpos, col=col_pal[x],pch=20)
    lines(birdLines, col=col_pal[x],lwd=2)
  })
  legend(20, 15, legend = ID_list[list_idx],lty = 1, lwd = 2:1, col = col_pal, cex=0.8, title = "Fall routes")
  axis(side = 1, at = seq(-60,60, 5))
  axis(side = 2)
  box()
  
  #Plot spring routes without wind on subplot bottom left
  srBirdPts<- function(x){
    birdposdf<-data.frame(sr_idx_NoWind[[x]][,1],sr_idx_NoWind[[x]][,2])
    birdpos<-SpatialPointsDataFrame(birdposdf,data = birdposdf, proj4string = WGS84)
  }
  
  srBirdLns<- function(x){
    birdLines<- Lines(list(Line(cbind(sr_idx_NoWind[[x]][,1],sr_idx_NoWind[[x]][,2]))),ID=as.character(1))
    birdSpatialLines<- SpatialLines(list(birdLines), proj4string =WGS84)
  }
  plot(birdMap, xlim=c(-20,20), ylim=c(0,50), xlab=xText, ylab=yText, main = "Spring routes no wind")
  lapply(list_idx, FUN = function(x){
    birdpos<- srBirdPts(x)
    birdLines<- srBirdLns(x)
    points(birdpos, col=col_pal[x],pch=20)
    lines(birdLines, col=col_pal[x],lwd=2)
  })
  legend(20, 15, legend = ID_list[list_idx],lty = 1, lwd = 2:1, col = col_pal, cex=0.8, title = "Spring routes")
  axis(side = 1, at = seq(-60,60, 5))
  axis(side = 2)
  box()
  
  #Plot fall routes without wind on subplot bottom right
  frBirdPts<- function(x){
    birdposdf<-data.frame(fr_idx_NoWind[[x]][,1],fr_idx_NoWind[[x]][,2])
    birdpos<-SpatialPointsDataFrame(birdposdf,data = birdposdf, proj4string = WGS84)
  }
  
  frBirdLns<- function(x){
    birdLines<- Lines(list(Line(cbind(fr_idx_NoWind[[x]][,1],fr_idx_NoWind[[x]][,2]))),ID=as.character(1))
    birdSpatialLines<- SpatialLines(list(birdLines), proj4string =WGS84)
  }
  
  
  plot(birdMap, xlim=c(-20,20), ylim=c(0,50), xlab=xText, ylab=yText, main = "Fall routes no wind")
  lapply(list_idx, FUN = function(x){
    birdpos<- frBirdPts(x)
    birdLines<- frBirdLns(x)
    points(birdpos, col=col_pal[x],pch=20)
    lines(birdLines, col=col_pal[x],lwd=2)
  })
  legend(20, 15, legend = ID_list[list_idx],lty = 1, lwd = 2:1, col = col_pal, cex=0.8, title = "Fall routes")
  axis(side = 1, at = seq(-60,60, 5))
  axis(side = 2)
  box()
  
  dev.off()
  
  
}

route_analyser<- function(ID_list,writeVar=F,results_path){
  #' With the function route_analyser it is possible to extract all relevant metrics of the inferred routes, which can be used to compare
  #' the different models applied.
  #' The metrics computed in this function are the standard deviaton, confidence intervall, harmonic mean and min, max, mean and median on speed data
  #' @param ID_list list of folders to work trough
  #' @retrun list of comparison metrics
  idx<- get_idx(ID_list)
  BI_dta<- read_dta(ID_list)  
  sd_dta<- get_sd_dta(idx, BI_dta)
  cv_dta<- get_cv_dta(idx,BI_dta)
  hm_dta<- get_hm_dta(idx, BI_dta)
  spd_dta<- get_spd_dta(idx,BI_dta)
  
  #get CV values
  CV_Vals_wind<- cv_dta[1:idx]
  CV_Vals_Nowind<- cv_dta[(idx+1):(2*idx)]  
  
  rm(BI_dta)
  rm(cv_dta)
  
  cvRange_lat_wind<-lapply(c(1:idx), FUN = function(x){
    CV_Vals_wind[[x]][,2] - CV_Vals_wind[[x]][,1]
  })
  
  cvRange_lon_wind<-lapply(c(1:idx), FUN = function(x){
    CV_Vals_wind[[x]][,4] - CV_Vals_wind[[x]][,3]
  })
  
  cvRange_lat_Nowind<-lapply(c(1:idx), FUN = function(x){
    CV_Vals_Nowind[[x]][,2] - CV_Vals_Nowind[[x]][,1]
  })
  
  cvRange_lon_Nowind<-lapply(c(1:idx), FUN = function(x){
    CV_Vals_Nowind[[x]][,4] - CV_Vals_Nowind[[x]][,3]
  })
  
  #min, max, mean cv
  cvDataNowind<- lapply(c(1:idx), FUN = function(x){
    mincvlat<-min(cvRange_lat_Nowind[[x]])
    maxcvlat<-max(cvRange_lat_Nowind[[x]])
    meancvlat<-mean(cvRange_lat_Nowind[[x]])
    meadiancvlat<-median(cvRange_lat_Nowind[[x]])
    cv_metrics_lat<- cbind(mincvlat,maxcvlat,meancvlat,meadiancvlat)
    
    mincvlon<-min(cvRange_lon_Nowind[[x]])
    maxcvlon<-max(cvRange_lon_Nowind[[x]])
    meancvlon<-mean(cvRange_lon_Nowind[[x]])
    mediancvlon<-median(cvRange_lon_Nowind[[x]])
    cv_metrics_lon<- cbind(mincvlon,maxcvlon,meancvlon,mediancvlon)
    cv_metrics<-rbind(cv_metrics_lat,cv_metrics_lon)
  })
  
  cvDatawind<- lapply(c(1:idx), FUN = function(x){
    mincvlat<-min(cvRange_lat_wind[[x]])
    maxcvlat<-max(cvRange_lat_wind[[x]])
    meancvlat<-mean(cvRange_lat_wind[[x]])
    meadiancvlat<-median(cvRange_lat_wind[[x]])
    cv_metrics_lat<- cbind(mincvlat,maxcvlat,meancvlat,meadiancvlat)
    
    mincvlon<-min(cvRange_lon_wind[[x]])
    maxcvlon<-max(cvRange_lon_wind[[x]])
    meancvlon<-mean(cvRange_lon_wind[[x]])
    mediancvlon<-median(cvRange_lon_wind[[x]])
    cv_metrics_lon<- cbind(mincvlon,maxcvlon,meancvlon,mediancvlon)
    cv_metrics<-rbind(cv_metrics_lat,cv_metrics_lon)
  })
  
  #sd 
  sd_vals_wind<- sd_dta[1:idx]
  sd_vals_nowind<- sd_dta[(idx+1):(2*idx)]
  rm(sd_dta)
  
  sd_vals_wind<- lapply(c(1:idx), FUN=function(x){
    sd_min<- min(sd_vals_wind[[x]])
    sd_max<- max(sd_vals_wind[[x]])
    sd_mean<- mean(sd_vals_wind[[x]])
    sd_median<- median(sd_vals_wind[[x]])
    sd_dta_wind<- cbind(sd_min, sd_max, sd_mean, sd_median)
    
  })
  
  sd_vals_nowind<- lapply(c(1:idx), FUN=function(x){
    sd_min<- min(sd_vals_nowind[[x]])
    sd_max<- max(sd_vals_nowind[[x]])
    sd_mean<- mean(sd_vals_nowind[[x]])
    sd_median<- median(sd_vals_nowind[[x]])
    sd_dta_nowind<- cbind(sd_min, sd_max, sd_mean, sd_median)
  })
  
  #grounds
  spd_vals_wind<- spd_dta[1:idx]
  spd_vals_nowind<- spd_dta[(idx+1):(2*idx)]
  rm(spd_dta)
  
  spd_metrics_wind<- lapply(c(1:idx), FUN=function(x){
    spd_min<- min(spd_vals_wind[[x]])
    spd_max<- max(spd_vals_wind[[x]])
    spd_mean<- mean(spd_vals_wind[[x]])
    spd_median<- median(spd_vals_wind[[x]])
    spd_dta_wind<- cbind(spd_min, spd_max, spd_mean, spd_median)
  })
  
  spd_metrics_nowind<- lapply(c(1:idx), FUN=function(x){
    spd_min<- min(spd_vals_nowind[[1]])
    spd_max<- max(spd_vals_nowind[[x]])
    spd_mean<- mean(spd_vals_nowind[[x]])
    spd_median<- median(spd_vals_nowind[[x]])
    spd_dta_wind<- cbind(spd_min, spd_max, spd_mean, spd_median)
  })
  
  spd_comp<- lapply(c(1:idx), FUN=function(x){
    spd_vals_wind[[x]] - spd_vals_nowind[[x]]
  })
  #BF
  hm_Vals_wind<- hm_dta[1:idx]
  hm_Vals_nowind<- hm_dta[(idx+1):(2*idx)] 
  rm(hm_dta)
  
  BF_vals<- lapply(c(1:idx), FUN=function(x){
    BF_p<-hm_Vals_wind[[x]][1]/hm_Vals_nowind[[x]][1]
    BF_b<-hm_Vals_wind[[x]][2]/hm_Vals_nowind[[x]][2]
    BF_i<-hm_Vals_wind[[x]][3]/hm_Vals_nowind[[x]][3]
    BF_RI<-cbind(BF_p, BF_b, BF_i)
  })
  
  
  
  if(writeVar==T){
    lapply(c(1:idx), FUN=function(x){
      dta<- rbind(cvDatawind[[x]],cvDataNowind[[x]],sd_vals_wind[[x]],sd_vals_nowind[[x]],spd_metrics_wind[[x]],spd_metrics_nowind[[x]])
      bf_dta<- rbind(BF_vals[[x]])
      
      write.csv(dta,file = paste0(results_path,ID_list[x],"_RouteAnalysis.csv"), row.names = F)
      write.csv(bf_dta,file = paste0(results_path,ID_list[x],"_BF.csv"),row.names=F)
    })
  }
  
  return(list(cv_dta_wind=cvDatawind,
              cv_dta_nowind=cvDataNowind,
              sd_dta_wind=sd_vals_wind,
              sd_dta_nowind=sd_vals_nowind,
              spd_dta_wind=spd_vals_wind,
              sped_dta_nowind=spd_vals_nowind,
              spd_comparison=spd_comp,
              spd_wind=spd_metrics_wind,
              spd_nowind=spd_metrics_nowind,
              BayesFactor=BF_vals))
}



read_wind_dta<- function(ID_list){
  #' function to read the wind data used to compute the mean migration route
  #' @param ID_list list of folders avaliable to work trough
  #' @return list of data read from the csv files
  wind_dta_wm<- lapply(ID_list, FUN=function(x){
    df<- read.csv(paste0("C:/Users/Mike_Werfeli/Documents/Arbeit/VoWa/Results/Wind/",x, "/",x,"_WM_windDTA.csv"), sep = ',')
    wspeed<- df$X
    wsup<- df$wind_support
    wcomp<- df$compensation
    wdta<- cbind(wspeed, wsup, wcomp)
  })
  
  wind_dta_nw<- lapply(ID_list, FUN=function(x){
    df<- read.csv(paste0("C:/Users/Mike_Werfeli/Documents/Arbeit/VoWa/Results/Nowind/",x, "/",x,"_NW_windDTA.csv"), sep = ',')
    wspeed<- df$X
    wsup<- df$wind_support
    wcomp<- df$compensation
    wdta<- cbind(wspeed, wsup, wcomp)
  })
  return(c(wind_dta_wm, wind_dta_nw))
}

get_wind_spd<- function(idx, WI_dta){
  #' getter fuction for wind speed as well for WM as for NW wind data
  #' @param idx parameter to indicate the number of routes to be analyzed. It also indicates the index position of the last 
  #' element to be extracted
  #' @param WI_dta data extracted form the wind summary files
  #' @return wind speeds along the migration route
  spd_wm<- lapply(c(1:idx), FUN=function(x){
    dta<- WI_dta[[x]][,1]
  })
    
  spd_nw<- lapply(c((idx+1):(2*idx)), FUN=function(x){
    dta<- WI_dta[[x]][,1]
  })
  
  return(c(spd_wm,spd_nw))
}

get_wind_sup<- function(idx, WI_dta){
  #' getter fuction for wind support as well for WM as for NW wind data
  #' @param idx parameter to indicate the number of routes to be analyzed. It also indicates the index position of the last 
  #' element to be extracted
  #' @param WI_dta data extracted form the wind summary files
  #' @return wind supports along the migration route
  sup_wm<- lapply(c(1:idx), FUN=function(x){
    dta<- WI_dta[[x]][,2]
  })
 
  sup_nw<- lapply(c((idx+1):(2*idx)), FUN=function(x){
    dta<- WI_dta[[x]][,2]
  })
  return(c(sup_wm, sup_nw))
}

get_wind_com<- function(idx, WI_dta){
  #' getter fuction for wind compensation as well for WM as for NW wind data
  #' @param idx parameter to indicate the number of routes to be analyzed. It also indicates the index position of the last 
  #' element to be extracted
  #' @param WI_dta data extracted form the wind summary files
  #' @return wind compensation along the migration route
  com_wm<- lapply(c(1:idx), FUN=function(x){
    dta<- WI_dta[[x]][,3]
  })

  com_nw<- lapply(c((idx+1):(2*idx)), FUN=function(x){
    dta<- WI_dta[[x]][,3]
  })

  return(c(com_wm, com_nw))
}

wind_analyser<- function(ID_list,writeVar = T,results_path){
  #' function which returns usefull comparison metrics of the different wind components found along each migration route
  #' @param ID_list list of folders to work trought
  #' @return list of metrics of the wind conditions found en route
  idx<- get_idx(ID_list)
  WI_dta<- read_wind_dta(ID_list)
  wind_spd<- get_wind_spd(idx, WI_dta)
  wspd_wm<- wind_spd[1:idx]
  wspd_nw<- wind_spd[(idx+1):(2*idx)]
  rm(wind_spd)
  
  wspd_metrics_wm<- lapply(c(1:idx), FUN=function(x){
    wspd_min<- min(wspd_wm[[x]])
    wspd_max<- max(wspd_wm[[x]])
    wspd_mean<- mean(wspd_wm[[x]])
    wspd_median<- median(wspd_wm[[x]])
    wspd_dta<- cbind(wspd_min, wspd_max, wspd_mean, wspd_median)  
  })
  
  wspd_metrics_nw<- lapply(c(1:idx), FUN=function(x){
    wspd_min<- min(wspd_nw[[x]])
    wspd_max<- max(wspd_nw[[x]])
    wspd_mean<- mean(wspd_nw[[x]])
    wspd_median<- median(wspd_nw[[x]])
    wspd_dta<- cbind(wspd_min, wspd_max, wspd_mean, wspd_median)  
  })
  
  
  wind_support<- get_wind_sup(idx, WI_dta)
  wsup_wm<- wind_support[1:idx]
  wsup_nw<- wind_support[(idx+1):(2*idx)]
  rm(wind_support)
  
  wsup_metrics_wm<- lapply(c(1:idx), FUN=function(x){
    wsup_min<- min(wsup_wm[[x]])
    wsup_max<- max(wsup_wm[[x]])
    wsup_mean<- mean(wsup_wm[[x]])
    wsup_median<- median(wsup_wm[[x]])
    wsup_dta<- cbind(wsup_min, wsup_max, wsup_mean, wsup_median)  
  })
  
  wsup_metrics_nw<- lapply(c(1:idx), FUN=function(x){
    wsup_min<- min(wsup_nw[[x]])
    wsup_max<- max(wsup_nw[[x]])
    wsup_mean<- mean(wsup_nw[[x]])
    wsup_median<- median(wsup_nw[[x]])
    wsup_dta<- cbind(wsup_min, wsup_max, wsup_mean, wsup_median)  
  })
  
  wind_compensate<- get_wind_com(idx, WI_dta)
  wcom_wm<- wind_compensate[1:idx]
  wcom_nw<- wind_compensate[(idx+1):(2*idx)]
  rm(wind_compensate)
  
  wcom_metrics_wm<- lapply(c(1:idx), FUN=function(x){
    wcom_min<- min(wcom_wm[[x]])
    wcom_max<- max(wcom_wm[[x]])
    wcom_mean<- mean(wcom_wm[[x]])
    wcom_median<- median(wcom_wm[[x]])
    wcom_dta<- cbind(wcom_min, wcom_max, wcom_mean, wcom_median)  
  })
  
  wcom_metrics_nw<- lapply(c(1:idx), FUN=function(x){
    wcom_min<- min(wcom_nw[[x]])
    wcom_max<- max(wcom_nw[[x]])
    wcom_mean<- mean(wcom_nw[[x]])
    wcom_median<- median(wcom_nw[[x]])
    wcom_dta<- cbind(wcom_min, wcom_max, wcom_mean, wcom_median)  
  })
  
  rm(WI_dta)
  
  if(writeVar==T){
    lapply(c(1:idx), FUN=function(x){
      dta<- rbind(wspd_metrics_wm[[x]],wspd_metrics_nw[[x]],wsup_metrics_wm[[x]],wsup_metrics_nw[[x]],wcom_metrics_wm[[x]],wcom_metrics_nw[[x]])
      
      write.csv(dta,file = paste0(results_path,ID_list[x],"_WindAnalysis.csv"), row.names = F)
      })
  }
  return(list(wind_spd_metrics_wm=wspd_metrics_wm,
              wind_spd_metrics_nw=wspd_metrics_nw,
              wind_support_metrics_wm=wsup_metrics_wm,
              wind_support_metrics_nw=wsup_metrics_nw,
              wind_compensation_metrics_wm=wcom_metrics_wm,
              wind_compensation_metrics_wm=wcom_metrics_wm
              ))
}

wind_plotter<- function(ID_list, result_storage_path){
  idx<- get_idx(ID_list)
  WI_dta<- read_wind_dta(ID_list)
  BI_dta<- read_dta(ID_list)
  gspd<- get_spd_dta(idx, BI_dta)
  wspd<- get_wind_spd(idx, WI_dta)
  wsup<- get_wind_sup(idx, WI_dta)
  wcom<- get_wind_com(idx, WI_dta)
  rm(BI_dta)
  rm(WI_dta)
  
  airspd_wm<- lapply(c(1:idx), FUN=function(x){
    (gspd[[x]][-length(gspd[[x]])]-wsup[[x]])/3.6
  })
  
  lapply(c(1:idx), FUN=function(x){
    yl<- c(min(c(min(wspd[[x]]), min(wsup[[x]]), min(wcom[[x]])))/2, max(c(max(wspd[[x]]), max(wsup[[x]]), max(wcom[[x]])))/2)
    
    result_path<-paste0(result_storage_path,ID_list[x],"WindPlot.png")
    
    png(filename = result_path,width = 700,height = 500,units = "px")
    plot(wspd[[x]]/3.6, type="l", col="green", ylim=yl, ylab="Speed [m/s]", xlab="Position index", main=paste0("Wind conditions along the route ",ID_list[x]))
      lines(wsup[[x]]/3.6, lwd=2, col="blue")
      lines(wcom[[x]]/3.6, col="orange")
      abline(h=0, lty=2)
      legend("topright", legend = c("Wind speed", "Wind support", "Wind compensation"), cex=0.8, lwd=c(1,2,1),col = c("green","blue","orange"), title = "Wind data")
    dev.off()
    
    result_path<-paste0(result_path,ID_list[x],"SpdPlot.png")
    gspd<- gspd[[x]][-length(gspd[[x]])]
    yl<- c(min(c(min(gspd/3.6), min(airspd_wm[[x]]/3.6))), max(c(max(gspd/3.6), max(airspd_wm[[x]]/3.6))))
    
    png(filename = result_path,width = 700,height = 500,units = "px")
    plot(gspd/3.6, type="l", col="green", ylim=yl, xlab="Position index", ylab="Speed [m/s]", main=paste0("Bird speeds on route", ID_list[x]))
    lines(airspd_wm[[x]]/3.6, col="blue")
    legend("topleft", legend = c("Ground speed", "Air speed"), cex=0.8, lwd=c(1,1), col=c("green","blue"), title = "Speeds")
    dev.off()
    })
}




#run the functions
result_storage_path<- paste0("C:/Users/Mike_Werfeli/Documents/Arbeit/VoWa/Paper/Plots/Comparison/")
results_path<-"C:/Users/Mike_Werfeli/Documents/Arbeit/VoWa/Paper/Plots/Comparison/"
allRoutesFilePath<-"C:/Users/Mike_Werfeli/Documents/Arbeit/VoWa/Paper/Plots/Comparison/AllRoutes.png"

plot_Routes<- routes_plotter(ID_list, results_path, allRoutesFilePath)
analyse_Routes<- route_analyser(ID_list,writeVar = T,results_path)
analyse_Wind<- wind_analyser(ID_list,writeVar = T,results_path)

result_storage_path<- paste0("C:/Users/Mike_Werfeli/Documents/Arbeit/VoWa/Paper/Plots/Comparison/")
plot_wind<- wind_plotter(ID_list, result_storage_path)



