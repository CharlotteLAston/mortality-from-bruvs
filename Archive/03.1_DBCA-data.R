###################################################

# Calculating proportions of different length
# classes from real data

###################################################
library(tidyverse)
library(dplyr)
library(ggplot2)
library(sf)
library(raster)
library(stringr)
library(forcats)
library(RColorBrewer)
library(geosphere)
library(abind)

# If CheckEM package has changed
# remove.packages("CheckEM")
# renv::install("GlobalArchiveManual/CheckEM")
# devtools::install_github("GlobalArchiveManual/CheckEM")
library(CheckEM)

rm(list = ls())

working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # to directory of current file - or type your own

data_dir <- paste(working.dir, "Data", sep="/")
fig_dir <- paste(working.dir, "Figures", sep="/")
pop_dir <- paste(working.dir, "Population_Files", sep="/")

googledrive_download_data <- function() {
  
  library(dplyr)
  
  options(gargle_oauth_cache = ".secrets")
  # check the value of the option, if you like
  gargle::gargle_oauth_cache()
  googlesheets4::gs4_auth(scopes="drive")
  2
  
  main.dir <- "/Users/22291142/Documents/PhD/R Code/fish-size-condition-categories/"
  # popup.dir <- "inst/app/www/popups"
  
  # Main folder with all marine parks
  drive.folder <- "https://drive.google.com/drive/folders/1GDcemLHxOvAjyecWfv83hhH-ZZh1hqfT"
  folder.id <- googledrive::drive_get(googledrive::as_id(drive.folder))
  2
  # find all folders in marine park folder
  files <- googledrive::drive_ls(folder.id, type = "folder")
  
  # files <- files %>% filter(name %in% "NMP")
  
  # loop through folders and create a directory
  for (parent in seq_along(files$name)) {
    
    # parent <- "NMP"
    
    # Print the current directory
    print(files$name[parent])
    
    # Make directory
    dir.create(file.path(main.dir, files$name[parent]), recursive = TRUE)
    
    # Find all children folders in the current folder
    current.folder <- files$id[parent]
    children.folders <- googledrive::drive_ls(current.folder, type = "folder")
    
    for (child in seq_along(children.folders$name)) {
      
      # Print the current directory
      print(children.folders$name[child])
      
      # Make directory
      dir.create(paste(main.dir, files$name[parent], children.folders$name[child], sep = "/"))
      
      if(nrow(children.folders) > 0){
        
        # Find all children folders in the current folder
        current.child.folder <- children.folders$id[child]
        
        print("view baby names")
        
        baby.folders <- googledrive::drive_ls(current.child.folder, type = "folder") %>%
          dplyr::glimpse()
        
        for (baby in seq_along(baby.folders$name)) {
          
          # Print the current directory
          print(baby.folders$name[baby])
          
          # Make directory
          dir.create(paste(main.dir, files$name[parent], children.folders$name[child], baby.folders$name[baby], sep = "/"))
          
          # list files
          i_dir <- googledrive::drive_ls(baby.folders[baby, ])
          
          #download files
          for (file_i in seq_along(i_dir$name)) {
            #fails if already exists
            try({
              
              print(i_dir$name[file_i])
              
              # print(getwd())
              
              googledrive::drive_download(googledrive::as_id(i_dir$id[file_i]),
                                          path = stringr::str_c(main.dir, files$name[parent], children.folders$name[child], baby.folders$name[baby], i_dir$name[file_i], sep = "/"))
            })
          }
        }
      }
    }
  }
  
  
  # TODO implement this in the metadata reading in - just change the codes rather than changing the folder names
  codes <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1OuOt80TvJBCMPLR6oy7YhfoSD4VjC73cuKovGobxiyI/edit?usp=sharing",
                                     sheet = "park_codes")
  
  setwd(main.dir)
  
  for(i in unique(codes$code)){
    code <- i
    new.name <- unique((codes %>% dplyr::filter(code %in% i))$full.name)
    
    print(new.name)
    
    file.rename(code, new.name)
  }
  
  # setwd("G:/mpaviewer") # NEED to find a way to fix this
  # 
  # # Folder with images
  # drive.folder <- "https://drive.google.com/drive/folders/1PeEcdENN0BhXpkzryqsBbq-0kFn_1N7z"
  # folder.id <- googledrive::drive_get(googledrive::as_id(drive.folder))
  # 
  # # find all folders in marine park folder
  # files <- googledrive::drive_ls(folder.id)
  # 
  # #download files
  # for (i in seq_along(files$name)) {
  #   try({
  #     googledrive::drive_download(googledrive::as_id(files$id[i]),
  #                                 path = stringr::str_c("inst/app/www/images", files$name[i], sep = "/"),
  #                                 overwrite = TRUE)
  #   })
  # }
  # 
  # # Folder with pop-ups
  # drive.folder <- "https://drive.google.com/drive/folders/1laVfBAmFlnrxGInyOEx-2Sp8-s5Dw4rT"
  # folder.id <- googledrive::drive_get(googledrive::as_id(drive.folder))
  # 
  # # find all folders in marine park folder
  # files <- googledrive::drive_ls(folder.id)
  # 
  # #download files
  # for (i in seq_along(files$name)) {
  #   try({
  #     googledrive::drive_download(googledrive::as_id(files$id[i]),
  #                                 path = stringr::str_c("inst/app/www/popups", files$name[i], sep = "/"),
  #                                 overwrite = TRUE)
  #   })
  # }
  # 
  # 
  # # Find all files in the pop-ups folder to turn into html
  # files <- dir("inst/app/www/popups")
  # 
  # for (i in unique(files)) {
  #   
  #   rmarkdown::render(paste0("inst/app/www/popups/", i))
  #   
  # }
  
}

read_dbca_files_csv <- function(flnm, data_dir = here::here("inst/data/raw")) {
  flnm %>%
    readr::read_csv(col_types = readr::cols(.default = "c")) %>%
    dplyr::mutate(folder.structure = stringr::str_replace_all(flnm, paste(data_dir, "/", sep = ""), "")) %>%
    tidyr::separate(folder.structure, into = c("marine_park","indicator", "method", "campaignid"), sep = "/", extra = "drop", fill = "right") %>%
    GlobalArchive::ga.clean.names()
}

read_dbca_files_txt <- function(flnm, data_dir = here::here("inst/data/raw")) {
  readr::read_tsv(flnm, col_types = readr::cols(.default = "c")) %>%
    dplyr::mutate(folder.structure = stringr::str_replace_all(flnm, paste(data_dir, "/", sep = ""), "")) %>%
    tidyr::separate(folder.structure, into = c("marine_park", "indicator", "method", "campaignid"), sep = "/", extra = "drop", fill = "right") %>%
    GlobalArchive::ga.clean.names()
}

### ► Life history sheet ----
## Will need to replace with DBCA's own version eventually but this will work for time being
# TODO use the caab common names here instead

synonyms <- CheckEM::aus_synonyms %>% dplyr::distinct()

common_names <- CheckEM::australia_life_history %>%
  dplyr::select(scientific_name, family, genus, species, australian_common_name)

dbca_googlesheet_url <- "https://docs.google.com/spreadsheets/d/1OuOt80TvJBCMPLR6oy7YhfoSD4VjC73cuKovGobxiyI/edit?usp=sharing"

life_history <- googlesheets4::read_sheet(dbca_googlesheet_url, sheet = "life_history") %>%
  CheckEM::clean_names() %>%
  dplyr::rename(code = region_code)

2

life_history <- googlesheets4::read_sheet(dbca_googlesheet_url, sheet = "functional_traits") %>%
  CheckEM::clean_names() %>%
  dplyr::select(!complex_functional_group) %>%
  dplyr::rename(trophic_group = simple_functional_group) %>%
  dplyr::glimpse()

# unique(life_history$trophic_group)

complete_sites <- googlesheets4::read_sheet(dbca_googlesheet_url, sheet = "temporal_years_sites") %>%
  CheckEM::clean_names() %>%
  dplyr::mutate(year = strsplit(as.character(include_years), split = ", "))%>%
  tidyr::unnest(year) %>%
  dplyr::mutate(site = strsplit(as.character(include_sites), split = ", "))%>%
  tidyr::unnest(site) %>%
  dplyr::distinct() %>%
  dplyr::mutate(year = as.numeric(year)) %>%
  dplyr::select(marine_park, method, year, site) %>%
  dplyr::mutate(complete = "Consistently sampled")

complete_needed_campaigns <- complete_sites %>%
  dplyr::distinct(marine_park, method) %>%
  dplyr::mutate(complete_needed = "Consistently sampled")

codes <- googlesheets4::read_sheet(dbca_googlesheet_url, sheet = "park_codes") %>%
  dplyr::rename(marine_park = full.name) %>%
  dplyr::select(marine_park, code)

fished_species <- googlesheets4::read_sheet(dbca_googlesheet_url, sheet = "target_species") %>%
  CheckEM::clean_names()

trophic_groups <- life_history %>%
  # dplyr::select(-c(target.code, feeding.guild, trophic.guild)) %>% # TODO use the regions for this
  # dplyr::mutate(trophic_group = stringr::str_replace_all(.$community, c("NANA" = "Unknown",
  #                                                                           "NA" = "Unknown",
  #                                                                           "planktivore" = "Planktivore",
  #                                                                           "Algal Feeder" = "Algal feeder"))) %>%
  # tidyr::replace_na(list(trophic_group = "Unknown")) %>%
  dplyr::distinct() %>%
  # dplyr::full_join(codes)%>%
  # dplyr::filter(!is.na(marine_park)) %>% # To get rid of ones that don't have data yet
  dplyr::filter(!is.na(genus))

test <- trophic_groups %>%
  dplyr::group_by(family, genus, species) %>%
  dplyr::summarise(n = dplyr::n())

# There are a few duplicate trophic groups that will cause errors
# TODO Use the region matching for trophic and fish!!!

# park_popups <- here::here("inst/data/parks.popups.csv") |> # BG TO DO -  CHANGE THIS
#   read.csv(na.strings = c("NA", "NaN", " ", "", NA)) |>
#   CheckEM::clean_names()

zoning <- googlesheets4::read_sheet(dbca_googlesheet_url, sheet = "park_gazettal") %>%
  CheckEM::clean_names()

foa_codes <- googlesheets4::read_sheet(dbca_googlesheet_url, sheet = "fishes_of_australia") %>%
  CheckEM::clean_names() %>%
  dplyr::select(-c(number)) %>%
  dplyr::mutate(scientific_name = paste(genus, species, sep = " ")) %>%
  dplyr::left_join(common_names) %>%
  dplyr::mutate(scientific_name = paste0(scientific_name, " (", australian_common_name, ")"))

interpretation_trends <- googlesheets4::read_sheet(dbca_googlesheet_url, sheet = "interpretation_trends") %>%
  CheckEM::clean_names()

# _______________________________________________________ ----
#                        READ IN DATA                     ----
# _______________________________________________________ ----

## ► Metadata (same for every method and data type) ----

folders <- list.files(path = paste0(data_dir), recursive = T, pattern = "_Metadata.csv", full.names = T) %>%
  as.data.frame() %>%
  dplyr::mutate(folder_structure = stringr::str_replace_all(., paste(data_dir, "/", sep = ""), "")) %>%
  tidyr::separate(folder_structure, into = c("marine_park","indicator", "method", "campaignid"), sep = "/", extra = "drop", fill = "right") %>%
  dplyr::mutate(read_method = forcats::fct_recode(method,
                                                  "point" = "BRUVs",
                                                  "point" = "BRUVS",
                                                  "transect" = "DOVs",
                                                  "transect" = "ROVs",
                                                  "transect" = "UVC_ROV")) %>%
  dplyr::distinct(marine_park, indicator, read_method, method)

metadata <- data.frame()

for(i in 1:nrow(folders)){
  
  
  folder <- folders[i,]
  path <- paste(data_dir, unique(folder$marine_park), unique(folder$indicator), unique(folder$method), sep = "/")
  
  message(path)
  
  read_method <- unique(folder$read_method)
  marine_park <- unique(folder$marine_park)
  method <- unique(folder$method)
  
  if(read_method %in% "point"){
    
    metadata_temp <- CheckEM::read_metadata(dir = path, method = "BRUVs") %>%
      dplyr::mutate(marine_park = marine_park) %>%
      dplyr::mutate(method = method)
    
  } else {
    
    metadata_temp <- CheckEM::read_metadata(dir = path, method = "DOVs")%>%
      dplyr::mutate(marine_park = marine_park) %>%
      dplyr::mutate(method = method)
    
  }
  
  metadata <- dplyr::bind_rows(metadata, metadata_temp)
  
}

metadata <- metadata %>%
  dplyr::mutate(latitude_dd = as.numeric(latitude_dd),
                longitude_dd = as.numeric(longitude_dd)) %>%
  dplyr::left_join(zoning) %>%
  
  dplyr::mutate(status = stringr::str_replace_all(.$status, c("Sanctuary" = "No-take",
                                                              "No Take" = "No-take",
                                                              "MPA" = "No-take",
                                                              "Reserve" = "No-take",
                                                              "No-Take" = "No-take",
                                                              "Protected" = "No-take"))) %>%
  dplyr::mutate(dbca_zone = stringr::str_replace_all(.$dbca_zone, c("Sanctuary Zone" = "Sanctuary"))) %>%
  dplyr::mutate(method = forcats::fct_recode(method,
                                             "stereo-BRUVs" = "BRUVs",
                                             "stereo-BRUVs" = "BRUVS",
                                             "stereo-DOVs" = "DOVs",
                                             "stereo-ROVs" = "ROVs",
                                             "stereo-ROVs+UVC" = "UVC_ROV")) %>%
  dplyr::mutate(year = as.numeric(substr(campaignid, 1, 4))) %>%
  dplyr::left_join(.,complete_sites) %>%
  dplyr::left_join(.,complete_needed_campaigns) %>%
  dplyr::mutate(complete = dplyr::if_else(is.na(complete_needed), "Consistently sampled", complete)) %>%
  tidyr::replace_na(list(complete = "Intermittently sampled")) %>%
  
  dplyr::mutate(dbca_zone = as.character(dbca_zone)) %>%
  
  dplyr::mutate(dbca_zone = stringr::str_replace_all(dbca_zone, c("Sanctaury" = "Sanctuary",
                                                                  "SP Benthic" = "SP Benthic Protection",
                                                                  "SP Benthic Protection Protection" = "SP Benthic Protection",
                                                                  "Marine Management Area" = "General Use"))) %>%
  
  dplyr::select(marine_park, method, campaignid, sample, latitude_dd, longitude_dd, date_time,
                location, status, site,
                successful_count, successful_length,
                depth_m,
                year,
                # month, day,
                gazetted, re_zoned, complete, dbca_zone, dbca_sanctuary) %>% # Trying to remove columns to save space/time to load the app
  dplyr::filter(!campaignid %in% c("2021-05_JurienBay.MP.Monitoring_UVC"))


unique(metadata$dbca_zone)

names(metadata) %>% sort()
unique(metadata$complete)

unique(metadata$year)
unique(metadata$campaignid) %>% sort()

test <- metadata %>%
  dplyr::filter(is.na(year))

unique(metadata$marine_park)

test_complete <- metadata %>%
  dplyr::filter(complete %in% "Consistently sampled")

testing <- metadata %>%
  dplyr::filter(marine_park %in% "Rottnest Island Marine Reserve")

unique(testing$status)
unique(metadata$marine_park) %>% sort()
unique(metadata$method) %>% sort()
unique(metadata$campaignid) %>% sort()


unique(metadata$dbca_sanctuary)
unique(metadata$status)

campaign_list <- metadata %>% dplyr::distinct(marine_park, method, campaignid, sample) # Want to create a list of every sample
# DOes it have maxn and length associated with it??

lats <- metadata %>%
  dplyr::group_by(marine_park) %>%
  dplyr::summarise(mean_lat = mean(latitude_dd)) %>% # biggest is the most north
  dplyr::arrange(desc(mean_lat)) # Could make this again on the server side

metadata$marine_park <- forcats::fct_relevel(metadata$marine_park, c(unique(lats$marine_park)))

# Testing
bruv_test <- metadata %>%
  dplyr::filter(method %in% "stereo-BRUVs") %>%
  dplyr::filter(marine_park %in% "Ningaloo Marine Park")

## ► Summarise to find sampling effort, this is used for the leaflet maps ----
sampling_effort <- metadata %>%
  dplyr::group_by(marine_park, method, sample, status, site, location, latitude_dd, longitude_dd, depth_m, complete) %>%
  dplyr::summarise(number_of_times_sampled = dplyr::n()) %>%
  dplyr::ungroup()

## ► Generic data (still using "sample") ----
## ► Count ----

count <- data.frame()

for(i in 1:nrow(folders)){
  
  folder <- folders[i,]
  path <- paste(data_dir, unique(folder$marine_park), unique(folder$indicator), unique(folder$method), sep = "/")
  
  message(path)
  
  read_method <- unique(folder$read_method)
  marine_park <- unique(folder$marine_park)
  method <- unique(folder$method)
  
  if(read_method %in% "point"){
    
    count_temp <- CheckEM::read_counts(dir = path, method = "BRUVs") %>%
      dplyr::mutate(marine_park = marine_park) %>%
      dplyr::mutate(method = method)
    
  } else {
    
    count_temp <- CheckEM::read_counts(dir = path, method = "DOVs")%>%
      dplyr::mutate(marine_park = marine_park) %>%
      dplyr::mutate(method = method)
    
  }
  
  count <- dplyr::bind_rows(count, count_temp)
  
}

count <- count %>%
  dplyr::mutate(count = as.numeric(count)) %>%
  dplyr::mutate(method = forcats::fct_recode(method,
                                             "stereo-BRUVs" = "BRUVs",
                                             "stereo-BRUVs" = "BRUVS",
                                             "stereo-DOVs" = "DOVs",
                                             "stereo-ROVs" = "ROVs",
                                             "stereo-ROVs+UVC" = "UVC_ROV")) %>%
  # Attempt to partially tidy the data ---
  dplyr::filter(!family %in% c("Unknown", NA)) %>%
  dplyr::mutate(species = dplyr::if_else(is.na(species), "spp", species)) %>%
  dplyr::mutate(genus = dplyr::if_else(is.na(genus), family, genus)) %>%
  dplyr::mutate(genus = dplyr::if_else(genus %in% "Unknown", family, genus))%>% dplyr::filter(!campaignid %in% c("2021-05_JurienBay.MP.Monitoring_UVC"))

unique(count$campaignid)

names(count)

### ► Length ----

length <- data.frame()

for(i in 1:nrow(folders)){
  
  folder <- folders[i,]
  path <- paste(data_dir, unique(folder$marine_park), unique(folder$indicator), unique(folder$method), sep = "/")
  
  message(path)
  
  read_method <- unique(folder$read_method)
  marine_park <- unique(folder$marine_park)
  method <- unique(folder$method)
  
  if(read_method %in% "point"){
    
    length_temp <- CheckEM::read_gen_length(dir = path, method = "BRUVs") %>%
      dplyr::mutate(marine_park = marine_park) %>%
      dplyr::mutate(method = method)
    
  } else {
    
    length_temp <- CheckEM::read_gen_length(dir = path, method = "DOVs")%>%
      dplyr::mutate(marine_park = marine_park) %>%
      dplyr::mutate(method = method)
    
  }
  
  length <- dplyr::bind_rows(length, length_temp)
  
}

length <- length %>%
  dplyr::mutate(number = as.numeric(count)) %>%
  dplyr::mutate(length = as.numeric(length_mm)) %>%
  dplyr::mutate(method = forcats::fct_recode(method,
                                             "stereo-BRUVs" = "BRUVs",
                                             "stereo-BRUVs" = "BRUVS",
                                             "stereo-DOVs" = "DOVs",
                                             "stereo-ROVs" = "ROVs",
                                             "stereo-ROVs+UVC" = "UVC_ROV")) %>%
  dplyr::filter(!is.na(length)) %>%
  # Attempt to partially tidy the data ---
  dplyr::filter(!family %in% c("Unknown", NA)) %>%
  dplyr::mutate(species = dplyr::if_else(is.na(species), "spp", species)) %>%
  dplyr::mutate(genus = dplyr::if_else(is.na(genus), family, genus)) %>%
  dplyr::mutate(genus = dplyr::if_else(genus %in% "Unknown", family, genus))%>% dplyr::filter(!campaignid %in% c("2021-05_JurienBay.MP.Monitoring_UVC"))


names(length)

### ► EventMeasure data ----
em_campaigns <- list.files(path = paste(data_dir, sep = "/"), recursive = T, pattern = "_Lengths.txt|_Lengths.TXT", full.names = T) %>%
  purrr::map_df(~ read_dbca_files_txt(.)) %>%
  dplyr::mutate(campaignid = stringr::str_replace_all(.$campaignid, c("_Lengths.txt" = "",
                                                                      "_Lengths.TXT" = ""))) %>%
  dplyr::distinct(campaignid) %>%
  dplyr::filter(!campaignid %in% c("2021-05_JurienBay.MP.Monitoring_UVC")) %>%
  dplyr::pull("campaignid")

# Read in points ----
points <- data.frame()

for(i in 1:nrow(folders)){
  
  folder <- folders[i,]
  path <- paste(data_dir, unique(folder$marine_park), unique(folder$indicator), unique(folder$method), sep = "/")
  
  message(path)
  
  read_method <- unique(folder$read_method)
  marine_park <- unique(folder$marine_park)
  method <- unique(folder$method)
  
  if(read_method %in% "point"){
    
    points_temp <- CheckEM::read_points(dir = path, method = "BRUVs") %>%
      dplyr::mutate(marine_park = marine_park) %>%
      dplyr::mutate(method = method)
    
  } else {
    
    points_temp <- CheckEM::read_points(dir = path, method = "DOVs")%>%
      dplyr::mutate(marine_park = marine_park) %>%
      dplyr::mutate(method = method)
    
  }
  
  points <- dplyr::bind_rows(points, points_temp)
  
}



points <- points %>%
  dplyr::mutate(method = forcats::fct_recode(method,
                                             "stereo-BRUVs" = "BRUVs",
                                             "stereo-BRUVs" = "BRUVS",
                                             "stereo-DOVs" = "DOVs",
                                             "stereo-ROVs" = "ROVs",
                                             "stereo-ROVs+UVC" = "UVC_ROV")) %>%
  
  # Attempt to partially tidy the data ---
  dplyr::filter(!family %in% c("Unknown", NA)) %>%
  dplyr::mutate(species = dplyr::if_else(is.na(species), "spp", species)) %>%
  dplyr::mutate(genus = dplyr::if_else(is.na(genus), family, genus)) %>%
  dplyr::mutate(genus = dplyr::if_else(genus %in% "Unknown", family, genus)) %>%
  dplyr::mutate(sample = stringr::str_replace_all(.$sample, "SIMP_20200323_PP_DOV_3.", "SIMP_20200323_PP_DOV_3"))%>% dplyr::filter(!campaignid %in% c("2021-05_JurienBay.MP.Monitoring_UVC")) # to fix mistake

unique(points$campaignid) %>% sort()

length_threed_points <- data.frame()

for(i in 1:nrow(folders)){
  
  folder <- folders[i,]
  path <- paste(data_dir, unique(folder$marine_park), unique(folder$indicator), unique(folder$method), sep = "/")
  
  message(path)
  
  read_method <- unique(folder$read_method)
  marine_park <- unique(folder$marine_park)
  method <- unique(folder$method)
  
  if(read_method %in% "point"){
    
    length_threed_points_temp <- CheckEM::read_em_length(dir = path, method = "BRUVs") %>%
      dplyr::mutate(marine_park = marine_park) %>%
      dplyr::mutate(method = method)
    
  } else {
    
    length_threed_points_temp <- CheckEM::read_em_length(dir = path, method = "DOVs")%>%
      dplyr::mutate(marine_park = marine_park) %>%
      dplyr::mutate(method = method)
    
  }
  
  length_threed_points <- dplyr::bind_rows(length_threed_points, length_threed_points_temp)
  
}

length_threed_points <- length_threed_points %>%
  dplyr::mutate(method = forcats::fct_recode(method,
                                             "stereo-BRUVs" = "BRUVs",
                                             "stereo-BRUVs" = "BRUVS",
                                             "stereo-DOVs" = "DOVs",
                                             "stereo-ROVs" = "ROVs",
                                             "stereo-ROVs+UVC" = "UVC_ROV")) %>%
  dplyr::mutate(sample = stringr::str_replace_all(.$sample, "SIMP_20200323_PP_DOV_3.", "SIMP_20200323_PP_DOV_3"))%>% dplyr::filter(!campaignid %in% c("2021-05_JurienBay.MP.Monitoring_UVC")) # to fix mistake

setwd(data_dir)
saveRDS(length, "EM_lengths_DBCA.RDS")
saveRDS(length_threed_points, "Point_lengths_DBCA.RDS")

#### Filter data to the method and species we are using ####
# setwd(data_dir)
# lengths <- readRDS("EM_lengths_DBCA.RDS")
# lenth_threed_points <-  readRDS("Point_lengths_DBCA.RDS")


## EM data
lengths1 <- length %>% 
  mutate(scientific = paste0(genus, sep=" ", species)) %>% 
  filter(scientific %in% c("Lethrinus nebulosus", "Chrysophrys auratus", "Coris auricularis", 
                           "Ophthalmolepis lineolatus", "Glaucosoma hebraicum", "Choerodon rubescens",
                           "Epinephelides armatus")) %>% 
  filter(campaignid %in% c("2019-02_NgariCapes_stereoBRUVs")) %>% 
  dplyr::select("campaignid", "sample", "length", "scientific")
  
table(lengths1$scientific)

lengths2 <- length_threed_points %>% 
  filter(!is.na(length_mm)) %>% 
  mutate(scientific = paste0(genus, sep=" ", species)) %>% 
  filter(scientific %in% c("Lethrinus nebulosus", "Chrysophrys auratus", "Coris auricularis", 
                           "Ophthalmolepis lineolatus", "Glaucosoma hebraicum", "Choerodon rubescens",
                           "Epinephelides armatus")) %>% 
  filter(campaignid %in% c("2021-05_Jurien.Bay.MP.Monitoring_stereoBRUVs", "2015-08_Ningaloo.deep.sanctuaries_stereoBRUVs",
                           "2019-08_Ningaloo.MP.Monitoring_stereoBRUVs", "2020-03_Sharkbay.MP.Monitoring_stereoDOVs")) %>% 
  dplyr::select("campaignid", "sample", "length_mm", "scientific") %>% 
  rename(length = "length_mm")

table(lengths2$scientific)

lengths_full <- rbind(lengths1, lengths2)

lengths_met <- lengths_full %>% 
  left_join(., metadata, by=c("campaignid", "sample")) %>% 
  dplyr::select("campaignid", "sample", "length", "scientific", "latitude_dd", "longitude_dd") %>% 
  rename(latitude = "latitude_dd",
         longitude = "longitude_dd") 
  

setwd(data_dir)
saveRDS(lengths_met, "DBCA_lengths.RDS")

