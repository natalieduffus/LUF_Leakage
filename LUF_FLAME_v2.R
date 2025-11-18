#LEAKAGE ANALYSIS
#AUTHOR: NATALIE ELIZABETH DUFFUS
#DATE CREATED: 20TH OCTOBER 2025
#LAST EDITED: 13TH NOVEMBER 2025

#REQUIRED PACKAGES ----
#need to double check - some of these might not actually be needed
library(sf)
library(dplyr)
library(units)
library(tibble)
library(readr)
library(here)
library(terra)
library(raster)
library(tidyr)
library(tidyverse)
library(exactextractr)
library(ggplot2)
library(countrycode)
library(scales)
library(boot)
library(multcompView)
library(forcats)
library(ggbeeswarm)
library(rnaturalearth)
library(rnaturalearthdata)
library(countrycode)


#define iterations for running code in HPC
#pass cluster array identifier to R
iteration <- commandArgs(trailingOnly = TRUE)[1]

#for local debugging, make iteration 1 if it is NA
if(is.na(iteration)){
  iteration <- 1
}

cat('iteration is', iteration, '\n')

#LOAD RAW DATA ----
#required raw datasets
#farm polygons (field level asset mapping in England (FLAME))
FLAME <- st_read(here('Raw_data', 'farm_base_with_geo', 'farm_base_with_geometry.shp'))

#set up the batch IDS for each array iteration
n = 1000 #the number of array iterations
x <- seq(1, nrow(FLAME))
batches <- split(FLAME$ID, sort(x%%n)) #a list 1000 items long, each consisting of a vector of IDs

cat('IDs being used in this iteration are', batches[[iteration]], '\n')

FLAME <- FLAME %>%
  filter(ID %in% batches[[iteration]])

#provisional agricultural land classification polygons - soil types
soil <- st_read(here('Raw_data', 'Provisional Agricultural Land Classification (ALC) (England)_1909723263035822565.gpkg'))

#local nature recovery strategy (LNRS) area polygons 
LNRS_shapefiles <- st_read(here('Raw_data', 'Local_Nature_Recovery_Strategy_Areas_England.shp', 'Local_Nature_Recovery_Strategy_Areas_England.shp'))

#crop datasets 2016-2023 - polygons
crop_2016 <- st_read('Raw_data/lccm-2016_6040999.gpkg')
crop_2017 <- st_read('Raw_data/lccm-2017_6041000.gpkg')
crop_2018 <- st_read('Raw_data/lccm-2018_6040275.gpkg')
crop_2019 <- st_read('Raw_data/lccm-2019_6040276.gpkg')
crop_2020 <- st_read('Raw_data/lccm-2020_6040277.gpkg')
crop_2021 <- st_read('Raw_data/lccm-2021_6038560.gpkg')
crop_2022 <- st_read('Raw_data/lccm-2022_6038561.gpkg')
crop_2023 <- st_read('Raw_data/lccm-2023_6038562.gpkg')

#land use framework (LUF) pixel type polygons (1:9)
LUF <- st_read(here('Raw_data', 'LUF_2010'))

#crop yield dataset (from John Nix Pocketbook of Farm Management)
yield_data <- read.csv(file = here('Raw_data', 'CropYields.csv'))

#extinction risk values from Tom Ball dataset
extinctions <- read.csv(file = here('Raw_data', 'GBR_leakage.csv'))

#life metric restoration layer
life <- rast('Raw_data/normalised_v6.tif')

#offsets polygons
offsets <- st_read('Raw_data/offsets_1808/offsetsBGS1808.shp')

print('finished reading in memory')

#DATA PREPARATION ----

#if needed to sort out the FLAME geometry
#FLAME <- sf::st_make_valid(FLAME)

#calculate total farm area in hectares and in metres-squared
FLAME  <- FLAME  %>% mutate(farm_area_m2 = as.numeric(st_area(geometry)),
                            farm_area_ha = farm_area_m2 / 10000)

#omit any farms bigger than 400ha - unlikely to be offsets and this cuts the
#run time down
FLAME$farm_area_ha <- as.numeric(FLAME$farm_area_ha)
FLAME <- FLAME %>% filter(!is.na(farm_area_ha) & farm_area_ha <= 400)

#fix known geometry issue with 2021 crop data
crop_2021 <- sf::st_make_valid(crop_2021)

#rename the extinction risk dataset columns and align crop group names
extinctions <- extinctions %>% 
  rename(Crop_group = ITEM_OF_INTEREST)

extinction_regroup <- c(
  "Wheat" = "wheat",
  "Barley"                 = "barley",
  "Oats"                 = "oats",
  "Grain; mixed"                 = "maize",
  "Potatoes"                  = "potatoes",
  "Rapeseed" = "oilseed rape",
  "Broad beans; horse beans; dry" = "beans",
  "Peas; green" = "peas"
  
)

extinctions <- extinctions %>%
  mutate(Crop_group = recode(Crop_group, !!!extinction_regroup))


#trim life layer to england boundary

uk <- geodata::gadm(country = "GBR", level = 1, path = tempdir())

levels(uk$NAME_1)
england <- na_rows

life_eng <- crop(life, england)
life_eng <- mask(life_eng, england)

#save it and load back in
writeRaster(life_eng, "Processed_data/life_raster_england.tif", overwrite = TRUE)
life <- rast('Processed_data/life_raster_england.tif')

#FLAME APPEND LUF VALUES ----

#append the LUF pixel type to the FLAME dataset to enable future analysis
#in this code we want to identify the dominant LUF pixel for each farm (can be a value of)
land_use_class <- "Type"

farm_area_m2 <- as.numeric(st_area(FLAME))

#identify candidate pairs via spatial index
hits <- st_intersects(FLAME, LUF, sparse = TRUE)

#pre-compute union per Type to avoid double-counting
#sorted vector of land use types: convert land use class types into integers and remove any NAs
types_present <- sort(unique(na.omit(as.integer(LUF[[land_use_class]]))))

#create empty list of the types present
type_unions <- vector("list", length(types_present))
names(type_unions) <- as.character(types_present)

#for each type subset it to the relevant rows, dissolve all geometries of that type
#into a single geometry
for (t in types_present) {
  type_unions[[as.character(t)]] <- st_union(LUF[which(as.integer(LUF[[land_use_class]]) == t), ])
}

#create a vector that will store our results
n <- nrow(FLAME)
dom_landuse <- integer(n)
dom_landuse_area_m2 <- numeric(n)
dom_landuse_pct <- numeric(n)

for (i in seq_len(n)) {
  idxs <- hits[[i]]
  if (length(idxs) == 0) {
    dom_landuse[i] <- NA_integer_
    dom_landuse_area_m2[i] <- 0
    dom_landuse_pct[i] <- 0
    next
  }
  
  #determine which LUF pixel types are present in the candidate polygons
  candidate_types <- sort(unique(as.integer(LUF[[land_use_class]][idxs])))
  
  #calculate the area per candidate type by intersecting farm with the precomputed union from above
  areas_t <- numeric(length(candidate_types))
  for (k in seq_along(candidate_types)) {
    tval <- candidate_types[k]
    ugeom <- type_unions[[as.character(tval)]]
    #if the union has no geometry or doesn't touch farm, area = 0
    if (is.null(ugeom) || length(ugeom) == 0 || !st_intersects(FLAME[i, , drop = FALSE], ugeom, sparse = FALSE)[1,1]) {
      areas_t[k] <- 0
      next
    }
    inter_sf <- st_intersection(FLAME[i, , drop = FALSE], ugeom)
    if (length(inter_sf) == 0 || all(st_is_empty(inter_sf))) {
      areas_t[k] <- 0
    } else {
      areas_t[k] <- sum(as.numeric(st_area(inter_sf)), na.rm = TRUE)
    }
  }
  
  #where they are more than 1 overlaps, pick type with max area as dominant type;
  #if tie (improbable), pick smaller number out of 1:9
  #add in % covered by dominant land use pixel type so we can double check this later
  if (all(areas_t == 0)) {
    dom_landuse[i] <- NA_integer_
    dom_landuse_area_m2[i] <- 0
    dom_landuse_pct[i] <- 0
    next
  }
  max_area <- max(areas_t, na.rm = TRUE)
  tied_types <- candidate_types[which(areas_t == max_area)]
  chosen_type <- min(tied_types, na.rm = TRUE)
  
  dom_landuse[i] <- chosen_type
  dom_landuse_area_m2[i] <- as.numeric(max_area)
  dom_landuse_pct[i] <- if (farm_area_m2[i] > 0) 100 * as.numeric(max_area) / farm_area_m2[i] else NA_real_
}

#append values to FLAME data
FLAME$dom_landuse <- dom_landuse
FLAME$dom_landuse_area_m2 <- dom_landuse_area_m2
FLAME$dom_landuse_pct <- dom_landuse_pct

print("appended land use pixel types")



#FLAME APPEND LNRS AREA ----

#here we want to append the LNRS area name to each farm in FLAME dataset

LNRS_area <- "name"
#again make spatial index of intersections
hits_LNRS <- st_intersects(FLAME, LNRS_shapefiles, sparse = TRUE)

#prepare output vectors
n <- nrow(FLAME)
county_name <- character(n)
county_pct  <- numeric(n)

for (i in seq_len(n)) {
  idxs <- hits_LNRS[[i]]
  if (length(idxs) == 0) {
    county_name[i] <- NA_character_
    county_pct[i]  <- NA_real_
    next
  }
  
  #intersect each farm with only candidate counties
  inter_sf <- st_intersection(FLAME[i, , drop = FALSE], LNRS_shapefiles[idxs, , drop = FALSE])
  if (nrow(inter_sf) == 0) {
    county_name[i] <- NA_character_
    county_pct[i]  <- NA_real_
    next
  }
  
  #compute area per county
  inter_sf$area_m2 <- as.numeric(st_area(inter_sf))
  names_vec <- as.character(inter_sf[[LNRS_area]])
  area_by_county <- tapply(inter_sf$area_m2, names_vec, sum, simplify = TRUE)
  
  #compute percent of farm area for each county
  farm_area <- farm_area_m2[i]
  pct_by_county <- 100 * as.numeric(area_by_county) / farm_area
  
  #choose county: if any >50% pick that
  over50 <- which(pct_by_county > 50)
  if (length(over50) > 0) {
    #tie-break (unlikely), but if so, pick alphabetically
    sel_idx <- over50[which.max(pct_by_county[over50])]
  } else {
    sel_idx <- which.max(pct_by_county)
  }
  
  chosen_name <- names(area_by_county)[sel_idx]
  chosen_pct  <- pct_by_county[sel_idx]
  
  county_name[i] <- chosen_name
  county_pct[i]  <- chosen_pct
}

#append to FLAME dataset
FLAME$county_name <- county_name
FLAME$county_pct  <- county_pct

print("appended lnrs area")

#FLAME CROP AND SOIL TYPE CALCULATION ----

#for each farm we want to calculate the intersection of crop type and soil type
#to get the total area of each crop type on each soil type on each farm

#quickly tidy up the dataset and remove columns not needed now
FLAME_tidy <- FLAME %>%
  dplyr::select(ID,
                farm_area_m2,
                Area.hec.,
                dom_landuse,
                county_name,
                geometry) 

#assume constant geometries
sf::sf_use_s2(FALSE)
sf::st_agr(FLAME_tidy) <- "constant"
sf::st_agr(LUF) <- "constant"
sf::st_agr(soil) <- "constant"

#calculate intersections between FLAME and soil dataset (these are the same for
#every year so only need computed once)
soil_farm_intersect <- st_intersection(FLAME_tidy, soil)
#compute area m2 of the intersections
soil_farm_intersect$area_m2 <- as.numeric(st_area(soil_farm_intersect))

#keep the non-spatial attributes table for FLAME to re-attach later (it gets lost in next step)
flame_attrs <- st_drop_geometry(FLAME_tidy) %>% distinct(ID, .keep_all = TRUE)

#function which calculates the intersection between farm, soil and crop type in
#a given year - with chunking to speed it up
process_crop <- function(crop_sf, crop_name_field = "crop_name", ps_obj, chunk_size = 5000L) {
  if (missing(ps_obj)) stop("please pass the soil_farm_intersect object as ps_obj")
  
  crop_sf <- dplyr::select(crop_sf, all_of(crop_name_field))
  hits <- sf::st_intersects(ps_obj, crop_sf, sparse = TRUE)
  
  pair_list <- lapply(seq_along(hits), function(i) {
    if (length(hits[[i]]) == 0) return(NULL)
    data.frame(ps_idx = i, crop_idx = hits[[i]], stringsAsFactors = FALSE)
  })
  
  pair_df <- if (length(pair_list) == 0) NULL else do.call(rbind, pair_list)
  if (is.null(pair_df) || nrow(pair_df) == 0) 
    return(tibble::tibble(ID = integer(0), ALC_GRADE = character(0), crop_name = character(0), area_ha = numeric(0)))
  
  n_pairs <- nrow(pair_df)
  results <- vector("list", ceiling(n_pairs / chunk_size))
  block_i <- 1
  
  for (s in seq(1, n_pairs, by = chunk_size)) {
    e <- min(n_pairs, s + chunk_size - 1L)
    block <- pair_df[s:e, , drop = FALSE]
    
    areas_block <- mapply(function(p_idx, c_idx) {
      inter <- sf::st_intersection(ps_obj[p_idx, , drop = FALSE],
                                   crop_sf[c_idx, , drop = FALSE])
      if (nrow(inter) == 0 || all(sf::st_is_empty(inter$geometry))) return(NULL)
      a_m2 <- sum(as.numeric(sf::st_area(inter)), na.rm = TRUE)
      tibble::tibble(
        ID = ps_obj$ID[p_idx],
        ALC_GRADE = ps_obj$ALC_GRADE[p_idx],
        crop_name = inter[[crop_name_field]][1],
        area_m2 = a_m2
      )
    }, block$ps_idx, block$crop_idx, SIMPLIFY = FALSE, USE.NAMES = FALSE)
    
    results[[block_i]] <- dplyr::bind_rows(areas_block[!vapply(areas_block, is.null, logical(1))])
    block_i <- block_i + 1
  }

  inter_tbl <- dplyr::bind_rows(results)
  if (nrow(inter_tbl) == 0)
    return(tibble::tibble(ID = integer(0), ALC_GRADE = character(0), crop_name = character(0), area_ha = numeric(0)))
  
  out <- inter_tbl %>%
    dplyr::group_by(ID, ALC_GRADE, crop_name) %>%
    dplyr::summarise(area_m2 = sum(area_m2, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(area_ha = area_m2 / 10000) %>%
    dplyr::arrange(ID, dplyr::desc(area_ha)) %>%
    dplyr::select(ID, ALC_GRADE, crop_name, area_ha)
  
  #reattach flame attributes to the dataset
  out_with_attrs <- dplyr::left_join(out, flame_attrs, by = "ID")
  
  return(out_with_attrs)
}


#run the function for all crop years - using soil_farm_intersect calculated earlier
FLAME_2016 <- process_crop(crop_2016, crop_name_field = "crop_name", ps_obj = soil_farm_intersect)
FLAME_2017 <- process_crop(crop_2017, crop_name_field = "crop_name", ps_obj = soil_farm_intersect)
FLAME_2018 <- process_crop(crop_2018, crop_name_field = "crop_name", ps_obj = soil_farm_intersect)
FLAME_2019 <- process_crop(crop_2019, crop_name_field = "crop_name", ps_obj = soil_farm_intersect)
FLAME_2020 <- process_crop(crop_2020, crop_name_field = "crop_name", ps_obj = soil_farm_intersect)
FLAME_2021 <- process_crop(crop_2021, crop_name_field = "crop_name", ps_obj = soil_farm_intersect)
FLAME_2022 <- process_crop(crop_2022, crop_name_field = "crop_name", ps_obj = soil_farm_intersect)
FLAME_2023 <- process_crop(crop_2023, crop_name_field = "crop_name", ps_obj = soil_farm_intersect)

#drop the geometry and make a df of the results
FLAME_2016 <- sf::st_drop_geometry(FLAME_2016)
FLAME_2017 <- sf::st_drop_geometry(FLAME_2017)
FLAME_2018 <- sf::st_drop_geometry(FLAME_2018)
FLAME_2019 <- sf::st_drop_geometry(FLAME_2019)
FLAME_2020 <- sf::st_drop_geometry(FLAME_2020)
FLAME_2021 <- sf::st_drop_geometry(FLAME_2021)
FLAME_2022 <- sf::st_drop_geometry(FLAME_2022)
FLAME_2023 <- sf::st_drop_geometry(FLAME_2023)

FLAME_crops <- bind_rows(FLAME_2016 %>% mutate(year = 2016),
                  FLAME_2017 %>% mutate(year = 2017),
                  FLAME_2018 %>% mutate(year = 2018),
                  FLAME_2019 %>% mutate(year = 2019),
                  FLAME_2020 %>% mutate(year = 2020),
                  FLAME_2021 %>% mutate(year = 2021),
                  FLAME_2022 %>% mutate(year = 2022),
                  FLAME_2023 %>% mutate(year = 2023),
                  .id = NULL) %>%
  dplyr::select(year, everything())

print("soil and crop calculations done")

#FLAME CROP YIELD CALCULATION ----


#exclude grass and solar panels from dataset
FLAME_crops <- FLAME_crops %>% filter(crop_name != "Grass")
FLAME_crops <- FLAME_crops %>% filter(crop_name != "Solar panels")

#remove small polygons which are likely spatial mismatches rather than actual fields
FLAME_crops <- FLAME_crops %>% filter(area_ha >= 0.05)

#convert areas into km2 to fit what the extinction dataset needs
FLAME_crops <- FLAME_crops %>%
  mutate(area_km2 = area_ha * 0.01)

#need to make crop names columns match between yield and FLAME dataset
yield_data$Crop <- tolower(yield_data$Crop)
FLAME_crops$crop_name <- tolower(FLAME_crops$crop_name)

FLAME_crops <- FLAME_crops %>% 
  rename(Crop = crop_name)

#change value of wheat to match
FLAME_crops <- FLAME_crops %>%
  mutate(Crop = gsub("winter wheat \\(includes winter oats\\)", "winter wheat", Crop))

#yield data is currently in tonnes per ha but the extinction dataset needs kg/km2
#convert yield data to kg/km2 from tonnes per ha
yield_data <- yield_data %>%
  mutate(LowProd_km2 = LowProd * 100000)

yield_data <- yield_data %>%
  mutate(AvgProd_km2 = AvgProd * 100000)

yield_data <- yield_data %>%
  mutate(HighProd_km2 = HighProd * 100000)

#multiply through yield values by FLAME dataset according to ALC type
#ALC grade 4 and urban gets low yield class, grade 2, 3, and non-agricultural gets average yield class,
#ALC grade 1 gets high yield class.
FLAME_yields <- FLAME_crops %>%
  left_join(yield_data, by = "Crop") %>%
  mutate(
    yield_per_km2 = case_when(
      ALC_GRADE %in% c("Grade 4", "Urban") ~ LowProd_km2,
      ALC_GRADE %in% c("Grade 2", "Grade 3", "Non Agricultural") ~ AvgProd_km2,
      ALC_GRADE == "Grade 1" ~ HighProd_km2
    ),
    total_yield = yield_per_km2 * area_km2
  )


#aggregate crop yields by farm, crop type, and LUF pixel type for further analysis
totals_crop_year <- FLAME_yields %>%
  group_by(year, Crop, dom_landuse, ID, county_name) %>%
  summarise(
    total_area_km2 = sum(area_km2, na.rm = TRUE),
    total_yield = sum(total_yield, na.rm = TRUE), 
    .groups = "drop") %>%
  arrange(year, Crop, dom_landuse, ID)

#also need to group crops of the same type
crop_map <- c(
  "spring wheat" = "wheat",
  "winter wheat" = "wheat",
  "winter barley"       = "barley",
  "spring barley"       = "barley",
  "field beans"       = "beans",
  "winter field beans"       = "beans",
  "spring field beans"       = "beans",
  "spring oats"       = "oats",
  "winter oats"       = "oats",
  "	beet (sugar beet / fodder beet)" = "	beet (sugar beet / fodder beet)",
  "potatoes" = "potatoes",
  "oilseed rape" = "oilseed rape",
  "maize" = "maize",
  "peas"         = "peas"
)

totals_crop_year <- totals_crop_year %>%
  mutate(Crop_group = recode(Crop, !!!crop_map)) %>%
  group_by(year, Crop_group, dom_landuse, ID, county_name) %>%
  summarise(
    total_area_km2 = sum(total_area_km2, na.rm = TRUE),
    total_yield = sum(total_yield),
    .groups = "drop"
  )

#extinction dataset cannot compute for other crops and for sugar beet - remove at this stage
totals_crop_year <- totals_crop_year %>% filter(Crop_group != "other crops")
totals_crop_year <- totals_crop_year %>% filter(Crop_group !="beet (sugar beet / fodder beet)")

print("yields calculated")


total_path <- here('Processed_data', 'FLAME_total_crop_yields_2016_23')
if(!dir.exists(total_path)){
  dir.create(total_path)
}

write.csv(x = totals_crop_year, 
          file = paste0(total_path, '/iteration_', iteration, '.csv'))

totals_crop_year <- read.csv(file= here('Results', 'england_crop_totals_0111.csv'))

#now want to calculate average yield across the years to get mean annual productivity
#this will enable us to compute extinction values later 
#for calculating averages, need to add in zeroes for the missing years - as this doesnt't happen in step above
totals_crop_year <- totals_crop_year %>%
  complete(ID, year, Crop_group, fill = list(total_area_km2 = 0, total_yield = 0))

#calculate average crop yield and crop area across the years that the dataset covers,
#average for each crop by each farm - both yield and area
avg_total_by_crop <- totals_crop_year %>%
  group_by(ID, Crop_group) %>%
  summarise(
    mean_crop_area_km2 = mean(total_area_km2, na.rm = TRUE),
    mean_yield    = mean(total_yield, na.rm = TRUE),
    .groups = "drop"
  )

#sum the crop areas and save them for later - this will be used as 'productive area'
crop_areas_per_farm <- avg_total_by_crop %>%
  group_by(ID) %>%
  summarise(
    mean_total_crop_area_km2 = sum(mean_crop_area_km2, na.rm = TRUE),
    .groups = "drop"
  )


#save crop yield dataset
avg_path <- here('Processed_data', 'FLAME_avg_crop_yields_2016_23')
if(!dir.exists(avg_path)){
  dir.create(avg_path)
}

write.csv(x = avg_total_by_crop, 
          file = paste0(avg_path, '/iteration_', iteration, '.csv'))

avg_total_by_crop <- read.csv('Results/england_crop_avg_0111.csv')
#FLAME LEAKAGE LOSS VALUES ----

#extinction risk dataset is not easily compatible with our average crop yields
#but if we calculate the leakage value per 1kg of displaced crop it can be
#multiplied through the above values
dummy_yield_df <- tibble(
  Crop_group = unname(extinction_regroup),
  yield = 1                         
)

#join the df to extinctions, multiply 1kg yield by leakage_per_kg_halted production to get
#leakage values - filter it so it only looks at Loss values (not restoration values which we will compute with LIFE layer)
extinctions_perkg <- extinctions %>%
  left_join(
    dummy_yield_df %>% dplyr::select(Crop_group, yield),
    by = c("Crop_group")
  ) %>%
  mutate(
    leakage_values = leakage_per_kg_halted_production * yield
  ) %>%
  group_by(Crop_group, yield) %>%
  ungroup() %>%
  filter(Loss_Restore == "Loss")

#computes leakage for each country the crop would be imported from, so need to sum them to get total leakage
leakage_sum <- extinctions_perkg %>%
  group_by(Crop_group) %>%
  summarise(total_leak = sum(leakage_values, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(total_leak) & total_leak != 0)

net_leak_perkg <- leakage_sum %>%
  group_by(Crop_group) %>%
  summarise(
    net_leak = sum(total_leak, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Crop_group, desc(net_leak))

#now we have values for the total leakage for 1kg of displaced crop production, so we can
#add this to our average crop yield for each farm and multiply it through to get total leakage per crop per farm
avg_total_by_crop <- avg_total_by_crop %>% 
  left_join(net_leak_perkg, by = "Crop_group") %>%
  mutate(
    total_leak = net_leak * mean_yield
  )

#now calculate total leakage per farm
per_farm_leakage <- avg_total_by_crop %>% group_by(ID) %>% 
  summarise(total_leakage = sum(total_leak, na.rm = TRUE), 
             .groups = "drop" )

#bring back in mean productive crop area per farm
per_farm_leakage <- per_farm_leakage %>%
  left_join(crop_areas_per_farm, by = "ID")

farm_crop_area <- avg_total_by_crop %>%
  group_by(ID) %>%
  summarise(
    mean_total_crop_area_km2 = sum(mean_total_area_km2, na.rm = TRUE),
    .groups = "drop"
  )

#add back in all the bits of FLAME metadata that has been lost in the above steps
farm_metadata <- FLAME_tidy %>% 
  dplyr::select(ID, dom_landuse, Area.hec., county_name) %>%  
  distinct(ID, .keep_all = TRUE)

farm_metadata <- farm_metadata %>%
  rename(total_farm_ha = Area.hec.)

#convert farm area to km2
farm_metadata <- farm_metadata %>%
  mutate(total_farm_km2 = total_farm_ha * 0.01)

#calculate the loss per km2 for farm and for productive crop area
#computing both leakage per km2 for total farm area and for productive crop area 
#in analysis should use leakage per km2 of productive crop area as this makes farms more
#comparable to offsets
per_farm_leakage <- per_farm_leakage %>%
  left_join(farm_metadata, by = "ID") %>%
  mutate(
    leak_per_farm_km2 = total_leakage / total_farm_km2,
    leak_per_croparea_km2 = total_leakage / mean_total_crop_area_km2
  )

print("leakage appended")

leakage_path <- here('Processed_data', 'FLAME_leakage')
if(!dir.exists(leakage_path)){
  dir.create(leakage_path)
}

write.csv(x = per_farm_leakage, 
          file = paste0(leakage_path, '/iteration_', iteration, '.csv'))

write.csv(per_farm_leakage, 'Results/england_leakage_0111.csv')


print("All done! Cheers :)")

#FLAME COMBINING DATASETS ----

#if running code in iterations on the high performance computing cluster, need to build datasets out of the
#1000 iterations

leakage_files <- list.files(path = here('Processed_data', 'FLAME_leakage'),
                            pattern = "\\.csv$",
                            full.names = TRUE)

leakage_list <- map(leakage_files, ~ read_csv(.x, col_types = cols(.default = "c"), na = c("", "NA", "n/a")) %>%
                      mutate(iteration = basename(.x)))

england_leakage_values <- bind_rows(leakage_list)

england_leakage_values <- england_leakage_values %>%
  mutate(
    ...1 = as.integer(...1),
    ID = as.character(ID),
    total_leakage = as.numeric(total_leakage),
    dom_landuse   = as.integer(dom_landuse),
    `Area.hec.`   = as.numeric(`Area.hec.`),
    area_km2      = as.numeric(area_km2),
    leak_per_km2  = as.numeric(leak_per_km2),
    county_name   = as.character(county_name)
  )

write.csv(england_leakage_values, here('Results', 'england_leakage_0111.csv'))

crop_t_files <- list.files(path = here('Processed_data', 'FLAME_total_crop_yields_2016_23'),
                           pattern = "\\.csv$",
                           full.names = TRUE)

crop_t_list <- map(crop_t_files, ~ read_csv(.x, col_types = cols(.default = "c"), na = c("", "NA", "n/a")) %>%
                     mutate(iteration = basename(.x)))

england_crop_t_values <- bind_rows(crop_t_list)

england_crop_t_values <- england_crop_t_values %>%
  mutate(
    ...1 = as.integer(...1),
    ID = as.character(ID),
    year = as.numeric(year),
    dom_landuse   = as.integer(dom_landuse),
    Crop_group   = as.character(Crop_group),
    total_yield     = as.numeric(total_yield),
    total_area_km2 = as.numeric(total_area_km2),
    county_name   = as.character(county_name)
  )

write.csv(england_crop_t_values, here('Results', 'england_crop_totals_0111.csv'))

crop_avg_files <- list.files(path = here('Processed_data', 'FLAME_avg_crop_yields_2016_23'),
                             pattern = "\\.csv$",
                             full.names = TRUE)

crop_avg_list <- map(crop_avg_files, ~ read_csv(.x, col_types = cols(.default = "c"), na = c("", "NA", "n/a")) %>%
                       mutate(iteration = basename(.x)))

england_crop_avg_values <- bind_rows(crop_avg_list)

england_crop_avg_values <- england_crop_avg_values %>%
  mutate(
    ...1 = as.integer(...1),
    ID = as.character(ID),
    dom_landuse   = as.integer(dom_landuse),
    Crop_group   = as.character(Crop_group),
    mean_total_yield     = as.numeric(mean_total_yield),
    sd_total_yield     = as.numeric(sd_total_yield),
    se_total_yield     = as.numeric(se_total_yield),
    n_years     = as.integer(n_years),
    county_name   = as.character(county_name),
    mean_total_area_km2 = as.numeric(mean_total_area_km2),
    sd_total_area_km2 = as.numeric(sd_total_area_km2),
    se_total_area_km2 = as.numeric(se_total_area_km2)
  )

write.csv(england_crop_avg_values, here('Results', 'england_crop_avg_0111.csv'))



#FLAME RESTORATION VALUES ----

#have calculated leakage due to displaced crop production, but now need to calculate avoided extinctions
#due to the restoration of farms
#to do so, we need to get the mean restoration value for each farm from the life layer
#need to make life layer compatible with flame layer
life_terra <- terra::rast('Processed_data/life_raster_england.tif')
print(life_terra)
print(terra::ext(life_terra))
print(terra::crs(life_terra))

#need life layer to to be in the bng CRS (british national grid) to match the FLAME (and later offsets) layer
life_bng <- terra::project(life_terra, "EPSG:27700", method = "bilinear")
print(life_bng)

#have problems with the conversion to alternative CRS, need to create temporary file with single layer raster
tmp <- tempfile(fileext = ".tif")
terra::writeRaster(life_bng, tmp, overwrite = TRUE)
life_fixed <- raster::raster(tmp) 
print(life_fixed)

#extract the mean and median raster values for each farm in the dataset and append ID to these values
flame_restore_stats <- exact_extract(
  life_fixed,
  FLAME,
  fun = c("mean", "median"),
  append_cols = "ID"
)

#FLAME NET LEAKAGE VALUES ----

#now need to add together the loss and restoration values to get the net values
#net change in extinction risk = increased extinction risk due to leakage - avoided extinctions due to restoration
england_leakage_values <- read.csv('Results/england_leakage_0111.csv')
england_leakage_values <- england_leakage_values %>% mutate(ID = as.character(ID))
flame_restore_stats <- flame_restore_stats %>% mutate(ID = as.character(ID))


#append restoration values to the leakage dataset
england_leakage_values <- england_leakage_values %>%
  left_join(
    flame_restore_stats %>% dplyr::select(ID, mean, median),
    by = "ID"
  )

#convert restoration values into a negative number as these are 'avoided extinctions'
#note: layer in raster is 'changes in extinction per lm2 of restored agriculture' so these values can be treated as restoration per km2 
england_leakage_values$median_restoration_value_km2 <- -abs(england_leakage_values$median)
england_leakage_values$mean_restoration_value_km2 <- -abs(england_leakage_values$mean)

#need a comparable value which gives leakage per productive crop area
england_leakage_values$leak_per_crop_km2 <- england_leakage_values$total_leakage/england_leakage_values$mean_total_crop_area_km2

#need to multiply the restoration values by the area to get the total value of restoration for each farm
#note doing it both per farm area and per productive crop area
england_leakage_values$restored_crop_area <- england_leakage_values$mean_restoration_value_km2 * england_leakage_values$mean_total_crop_area_km2
england_leakage_values$restored_farm_area <- england_leakage_values$mean_restoration_value_km2 * england_leakage_values$total_farm_km2

#get net change in extinction risk by adding together total leakage and restoration values for each farm
england_leakage_values$net_extinctions_crop_area <- england_leakage_values$total_leakage + england_leakage_values$restored_crop_area
england_leakage_values$net_extinctions_farm_area <- england_leakage_values$total_leakage + england_leakage_values$restored_farm_area

#calculate net change in extinction risk per km2 for each farm
england_leakage_values$change_per_farm_area <- england_leakage_values$net_extinctions_farm_area/england_leakage_values$total_farm_km2
england_leakage_values$change_per_crop_area <- england_leakage_values$net_extinctions_crop_area/england_leakage_values$mean_total_crop_area_km2


write.csv(england_leakage_values, 'Results/england_net_leakage_1211.csv')

#FLAME VISUALISATION ----

england_leakage_values <- read.csv(file = here('Results', 'england_net_leakage_1211.csv'))
#few random 0 values - there should not be zeroes in this dataset. it's a very small number so will remove for now and inspect them later
england_leakage_values <- england_leakage_values[england_leakage_values$net_extinctions_crop_area != 0, ]

summary(england_leakage_values$change_per_crop_area)
hist(england_leakage_values$change_per_crop_area, breaks = 50, main = "Distribution", xlab = "change in extinctions per km2")


#OFFSET SOIL AND CROP TYPE CALCULATION ----

#done all the calculations on the flame dataset, now need to do the same for the offsets
#note: aware that this is not as clean as the FLAME code - its from a separate block, needs tidying

#assign ID to each offset
offsets$id <- 1:nrow(offsets)

#calculate offset polygon areas
offsets <- offsets %>%
  mutate(
    total_offset_area_m2  = as.numeric(st_area(geometry)),
    total_offset_area_ha  = total_offset_area_m2 / 10000,                
    total_offset_area_km2 = total_offset_area_m2 / 1e6   
  )

#ALC and offset intersection
offset_soil_intersect <- st_intersection(
  offsets |> dplyr::select(id),
  soil |> dplyr::select(ALC_GRADE)
)

#2016 crop types
psp <- st_intersection(
  offset_soil_intersect |> dplyr::select(id, ALC_GRADE),
  crop_2016 |> dplyr::select(crop_name)
)

result_2016 <- psp |>
  mutate(area_ha = set_units(st_area(geometry), ha)) |>
  mutate(area_ha = as.numeric(area_ha)) |>
  filter(area_ha > 0) |>
  group_by(id, ALC_GRADE, crop_name) |>
  summarise(
    area_ha = sum(area_ha),
    .groups = "drop"
  ) |>
  arrange(id, desc(area_ha))

#2017 crop types
psp <- st_intersection(
  offset_soil_intersect |> dplyr::select(id, ALC_GRADE),
  crop_2017 |> dplyr::select(crop_name)
)

result_2017 <- psp |>
  mutate(area_ha = set_units(st_area(geometry), ha)) |>
  mutate(area_ha = as.numeric(area_ha)) |>
  filter(area_ha > 0) |>
  group_by(id, ALC_GRADE, crop_name) |>
  summarise(
    area_ha = sum(area_ha),
    .groups = "drop"
  ) |>
  arrange(id, desc(area_ha))

#2018 crop types
psp <- st_intersection(
  offset_soil_intersect |> dplyr::select(id, ALC_GRADE),
  crop_2018 |> dplyr::select(crop_name)
)

result_2018 <- psp |>
  mutate(area_ha = set_units(st_area(geometry), ha)) |>
  mutate(area_ha = as.numeric(area_ha)) |>
  filter(area_ha > 0) |>
  group_by(id, ALC_GRADE, crop_name) |>
  summarise(
    area_ha = sum(area_ha),
    .groups = "drop"
  ) |>
  arrange(id, desc(area_ha))

#2019 crop types
psp <- st_intersection(
  offset_soil_intersect |> dplyr::select(id, ALC_GRADE),
  crop_2019 |> dplyr::select(crop_name)
)

result_2019 <- psp |>
  mutate(area_ha = set_units(st_area(geometry), ha)) |>
  mutate(area_ha = as.numeric(area_ha)) |>
  filter(area_ha > 0) |>
  group_by(id, ALC_GRADE, crop_name) |>
  summarise(
    area_ha = sum(area_ha),
    .groups = "drop"
  ) |>
  arrange(id, desc(area_ha))

#2020 crop types
psp <- st_intersection(
  offset_soil_intersect |> dplyr::select(id, ALC_GRADE),
  crop_2020 |> dplyr::select(crop_name)
)
result_2020 <- psp |>
  mutate(area_ha = set_units(st_area(geometry), ha)) |>
  mutate(area_ha = as.numeric(area_ha)) |>
  filter(area_ha > 0) |>
  group_by(id, ALC_GRADE, crop_name) |>
  summarise(
    area_ha = sum(area_ha),
    .groups = "drop"
  ) |>
  arrange(id, desc(area_ha))

#2021 crop types
psp <- st_intersection(
  offset_soil_intersect |> dplyr::select(id, ALC_GRADE),
  crop_2021 |> dplyr::select(crop_name)
)

result_2021 <- psp |>
  mutate(area_ha = set_units(st_area(geometry), ha)) |>
  mutate(area_ha = as.numeric(area_ha)) |>
  filter(area_ha > 0) |>
  group_by(id, ALC_GRADE, crop_name) |>
  summarise(
    area_ha = sum(area_ha),
    .groups = "drop"
  ) |>
  arrange(id, desc(area_ha))

#2022 crop types
psp <- st_intersection(
  offset_soil_intersect |> dplyr::select(id, ALC_GRADE),
  crop_2022 |> dplyr::select(crop_name)
)

result_2022 <- psp |>
  mutate(area_ha = set_units(st_area(geometry), ha)) |>
  mutate(area_ha = as.numeric(area_ha)) |>
  filter(area_ha > 0) |>
  group_by(id, ALC_GRADE, crop_name) |>
  summarise(
    area_ha = sum(area_ha),
    .groups = "drop"
  ) |>
  arrange(id, desc(area_ha))

#2023 crop types
psp <- st_intersection(
  offset_soil_intersect |> dplyr::select(id, ALC_GRADE),
  crop_2023 |> dplyr::select(crop_name)
)

result_2023 <- psp |>
  mutate(area_ha = set_units(st_area(geometry), ha)) |>
  mutate(area_ha = as.numeric(area_ha)) |>
  filter(area_ha > 0) |>
  group_by(id, ALC_GRADE, crop_name) |>
  summarise(
    area_ha = sum(area_ha),
    .groups = "drop"
  ) |>
  arrange(id, desc(area_ha))

#drop geometries and bind results together

result_2016 <- sf::st_drop_geometry(result_2016)
result_2017 <- sf::st_drop_geometry(result_2017)
result_2018 <- sf::st_drop_geometry(result_2018)
result_2019 <- sf::st_drop_geometry(result_2019)
result_2020 <- sf::st_drop_geometry(result_2020)
result_2021 <- sf::st_drop_geometry(result_2021)
result_2022 <- sf::st_drop_geometry(result_2022)
result_2023 <- sf::st_drop_geometry(result_2023)

years <- 2016:2023
long <- bind_rows(mget(paste0("result_", years)), .id = "src") %>%
  mutate(year = as.integer(sub("^result_", "", src))) %>%
  dplyr::select(year, everything(), -src)

write.csv(long, "crop_data_2016_2023.csv")
offset_yields <- read.csv('Processed_data/crop_data_2016_2023.csv')

#OFFSET YIELD CALCULATIONS ----

#again fix all the names so they match and remove things that can't be calculated later
offset_yields$crop_name <- tolower(offset_yields$crop_name)

offset_yields <- offset_yields %>% 
  rename(Crop = crop_name)

offset_yields <- offset_yields %>% filter(Crop != "grass")
offset_yields <- offset_yields %>% filter(Crop != "solar panels")
offset_yields <- offset_yields %>% filter(Crop != "beet (sugar beet / fodder beet)")
offset_yields <- offset_yields %>% filter(Crop != "other crops")

#remove very small areas where there have obviously been snipping errors
offset_yields <- offset_yields %>% filter(area_ha >= 0.05)

#convert area of offsets to km2
offset_yields <- offset_yields %>%
  mutate(area_km2 = area_ha * 0.01)

offset_yields <- offset_yields %>%
  mutate(Crop = gsub("winter wheat \\(includes winter oats\\)", "winter wheat", Crop))

yield_data <- yield_data %>%
  mutate(LowProd_km2 = LowProd * 100000)

yield_data <- yield_data %>%
  mutate(AvgProd_km2 = AvgProd * 100000)

yield_data <- yield_data %>%
  mutate(HighProd_km2 = HighProd * 100000)

offset_data_w_yields <- offset_yields %>%
  left_join(yield_data, by = "Crop") %>%
  mutate(
    yield_per_km2 = case_when(
      ALC_GRADE %in% c("Grade 4", "Urban") ~ LowProd_km2,
      ALC_GRADE %in% c("Grade 2", "Grade 3", "Non Agricultural") ~ AvgProd_km2,
      ALC_GRADE == "Grade 1" ~ HighProd_km2
    ),
    total_yield = yield_per_km2 * area_km2
  )

#need to add in zeroes for years that have crop production
offset_data_w_yields <- offset_data_w_yields %>%
  complete(id, year, Crop, fill = list(area_km2 = 0, area_ha = 0, total_yield = 0))

#group crops of the same type
offset_data_w_yields <- offset_data_w_yields %>%
  mutate(Crop_group = recode(Crop, !!!crop_map)) %>%
  group_by(id, year, Crop_group) %>%
  summarise(
    total_yield = sum(total_yield),
    area_ha = sum(area_ha),
    area_km2 = sum(area_km2),
    .groups = "drop"
  )

#calculate the average offset yields per offset and per crop
avg_offset_yields <- offset_data_w_yields %>%
  group_by(id, Crop_group) %>%
  summarise(
    mean_crop_area_km2 = mean(area_km2, na.rm = TRUE),
    mean_yield    = mean(total_yield, na.rm = TRUE),
    .groups = "drop"
  )

#sum the crop areas and save them for later: again, productive crop area important
crop_areas_per_offset <- avg_offset_yields %>%
  group_by(id) %>%
  summarise(
    mean_total_crop_area_km2 = sum(mean_crop_area_km2, na.rm = TRUE),
    .groups = "drop"
  )


#OFFSET LEAKAGE LOSS VALUES ----

#calculate the leakage caused by displaced production using dummy_df created earlier
offset_avg_total_by_crop <- avg_offset_yields %>% 
  left_join(net_leak_perkg, by = "Crop_group") %>%
  mutate(
    total_leak = net_leak * mean_yield
  )

#calculate total leakage per offset
per_offset_leakage <- offset_avg_total_by_crop %>%
  group_by(id) %>%
  summarise(
    total_leakage = sum(total_leak, na.rm = TRUE),
    .groups = "drop"
  )

#add back in the important metadata but sans geometry
offsets <- sf::st_drop_geometry(offsets)

offset_metadata <- offsets %>%
  dplyr::select(id, GainSite, total_offset_area_km2, total_offset_area_ha) %>%  
  distinct(id, .keep_all = TRUE)

per_offset_leakage <- per_offset_leakage %>%
  left_join(crop_areas_per_offset, by = "id")

#calculate leak per crop area and per offset area 
per_offset_leakage <- per_offset_leakage %>%
  left_join(offset_metadata, by = "id") %>%
  mutate(
    leak_per_crop_area = total_leakage / mean_total_crop_area_km2,
    leak_per_offset_area = total_leakage / total_offset_area_km2
  )


#OFFSET RESTORATION VALUES ----

#need to get the mean restoration value for each farm from the life layer
offset_restore_stats <- exact_extract(
  life_fixed,
  offsets,
  fun = c("mean", "median"),
  append_cols = "id"
)

per_offset_leakage <- per_offset_leakage %>% mutate(id = as.character(id))
offset_restore_stats <- offset_restore_stats %>% mutate(id = as.character(id))

#append restoration values to the leakage dataset
offset_leakage_values <- per_offset_leakage %>%
  left_join(
    offset_restore_stats %>% dplyr::select(id, mean, median),
    by = "id"
  )

#store restoration values per km2 (as above)
offset_leakage_values$mean_restoration_value_km2 <- -abs(offset_leakage_values$mean)
offset_leakage_values$median_restoration_value_km2 <- -abs(offset_leakage_values$median)

#need to multiply the restoration values by the area to get total restoration value
offset_leakage_values$restored_crop_area <- offset_leakage_values$mean_restoration_value_km2 * offset_leakage_values$mean_total_crop_area_km2
offset_leakage_values$restored_offset_area <- offset_leakage_values$mean_restoration_value_km2 * offset_leakage_values$total_offset_area_km2

#leak per area
offset_leakage_values$leak_per_offset_area <- offset_leakage_values$total_leakage/offset_leakage_values$total_offset_area_km2
offset_leakage_values$leak_per_crop_area <- offset_leakage_values$total_leakage/offset_leakage_values$mean_total_crop_area_km2


#get net extinctions by summing
offset_leakage_values$net_extinctions_crop_area <- offset_leakage_values$total_leakage + offset_leakage_values$restored_crop_area
offset_leakage_values$net_extinctions_offset_area <- offset_leakage_values$total_leakage + offset_leakage_values$restored_offset_area

#calculate net change per area 
offset_leakage_values$change_per_offset_area <- offset_leakage_values$net_extinctions_offset_area/offset_leakage_values$total_offset_area_km2
offset_leakage_values$change_per_crop_area <- offset_leakage_values$net_extinctions_crop_area/offset_leakage_values$mean_total_crop_area_km2



#OFFSET APPEND LUF VALUES ----

#add LUF values to offsets - note: need to move to earlier in the code.

land_use_class <- "Type"

total_offset_area_m2 <- as.numeric(st_area(offsets))

#identify candidate pairs via spatial index
hits <- st_intersects(offsets, LUF, sparse = TRUE)

#pre-compute union per Type to avoid double-counting
#sorted vector of land use types: convert land use class types into integers and remove any NAs
types_present <- sort(unique(na.omit(as.integer(LUF[[land_use_class]]))))

#create empty list of the types present
type_unions <- vector("list", length(types_present))
names(type_unions) <- as.character(types_present)

#for each type subset it to the relevant rows, dissolve all geometries of that type
#into a single geometry
for (t in types_present) {
  type_unions[[as.character(t)]] <- st_union(LUF[which(as.integer(LUF[[land_use_class]]) == t), ])
}

#create a vector that will store our results
n <- nrow(offsets)
dom_landuse <- integer(n)
dom_landuse_area_m2 <- numeric(n)
dom_landuse_pct <- numeric(n)

for (i in seq_len(n)) {
  idxs <- hits[[i]]
  if (length(idxs) == 0) {
    dom_landuse[i] <- NA_integer_
    dom_landuse_area_m2[i] <- 0
    dom_landuse_pct[i] <- 0
    next
  }
  
  #determine which LUF pixel types are present in the candidate polygons
  candidate_types <- sort(unique(as.integer(LUF[[land_use_class]][idxs])))
  
  #ccalculate the area per candidate type by intersecting offset with the precomputed union
  areas_t <- numeric(length(candidate_types))
  for (k in seq_along(candidate_types)) {
    tval <- candidate_types[k]
    ugeom <- type_unions[[as.character(tval)]]
    #if the union has no geometry or doesn't touch offset, area = 0
    if (is.null(ugeom) || length(ugeom) == 0 || !st_intersects(offsets[i, , drop = FALSE], ugeom, sparse = FALSE)[1,1]) {
      areas_t[k] <- 0
      next
    }
    inter_sf <- st_intersection(offsets[i, , drop = FALSE], ugeom)
    if (length(inter_sf) == 0 || all(st_is_empty(inter_sf))) {
      areas_t[k] <- 0
    } else {
      areas_t[k] <- sum(as.numeric(st_area(inter_sf)), na.rm = TRUE)
    }
  }
  
  #where they are more than 1 overlaps, pick type with max area as dominant type;
  #if tie (improbable), pick smaller number out of 1:9
  #add in % covered by dominant land use pixel type
  if (all(areas_t == 0)) {
    dom_landuse[i] <- NA_integer_
    dom_landuse_area_m2[i] <- 0
    dom_landuse_pct[i] <- 0
    next
  }
  max_area <- max(areas_t, na.rm = TRUE)
  tied_types <- candidate_types[which(areas_t == max_area)]
  chosen_type <- min(tied_types, na.rm = TRUE)
  
  dom_landuse[i] <- chosen_type
  dom_landuse_area_m2[i] <- as.numeric(max_area)
  dom_landuse_pct[i] <- if (total_offset_area_m2[i] > 0) 100 * as.numeric(max_area) / total_offset_area_m2[i] else NA_real_
}

#append values to offsets data
offsets$dom_landuse <- dom_landuse
offsets$dom_landuse_area_m2 <- dom_landuse_area_m2
offsets$dom_landuse_pct <- dom_landuse_pct

print("appended land use pixel types")

#add this to leakage values
offset_leakage_values <- offset_leakage_values %>%
  left_join(offsets %>% dplyr::select(id, dom_landuse), by = "id")

offset_leakage_values <- sf::st_drop_geometry(offset_leakage_values)
offset_leakage_values <- dplyr::select(offset_leakage_values, -geometry)

write.csv(offset_leakage_values, 'Results/offsets_net_leakage_1211.csv')

#OFFSET VISUALISATION ----
offset_leakage_values <- read.csv(file = here('Results', 'offsets_net_leakage_1211.csv'))

summary(offset_leakage_values$change_per_crop_area)
hist(offset_leakage_values$change_per_crop_area, breaks = 50, main = "Distribution", xlab = "change in extinctions per km2")


#remove really strange value for project 34 for timebeing
offset_leakage_values <- offset_leakage_values[offset_leakage_values$id != "34", ]

#order offsets by total value and make df long
offset_leakage_values_long <- offset_leakage_values %>%
  pivot_longer(cols = c(leak_per_crop_area, median_restoration_value_km2),
               names_to = "type", values_to = "value") %>%
  mutate(
    id = fct_reorder(
      as.character(id),
      change_per_crop_area,
      .fun = function(x) min(x, na.rm = TRUE),
      .na_rm = TRUE,
      .desc = FALSE
    )
  )

#visualise distribution of leakage and restoration values
ggplot(offset_leakage_values_long, aes(x = id, y = value, fill = type)) +
  geom_col() +
  geom_hline(yintercept = 0, color = "black") +
  scale_fill_manual(
    values = c(leak_per_crop_area = "steelblue", median_restoration_value_km2 = "firebrick"),
    labels = c(leak_per_crop_area = "Increased extinction risk (leakage)",
               median_restoration_value_km2 = "Avoided extinction risk (restoration)")
  ) +
  labs(x = "Offset ID",
       y = "Change in extinction risk per km2",
       fill = NULL,
       title = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



#TOTAL OFFSET LEAKAGE VALUES BY COUNTRY ----

#we want to look at the leakage values per country too for offsets to map where the
#displaced biodiversity impacts are going to 

#sum the impact on each crop type
total_offset_impacts <- avg_offset_yields %>%
  group_by(Crop_group) %>%
  summarise(
    mean_total_crop_area_km2 = sum(mean_crop_area_km2, na.rm = TRUE),
    total_crop_yield = sum(mean_yield, na.rm = TRUE),
    .groups = "drop"
  )

#get the impacts for each crop
#attaching the yields to the leakage dataset to do so 
total_offset_impacts <- extinctions %>%
  left_join(
    total_offset_impacts %>% dplyr::select(Crop_group, total_crop_yield, mean_total_crop_area_km2),
    by = c("Crop_group")
  ) %>%
  mutate(
    leakage_values = leakage_per_kg_halted_production * total_crop_yield
  ) %>%
  group_by(Crop_group)

#add country names to the dataset for plotting impacts later
total_offset_impacts <- total_offset_impacts %>%
  mutate(
    continent   = countrycode(Country_ISO, origin = "iso3c", destination = "continent"),
    wb_region   = countrycode(Country_ISO, origin = "iso3c", destination = "region")
  )

#check for missing ones
total_offset_impacts %>% 
  filter(is.na(continent)) %>% 
  distinct(Country_ISO)

total_offset_impacts %>% 
  filter(is.na(wb_region)) %>% 
  distinct(Country_ISO)

#now need to sum the leakage values and restoration values to get the net impact 
#nice save point
write.csv(total_offset_impacts, here('Results', 'offsets_crop_summed_impacts1211.csv'))
total_offset_impacts <- read.csv(file = here('Results', 'offsets_crop_summed_impacts1211.csv'))

net_impact_total_offsets <- total_offset_impacts %>%
  group_by(Crop_group, Loss_Restore, continent) %>%
  summarise(total_leak = sum(leakage_values, na.rm = TRUE), .groups = "drop")

net_impact_total_offsets <- net_impact_total_offsets %>%
  filter(!is.na(total_leak)) %>%
  filter(total_leak != 0)

#all restoration happens in England, so for purposes of plotting - put England as contintental value for restoration
net_impact_total_offsets <- net_impact_total_offsets %>%
  mutate(
    fill_group = ifelse(Loss_Restore == "Restore", "England", continent)
  )

impact_total_offset_losses <- net_impact_total_offsets %>%
  filter(Loss_Restore == "Loss") %>%
  group_by(Crop_group) %>%
  summarise(loss_total = sum(total_leak, na.rm = TRUE), .groups = "drop")


#plot offset impacts by continent
ggplot(net_impact_total_offsets,
       aes(x = Crop_group, y = total_leak, fill = fill_group)) +
  geom_col(position = "stack") +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.2) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "England" = "#009E73",
      "Africa"        = "#0072B2",
      "Americas"      = "#F0E442",
      "Asia"          = "#56B4E9",
      "Europe"        = "#CC79A7",
      "Oceania"       = "#E69F00"
    ),
    name = "Region"
  ) +
  labs(
    x = "Crop",
    y = "Change in expected extinctions",
    title = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14, colour ="black", angle=45, hjust = 1),
    axis.text.y = element_text(size = 14, colour ="black"),
    axis.title.x = element_text(size = 16, colour ="black"),
    axis.title.y = element_text(size = 16, colour ="black")
  )

#net impact per crop type
net_by_crop_offsets <- net_impact_total_offsets %>%
  group_by(Crop_group) %>%
  summarise(
    restore = sum(total_leak[Loss_Restore == "Restore"], na.rm = TRUE),
    losses  = sum(total_leak[grepl("Loss", Loss_Restore, ignore.case = TRUE)], na.rm = TRUE),
    net     = restore + losses
  )

#going to create a map that plots where all the impacts are

total_offset_impacts <- total_offset_impacts %>%
  filter(!is.na(leakage_values))

total_offset_impacts <- total_offset_impacts %>%
  mutate(
    country_name = countrycode(Country_ISO, origin = "iso3c", destination = "country.name")
  )

total_offset_impacts_iso <- total_offset_impacts %>%
  mutate(
    iso3 = countrycode(country_name,
                       origin = "country.name",
                       destination = "iso3c",
                       warn = FALSE)
  ) 

#summarising the impacts by country - using net impact, but for every country apart from GBR
#there will only be leakage
impact_by_country <- total_offset_impacts_iso %>%
  group_by(iso3, Loss_Restore) %>%
  summarise(total = sum(leakage_values, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Loss_Restore, values_from = total, values_fill = 0) %>%
  mutate(
    net_impact = Loss - Restore )

world_raw <- ne_countries(scale = "medium", returnclass = "sf")

world <- world_raw %>%
  mutate(
    iso3_ne = iso_a3,
    iso3_ne = ifelse(iso3_ne == "-99" | is.na(iso3_ne), wb_a3,   iso3_ne),
    iso3_ne = ifelse(iso3_ne == "-99" | is.na(iso3_ne), adm0_a3, iso3_ne),
    iso3_ne = ifelse(iso3_ne == "-99" | is.na(iso3_ne), sov_a3,  iso3_ne),
    iso3_ne = ifelse(iso3_ne == "-99" | is.na(iso3_ne),
                     countrycode(name_long, "country.name", "iso3c", warn = FALSE),
                     iso3_ne),
    iso3_ne = toupper(iso3_ne)
  ) %>%
  mutate(
    iso3_ne = case_when(
      name_long == "Somaliland"                   ~ "SOM",   
      name_long == "Northern Cyprus"              ~ "CYP",   
      name_long == "Kosovo"                       ~ "KOS",   
      name_long == "Indian Ocean Territories"     ~ "IOT",   #
      name_long == "Ashmore and Cartier Islands"  ~ "AUS",   
      name_long == "Siachen Glacier"              ~ NA_character_, 
      name_long == "France"                       ~ "FRA",
      name_long == "Norway"                       ~ "NOR",
      TRUE                                        ~ iso3_ne
    )
  ) %>%
  transmute(iso3_ne, country_label = name_long, geometry)

world <- world %>%
  filter(country_label != "Ashmore and Cartier Islands")

impact_by_country <- impact_by_country %>% mutate(iso3 = toupper(iso3))

world_biodiversity <- world %>%
  left_join(impact_by_country, by = c("iso3_ne" = "iso3"))

#plot - using a log10 so that the values can be differentiated from each other
#on regular scale the values are too close together
ggplot(world_biodiversity) +
  geom_sf(aes(fill = Loss), colour = "grey40", linewidth = 0.1) +
  scale_fill_gradient(
    low = "lightyellow", high = "firebrick",
    na.value = "grey90",
    trans = "log10",
    labels = scales::label_scientific(digits =2),
    name = "Increased extinction risk (log10)"
  ) +
  labs(
    title = "",
    subtitle = ""
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.title = element_text(face = "bold", size = "14"),
    legend.text = element_text(size = "12"),
    plot.title   = element_text(face = "bold")
  )


#OFFSET AND FLAME SIDE BY SIDE ----

#want to compare the distribution of the leakage values from FLAME with the leakage values from offsets
#going to start off with change_per_crop_area

combined_change_per_area <- bind_rows(
  farms  = tibble(value = england_leakage_values$change_per_crop_area),
  offsets = tibble(value = offset_leakage_values$change_per_crop_area),
  .id = "source"
)

#plot density of values
combined_change_per_area_density <- ggplot(combined_change_per_area, aes(x = value)) +
  geom_density(fill = "pink", alpha = 0.6) +
  facet_wrap(~ source, ncol = 1, scales = "fixed") +
  labs(x = "Change in extinction risk per km2", y = "Density") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14, colour ="black"),
    axis.text.y = element_text(size = 14, colour ="black"),
    axis.title.x = element_text(size = 16, colour ="black"),
    axis.title.y = element_text(size = 16, colour ="black")
  )
print(combined_change_per_area_density)

#there is a really long tail to the negative values, will inspect further
summary(combined_change_per_area$value)
quantile(combined_change_per_area$value, probs = c(0.001, 0.01, 0.05, 0.5, 0.95), na.rm = TRUE)

#also plot a histogram
ggplot(combined_change_per_area, aes(x = value)) +
  geom_histogram(bins = 500) +
  scale_x_continuous(limits = c(-0.0005, 0.0005)) +
  theme_minimal()

#there are a few rare extreme negative values which comprise less than 1% of observations
#they are of a much more extreme magnitude - should treat these as outliers or analyse them differently?
#otherwise they will potentially dominante any analysis of variance
#note: ask about LIFE extreme values / check which areas of England get these extreme LIFE values
#e.g., if they are clustered in a particularly geographic area

#but for now, focus on values within the 0.001 and 0.999 quantiles
combined_q_lo <- quantile(combined_change_per_area$value, probs = 0.001, na.rm = TRUE)
combined_q_hi <- quantile(combined_change_per_area$value, probs = 0.999, na.rm = TRUE)

trimmed_combined_change_per_area <- combined_change_per_area %>%
  filter(!is.na(value)) %>%
  filter(value > combined_q_lo, value < combined_q_hi)

#plot again to see whats going on
trimmed_combined_change_per_area_density <- ggplot(trimmed_combined_change_per_area, aes(x = value)) +
  geom_density(fill = "pink", alpha = 0.6) +
  facet_wrap(~ source, ncol = 1, scales = "fixed") +
  labs(x = "Change in extinction risk per km2", y = "Density") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14, colour ="black"),
    axis.text.y = element_text(size = 14, colour ="black"),
    axis.title.x = element_text(size = 16, colour ="black"),
    axis.title.y = element_text(size = 16, colour ="black")
  )
print(trimmed_combined_change_per_area_density)

ggplot(trimmed_combined_change_per_area, aes(x = value)) +
  geom_histogram(bins = 500) +
  scale_x_continuous(limits = c(-0.0005, 0.0005)) +
  theme_minimal()

#still quite long tails but will proceed
#now want to identify if the medians are different between offsets and farms 

#nonparametric test to identify differences
wilcox.test(value ~ source, data = trimmed_combined_change_per_area)

#statistically significant differences between medians, but should do median
#bootstrapping so can get estimates of variance and the 95% confidence intervals
farms_vec  <- trimmed_combined_change_per_area$value[trimmed_combined_change_per_area$source == "farms"]
offsets_vec <- trimmed_combined_change_per_area$value[trimmed_combined_change_per_area$source == "offsets"]
farms_vec  <- na.omit(farms_vec)
offsets_vec <- na.omit(offsets_vec)

boot_median <- function(data, idx) median(data[idx])

boot_farms <- boot(farms_vec, statistic = boot_median, R = 10000)
boot_offsets <- boot(offsets_vec, statistic = boot_median, R = 10000)

ci_farms <- boot.ci(boot_farms, type = "perc")$percent[4:5]
ci_offsets <- boot.ci(boot_offsets, type = "perc")$percent[4:5]

# calculate the observed difference in bootstrapped medians
obs_diff <- median(farms_vec) - median(offsets_vec)
obs_diff

#bootstrap the difference between medians to get 95% CIs
diff_fun <- function(data, i1, i2) median(data[i1]) - median(data[i2])
combined <- c(farms_vec, offsets_vec)
group <- c(rep(1, length(farms_vec)), rep(2, length(offsets_vec)))
boot_diff <- boot(
  data = combined,
  statistic = function(data, idx) {
    x <- data[idx]
    g <- group[idx]
    median(x[g==1]) - median(x[g==2])
  },
  R = 10000
)

boot.ci(boot_diff, type = "perc")

#histogram of bootstrapped median differences
hist(boot_diff$t, 
     breaks = 50, 
     main = "Distribution of bootstrapped median difference (x10,000)",
     xlab = "Median difference (Farms - Offsets)",
     col = "lightblue", border = "white")
abline(v = median(boot_diff$t), col = "red", lwd = 2)
abline(v = quantile(boot_diff$t, c(0.025, 0.975)), col = "darkgray", lty = 2)

plot(density(boot_diff$t), 
     main = "Distribution of bootstrapped median difference",
     xlab = "Median difference", 
     lwd = 2)
abline(v = 0, col = "black", lty = 3)
abline(v = quantile(boot_diff$t, c(0.025, 0.975)), col = "gray", lty = 2)
abline(v = median(boot_diff$t), col = "red", lwd = 2)

#the difference between medians for farms and offsets is really tiny, offsets have
#median leakage per km2 0.0000021583 greater than for farms.


#FLAME LUF ANALYSIS NET CHANGE ----

#now want to look at how farms in different parts of the land use framework compare in terms of
#their change in extinctions per km2

#plot cube root
england_leakage_values$cbrt_change_per_crop_area <- england_leakage_values$change_per_crop_area^(1/3)
offset_leakage_values$cbrt_change_per_crop_area <- offset_leakage_values$change_per_crop_area^(1/3)

#set colours (taking from defra land use framework document for LUF)
LUF_colours <- c("#3A7050", "#4D9C59", "#5FCB63","#AFDDB7", "#F0EBFF", "#C2B3E8", "#967FD3", "#6A7898", "#8BABA9")
offset_colour <- "darkblue" 

#treat the land_use categories as a character in the farm df so that we can append offsets to them
df_farms <- england_leakage_values %>%
  filter(!is.na(dom_landuse) & !is.na(cbrt_change_per_crop_area)) %>%
  mutate(dom_landuse = as.character(dom_landuse)) %>%
  dplyr::select(dom_landuse, cbrt_change_per_crop_area)

#create an offsets df where the land_use is called 'offsets' to make it easier to plot altogether
df_offsets <- offset_leakage_values %>%
  filter(!is.na(cbrt_change_per_crop_area)) %>%
  mutate(dom_landuse = "offsets") %>%
  dplyr::select(dom_landuse, cbrt_change_per_crop_area)

all_leakage_values <- bind_rows(df_farms, df_offsets)

#for plotting, i want farm medians lowest - highest, so will sort that order here
median_farm_order <- df_farms %>%
  group_by(dom_landuse) %>%
  summarise(med = stats::median(cbrt_change_per_crop_area, na.rm = TRUE), .groups = "drop") %>%
  arrange(med) %>%
  pull(dom_landuse) %>%
  as.character()

#plotting order
plot_levels <- c(median_farm_order, "offsets")

#trying to keep correct colours for each LUF category
canonical_labels <- as.character(1:9)
canonical_map <- setNames(LUF_colours, canonical_labels)

#add offsets
canonical_map <- c(canonical_map, "offsets" = offset_colour)

all_leakage_values <- all_leakage_values %>%
  mutate(dom_landuse = factor(as.character(dom_landuse), levels = plot_levels))

plot_levels_now <- levels(all_leakage_values$dom_landuse)  
cols_for_plot <- canonical_map[plot_levels_now]          

#create a colour vector that makes sure the order of LUF colours matches
named_cols_plot <- setNames(as.character(cols_for_plot), plot_levels_now)


ggplot(all_leakage_values, aes(x = dom_landuse, y = cbrt_change_per_crop_area, fill = dom_landuse)) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.03, size = 0.35, colour = "black") +
  stat_summary(fun = stats::median, geom = "point", shape = 21, size = 2, colour = "black", fill = "white") +
  scale_fill_manual(values = named_cols_plot, guide = "none") +
  coord_flip() +
  labs(x = "Land use framework category", y = "cube-root (change in extinctions per km2)") +
  theme_minimal(base_size = 13)

#boxplot is messy and unintuitive - try dot plot that maybe shows distribution of data better
ggplot(all_leakage_values, aes(x = dom_landuse, y = cbrt_change_per_crop_area)) +
  geom_quasirandom(aes(color = dom_landuse),
                   size = 0.7,         
                   alpha = 0.25,           
                   width = 0.3,           
                   groupOnX = TRUE) +      
  stat_summary(fun = stats::median, geom = "point", shape = 21, size = 2.2,
               colour = "black", fill = "white") +
  scale_color_manual(values = named_cols_plot, guide = "none") +
  coord_flip() +
  labs(x = "Land use framework category", y = "cube-root (change in extinctions per km2)") +
  theme_minimal()


#non-parametric test to identify differences between groups
kruskal.test(cbrt_change_per_crop_area ~ dom_landuse, data = england_leakage_values)

#pairwise comparison
pairwise.wilcox.test(england_leakage_values$cbrt_change_per_crop_area, england_leakage_values$dom_landuse,
                     p.adjust.method = "BH")

#now bootstrap median differences between groups
set.seed(42)
R_boot <- 10000

#boostrapping: compute median and 95% confidence intervals
boot_ci_by_group <- all_leakage_values %>%
  group_by(dom_landuse) %>%
  summarise(
    median = stats::median(cbrt_change_per_crop_area, na.rm = TRUE),
    n = sum(!is.na(cbrt_change_per_crop_area)),
    ci = list({
      x <- na.omit(cbrt_change_per_crop_area)
      if (length(x) < 2) {
        c(NA_real_, NA_real_)
      } else {
        med_boots <- replicate(R_boot, median(sample(x, size = length(x), replace = TRUE)))
        quantile(med_boots, probs = c(0.025, 0.975), na.rm = TRUE)
      }
    }),
    .groups = "drop"
  ) %>%
  mutate(ci_lower = sapply(ci, `[`, 1),
         ci_upper = sapply(ci, `[`, 2))

#plot
ggplot(all_leakage_values, aes(x = dom_landuse, y = cbrt_change_per_crop_area, fill = dom_landuse)) +
  geom_errorbar(data = boot_ci_by_group, inherit.aes = FALSE,
                mapping = aes(x = dom_landuse, ymin = ci_lower, ymax = ci_upper),
                width = 0.3, linewidth = 1, colour = "black") +
  geom_point(data = boot_ci_by_group, inherit.aes = FALSE,
                 mapping = aes(x = dom_landuse, y = median, fill = dom_landuse),
                shape = 21, size = 5, colour = "black") +
  scale_fill_manual(values = named_cols_plot, guide = "none")  +  coord_flip() +
  labs(x = "Land use framework category", y = "cube-root (change in extinctions per km2)",
                  title = "") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(size = 14, colour ="black", angle=45, hjust = 1),
    axis.text.y = element_text(size = 14, colour ="black"),
    axis.title.x = element_text(size = 16, colour ="black"),
    axis.title.y = element_text(size = 16, colour ="black"),
    legend.position = "none"
  )



#FLAME LUF LEAKAGE ----
#now want to look at how farms in different parts of the land use framework compare in terms of
#their leakage - increase in extinctions per km2 due to displaced crop production

#plot cube root
england_leakage_values$cbrt_leak_per_crop_area <- england_leakage_values$leak_per_crop_km2^(1/3)
offset_leakage_values$cbrt_leak_per_crop_area <- offset_leakage_values$leak_per_crop_km2^(1/3)

#set colours (taking from defra land use framework document for LUF)
LUF_colours <- c("#3A7050", "#4D9C59", "#5FCB63","#AFDDB7", "#F0EBFF", "#C2B3E8", "#967FD3", "#6A7898", "#8BABA9")
offset_colour <- "darkblue" 

#treat the land_use categories as a character in the farm df so that we can append offsets to them
df_farms <- england_leakage_values %>%
  filter(!is.na(dom_landuse) & !is.na(cbrt_leak_per_crop_area)) %>%
  mutate(dom_landuse = as.character(dom_landuse)) %>%
  dplyr::select(dom_landuse, cbrt_leak_per_crop_area)

#create an offsets df where the land_use is called 'offsets' to make it easier to plot altogether
df_offsets <- offset_leakage_values %>%
  filter(!is.na(cbrt_leak_per_crop_area)) %>%
  mutate(dom_landuse = "offsets") %>%
  dplyr::select(dom_landuse, cbrt_leak_per_crop_area)

all_leakage_values_leak <- bind_rows(df_farms, df_offsets)

#for plotting, i want farm medians lowest - highest, so will sort that order here
median_farm_order_leak <- df_farms %>%
  group_by(dom_landuse) %>%
  summarise(med = stats::median(cbrt_leak_per_crop_area, na.rm = TRUE), .groups = "drop") %>%
  arrange(med) %>%
  pull(dom_landuse) %>%
  as.character()

#plotting order
plot_levels <- c(median_farm_order_leak, "offsets")

#trying to keep correct colours for each LUF category
canonical_labels <- as.character(1:9)
canonical_map <- setNames(LUF_colours, canonical_labels)

#add offsets
canonical_map <- c(canonical_map, "offsets" = offset_colour)

all_leakage_values_leak <- all_leakage_values_leak %>%
  mutate(dom_landuse = factor(as.character(dom_landuse), levels = plot_levels))

plot_levels_now <- levels(all_leakage_values_leak$dom_landuse)  
cols_for_plot <- canonical_map[plot_levels_now]          

#create a colour vector that makes sure the order of LUF colours matches
named_cols_plot <- setNames(as.character(cols_for_plot), plot_levels_now)


ggplot(all_leakage_values_leak, aes(x = dom_landuse, y = cbrt_leak_per_crop_area, fill = dom_landuse)) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.03, size = 0.35, colour = "black") +
  stat_summary(fun = stats::median, geom = "point", shape = 21, size = 2, colour = "black", fill = "white") +
  scale_fill_manual(values = named_cols_plot, guide = "none") +
  coord_flip() +
  labs(x = "Land use framework category", y = "cube-root (change in extinctions per km2)") +
  theme_minimal(base_size = 13)

#boxplot is messy and unintuitive - try dot plot that maybe shows distribution of data better
ggplot(all_leakage_values_leak, aes(x = dom_landuse, y = cbrt_leak_per_crop_area)) +
  geom_quasirandom(aes(color = dom_landuse),
                   size = 0.7,         
                   alpha = 0.25,           
                   width = 0.3,           
                   groupOnX = TRUE) +      
  stat_summary(fun = stats::median, geom = "point", shape = 21, size = 2.2,
               colour = "black", fill = "white") +
  scale_color_manual(values = named_cols_plot, guide = "none") +
  coord_flip() +
  labs(x = "Land use framework category", y = "cube-root (change in extinctions per km2)") +
  theme_minimal()

kruskal.test(cbrt_leak_per_crop_area ~ dom_landuse, data = england_leakage_values)

#pairwise comparison
pairwise.wilcox.test(england_leakage_values$cbrt_leak_per_crop_area, england_leakage_values$dom_landuse,
                     p.adjust.method = "BH")

#now bootstrap median differences between groups
set.seed(42)
R_boot <- 10000

#boostrapping: compute median and 95% confidence intervals
leak_boot_ci_by_group <- all_leakage_values_leak %>%
  group_by(dom_landuse) %>%
  summarise(
    median = stats::median(cbrt_leak_per_crop_area, na.rm = TRUE),
    n = sum(!is.na(cbrt_leak_per_crop_area)),
    ci = list({
      x <- na.omit(cbrt_leak_per_crop_area)
      if (length(x) < 2) {
        c(NA_real_, NA_real_)
      } else {
        med_boots <- replicate(R_boot, median(sample(x, size = length(x), replace = TRUE)))
        quantile(med_boots, probs = c(0.025, 0.975), na.rm = TRUE)
      }
    }),
    .groups = "drop"
  ) %>%
  mutate(ci_lower = sapply(ci, `[`, 1),
         ci_upper = sapply(ci, `[`, 2))

#plot
ggplot(all_leakage_values_leak, aes(x = dom_landuse, y = cbrt_leak_per_crop_area, fill = dom_landuse)) +
  geom_errorbar(data = leak_boot_ci_by_group, inherit.aes = FALSE,
                mapping = aes(x = dom_landuse, ymin = ci_lower, ymax = ci_upper),
                width = 0.3, linewidth = 1, colour = "black") +
  geom_point(data = leak_boot_ci_by_group, inherit.aes = FALSE,
             mapping = aes(x = dom_landuse, y = median, fill = dom_landuse),
             shape = 21, size = 5, colour = "black") +
  scale_fill_manual(values = named_cols_plot, guide = "none")  +  coord_flip() +
  labs(x = "Land use framework category", y = "cube-root (change in extinctions per km2)",
       title = "") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(size = 14, colour ="black", angle=45, hjust = 1),
    axis.text.y = element_text(size = 14, colour ="black"),
    axis.title.x = element_text(size = 16, colour ="black"),
    axis.title.y = element_text(size = 16, colour ="black"),
    legend.position = "none"
  )

#FLAME LUF RESTORATION ----
#now want to look at how farms in different parts of the land use framework compare in terms of
#their restoration potential - decrease in extinctions per km2 due to habitat restoration

#all restoration values are negative so need sign preserving transformation

#plot cube root
england_leakage_values <- england_leakage_values %>%
  mutate(cbrt_restoration_km2 = sign(median_restoration_value_km2) * abs(median_restoration_value_km2)^(1/3))
offset_leakage_values <- offset_leakage_values %>%
  mutate(cbrt_restoration_km2 = sign(median_restoration_value_km2) * abs(median_restoration_value_km2)^(1/3))

#set colours (taking from defra land use framework document for LUF)
LUF_colours <- c("#3A7050", "#4D9C59", "#5FCB63","#AFDDB7", "#F0EBFF", "#C2B3E8", "#967FD3", "#6A7898", "#8BABA9")
offset_colour <- "darkblue" 

#treat the land_use categories as a character in the farm df so that we can append offsets to them
df_farms <- england_leakage_values %>%
  filter(!is.na(dom_landuse) & !is.na(cbrt_restoration_km2)) %>%
  mutate(dom_landuse = as.character(dom_landuse)) %>%
  dplyr::select(dom_landuse, cbrt_restoration_km2)

#create an offsets df where the land_use is called 'offsets' to make it easier to plot altogether
df_offsets <- offset_leakage_values %>%
  filter(!is.na(cbrt_restoration_km2)) %>%
  mutate(dom_landuse = "offsets") %>%
  dplyr::select(dom_landuse, cbrt_restoration_km2)

all_leakage_values_restore <- bind_rows(df_farms, df_offsets)

#for plotting, i want farm medians lowest - highest, so will sort that order here
median_farm_order_restore <- df_farms %>%
  group_by(dom_landuse) %>%
  summarise(med = stats::median(cbrt_restoration_km2, na.rm = TRUE), .groups = "drop") %>%
  arrange(med) %>%
  pull(dom_landuse) %>%
  as.character()

#plotting order
plot_levels <- c(median_farm_order_restore, "offsets")

#trying to keep correct colours for each LUF category
canonical_labels <- as.character(1:9)
canonical_map <- setNames(LUF_colours, canonical_labels)

#add offsets
canonical_map <- c(canonical_map, "offsets" = offset_colour)

all_leakage_values_restore <- all_leakage_values_restore %>%
  mutate(dom_landuse = factor(as.character(dom_landuse), levels = plot_levels))

plot_levels_now <- levels(all_leakage_values_restore$dom_landuse)  
cols_for_plot <- canonical_map[plot_levels_now]          

#create a colour vector that makes sure the order of LUF colours matches
named_cols_plot <- setNames(as.character(cols_for_plot), plot_levels_now)


ggplot(all_leakage_values_restore, aes(x = dom_landuse, y = cbrt_restoration_km2, fill = dom_landuse)) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.03, size = 0.35, colour = "black") +
  stat_summary(fun = stats::median, geom = "point", shape = 21, size = 2, colour = "black", fill = "white") +
  scale_fill_manual(values = named_cols_plot, guide = "none") +
  coord_flip() +
  labs(x = "Land use framework category", y = "cube-root (change in extinctions per km2)") +
  theme_minimal(base_size = 13)

#boxplot is messy and unintuitive - try dot plot that maybe shows distribution of data better
ggplot(all_leakage_values_restore, aes(x = dom_landuse, y = cbrt_restoration_km2)) +
  geom_quasirandom(aes(color = dom_landuse),
                   size = 0.7,         
                   alpha = 0.25,           
                   width = 0.3,           
                   groupOnX = TRUE) +      
  stat_summary(fun = stats::median, geom = "point", shape = 21, size = 2.2,
               colour = "black", fill = "white") +
  scale_color_manual(values = named_cols_plot, guide = "none") +
  coord_flip() +
  labs(x = "Land use framework category", y = "cube-root (change in extinctions per km2)") +
  theme_minimal()

kruskal.test(cbrt_restoration_km2 ~ dom_landuse, data = england_leakage_values)

#pairwise comparison
pairwise.wilcox.test(england_leakage_values$cbrt_restoration_km2, england_leakage_values$dom_landuse,
                     p.adjust.method = "BH")

#now bootstrap median differences between groups
set.seed(42)
R_boot <- 10000

#boostrapping: compute median and 95% confidence intervals
restore_boot_ci_by_group <- all_leakage_values_restore %>%
  group_by(dom_landuse) %>%
  summarise(
    median = stats::median(cbrt_restoration_km2, na.rm = TRUE),
    n = sum(!is.na(cbrt_restoration_km2)),
    ci = list({
      x <- na.omit(cbrt_restoration_km2)
      if (length(x) < 2) {
        c(NA_real_, NA_real_)
      } else {
        med_boots <- replicate(R_boot, median(sample(x, size = length(x), replace = TRUE)))
        quantile(med_boots, probs = c(0.025, 0.975), na.rm = TRUE)
      }
    }),
    .groups = "drop"
  ) %>%
  mutate(ci_lower = sapply(ci, `[`, 1),
         ci_upper = sapply(ci, `[`, 2))

#plot
ggplot(all_leakage_values_restore, aes(x = dom_landuse, y = cbrt_restoration_km2, fill = dom_landuse)) +
  geom_errorbar(data = restore_boot_ci_by_group, inherit.aes = FALSE,
                mapping = aes(x = dom_landuse, ymin = ci_lower, ymax = ci_upper),
                width = 0.3, linewidth = 1, colour = "black") +
  geom_point(data = restore_boot_ci_by_group, inherit.aes = FALSE,
             mapping = aes(x = dom_landuse, y = median, fill = dom_landuse),
             shape = 21, size = 5, colour = "black") +
  scale_fill_manual(values = named_cols_plot, guide = "none")  +  coord_flip() +
  labs(x = "Land use framework category", y = "cube-root (change in extinctions per km2)",
       title = "") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(size = 14, colour ="black", angle=45, hjust = 1),
    axis.text.y = element_text(size = 14, colour ="black"),
    axis.title.x = element_text(size = 16, colour ="black"),
    axis.title.y = element_text(size = 16, colour ="black"),
    legend.position = "none"
  )
