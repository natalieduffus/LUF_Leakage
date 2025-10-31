library(sf)
library(dplyr)
library(units)
library(tibble)
library(readr)
library(here)

#define iterations for running code in HPC
# pass cluster array identifier to R
iteration <- commandArgs(trailingOnly = TRUE)[1]


# for local debugging, make iteration 1 if it is NA
if(is.na(iteration)){
  iteration <- 1
}

cat('iteration is', iteration, '\n')

#LOAD RAW DATA ----
#required raw datasets
#farm polygons
FLAME <- st_read(here('Raw_data', 'farm_base_with_geo', 'farm_base_with_geometry.shp'))

# set up the batch IDS for each array iteration

n = 1000 # the number of array iterations
x <- seq(1, nrow(FLAME))
batches <- split(FLAME$ID, sort(x%%n)) # a list 1000 items long, each consisting of a vector of IDs

cat('IDs being used in this iteration are', batches[[iteration]], '\n')

FLAME <- FLAME %>%
  filter(ID %in% batches[[iteration]])
#provisional agricultural land classification polygons
soil <- st_read(here('Raw_data', 'Provisional Agricultural Land Classification (ALC) (England)_1909723263035822565.gpkg'))
#local nature recovery strategy area polygons
LNRS_shapefiles <- st_read(here('Raw_data', 'Local_Nature_Recovery_Strategy_Areas_England.shp', 'Local_Nature_Recovery_Strategy_Areas_England.shp'))
#crop datasets 2016-2023 - polygons
crop_2016 <- st_read(here( 'Raw_data', 'lccm-2016_6040999.gpkg'))
crop_2017 <- st_read(here('Raw_data', 'lccm-2017_6041000.gpkg'))
crop_2018 <- st_read(here('Raw_data', 'lccm-2018_6040275.gpkg'))
crop_2019 <- st_read(here('Raw_data', 'lccm-2019_6040276.gpkg'))
crop_2020 <- st_read(here('Raw_data', 'lccm-2020_6040277.gpkg'))
crop_2021 <- st_read(here('Raw_data', 'lccm-2021_6038560.gpkg'))
crop_2022 <- st_read(here('Raw_data', 'lccm-2022_6038561.gpkg'))
crop_2023 <- st_read(here('Raw_data', 'lccm-2023_6038562.gpkg'))
#land use framework pixel type polygons (1:9)
LUF <- st_read(here('Raw_data', 'LUF_2010'))
#john nix crop yields dataset
yield_data <- read.csv(file = here('Raw_data', 'CropYields.csv'))
#leakage values from ball dataset
extinctions <- read.csv(file = here('Raw_data', 'GBR_leakage.csv'))


#calculate farm area in hectares and in metres-squared
FLAME  <- FLAME  %>% mutate(farm_area_m2 = as.numeric(st_area(geometry)),
                            farm_area_ha = farm_area_m2 / 10000)

#omit any farms bigger than 400ha - unlikely to be offsets and this cuts the
#run time down
FLAME$farm_area_ha <- as.numeric(FLAME$farm_area_ha)
FLAME <- FLAME %>% filter(!is.na(farm_area_ha) & farm_area_ha <= 400)

#fix known geometry issue with 2021 crop data
crop_2021 <- sf::st_make_valid(crop_2021)

#need to tidy the extinctions dataset
extinctions <- extinctions %>% 
  rename(Crop_group = ITEM_OF_INTEREST)

second_map <- c(
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
  mutate(Crop_group = recode(Crop_group, !!!second_map))

#FLAME APPEND LUF VALUES ----
#append the LUF pixel type to the FLAME dataset to enable future analysis
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
  
  #ccalculate the area per candidate type by intersecting farm with the precomputed union
  areas_t <- numeric(length(candidate_types))
  for (k in seq_along(candidate_types)) {
    tval <- candidate_types[k]
    ugeom <- type_unions[[as.character(tval)]]
    #skip: if the union has no geometry or doesn't touch farm, area = 0
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
  dom_landuse_pct[i] <- if (farm_area_m2[i] > 0) 100 * as.numeric(max_area) / farm_area_m2[i] else NA_real_
}

#append values to FLAME data
FLAME$dom_landuse <- dom_landuse
FLAME$dom_landuse_area_m2 <- dom_landuse_area_m2
FLAME$dom_landuse_pct <- dom_landuse_pct

print("appended land use pixel types")



#FLAME APPEND LNRS AREA ----

LNRS_area <- "name"
#intersections
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
  
  #intersect this farm with only candidate counties
  inter_sf <- st_intersection(FLAME[i, , drop = FALSE], LNRS_shapefiles[idxs, , drop = FALSE])
  if (nrow(inter_sf) == 0) {
    county_name[i] <- NA_character_
    county_pct[i]  <- NA_real_
    next
  }
  
  #compute area per county name
  inter_sf$area_m2 <- as.numeric(st_area(inter_sf))
  names_vec <- as.character(inter_sf[[LNRS_area]])
  area_by_county <- tapply(inter_sf$area_m2, names_vec, sum, simplify = TRUE)
  
  #compute percent of farm area for each county
  farm_area <- farm_area_m2[i]
  pct_by_county <- 100 * as.numeric(area_by_county) / farm_area
  
  #choose county: if any >50%? pick that
  over50 <- which(pct_by_county > 50)
  if (length(over50) > 0) {
    #tie-break alphabetically
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

#tidy up the dataset and remove useless columns
FLAME_tidy <- FLAME %>%
  dplyr::select(ID,
                farm_area_m2,
                Area.hec.,
                dom_landuse,
                county_name,
                geometry) 

sf::sf_use_s2(FALSE)
sf::st_agr(FLAME_tidy) <- "constant"
sf::st_agr(LUF) <- "constant"
sf::st_agr(soil) <- "constant"

#calculate intersections between FLAME and soil dataset (these are the same for
# every year)
ps <- st_intersection(FLAME_tidy, soil)
#compute area m2 of the intersections
ps$area_m2 <- as.numeric(st_area(ps))

#keep the non-spatial attributes table for FLAME to re-attach later 
flame_attrs <- st_drop_geometry(FLAME_tidy) %>% distinct(ID, .keep_all = TRUE)

#function which calculates the intersection between farm, soil and crop type in
#a given year - with chunking to speed it up
process_crop <- function(crop_sf, crop_name_field = "crop_name", ps_obj, chunk_size = 5000L) {
  if (missing(ps_obj)) stop("please pass the ps object as ps_obj")
  
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
  
  out_with_attrs <- dplyr::left_join(out, flame_attrs, by = "ID")
  
  return(out_with_attrs)
}

#run it for all crop years
FLAME_2016 <- process_crop(crop_2016, crop_name_field = "crop_name", ps_obj = ps)
FLAME_2017 <- process_crop(crop_2017, crop_name_field = "crop_name", ps_obj = ps)
FLAME_2018 <- process_crop(crop_2018, crop_name_field = "crop_name", ps_obj = ps)
FLAME_2019 <- process_crop(crop_2019, crop_name_field = "crop_name", ps_obj = ps)
FLAME_2020 <- process_crop(crop_2020, crop_name_field = "crop_name", ps_obj = ps)
FLAME_2021 <- process_crop(crop_2021, crop_name_field = "crop_name", ps_obj = ps)
FLAME_2022 <- process_crop(crop_2022, crop_name_field = "crop_name", ps_obj = ps)
FLAME_2023 <- process_crop(crop_2023, crop_name_field = "crop_name", ps_obj = ps)

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

#tidy up new FLAME dataset
#exclude grass and solar panels
FLAME_crops <- FLAME_crops %>% filter(crop_name != "Grass")
FLAME_crops <- FLAME_crops %>% filter(crop_name != "Solar panels")

#remove small polygons which are likely spatial mismatches
FLAME_crops <- FLAME_crops %>% filter(area_ha >= 0.05)

#convert areas into km2 to fit what the other datasets need
FLAME_crops <- FLAME_crops %>%
  mutate(area_km2 = area_ha * 0.01)

#need to make crop names match between yield and FLAME dataset
yield_data$Crop <- tolower(yield_data$Crop)
FLAME_crops$crop_name <- tolower(FLAME_crops$crop_name)

FLAME_crops <- FLAME_crops %>% 
  rename(Crop = crop_name)

FLAME_crops <- FLAME_crops %>%
  mutate(Crop = gsub("winter wheat \\(includes winter oats\\)", "winter wheat", Crop))

#convert yield data to kg/km2 from tonnes per ha
yield_data <- yield_data %>%
  mutate(LowProd_km2 = LowProd * 100000)

yield_data <- yield_data %>%
  mutate(AvgProd_km2 = AvgProd * 100000)

yield_data <- yield_data %>%
  mutate(HighProd_km2 = HighProd * 100000)

#multiply through yield values by FLAME dataset according to ALC type
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

#aggregate crop yields by farm, crop type, and LUF pixel type
totals_crop_year <- FLAME_yields %>%
  group_by(year, Crop, dom_landuse, ID, county_name) %>%
  summarise(
    total_area_km2 = sum(area_km2, na.rm = TRUE),
    total_yield = sum(total_yield, na.rm = TRUE), 
    .groups = "drop") %>%
  arrange(year, Crop, dom_landuse, ID)

#group crops of the same type
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

print("yields calculated")

total_path <- here('Processed_data', 'FLAME_total_crop_yields_2016_23')
if(!dir.exists(total_path)){
  dir.create(total_path)
}

write.csv(x = totals_crop_year, 
          file = paste0(total_path, '/iteration_', iteration, '.csv'))

#calculate average crop yield across the years
avg_total_by_crop <- totals_crop_year %>%
  group_by(Crop_group, dom_landuse, ID, county_name) %>%
  summarise(
    mean_total_area_km2 = mean(total_area_km2, na.rm = TRUE),   
    sd_total_area_km2   = sd(total_area_km2, na.rm = TRUE),    
    se_total_area_km2   = sd_total_area_km2 / sqrt(n_distinct(year)), 
    mean_total_yield = mean(total_yield, na.rm = TRUE),
    sd_total_yield   = sd(total_yield, na.rm = TRUE),
    n_years          = n_distinct(year),
    se_total_yield   = sd_total_yield / sqrt(n_years)
  ) %>%
  arrange(desc(mean_total_yield))

#save crop yield dataset

avg_path <- here('Processed_data', 'FLAME_avg_crop_yields_2016_23')
if(!dir.exists(avg_path)){
  dir.create(avg_path)
}

write.csv(x = avg_total_by_crop, 
          file = paste0(avg_path, '/iteration_', iteration, '.csv'))

#FLAME LEAKAGE LOSS VALUES ----

#calculate the leakage value per 1kg of displaced crop so that it can be
#multiplied through the above values


dummy_yield_df <- tibble(
  Crop_group = unname(second_map),
  yield = 1                         
)

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

avg_total_by_crop <- avg_total_by_crop %>% 
  left_join(net_leak_perkg, by = "Crop_group") %>%
  mutate(
    total_leak = net_leak * mean_total_yield,
    leak_per_crop_area = dplyr::if_else(!is.na(total_leak) & !is.na(mean_total_area_km2) & mean_total_area_km2 > 0,
                                        total_leak / mean_total_area_km2,
                                        NA_real_)
  )

avg_total_by_crop <- avg_total_by_crop %>% filter(Crop_group != "other crops")
avg_total_by_crop <- avg_total_by_crop %>% filter(Crop_group !="beet (sugar beet / fodder beet)")

# calculate total leakage per farm
per_farm_leakage <- avg_total_by_crop %>%
  group_by(ID) %>%
  summarise(
    total_leakage = sum(total_leak, na.rm = TRUE),
    .groups = "drop"
  )

farm_crop_area <- avg_total_by_crop %>%
  group_by(ID) %>%
  summarise(
    total_crop_area_km2 = sum(mean_total_area_km2, na.rm = TRUE),
    .groups = "drop"
  )

per_farm_leakage <- per_farm_leakage %>%
  left_join(farm_crop_area, by = "ID") %>%
  mutate(
    leak_per_crop_area = dplyr::if_else(total_crop_area_km2 > 0,
                                        total_leakage / total_crop_area_km2,
                                        NA_real_)
  )


# add back in the farm area and land use pixel metadata
farm_metadata <- england_leakage_values %>%
  dplyr::select(ID, dom_landuse, Area.hec., county_name) %>%  
  distinct(ID, .keep_all = TRUE)

# convert area to km2
farm_metadata <- farm_metadata %>%
  mutate(area_km2 = Area.hec. * 0.01)

# calculate the leakage per km2 for each farm
per_farm_leakage <- per_farm_leakage %>%
  left_join(farm_metadata, by = "ID") %>%
  mutate(
    leak_per_km2 = total_leakage / area_km2
  )

print("leakage appended")

leakage_path <- here('Processed_data', 'FLAME_leakage')
if(!dir.exists(leakage_path)){
  dir.create(leakage_path)
}

write.csv(x = per_farm_leakage, 
          file = paste0(leakage_path, '/iteration_', iteration, '.csv'))

print("All done! Cheers :)")
