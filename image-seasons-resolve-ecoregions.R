# Earth Engine output.
# Use the VNP22Q2: Land Surface Phenology Yearly L3 Global 500m data 
# (https://doi.org/10.5067/VIIRS/VNP22Q2.001) and calculate the per-pixel mean 
# day of Date_Mid_Greenup_Phase_1 across the 2013 to 2022 period. This gives us
# a per-pixel "image season start date". Calculate the per-pixel mean day of
# Date_Mid_Senescence_Phase_1 across the 2013 to 2022 period. This gives us a
# per-pixel "image season end date". For each Resolve ecoregion, calculate the
# spatial mean of the image season start date and image season end date and 
# export.

ee_output = read.csv("data/resolve-ecoregion-image-seasons.csv") |> 
  dplyr::select(-.geo) |> 
  dplyr::rename(img_season_start = img_s_strt, 
                img_season_end = img_s_end) |> 
  dplyr::rename_with(.fn = tolower) |> 
  dplyr::select(biome_name, biome_num, eco_name, eco_id, img_season_start, img_season_end) |> 
  dplyr::mutate(img_season_start = round(img_season_start),
                img_season_end = round(img_season_end))

# Download Resolve ecoregions from https://hub.arcgis.com/datasets/37ea320eebb647c6838c23f72abae5ef/explore
resolve = sf::st_read("data/RESOLVE_Ecoregions_and_Biomes.geojson")
names(resolve) = tolower(names(resolve))

out = dplyr::inner_join(x = resolve, y = ee_output)

plot(out[, "img_season_start"], pal = viridis::viridis, main = "Image season start DOY")
plot(out[, "img_season_end"], pal = viridis::viridis, main = "Image season end DOY")

sf::st_write(obj = out, dsn = "data/image-seasons-resolve-ecoregions.gpkg")
