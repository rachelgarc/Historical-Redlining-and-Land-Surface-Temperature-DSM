# Redlining and Heat Risk in Des Moines, Iowa #
###############################################

## Loading in Map from Inequality Mapping for DSM Redlining
# install.packages("rjson")
library(dplyr)

## Loading in Map from Google Earth Engine
# install.packages("rgee")
library(rgee)
# ee_install() 
ee_Initialize()
# Des Moines bounding box
library(sf)
des_moines <- ee$Geometry$Rectangle(c(-93.85, 41.45, -93.45, 41.70))
print(des_moines)
# Loading the DSM LST File into Google Cloud
lst_collection <- ee$ImageCollection("LANDSAT/LC08/C02/T1_L2")$
  filterBounds(des_moines)$
  filterDate("2023-07-01", "2023-08-31")$
  filter(ee$Filter$lt("CLOUD_COVER", 10))

lst_median <- lst_collection$median()

lst_celsius <- lst_median$select("ST_B10")$
  multiply(0.00341802)$
  add(149.0)$
  subtract(273.15)$
  rename("LST_C")

task <- ee_image_to_drive(
  image = lst_celsius,
  region = des_moines,
  description = "DesMoines_LST_Summer2023",
  folder = "EarthEngine_exports",
  scale = 30,
  fileFormat = "GeoTIFF"
)
task$start()
ee_monitoring(task)

## Manually downloaded from G Cloud and uploading from personal file
library(terra)
LST_raster <- rast("~/Downloads/DesMoines_LST_Summer2023_2026_07_09_19_27_29.tif")

# Visualizing DSM Land Surface Temperature - very distinctly can see the river!
plot(LST_raster)
summary(values(LST_raster))

## Downloaded the Mapping Inequality .json File of HOLC Redlining Letter Grades
library(exactextractr)
HOLC <- st_read("~/Downloads/mappinginequality.json")
print(HOLC)
plot(st_geometry(HOLC))
table(HOLC$grade)

# Closing in on Des Moines
library(tidyterra)
HOLC_dsm <- HOLC %>%
  filter(city == "Des Moines")

# Combining for comparison 
st_crs(HOLC_dsm)
crs(LST_raster) 

HOLC_dsm <- st_transform(HOLC_dsm, crs = crs(LST_raster))
HOLC_dsm$mean_LST <- exact_extract(LST_raster, HOLC_dsm, "mean")
HOLC_dsm %>% st_drop_geometry() %>% select(grade, mean_LST) %>% head(10)

# Cleaning for grade replicas
table(HOLC_dsm$grade, useNA = "always")
unique(HOLC_dsm$grade)

# Looking at mean & median simply
grade_summary <- HOLC_dsm %>%
  st_drop_geometry() %>%
  filter(grade %in% c("A","B","C","D")) %>%
  group_by(grade) %>%
  summarize(
    mean_temp = mean(mean_LST, na.rm = TRUE),
    sd_temp = sd(mean_LST, na.rm = TRUE),
    median_temp = median(mean_LST, na.rm = TRUE),
    n = n()
  )
print(grade_summary)

# Visualizing with a boxplot the mean LST for each HOLC grade
library(ggplot2)

HOLC_dsm <- HOLC_dsm %>%
  mutate(grade = trimws(grade)) %>%
  filter(grade %in% c("A", "B", "C", "D")) # drop NA, E, F, and any non-residential codes

ggplot(HOLC_dsm, aes(x = grade, y = mean_LST, fill = grade)) +
  geom_boxplot() +
  scale_fill_manual(values = c("A" = "#2ca02c", "B" = "#1f77b4", 
                               "C" = "#e6b800", "D" = "#d62728")) +
  labs(title = "Land Surface Temperature by HOLC Redlining Grade",
       subtitle = "Des Moines, IA: Summer 2023 Landsat 8/9",
       x = "HOLC Grade", y = "Mean LST (°C)") +
  theme_minimal()

# Understanding the differences at a pixel level for stronger analysis
pixel_vals <- exact_extract(LST_raster, HOLC_dsm, fun = NULL, include_cols = "grade")
pixel_df <- do.call(rbind, pixel_vals)

kruskal.test(value ~ grade, data = pixel_df)

# Visualizing the LST as categorized by HOLC grades
HOLC_dsm <- st_transform(HOLC_dsm, crs = crs(LST_raster))
HOLC_bbox <- st_bbox(HOLC_dsm)
LST_cropped <- crop(LST_raster, HOLC_bbox)

overlay_map <- ggplot() +
  geom_spatraster(data = LST_cropped) +
  scale_fill_viridis_c(name = "LST (°C)", option = "magma", na.value = "transparent") +
  geom_sf(data = HOLC_dsm, aes(color = grade), fill = NA, linewidth = 0.8) +
  scale_color_manual(name = "HOLC Grade", 
                     values = c("A" = "#2ca02c", "B" = "#1f77b4", 
                                "C" = "#e6b800", "D" = "#d62728")) +
  labs(
    title = "Land Surface Temperature and Historical HOLC Redlining",
    subtitle = "Des Moines, IA — Summer 2023 Landsat 8/9 Composite",
    caption = "Sources: Mapping Inequality (DSL, U. Richmond); USGS Landsat 8-9 Collection 2"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "gray30"),
    axis.title = element_blank(),
    legend.position = "right"
  )

overlay_map

# Comparison of mean LST across HOLC Grade Pairs
library(dplyr)
library(tidyr)

# Comparison across each pair
grade_pairs <- combn(c("A","B","C","D"), 2, simplify = FALSE)

pairwise_results <- lapply(grade_pairs, function(pair) {
  g1 <- pixel_df$value[pixel_df$grade == pair[1]]
  g2 <- pixel_df$value[pixel_df$grade == pair[2]]
  
  mean_diff <- mean(g2, na.rm = TRUE) - mean(g1, na.rm = TRUE)
  median_diff <- median(g2, na.rm = TRUE) - median(g1, na.rm = TRUE)
  
  wtest <- wilcox.test(g2, g1, conf.int = TRUE)
  
  data.frame(
    comparison = paste(pair[2], "vs", pair[1]),
    mean_diff_C = round(mean_diff, 2),
    median_diff_C = round(median_diff, 2),
    p_value = wtest$p.value,
    ci_lower = round(wtest$conf.int[1], 2),
    ci_upper = round(wtest$conf.int[2], 2)
  )
})

pairwise_table <- bind_rows(pairwise_results)
print(pairwise_table)

# Rank-biserial effect size (magnitude comparison)
library(rcompanion)  

effect_sizes <- lapply(grade_pairs, function(pair) {
  g1 <- pixel_df$value[pixel_df$grade == pair[1]]
  g2 <- pixel_df$value[pixel_df$grade == pair[2]]
  
  wtest <- wilcox.test(g2, g1)
  n1 <- length(g1); n2 <- length(g2)
  
  r <- 1 - (2 * as.numeric(wtest$statistic)) / (as.numeric(n1) * n2)
  
  data.frame(comparison = paste(pair[2], "vs", pair[1]), rank_biserial_r = round(r, 3))
})

bind_rows(effect_sizes)

# AB (green lined) vs CD (red lined) binary grouping
pixel_df <- pixel_df %>%
  mutate(group_AB_CD = case_when(
    grade %in% c("A","B") ~ "AB (Greenlined)",
    grade %in% c("C","D") ~ "CD (Redlined)"
  ))

ab_vals <- pixel_df$value[pixel_df$group_AB_CD == "AB (Greenlined)"]
cd_vals <- pixel_df$value[pixel_df$group_AB_CD == "CD (Redlined)"]

# mean/median difference
mean(cd_vals, na.rm = TRUE) - mean(ab_vals, na.rm = TRUE)
median(cd_vals, na.rm = TRUE) - median(ab_vals, na.rm = TRUE)

# test + effect size
wtest_abcd <- wilcox.test(cd_vals, ab_vals, conf.int = TRUE)
wtest_abcd

n1 <- length(ab_vals)
n2 <- length(cd_vals)
r_abcd <- 1 - (2 * wtest_abcd$statistic) / (as.numeric(n1) * n2)
r_abcd
