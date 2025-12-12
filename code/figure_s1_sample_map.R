# Clear memory
cat("\014")
rm(list = ls())
gc()

library(dplyr)
library(tidyr)
library(geosphere)
library(patchwork)
library(ggspatial)
library(ggplot2)
library(sf)
library(geodata)
library(ggrepel)

setwd("C:/Users/qianh/Desktop/R/chapter3")
getwd()

# Original sampling data
df_raw <- data.frame(
  stringsAsFactors = FALSE,
  ID = c("G1", "G2", "G3", "G4", "G5", "G6"),
  Genotype = c("LS-5", "HZ-13", "NB-10", "TZ-8", "JX-22", "WZ-2"),
  Location = c("Lishui", "Hangzhou", "Ningbo", "Taizhou", "Jiaxing", "Wenzhou"),
  Latitude = c("28°23′11″N", "30°18′48″N", "29°53′50″N", "28°40′22″N", "30°41′55″N", "27°58′29″N"),
  Longitude = c("119°49′9″E", "120°23′25″E", "121°33′32″E", "121°25′37″E", "120°46′28″E", "120°45′38″E")
)

# Convert DMS to decimal degrees
df_dd <- df_raw %>%
  separate(Latitude, into = c("Lat_D", "Lat_M", "Lat_S"), sep = "[°′″N]", convert = TRUE) %>%
  separate(Longitude, into = c("Lon_D", "Lon_M", "Lon_S"), sep = "[°′″E]", convert = TRUE) %>%
  mutate(
    Lat = Lat_D + Lat_M / 60 + Lat_S / 3600,
    Lon = Lon_D + Lon_M / 60 + Lon_S / 3600
  ) %>%
  select(ID, Genotype, Location, Lat, Lon)

print(df_dd)

# Calculate the average distance
coords <- as.matrix(df_dd[, c("Lon", "Lat")])
dmat_m <- geosphere::distm(coords, fun = geosphere::distHaversine)  # meter
mean_km <- mean(dmat_m[upper.tri(dmat_m)]) / 1000
mean_km
# [1] 184.4621

# Load GADM data (China provinces and prefectures)
cache_dir <- "C:/Users/qianh/Desktop/R/GIS/gadm"
cn_l1 <- gadm(country = "CHN", level = 1, path = cache_dir) |> st_as_sf()
cn_l2 <- gadm(country = "CHN", level = 2, path = cache_dir) |> st_as_sf()

# Simplify provincial boundaries
cn_l1_proj   <- st_transform(cn_l1, 3857)
cn_l1_simplp <- st_simplify(cn_l1_proj, dTolerance = 3000, preserveTopology = TRUE)
cn_l1_simpl  <- st_transform(cn_l1_simplp, 4326)

# Load and transform South China Sea dashed line
south_sea <- st_read("./data/China_map/south_sea.shp", quiet = TRUE) |>
  st_make_valid() |>
  st_transform(4326)

# Extract Zhejiang province and prefectures
zhejiang_prov <- cn_l1 |> filter(NAME_1 == "Zhejiang")
zj_city <- cn_l2 |> filter(NAME_1 == "Zhejiang")
zj_city_proj   <- st_transform(zj_city, 3857)
zj_city_simplp <- st_simplify(zj_city_proj, dTolerance = 800, preserveTopology = TRUE)
zj_city_simpl  <- st_transform(zj_city_simplp, 4326)

# National map with South China Sea line
china.map <- ggplot() +
  geom_sf(data = cn_l1_simpl, fill = "#DBDBDB", color = "#868686", linewidth = 0.2) +
  geom_sf(data = zhejiang_prov, fill = "#737373", color = "black", linewidth = 0.4) +
  geom_sf(data = south_sea, color = "black", linewidth = 0.2, linetype = "dashed") +
  theme_void() +
  theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 1))

# Keep target prefectures
keep_cities <- c("Hangzhou","Ningbo","Wenzhou","Jiaxing","Huzhou",
                 "Shaoxing","Jinhua","Quzhou","Zhoushan","Taizhou","Lishui")
zj_keep <- zj_city_simpl |> filter(NAME_2 %in% keep_cities)

# Compute label positions (largest polygon centroid per prefecture)
zj_label_pts <- zj_keep %>%
  st_make_valid() %>%
  st_transform(3857) %>%
  st_cast("MULTIPOLYGON", warn = FALSE) %>%
  st_cast("POLYGON", warn = FALSE) %>%
  mutate(.area = st_area(geometry)) %>%
  group_by(NAME_2) %>%
  slice_max(.area, n = 1, with_ties = FALSE) %>%
  st_centroid() %>%
  st_transform(4326) %>%
  ungroup()

coords <- st_coordinates(zj_label_pts)
zj_labels <- zj_label_pts %>% mutate(X = coords[,1], Y = coords[,2])

# Subsets for each genotype
df_jx   <- filter(df_dd, Genotype == "JX-22")
df_hz   <- filter(df_dd, Genotype == "HZ-13")
df_ls   <- filter(df_dd, Genotype == "LS-5")
df_nb   <- filter(df_dd, Genotype == "NB-10")
df_tz   <- filter(df_dd, Genotype == "TZ-8")
df_wz   <- filter(df_dd, Genotype == "WZ-2")

# Map of Zhejiang prefectures and sampling sites
zj.map <- ggplot() +
  geom_sf(data = zj_keep, fill = "grey", color = "blue", linewidth = 0.3) +
  geom_point(data = df_dd, aes(x = Lon, y = Lat), color = "red", size = 3, alpha = 0.7) +
  geom_text_repel(data = df_jx, aes(x = Lon, y = Lat, label = ID),
                  size = 4, color = "black", 
                  nudge_x = 0.5, nudge_y = 0.3,
                  box.padding = 0.3, point.padding = 0.2,
                  min.segment.length = 0, segment.color = "black", segment.size  = 0.3) +
  geom_text_repel(data = df_hz, aes(x = Lon, y = Lat, label = ID),
                  size = 4, color = "black", 
                  nudge_x = -0.1, nudge_y = 0,
                  box.padding = 0.3, point.padding = 0.2) +
  geom_text_repel(data = df_ls, aes(x = Lon, y = Lat, label = ID),
                  size = 4, color = "black", 
                  nudge_x = 0.1, nudge_y = 0.1,
                  box.padding = 0.3, point.padding = 0.2) +
  geom_text_repel(data = df_nb, aes(x = Lon, y = Lat, label = ID),
                  size = 4, color = "black", 
                  nudge_x = -0.2, nudge_y = 0.1,
                  box.padding = 0.3, point.padding = 0.2) +
  geom_text_repel(data = df_tz, aes(x = Lon, y = Lat, label = ID),
                  size = 4, color = "black", 
                  nudge_x = -0.05, nudge_y = 0.25,
                  box.padding = 0.3, point.padding = 0.2,
                  min.segment.length = 0, segment.color = "black", segment.size  = 0.3) +
  geom_text_repel(data = df_wz, aes(x = Lon, y = Lat, label = ID),
                  size = 4, color = "black", 
                  nudge_x = 0, nudge_y = 0.1,
                  box.padding = 0.3, point.padding = 0.2) +
  coord_sf(xlim = c(117, 124), ylim = c(26.5, 31.5), expand = FALSE, clip = "off") +
  theme_classic(base_family = "serif") +
  theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
        axis.title   = element_text(size = 14),
        axis.text    = element_text(size = 14),
        plot.margin  = margin(8, 20, 8, 8)) +
  labs(x = "Longitude", y = "Latitude") +
  annotation_scale(location = "br", width_hint = 0.4, style = "bar",
                   bar_cols = c("grey20","white"),
                   text_cex = 1.2) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         style = north_arrow_fancy_orienteering, text_size = 10,
                         height = unit(2, "cm"), width = unit(2, "cm"))

# Labels for other prefectures vs. Zhoushan
lab_zs     <- filter(zj_labels, NAME_2 == "Zhoushan")
lab_others <- filter(zj_labels, NAME_2 != "Zhoushan")

zj.map1 <- zj.map +
  geom_text(data = lab_others, aes(x = X, y = Y, label = NAME_2),
            size = 4, family = "serif", color = "black", vjust = 0.8, alpha = 0.6) +
  geom_text_repel(data = lab_zs, aes(x = X, y = Y, label = NAME_2),
                  size = 4, family = "serif", color = "black",
                  nudge_x = 1, nudge_y = 0.15,
                  box.padding = 0.3, point.padding = 0.2, alpha = 0.6,
                  min.segment.length = 0, segment.color = "black", segment.size  = 0.3)

# Combine Zhejiang map with national inset
final_composite_map <- zj.map1 +
  patchwork::inset_element(china.map, left = 0.70, bottom = 0.05, right = 1.00, top = 0.35)

print(final_composite_map)

# Export final map
ggsave(
  filename = "./results/zhejiang_map_final.png",
  plot = final_composite_map,
  width = 6, height = 6, units = "in", dpi = 300
)
