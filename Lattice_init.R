#====Load Packages====
library(tidyverse)
library(ggspatial)
library(tidyterra)
library(terra)
library(cowplot)
library(scatterpie)
library(lubridate)
library(brms)

#====Set aesthetics====
theme_set(  theme_bw()+
              theme(panel.grid = element_line(color = "white"),
                    axis.text = element_blank(),
                    axis.title = element_blank(),
                    axis.ticks = element_blank(),
                    legend.position = "top",
                    legend.key.height = unit(0.5, "lines"),
                    legend.title = element_text(vjust = 1),
                    plot.title = element_text(hjust = 0.5)))

ggpreview <- function(...) {
  fname <- tempfile(fileext = ".png")
  ggsave(filename = fname, ...)
  system2("open", fname)
  invisible(NULL)
}

#====import data====


#Read Iceland from data downloaded from here: https://atlas.lmi.is/LmiData/index.php?id=3399612930 
#Downloaded IS 50V 17/06 2018 Strandlina. GDB and SHP. ISN 2004
# icelandshp_sf <- st_read("map/IS50V_STRANDLINA_17062019_ISN2004/IS50V_STRANDLINA_SHP/is50v_strandlina_linur_17062019.shp")

#sf seems to produce better maps than using rgdal imported data, 
#but to get the CRS in a comparable format, use readOGR
#same data
icelandshp <- vect("map/IS50V_STRANDLINA_17062019_ISN2004/IS50V_STRANDLINA_SHP/is50v_strandlina_linur_17062019.shp")

#extract Coordinate Reference System (CRS)
ISCRS <- crs(icelandshp)

#read in Myvatn shapefile (Lacks CRS in metadata)
myvatnshp <- vect("map/Myvatn fixed/Myvatn_fixed.shp") 

#set CRS with known epsg for the data
crs(myvatnshp) <- crs("+init=epsg:32628")

# #convert the CRS to match the Iceland file to plot together
# new <- spTransform(myvatnshp, CRS(ISCRS))



#=====Load data====
file <- "Lattice_MASTER_19Jan24.xlsx"

site_orig <- read_csv("1km_grid.csv")

grid_sites <- readxl::read_xlsx(file, sheet ="Coordinates", na = "NA") %>%  mutate(year = as.character(year(sampledate)),
                                                                              spot2021 = ifelse(year == "2021", spot,
                                                                                                ifelse(spot>=27+31, spot-30, spot-31)))
profiles <- readxl::read_xlsx(file, sheet = "profiles", na = "NA")
turner <- readxl::read_xlsx(file, sheet = "Turner", na = "NA")
cc <- readxl::read_xlsx(file, sheet = "chiro_counts", na = "NA")
cm <- readxl::read_xlsx(file, sheet = "chiro_measures", na = "NA")
bt <- readxl::read_xlsx(file, sheet = "benthotorch", na = "NA")
zoops <-  readxl::read_xlsx(file, sheet = "zooplankton", na = "NA")
esites <- read_csv("site_meta.csv", col_types =  ("cccddd"))
nutrients <- readxl::read_xlsx(file, sheet = "nutrients", na = "NA")

#=====combine data=====
# pull midges
midges <- cc %>% 
  mutate(chironominae = (chiro + tanyt)/fract_count) %>% 
  select(spot, chironominae)

# pull probes
pel_algae <- turner %>% 
  group_by(spot) %>% 
  summarise(chl = mean(chl),
            phyc = mean(phyc))

# pull data from profiles
site_depth <- profiles %>% 
  select(spot, depth) %>% 
  unique()

# prep nutrients

bothyears <- c("nh4", "no3", "po4") # select nutrients in both years

nut1 <- nutrients %>% 
  mutate(analyte2 = ifelse(analyte == "nh4_nh3", "nh4", analyte),
         layer_analyte = paste(layer, analyte2, sep = "_")) %>% 
  filter(analyte2 %in% bothyears) %>%                  
  select(spot, layer_analyte, mgl) %>%
  spread(layer_analyte, mgl)


# prep site information and join
sites <- grid_sites %>% 
  mutate(east.std = easting/1000 - min(easting/1000), #km distance from most western site
         north.std = northing/1000 - min(northing/1000)) %>% # km distance from most southern site
  select(site = spot2021, spot, year, sed_type, east.std, north.std, easting, northing)

nut <- sites %>% 
  full_join(site_depth) %>% 
  full_join(nut1) %>% 
  full_join(midges) %>% 
  full_join(pel_algae)

#====create adjacency matrix====
dist_mat <- as.matrix(dist(cbind(nut$northing, nut$easting), diag = TRUE, upper = TRUE)/1000)

row.names(dist_mat) <- nut$spot
colnames(dist_mat) <- nut$spot

adj_mat <- dist_mat
# plot to check

gs_check <- grid_sites %>% 
  group_by(spot2021) %>% 
  summarise(easting = mean(easting),
            northing = mean(northing))

thresh = 1.28 # select sites <thresh as being adjacent.

adj <- as.data.frame(as.matrix(dist_mat)) %>% 
  mutate(spot = row.names(.)) %>% 
  gather(spot2, dist.km, -spot) %>% 
  filter(as.numeric(spot)<=31,
         as.numeric(spot2)<=31) %>%
  mutate(adjacent = dist.km<thresh) 



adj_check <- adj %>% 
  filter(adjacent) %>% 
  mutate(spot2021 = as.numeric(spot),
         spot2 = as.numeric(spot2)) %>% 
  left_join(gs_check) %>% 
  left_join(gs_check %>% 
              rename(spot2 = spot2021, east2 = easting, north2 = northing))

# plot
adj_check %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", data = myvatnshp)+
  # geom_point(aes(easting, northing))+
  geom_text(aes(easting, northing, label = spot2021))+
  geom_segment(aes(x = easting, y = northing, xend = east2, yend = north2))+
  coord_sf()+
  guides(fill = "none")

# check for many points
adj_check %>% 
  ggplot()+
  facet_wrap(~spot) +
  geom_spatvector(fill = "dodgerblue", data = myvatnshp)+
  geom_point(aes(easting, northing))+
  # geom_text(aes(easting, northing, label = spot2021))+
  geom_segment(aes(x = easting, y = northing, xend = east2, yend = north2))+
  # coord_sf()+
  guides(fill = "none")

# create adjacency matrix
adj_mat[dist_mat<thresh] <-1
adj_mat[dist_mat>thresh] <- 0
adj_mat


#====Fit model====

start_time <- Sys.time()
nh4mod <- brm(mvbind(ben_nh4, pel_nh4)~year + depth + east.std + north.std + sed_type + chironominae + phyc +
                car(adj_mat, gr = site),
              data = nut,
              data2 = list(adj_mat = adj_mat, spot2021 = nut$site),
              iter = 500,
              chains = 4,
              cores = 4,
              family = "hurdle_gamma")
Sys.time() - start_time

start_time <- Sys.time()
po4mod <- brm(mvbind(ben_po4, pel_po4)~year + depth + east.std + north.std + sed_type + chironominae + phyc +
                car(adj_mat, gr = site),
              data = nut,
              data2 = list(adj_mat = adj_mat, spot2021 = nut$site),
              iter = 500,
              chains = 4,
              cores = 4,
              family = "hurdle_gamma")
Sys.time() - start_time


start_time <- Sys.time()
no3mod <- brm(mvbind(ben_no3, pel_no3)~year + depth + east.std + north.std + sed_type + chironominae + phyc +
                car(adj_mat, gr = site),
              data = nut,
              data2 = list(adj_mat = adj_mat, spot2021 = nut$site),
              iter = 500,
              chains = 4,
              cores = 4,
              family = "hurdle_gamma")
Sys.time() - start_time


#######################################################
# OLD (as of 5 August)
########################################################


#====Plot spot locations====

#site names
grid_sites %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", data = myvatnshp)+
  # geom_point(aes(easting, northing))+
  geom_text(aes(easting, northing, label = spot2021))+
  coord_sf()+
  guides(fill = "none")+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "left",
        plot.title = element_text(hjust = 0.5))+
  scale_color_viridis_c()+
  NULL

#sampling date
grid_sites %>% 
  ggplot()+
  geom_spatvector(data = myvatnshp)+
  # geom_point(aes(easting, northing))+
  geom_point(aes(easting, northing, color = sampledate))+
  facet_wrap(~year)+
  coord_sf()+
  guides(fill = "none")+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))+

  NULL



grid_sites %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", data = myvatnshp)+
  geom_point(aes(easting, northing))+
  geom_text(aes(lat, lon, label = site, color = site), size = 5, data = esites %>% mutate(site = str_remove(site, "e"), site=  str_remove(site, "st")) %>% filter(!is.na(as.numeric(site))))+
  # geom_text(aes(easting, northing, label = spot))+
  coord_sf()+
  guides(fill = "none")+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))+

  scale_color_brewer(palette= "Dark2")+
  NULL


ggpreview(last_plot(), dpi = 650, width = 3, height = 4, units = "in")


#Qualitative Sediment Characteristics
grid_sites %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  # geom_point(aes(easting, northing))+
  facet_wrap(~lubridate::year(sampledate))+
  geom_point(aes(easting, northing, color = sed_type), size = 3)+
  coord_sf()+
  guides(fill = "none")+
  labs(color = "Sediment Type")+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "left",
        plot.title = element_text(hjust = 0.5))+

  NULL
ggpreview(last_plot(), dpi = 650, width = 5, height = 4, units = "in")

grid_sites %>% 
  group_by(year, sed_type) %>% 
  count()
#how close were sites across years?

yeardists <- grid_sites %>% 
  select(year, easting, northing, spot2021) %>% 
  gather(var, val, -year, -spot2021) %>% 
  mutate(var = paste0(var, year)) %>% 
  select(-year) %>% 
  spread(var, val) %>% 
  mutate(distance = sqrt((easting2022- easting2021)^2 + (northing2022- northing2021)^2))

hist(yeardists$distance)

mean(yeardists$distance)
median(yeardists$distance)

sd(yeardists$distance)

yeardists %>% 
  # filter(distance<100) %>% 
  arrange(desc(distance)) %>% 
  pull(distance)

#pairwise distances
site_dist <- grid_sites %>% select(easting, northing) %>% dist() %>% as.matrix()


row.names(site_dist) <- grid_sites$spot
colnames(site_dist) <- grid_sites$spot



site_dist_long <- site_dist %>% 
  as.data.frame() %>% 
  mutate(spot1 = row.names(.)) %>% 
  gather(spot2, distance, -spot1) %>% 
  mutate(spot1 = as.numeric(spot1),
         spot2 = as.numeric(spot2)) %>% 
  filter(spot1 != spot2)


min(site_dist_long$distance)

mean(site_dist_long$distance)
median(site_dist_long$distance)
max(site_dist_long$distance)

site_dist_long %>% 
  mutate(spot1 = ifelse(spot1>27, spot1-1, spot1),
         spot2 = ifelse(spot2>27, spot2-1, spot2)) %>% 
  ggplot(aes(x = spot1, y = spot2, fill = distance/1000))+
  geom_tile()+
  scale_y_reverse()+
  scale_fill_viridis_c()+
  theme(axis.title = element_blank())+
  labs(fill = "Distance (km)")


#=====Limnological Measurements====
limno <- full_join(grid_sites, profiles) 

#surface temperatures
limno %>% 
  filter(sampledepth == 0.5) %>% 
  group_by(year) %>% 
  summarise(mean= mean(wtemp),
            median = median(wtemp),
            sd = sd(wtemp))


#depth
limno %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  # geom_point(aes(easting, northing))+
  geom_point(aes(easting, northing, color = depth))+
  facet_wrap(~lubridate::year(sampledate))+
  coord_sf()+
  guides(fill = "none")+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "left",
        plot.title = element_text(hjust = 0.5))+

  scale_color_gradient(low = "white", high = "black")+
  NULL
ggpreview(last_plot(), dpi = 650, width = 5, height = 4, units = "in")



limno %>% 
  # group_by(year) %>% 
  summarise(min = min(depth),
            median = median(depth),
            mean = mean_se(depth),
            max = max(depth))


depthcomp <- limno %>% 
  select(spot2021, year, depth) %>% 
  unique() %>% 
  mutate(year = paste0("depth", year)) %>% 
  spread(year, depth)

cor(depthcomp$depth2021, depthcomp$depth2022)^2

lm(depth2022~depth2021, depthcomp) %>% summary()

depthcomp %>% 
  ggplot(aes(x = depth2021, y = depth2022))+
  geom_point()

limno %>% 
  ggplot(aes(depth))+
  geom_histogram(center = TRUE)




#secchi
limno %>% 
  mutate(secchi_dn = as.numeric(secchi_dn),
         secchi_dn = ifelse(is.na(secchi_dn), 4, secchi_dn)) %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  # geom_point(aes(easting, northing))+
  geom_point(aes(easting, northing, color = secchi_dn))+
  facet_wrap(~year(sampledate))+
  coord_sf()+
  guides(fill = "none")+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "left",
        plot.title = element_text(hjust = 0.5))+

  scale_color_viridis_c(direction = -1, na.value = "black")

ggpreview(last_plot(), dpi = 650, width = 5, height = 4, units = "in")

limno %>% 
  mutate(secchi_dn = as.numeric(secchi_dn),
         secchi_dn = ifelse(is.na(secchi_dn), 4, secchi_dn),
         secchi_dn = ifelse(is.infinite(secchi_dn), 4, secchi_dn)) %>% 
  select(year, spot2021, secchi_dn) %>% 
  unique() %>% 
  mutate(secchi = paste0("secchi", year)) %>% 
  select(-year) %>% 
  spread(secchi, secchi_dn) %>% 
  ggplot(aes(x = secchi2021, secchi2022))+
  geom_point()


limno %>% 
  select(spot, secchi_dn) %>% 
  unique() %>% 
  filter(secchi_dn != "INF") %>% 
  mutate(secchi_dn = as.numeric(secchi_dn)) %>% 
  summarise(min = min(secchi_dn),
            median = median(secchi_dn),
            max = max(secchi_dn))


limno %>% 
  select(spot, secchi_dn) %>% 
  unique() %>% 
  filter(secchi_dn != "INF") %>% 
  mutate(secchi_dn = as.numeric(secchi_dn)) %>% 
  ggplot(aes(secchi_dn))+
  geom_histogram()
ggpreview(last_plot(), dpi = 650, width = 2, height = 2, units = "in")


#oxygen at benthos and surface
limno %>% 
  filter(sampledepth == 0) %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  # geom_point(aes(easting, northing))+
  geom_point(aes(easting, northing, color = do), size = 3)+
  coord_sf()+
  guides(fill = "none")+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "left",
        plot.title = element_text(hjust = 0.5))+

  scale_color_viridis_c()

limno %>% 
  group_by(spot) %>% 
  filter(sampledepth == max(sampledepth)) %>% 
  select(spot, sampledate, easting, northing, do) %>% 
  rename(benthic = do) %>% 
  left_join(limno %>% 
              filter(sampledepth == 0.5) %>% 
              select(spot, do)) %>% 
  rename(surface = do) %>% 
  gather(layer, do, benthic, surface) %>% 
  mutate(layer = fct_reorder(layer, c("surface", "benthic"))) %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  # geom_point(aes(easting, northing))+
  geom_point(aes(easting, northing, color = do), size = 3)+
  facet_grid(year(sampledate)~layer)+
  coord_sf()+
  guides(fill = "none")+
  labs(color = "DO mg/L")+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "left",
        plot.title = element_text(hjust = 0.5))+

  scale_color_viridis_c()

ggpreview(last_plot(), dpi = 650, width = 5, height = 4, units = "in")


limno %>% 
  group_by(spot) %>% 
  filter(sampledepth == max(sampledepth)) %>% 
  select(spot, sampledate, easting, northing, do) %>% 
  rename(benthic = do) %>% 
  left_join(limno %>% 
              filter(sampledepth == 0.5) %>% 
              select(spot, do)) %>% 
  rename(surface = do) %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  # geom_point(aes(easting, northing))+
  geom_point(aes(easting, northing, color = surface/benthic), size = 3)+
  facet_grid(~year(sampledate))+
  coord_sf()+
  guides(fill = "none")+
  labs(color = "surface DO/ benthic DO")+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "left",
        plot.title = element_text(hjust = 0.5))+

  scale_color_gradient2(low = "darkgoldenrod", high = "darkgreen", midpoint = 1)
ggpreview(last_plot(), dpi = 650, width = 5, height = 4, units = "in")



limno %>% 
  group_by(spot) %>% 
  filter(sampledepth == max(sampledepth)|
           sampledepth == 0) %>% 
  ungroup %>% 
  mutate(layer = ifelse(sampledepth == 0, "pel", "ben")) %>% 
  select(easting, northing, layer, do) %>% 
  spread(layer, do) %>% 
  mutate(auto_struc = pel/ben) %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  # geom_point(aes(easting, northing))+
  geom_point(aes(easting, northing, color = auto_struc), size = 3)+
  coord_sf()+
  facet_wrap(~year)+
  guides(fill = "none")+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "left",
        plot.title = element_text(hjust = 0.5))+

  scale_color_gradient2(low = "darkgoldenrod", high = "darkgreen", midpoint = 1)+
  labs(color = "surface DO/ benthic DO")


limno %>% 
  group_by(spot) %>% 
  filter(sampledepth == max(sampledepth)|
           sampledepth == 0) %>% 
  ungroup %>% 
  mutate(layer = ifelse(sampledepth == 0, "pel", "ben")) %>% 
  select(easting, northing, layer, do) %>% 
  spread(layer, do) %>% 
  mutate(auto_struc = (pel-ben)/((pel+ben)/2)) %>%
  select(- ben, -pel) %>% 
  mutate(layer = "Concentration") %>% 
  bind_rows(limno %>% 
              group_by(spot) %>% 
              filter(sampledepth == max(sampledepth)|
                       sampledepth == 0) %>% 
              ungroup %>% 
              mutate(layer = ifelse(sampledepth == 0, "pel", "ben")) %>% 
              select(easting, northing, layer, o2sat) %>% 
              spread(layer, o2sat) %>% 
              mutate(auto_struc = (pel-ben)/((pel+ben)/2)) %>%
              select(- ben, -pel) %>% 
              mutate(layer = "% Saturation")) %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  # geom_point(aes(easting, northing))+
  geom_point(aes(easting, northing, color = auto_struc), size = 3)+
  coord_sf()+
  guides(fill = "none")+
  facet_wrap(~layer)+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "left",
        plot.title = element_text(hjust = 0.5))+

  scale_color_gradient2(low = "darkgoldenrod", high = "darkgreen", midpoint = 0)+
  labs(color = "Autotrophic Structure")

ggpreview(last_plot(), dpi = 650, width = 8, height = 4, units = "in")

as2 <- limno %>% 
  group_by(spot) %>% 
  filter(sampledepth == max(sampledepth)|
           sampledepth == 0) %>% 
  ungroup %>% 
  mutate(layer = ifelse(sampledepth == 0, "pel", "ben")) %>% 
  select(easting, northing, layer, do, secchi_dn) %>% 
  spread(layer, do) %>% 
  mutate(auto_struc = (pel-ben)/ben) %>%
  select(- ben, -pel) %>% 
  mutate(layer = "Concentration") %>% 
  bind_rows(limno %>% 
              group_by(spot) %>% 
              filter(sampledepth == max(sampledepth)|
                       sampledepth == 0) %>% 
              ungroup %>% 
              mutate(layer = ifelse(sampledepth == 0, "pel", "ben")) %>% 
              select(easting, northing, layer, o2sat, secchi_dn) %>% 
              spread(layer, o2sat) %>% 
              mutate(auto_struc = (pel-ben)/ben) %>%
              select(- ben, -pel) %>% 
              mutate(layer = "% Saturation")) %>% 
  mutate(secchi_dn =  as.numeric(ifelse(secchi_dn == "INF", 4, secchi_dn))) %>% 
  ggplot(aes(x= secchi_dn, y= auto_struc))+
  facet_wrap(~layer)+
  geom_hline(yintercept = 0)+
  geom_point()+
  theme_bw()+
  labs(x = "Secchi Depth (m)",
       y = "Autotrophic\nStructure")+
  theme(panel.grid = element_line(color = "white"),
        legend.position = "left",
        plot.title = element_text(hjust = 0.5))

as2_trimmed <- limno %>% 
  group_by(spot) %>% 
  filter(sampledepth == max(sampledepth)|
           sampledepth == 0) %>% 
  ungroup %>% 
  mutate(layer = ifelse(sampledepth == 0, "pel", "ben")) %>% 
  select(easting, northing, layer, do, secchi_dn) %>% 
  spread(layer, do) %>% 
  mutate(auto_struc = (pel-ben)/ben) %>%
  select(- ben, -pel) %>% 
  mutate(layer = "Concentration") %>% 
  bind_rows(limno %>% 
              group_by(spot) %>% 
              filter(sampledepth == max(sampledepth)|
                       sampledepth == 0) %>% 
              ungroup %>% 
              mutate(layer = ifelse(sampledepth == 0, "pel", "ben")) %>% 
              select(easting, northing, layer, o2sat, secchi_dn) %>% 
              spread(layer, o2sat) %>% 
              mutate(auto_struc = (pel-ben)/ben) %>%
              select(- ben, -pel) %>% 
              mutate(layer = "% Saturation")) %>% 
  mutate(secchi_dn =  as.numeric(ifelse(secchi_dn == "INF", 4, secchi_dn))) %>% 
  filter(secchi_dn>1, secchi_dn<3) %>% 
  ggplot(aes(x= secchi_dn, y= auto_struc))+
  geom_hline(yintercept = 0)+
  facet_wrap(~layer)+
  geom_point()+
  theme_bw()+
  labs(x = "Secchi Depth (m)",
       y = "Autotrophic\nStructure")+
  theme(panel.grid = element_line(color = "white"),
        legend.position = "left",
        plot.title = element_text(hjust = 0.5))


plot_grid(as2, as2_trimmed, ncol = 1)
ggpreview(last_plot(), dpi = 650, width = 4, height = 4, units = "in")


limno %>% 
  group_by(spot) %>% 
  filter(sampledepth == max(sampledepth)|
           sampledepth == 0) %>% 
  ungroup %>% 
  mutate(layer = ifelse(sampledepth == 0, "pel", "ben")) %>% 
  select(easting, northing, layer, do, secchi_dn) %>% 
  spread(layer, do) %>% 
  mutate(auto_struc = (pel-ben)/ben) %>%
  select(- ben, -pel) %>% 
  mutate(layer = "Concentration") %>% 
  bind_rows(limno %>% 
              group_by(spot) %>% 
              filter(sampledepth == max(sampledepth)|
                       sampledepth == 0) %>% 
              ungroup %>% 
              mutate(layer = ifelse(sampledepth == 0, "pel", "ben")) %>% 
              select(easting, northing, layer, o2sat, secchi_dn) %>% 
              spread(layer, o2sat) %>% 
              mutate(auto_struc = (pel-ben)/ben) %>%
              select(- ben, -pel) %>% 
              mutate(layer = "% Saturation")) %>% 
  mutate(secchi_dn =  as.numeric(ifelse(secchi_dn == "INF", 4, secchi_dn))) %>% 
  filter(secchi_dn>1, secchi_dn<3) %>% 
  lm(auto_struc~secchi_dn, data = .) %>% summary()

limno %>% 
  group_by(spot) %>% 
  filter(sampledepth == max(sampledepth)|
           sampledepth == 0.5) %>% 
  ungroup %>% 
  mutate(layer = ifelse(sampledepth == 0.5, "pel", "ben")) %>% 
  select(year, easting, northing, layer, o2sat) %>% 
  spread(layer, o2sat) %>% 
  mutate(auto_struc = pel/ben) %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  # geom_point(aes(easting, northing))+
  geom_point(aes(easting, northing, color = auto_struc), size = 3)+
  coord_sf()+
  facet_wrap(~year)+
  guides(fill = "none")+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "left",
        plot.title = element_text(hjust = 0.5))+

  scale_color_gradient2(low = "darkgoldenrod", high = "darkgreen", midpoint = 1)+
  labs(color = "surface %DO/ benthic %DO")


limno %>% 
  gather(var, val, wtemp, do, o2sat) %>% 
  mutate(var = ifelse(var == "wtemp", "Temperature °C",
                      ifelse(var == "do", "DO mg/L", "% O2 Saturation"))) %>% 
  ggplot(aes(y = sampledepth, x = val, group = spot))+
  facet_grid(year~var, scales = "free_x")+
  geom_path(alpha = 0.7)+
  labs(y = "Depth (m)")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = "left",
        plot.title = element_text(hjust = 0.5),
        strip.placement = "outside",
        strip.background = element_blank(),
        axis.title.x = element_blank())+
  scale_y_reverse()

ggpreview(last_plot(), dpi = 650, width = 8, height = 4, units = "in")


  ggplot(aes(y = sampledepth, x = wtemp, group = spot))+
  geom_path(alpha = 0.7)+
  scale_y_reverse()

limno %>% 
  ggplot(aes(y = sampledepth, x = do, group = spot))+
  geom_path(alpha = 0.7)+
  scale_y_reverse()

limno %>% 
  ggplot(aes(y = sampledepth, x = o2sat, group = spot))+
  geom_path(alpha = 0.7)+
  geom_vline(xintercept = 100)+
  scale_y_reverse()

#====Plot nutrient data==== 

nut <- full_join(nutrients %>% 
                   mutate(analyte2 = ifelse(analyte == "nh4_nh3", "nh4", analyte)) %>% 
                   group_by(method, analyte2, layer) %>% 
                   mutate(z = scale(mgl)[,1],
                          mgLlog10 = log10(mgl + min(mgl)/2)) %>% 
                   ungroup,
                 grid_sites %>% mutate(year = as.factor(year))) 

bothyears <- c("nh4", "no3", "po4")

# summaries
summary(nut %>% filter(analyte2 %in% bothyears) %>% select(spot, analyte, year, layer, mgl) %>% 
          mutate(year_analyte = paste(analyte, year, layer)) %>% 
          select(spot, year_analyte, mgl) %>% 
          spread(year_analyte, mgl)) 

nut %>% 
  filter(analyte2 %in% bothyears) %>% 
  ggplot(aes(x = mgl)) + facet_wrap(layer~analyte, scales = "free") + geom_histogram() + theme_bw()

# how do years differ in analyte concentrations

lm(z~analyte2*year*layer, 
   data = nut %>%  filter(analyte2 %in% bothyears,
                                             paste(analyte2, layer)!="nh4 pel")) %>% 
  plot()

nut %>% 
  filter(analyte2 %in% bothyears,
         paste(analyte2, layer)!="nh4 pel") %>% 
  group_by(analyte2, layer) %>% 
  nest() %>% 
  mutate(fit = map(data, ~lm(mgl~year, data = .x)),
         out = map(fit, ~broom::tidy(car::Anova(.x)))) %>% 
  unnest(out)

nut %>% 
  filter(analyte2 %in% bothyears) %>% 
  select(analyte2, spot2021, layer, year, z) %>% 
  spread(year, z) %>% 
  ggplot(aes(x = `2021`, y = `2022`))+
  facet_wrap(~paste(layer,analyte2), scales = "free") + 
  geom_point()+
  theme_bw()

nut %>% 
  filter(analyte2 %in% bothyears) %>% 
  select(analyte2, spot2021, layer, year, mgLlog10) %>% 
  spread(year, mgLlog10) %>% 
  ggplot(aes(x = `2021`, y = `2022`))+
  facet_wrap(~paste(layer,analyte2), scales = "free") + 
  geom_point()+
  theme_bw()

# how correlated are nutrient concentratoins
nut %>% 
  filter(analyte2 %in% bothyears) %>% 
  select(analyte2, spot2021, layer, year, mgl) %>% 
  spread(layer, mgl) %>% 
  ggplot(aes(x = ben, y = pel, color = factor(year)))+
  facet_wrap(~analyte2, scales = "free") + 
  geom_point()+
  geom_abline(slope = 1) + 
  theme_bw()

# plot concentrations
nut %>% 
  filter(analyte2 %in% bothyears) %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  geom_point(aes(easting, northing, color = mgl), size = 2)+
  coord_sf()+
  facet_grid(paste(layer, year)~analyte2)+
  guides(fill = "none")+
  labs(color = "mg/L")+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "left",
        plot.title = element_text(hjust = 0.5)) +
  scale_color_viridis_c(trans = "log10", na.value = "gray20")

# plot scaled concentrations
nut %>% 
  filter(analyte2 %in% bothyears) %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  geom_point(aes(easting, northing, color = z), size = 2)+
  coord_sf()+
  facet_grid(paste(layer, year)~analyte2)+
  guides(fill = "none")+
  labs(color = "concentration\nZ-scored")+
  scale_color_viridis_c()

# plot wider
nz21 <- nut %>% 
  filter(analyte2 %in% bothyears) %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  geom_point(aes(easting, northing, color = z),  data = . %>% filter(year == 2021))+
  coord_sf()+
  facet_grid(fct_rev(layer)~analyte2, switch = "y")+
  guides(fill = "none")+
  labs(color = "concentration\nZ-scored",
       title = "2021")+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "top",
        legend.key.height = unit(0.5, "lines"),
        legend.title = element_text(vjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  scale_color_viridis_c(limits = c(min(nut$z, na.rm = T), max(nut$z, na.rm = T)))

nz22 <- nut %>% 
  filter(analyte2 %in% bothyears) %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  geom_point(aes(easting, northing, color = z),  data = . %>% filter(year == 2022))+
  coord_sf()+
  facet_grid(fct_rev(layer)~analyte2)+
  guides(fill = "none")+
  labs(color = "concentration\nZ-scored",
       title = "2022")+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom",
        legend.key.height = unit(0.5, "lines"),
        legend.title = element_text(vjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  scale_color_viridis_c(limits = c(min(nut$z, na.rm = T), max(nut$z, na.rm = T)))

nz22.2 <- nut %>% 
  filter(analyte2 %in% bothyears) %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  geom_point(aes(easting, northing, color = z),  data = . %>% filter(year == 2022))+
  coord_sf()+
  facet_grid(fct_rev(layer)~analyte2, switch = "y")+
  guides(fill = "none")+
  labs(color = "concentration\nZ-scored",
       title = "2022")+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom",
        legend.key.height = unit(0.5, "lines"),
        legend.title = element_text(vjust = 1),        plot.title = element_text(hjust = 0.5)) +
  scale_color_viridis_c(limits = c(min(nut$z, na.rm = T), max(nut$z, na.rm = T)))


plot_grid(nz21, nz22)

long <- plot_grid(nz21, nz22.2, ncol = 1)

ggsave(plot = long, width = 84, height= 200, units = "mm", dpi = 650, filename = "spat_nuts_fig2.jpeg")

plot_grid(long, nzleg, rel_heights = c(4, 1), ncol = 1)

nut %>% 
  select(-vial, -analyte) %>% 
  spread(analyte2, mgl) %>% 
  filter(layer == "ben") %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  # geom_point(aes(easting, northing))+
  geom_point(aes(easting, northing, color = nh4), size = 2)+
  facet_wrap(~year) +
  coord_sf()+
  guides(fill = "none")+
  labs(color = "NH4 NH3 mg/L",
       shape = "benthos class")+
  scale_shape_manual(values = c(21, 24, 16))+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "top",
        plot.title = element_text(hjust = 0.5))+

  scale_color_viridis_c()+
  ggspatial::annotation_scale()+
  ggspatial::annotation_north_arrow(style = ggspatial::north_arrow_orienteering, pad_y = unit(1.2, "lines"), height = unit(1, "lines"))

ggpreview(plot = last_plot(), width = 3, height = 5, units = "in")

nut %>% 
  group_by(analyte, layer) %>% 
  mutate(z = scale(mgl)[,1]) %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  # geom_point(aes(easting, northing))+
  geom_point(aes(easting, northing, color = z), size = 2)+
  coord_sf()+
  facet_grid(layer~analyte)+
  guides(fill = "none")+
  labs(color = "Z Score")+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "left",
        plot.title = element_text(hjust = 0.5))+

  scale_color_viridis_c()
# ggpreview(plot = last_plot(), width = 8, height = 6, units= "in", dpi = 650)


nutwide <- nut %>% 
  filter(analyte2 %in% bothyears) %>% 
  mutate(eaststd = (easting-min(easting))/1000,
        northstd = (northing-min(northing))/1000)%>% 
  select(year, spot,northing, easting, eaststd, northstd, layer, sed_type, analyte2, mgl) %>% 
  spread(analyte2, mgl) 


nutwide2 <- nut %>% 
  filter(analyte2 %in% bothyears) %>% 
  mutate(eaststd = (easting-min(easting))/1000,
         northstd = (northing-min(northing))/1000)%>% 
  select(year, spot2021, northing, easting, eaststd, northstd, layer, sed_type, analyte2, mgl) %>% 
  spread(analyte2, mgl) 


nutwide %>% 
  gather(chem, mgl, c(nh4:po4)) %>% 
  ggplot(aes(mgl, fill = factor(year)))+
  facet_wrap(layer~chem, scales = "free", nrow = 2)+
  geom_histogram()



nut %>% 
  filter(analyte2 %in% bothyears) %>% 
  left_join(limno %>% select(spot, secchi_dn) %>% mutate(secchi_dn = ifelse(secchi_dn == "INF", Inf, as.numeric(secchi_dn))) %>% unique()) %>% 
  ggplot(aes(x = secchi_dn, y = mgl, col = factor(year))) + 
  geom_point(size = 3) + 
  facet_grid(layer~analyte2, scales = "free_y") +
  theme_bw()


####Newer stats attempts####
# set up distance matrix
site.dist <- dist(cbind(nutwide2$northing, nutwide2$easting), diag = TRUE, upper = TRUE)/1000

brm(mvbind(nh4, no3, po4)~year*layer,
    family = "gaussian",
    data = nutwide2, 
    iter = 10,
    chains = 4,
    cores = 4)

# question :ben pel coupling

# ben_sub <- bf(ben_nh4~midges+depth+sed_type + year , autocorrelation) #if wanting a different set of predictors

# create adjacency matrix
# car(Adjacency, gr = site)

# maybe try 1 or 10 to check whether there is an error
# start with 500 to see if its working if it will converge
# 2000 iters at least for the full model

# this is the model
brm(mvn(ben_nh4, pel_nh4)~year + phc + year + depth + eaststd + nothstd + chironominae) # gamma, hurdle_gamma (basically 0inflated gamma), or zinfbinomial




#####OLD STATS####

library(nlme)

#===

nocor <- gls(nh4~1,
             data = nutwide %>% filter(layer == "ben"),
             method = "REML")


#===preliminary analyses of nutrient concentrations==
#pull benthic
nwb <- nutwide %>% 
  filter(layer == "ben")

#pull pelagic
nwp <- nutwide %>% 
  filter(layer == "pel")

bothyears


nutrient_resp_models <- function(data, rhs_formula){
  #effects of sediment type
  analyte_cols <- which(colnames(data) %in% bothyears)
  
  results <- list()
  
  for(i in analyte_cols){
    response <- colnames(data[,i])
    predictor <- rhs_formula
    formula <- formula(paste(response, predictor, sep = "~"))
    results[[i]] <- gls(formula, 
                                correlation = corExp(form = ~northing+easting, nugget = TRUE),
                                method = "REML",
                                data = data)
    names(results)[i] <- response
  }
  
  return(results[lapply(results, length)>0])
}
sedtypeb <- nutrient_resp_models(nwb, "sed_type * factor(year)")

sedtypeb

lapply(sedtypeb, summary)
lapply(sedtypeb, plot)


lapply(sedtypeb, function(x) anova(x))
nut

lapply(nutrient_resp_models(nwb %>% mutate(sed_type2 = ifelse(sed_type == "cladophorales", "cladophora", "sediment")), "sed_type2"), summary)


lm(nh4~(eaststd + northstd + factor(year))^2 + sed_type, data = nwb) %>% summary()
lm(nh4~(eaststd + northstd + factor(year))^2 + sed_type, data = nwp) %>% summary()

lm(po4~(eaststd + northstd + factor(year))^2 + sed_type, data = nwb) %>% summary()
lm(po4~(eaststd + northstd + factor(year))^2 + sed_type, data = nwp) %>% summary()


pelcheck <- nwp %>% 
  # select(spot, easting, northing, nh4_nh3, tn, po4) %>% 
  left_join(depthcomp %>% 
              select(spot=spot2021, depth2021)) %>% 
  mutate(depth = floor(depth2021-0.05),
         depth_unadj = floor(depth2021),
         depth_diff = depth2021 - floor(depth2021),
         depth_diff2 = depth2021 - depth)


pelcheck %>% 
  gather(var, val, nh4_nh3:tp) %>% 
  ggplot(aes(x = depth, y = val, color = factor(depth_unadj)))+
  facet_wrap(~var, scales = "free_y")+
  geom_jitter(width = 0.2, alpha = 0.5)

pelcheck %>% 
  gather(var, val, nh4_nh3:tp) %>% 
  ggplot(aes(x = depth2021, y = val))+
  facet_wrap(~var, scales = "free_y")+
  geom_point(alpha = 0.5)

pelcheck %>% 
  gather(var, val, nh4_nh3:tp) %>% 
  ggplot(aes(x = depth_diff, y = val, color = factor(depth_unadj)))+
  facet_wrap(~var, scales = "free_y")+
  geom_point(alpha = 0.5)

pelcheck %>% 
  gather(var, val, nh4_nh3:tp) %>% 
  ggplot(aes(x = depth_diff2, y = val, color = factor(depth_unadj)))+
  facet_wrap(~var, scales = "free_y")+
  geom_point(alpha = 0.5)


# lapply(nutrient_resp_models(pelcheck, "depth"), summary)
    
gls(nh4_nh3 ~ depth,

    data = pelcheck, 
    cor = corExp(form = ~easting + northing, nugget= TRUE)) %>% summary()

gls(tn ~ depth, 
    data = pelcheck, 
    cor = corExp(form = ~easting + northing, nugget = TRUE)) %>% summary()

gls(no2 ~ depth, 
    data = pelcheck, 
    cor = corExp(form = ~easting + northing, nugget = TRUE)) %>% summary()


#getting false convergence with nugget (nugget probably really small)
gls(po4 ~ depth2, 
    data = pelcheck, 
    cor = corExp(form = ~easting + northing)) %>% summary()


gls(tp ~ depth, 
    data = pelcheck, 
    cor = corExp(form = ~easting + northing)) %>% summary()


#### correlations in nutrients
cor(nutwide[fulllist])

nutcor1 <- data.frame(cor(nwb[which(colnames(nwb) %in% c(nutnames$chem, "np"))], method = "spearman"), layer = "Interstitial")

nutcor1$v1 = row.names(nutcor1)

nutcor2 <- data.frame(cor(nwp[which(colnames(nwp) %in% c(nutnames$chem, "np"))], method = "spearman"), layer = "Pelagic")

nutcor2$v1 = row.names(nutcor1)

full_join(nutcor1, nutcor2) %>% 
  gather(v2, cor, -v1, -layer) %>% 
  filter(v1!=v2) %>% 
  mutate(v2 = fct_rev(v2)) %>% 
  ggplot(aes(x = v1, y = v2, fill = cor))+
  facet_wrap(~layer)+
  geom_tile()+
  geom_text(aes(label = round(cor, 2)))+
  labs(x = element_blank(),
       y = element_blank(),
       fill  = "Spearman\nCorrelation")+
  coord_equal()+
  scale_fill_gradient2(low = "dodgerblue", high = "firebrick4", mid = "white", midpoint= 0, limits = c(-1,1) )

ggpreview(plot = last_plot(), width = 7, height= 4, dpi = 650, units = "in")

####relating benthic and pelagic nutrients


nutwide %>%
  gather(chem, mgl, nh4_nh3:dop) %>% 
  spread(layer, mgl) %>% 
  ggplot(aes(x = Ben, y = Pel, color = sed_type))+
  facet_wrap(~chem, scales = "free")+
  geom_point()

fulllist


pelbennut <- nutwide %>%
  gather(chem, mgl, nh4_nh3:np) %>% 
  spread(layer, mgl) %>% 
  mutate(pb = Pel/Ben,
         chem = fct_relevel(chem, levels = fulllist)) 


pelbennut%>% 
  ggplot(aes(pb))+
  geom_histogram()+
  facet_wrap(~chem, scales = "free", nrow = 2)+
  geom_vline(xintercept = 1)


plotlist <- list()
for(i in 1:length(fulllist)){
  plotlist[[i]] <- nutwide %>%
    gather(chem, mgl, nh4_nh3:np) %>% 
    spread(layer, mgl) %>% 
    mutate(pb = Pel/Ben) %>% 
    filter(chem == fulllist[i]) %>% 
    ggplot()+
    geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
    # geom_point(aes(easting, northing))+
    geom_point(aes(easting, northing, color = pb), size = 2)+
    coord_sf()+
    facet_wrap(~chem, nrow = 2)+
    guides(fill = "none")+
    labs(color = element_blank())+
    theme_bw()+
    theme(panel.grid = element_line(color = "white"),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "left",
          plot.title = element_text(hjust = 0.5))+
  
    scale_color_gradient2(low = "darkgoldenrod", high = "darkgreen", midpoint = 1, n.breaks = 4)+
    theme(legend.position = c(0.35, 0.9),
          legend.direction = "horizontal",
          legend.background = element_blank(),
          legend.key.height = unit(0.5, "lines"),
          legend.key.width = unit(0.8, "lines"))
  
}

plot_grid(plotlist = plotlist, nrow =2)

ggpreview(plot = last_plot(), width = 8, height = 4, dpi = 650, units = "in")

plotlist <- list()
for(i in 1:length(fulllist)){
  plotlist[[i]] <- nutwide %>%
    gather(chem, mgl, nh4_nh3:np) %>% 
    spread(layer, mgl) %>% 
    mutate(pb = Pel/Ben) %>% 
    filter(chem == fulllist[i]) %>% 
    ggplot()+
    geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
    # geom_point(aes(easting, northing))+
    geom_point(aes(easting, northing, color = pb, shape = pb>1), size = 2)+
    coord_sf()+
    facet_wrap(~chem, nrow = 2)+
    guides(fill = "none",
           shape = "none")+
    labs(color = element_blank())+
    theme_bw()+
    theme(panel.grid = element_line(color = "white"),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "left",
          plot.title = element_text(hjust = 0.5))+
  
    scale_color_viridis_c(option = "magma", n.breaks = 3)+
    theme(legend.position = c(0.3, 0.9),
          legend.direction = "horizontal",
          legend.background = element_blank(),
          legend.key.height = unit(0.5, "lines"),
          legend.key.width = unit(0.7, "lines"))
  
}

plot_grid(plotlist = plotlist, nrow =2)

ggpreview(plot = last_plot(), width = 8, height = 4, dpi = 650, units = "in")


####looking at phytoplankton composition

nutwide %>% 
  filter(layer == "Pel") %>% 
  gather(analyte, concentration, nh4_nh3:np) %>% 
  left_join(turner %>% 
              group_by(spot) %>% 
              summarise(chl = mean(chl),
                        phyc = mean(phyc))) %>% 
  mutate(chlphyc = chl/phyc) %>%
  ggplot(aes(x = concentration, y = chlphyc))+
  labs(x = "Concentration (mg/L)",
       y = "Pelagic Chlorophyll:Phycocyanin")+
  facet_wrap(~analyte, scales = "free", nrow = 2)+
  geom_point()+
  theme_bw()

nutwide %>% 
  filter(layer == "Ben") %>% 
  gather(analyte, concentration, nh4_nh3:np) %>% 
  left_join(turner %>% 
              group_by(spot) %>% 
              summarise(chl = mean(chl),
                        phyc = mean(phyc))) %>% 
  mutate(chlphyc = chl/phyc) %>%
  ggplot(aes(x = concentration, y = chlphyc))+
  labs(x = "Concentration (mg/L)",
       y = "Pelagic Chlorophyll:Phycocyanin")+
  facet_wrap(~analyte, scales = "free", nrow = 2)+
  geom_point()+
  theme_bw()




pelcomp1 <- nwp %>% 
  left_join(turner %>% 
                    group_by(spot) %>% 
                    summarise(chl = mean(chl),
                              phyc = mean(phyc))) %>% 
                    mutate(chlphyc = chl/phyc)


nutrient_exp_model <- function(data, lhs){
  analyte_cols <- which(colnames(data) %in% fulllist)
  
  results <- list()
  
  for(i in analyte_cols){
    predictor <- colnames(data[,i])
    response <- lhs
    formula <- formula(paste(response, predictor, sep = "~"))
    results[[i]] <- gls(formula, 
                                correlation = corExp(form = ~northing+easting, nugget = TRUE),
                                method = "REML",
                                data = data)
    names(results)[i] <- predictor
  }
  return(results[lapply(results, length)>0])
}

pelcomp <- nutrient_exp_model(pelcomp1, "chlphyc")

lapply(pelcomp, summary)

#####looking at epipelic algal composition

bencomp1 <- nwb %>% 
  filter(layer == "Ben") %>% 
  left_join(bt %>% filter(year == 2021)) %>% 
  mutate(cdr = sed_cyano/sed_diatoms)

bencomp1 %>% 
  gather(analyte, concentration, nh4_nh3:np) %>% 
  ggplot(aes(x = concentration, y = cdr))+
  labs(x = "Interstitial Concentration (mg/L)",
       y = "Benthic Cyanobacteria:Diatoms")+
  facet_wrap(~analyte, scales = "free", nrow = 2)+
  geom_point()+
  theme_bw()
nutrient_exp_model(data = bencomp1, lhs = "cdr")


#####looking at spatial gradients

#KEY TO EAST AND NORTH
nwb %>% 
  select(northing, easting, eaststd, northstd) %>% 
  gather(var, val, eaststd, northstd) %>% 
  ggplot()+
  facet_wrap(~var)+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  geom_point(aes(easting, northing, color = val), size = 2)+
  coord_sf()+
  guides(fill = "none",
         shape = "none")+
  labs(color = element_blank())+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "left",
        plot.title = element_text(hjust = 0.5))+

  scale_color_gradient(low = "black", high = "white")


gradb <- nutrient_resp_models(nwb,"eaststd+northstd")
gradb

extract_coefficients_gls <- function(outlist, name = NA){
  modcoef <- bind_rows(lapply(outlist, function(x) data.frame(coef(summary(x)))), .id = "analyte") %>% 
    mutate(term = row.names(.), 
           term = str_remove(term, "\\("), 
           term = str_remove(term, "\\)")) %>% 
    separate(term, c("term", "mod"))
  

  
  range = t(bind_rows(lapply(outlist, function(x) exp(x$modelStruc$corStruct[1])), .id = "analyte"))
  nugget = t(bind_rows(lapply(outlist, function(x) exp(x$modelStruc$corStruct[2])), .id = "analyte"))
  
  r <- data.frame(range)
  n <- data.frame(nugget)
  
  r$analyte = row.names(r)
  n$analyte = row.names(n)
  
  out <- modcoef %>% 
    janitor::clean_names() %>%
    left_join(r) %>% 
    left_join(n) %>% 
    mutate(response = name) %>% 
    select(-mod) %>% 
    select(response, 1, range, nugget, term, everything())
  
  return(out)
}

gradbcoef <- extract_coefficients_gls(gradb, "benthic")

gradp <- nutrient_resp_models(nwp,"eaststd+northstd")
gradp

gradpcoef <- extract_coefficients_gls(gradp, "pelagic")

pbready <- pelbennut %>% 
  filter(!chem %in% c("don", "dop", "no3")) %>% 
  select(-Ben, -Pel) %>% 
  spread(chem, pb)

pbcomp <- nutrient_resp_models(pbready, "eaststd+northstd")

gradpbcoef <- extract_coefficients_gls(pbcomp, "pelben")


full_join(gradbcoef, gradpcoef) %>% 
  full_join(gradpbcoef) %>% 
  filter(term!="Intercept",
         p_value<0.05)

ranplot <- full_join(gradbcoef, gradpcoef) %>% 
  full_join(gradpbcoef) %>% 
  select(analyte, response, range) %>% 
  unique() %>% 
  ggplot(aes(range))+
  facet_wrap(~response, scales = "free")+
  geom_histogram(aes(fill = analyte))+
  geom_vline(xintercept = 1000)+
  scale_x_log10()

nugplot <- full_join(gradbcoef, gradpcoef) %>% 
  full_join(gradpbcoef) %>% 
  select(analyte, response, nugget) %>% 
  unique() %>% 
  ggplot(aes(nugget))+
  facet_wrap(~response, scales = "free")+
  geom_histogram(aes(fill = analyte))+
  scale_x_log10()

plot_grid(nugplot, ranplot)


gb <- gradbcoef %>% 
  select(-response) %>% 
  select(analyte, term, everything())

colnames(gb)[3:8] <- paste("benthic", colnames(gb[,3:8]), sep = "_")

gp <- gradpcoef %>% 
  select(-response) %>% 
  select(analyte, term, everything())

colnames(gp)[3:8] <- paste("pelagic", colnames(gp[,3:8]), sep = "_")

gpb <- gradpbcoef %>% 
  select(-response) %>% 
  select(analyte, term, everything())


colnames(gpb)[3:8] <- paste("pelben", colnames(gpb[,3:8]), sep = "_")


# full_join(gb, gp) %>% 
#   full_join(gpb) %>% 
#   write_csv(., "gradient_gls_results.csv")



nwb %>% 
  gather(var, val, nh4_nh3:np) %>% 
  ggplot(aes(x = eaststd, y = val, color = northstd))+
  geom_point()+
  facet_wrap(~var, scale = "free_y")+
  scale_color_viridis_c(option = "plasma")

nwp %>% 
  gather(var, val, nh4_nh3:np) %>% 
  ggplot(aes(x = eaststd, y = val, color = northstd))+
  geom_point()+
  facet_wrap(~var, scale = "free_y")+
  scale_color_viridis_c(option = "plasma")


#####looking at midges and nutrients
midgenut <- nwb %>% left_join(cc)


midgenut %>% 
  gather(var, val, nh4_nh3:np) %>% 
  ggplot(aes(x = chiro, y = val))+
  facet_wrap(~var, scales = "free_y")+
  geom_point(alpha= 0.5, size = 2)

chironut <- nutrient_resp_models(midgenut, "chiro")

lapply(chironut, summary)
lapply(chrionut, function(x) plot(x, main = x))
extract_coefficients_gls(chironut, "midge") %>% 
  filter(term!="Intercept") %>% 
  tibble() %>% 
  ggplot(aes(x = analyte, y = value, ymin = value-std_error, ymax = value+std_error, color = p_value<0.05))+
  geom_hline(yintercept = 0)+
  geom_pointrange()+
  labs(x = element_blank(),
       y = "Effect of Chironomini")+
  theme(legend.position = "bottom")

#####older analysis#
#don
donut <- nut %>%
  filter(spot<32) %>% 
  select(spot, layer, chem, mgl, sed_type, northing, easting) %>% 
  mutate(mgl = ifelse(is.na(mgl), 0, mgl)) %>% 
  spread(chem, mgl) %>% 
  mutate(tin = nh4_nh3 + no2 + no3,
         don = tn - nh4_nh3 - no2 - no3,
         dop = tp - po4)


donut %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  # geom_point(aes(easting, northing))+
  geom_point(aes(easting, northing, color = don), size = 2)+
  coord_sf()+
  facet_grid(~layer)+
  guides(fill = "none")+
  labs(color = "DON mg/L")+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "left",
        plot.title = element_text(hjust = 0.5))+

  scale_color_viridis_c()


donut %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  # geom_point(aes(easting, northing))+
  geom_point(aes(easting, northing, color = tin/(po4+0.01)), size = 2)+
  coord_sf()+
  facet_grid(~layer)+
  guides(fill = "none")+
  labs(color = "tin:tip")+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "left",
        plot.title = element_text(hjust = 0.5))+

  scale_color_viridis_c()

nprat <- donut %>% 
  filter(layer == "Pel") %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  # geom_point(aes(easting, northing))+
  geom_point(aes(easting, northing, color = nh4_nh3/(po4)), size = 2)+
  coord_sf()+
  guides(fill = "none")+
  labs(color = "nh4:po4")+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "left",
        plot.title = element_text(hjust = 0.5))+

  scale_color_viridis_c(trans = "log10")


donut %>% 
  filter(layer == "Pel") %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  # geom_point(aes(easting, northing))+
  geom_point(aes(easting, northing, color = dop), size = 2)+
  coord_sf()+
  facet_grid(~layer)+
  guides(fill = "none")+
  labs(color = "DOP mg/L")+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "left",
        plot.title = element_text(hjust = 0.5))+

  scale_color_viridis_c()


peldon <- donut %>% 
  filter(layer == "Pel") %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  # geom_point(aes(easting, northing))+
  geom_point(aes(easting, northing, color = don), size = 2)+
  coord_sf()+
  guides(fill = "none")+
  labs(color = "Pel. DON")+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "left",
        plot.title = element_text(hjust = 0.5))+

  scale_color_viridis_c()


peldop <- donut %>% 
  filter(layer == "Pel") %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  # geom_point(aes(easting, northing))+
  geom_point(aes(easting, northing, color = dop), size = 2)+
  coord_sf()+
  guides(fill = "none")+
  labs(color = "Pel. DOP")+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "left",
        plot.title = element_text(hjust = 0.5))+

  scale_color_viridis_c()

plot_grid(peldon, peldop)

pelphyc <- turner %>% 
  group_by(spot) %>% 
  summarise(chl=  mean(chl), 
            phyc =mean(phyc))%>% 
  left_join(grid_sites) %>% 
  filter(year(sampledate)==2021) %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  # geom_point(aes(easting, northing))+
  geom_point(aes(easting, northing, color = phyc), size = 2)+
  coord_sf()+
  guides(fill = "none")+
  labs(color = "phyc")+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "left",
        plot.title = element_text(hjust = 0.5))+

  scale_color_viridis_c()


plot_grid(peldon, pelphyc)

plot_grid(nprat, pelphyc)



donut %>% 
  left_join(turner %>% group_by(spot) %>% summarise(chl=  mean(chl), phyc =mean(phyc))) %>% 
  ggplot(aes(x = phyc, y = don))+
  facet_wrap(~layer, scales = "free")+
  geom_point()+
  theme_bw()


donut %>% 
  left_join(turner %>% group_by(spot) %>% summarise(chl=  mean(chl), phyc =mean(phyc))) %>% 
  ggplot(aes(x = phyc, y = nh4_nh3/po4))+
  facet_wrap(~layer, scales = "free")+
  geom_point()+
  theme_bw()


donut %>% 
  left_join(turner %>% group_by(spot) %>% summarise(chl=  mean(chl), phyc =mean(phyc))) %>% 
  ggplot(aes(x = phyc, y = (nh4_nh3+no3+no2)/po4))+
  facet_wrap(~layer, scales = "free")+
  geom_point()+
  theme_bw()

donut %>% 
  filter(layer == "Ben") %>% 
  left_join(turner %>% group_by(spot) %>% summarise(chl=  mean(chl), phyc =mean(phyc))) %>% 
  ggplot(aes(x = po4, y = phyc))+
  geom_point(size = 3)

donut %>% 
  filter(layer == "Ben") %>% 
  left_join(turner %>% group_by(spot) %>% summarise(chl=  mean(chl), phyc =mean(phyc))) %>% 
  ggplot(aes(x = don, y = phyc))+
  geom_point(size = 3)


donut %>% 
  select(spot, northing, easting, sed_type, layer, tp) %>% 
  spread(layer, tp) %>% 
  ggplot(aes(x = Ben, y = Pel))+
  geom_point()

donut %>% 
  filter(layer == "Pel") %>% 
  lm(don~dop, data = .) %>% summary()

leakyp1 <-  
  gls(don~phyc,
      data = donut %>% 
        filter(layer == "Pel") %>% 
        left_join(turner %>% group_by(spot) %>% summarise(chl=  mean(chl), phyc =mean(phyc))),
      method = "REML")

donut %>% 
  filter(layer == "Pel") %>% 
  left_join(turner %>% group_by(spot) %>% summarise(chl=  mean(chl), phyc =mean(phyc))) %>%
  ggplot(aes(x = phyc, y = don))+
  geom_point()+
  labs(x = "Phycocyanin (Turner Units)",
       y = "DON (mg/L)")+
  geom_abline(intercept = leakyp1$coefficients[1], slope= leakyp1$coefficients[2])

ggpreview(plot = last_plot(), width = 3, height = 3, units = "in", dpi = 650)

lm(don~phyc,
    data = donut %>% 
      filter(layer == "Pel") %>% 
      left_join(turner %>% group_by(spot) %>% summarise(chl=  mean(chl), phyc =mean(phyc)))) %>% summary()

summary(leakyp1)

Variogram(leakyp1, form = ~northing+easting) %>% 
  plot()


leakyp2 <-  
  gls(don~phyc,
      correlation = corExp(form = ~northing+easting, nugget = TRUE),
      data = donut %>% 
        filter(layer == "Pel") %>% 
        left_join(turner %>% group_by(spot) %>% summarise(chl=  mean(chl), phyc =mean(phyc))),
      method = "REML")

summary(leakyp2)

donut %>% 
  left_join(turner %>% group_by(spot) %>% summarise(chl=  mean(chl), phyc =mean(phyc))) %>% 
  ggplot(aes(x = chl, y = don))+
  facet_wrap(~layer, scales = "free")+
  geom_point()+
  theme_bw()



#====Plot midge taxa====

cc <- full_join(cc, grid_sites) %>% filter(!is.na(tanyt))

tanytplot <- cc %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  # geom_point(aes(easting, northing))+
  geom_point(aes(easting, northing, size = tanyt), color = "#e7298a")+
  coord_sf()+
  facet_wrap(~year(sampledate))+
  guides(fill = "none")+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "left",
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank())+

  scale_color_viridis_c()+
  labs(title = "Tanytarsini",
       color = element_blank())


chiroplot <- cc %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  # geom_point(aes(easting, northing))+
  geom_point(aes(easting, northing, size = chiro), color = "#1b9e77")+
  coord_sf()+
  facet_wrap(~year(sampledate))+
  guides(fill = "none")+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "left",
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank())+

  scale_color_viridis_c()+
  labs(title = "Chironomini",
       color = element_blank())



orthoplot <- cc %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  # geom_point(aes(easting, northing))+
  geom_point(aes(easting, northing, size = ortho), color = "#d95f02")+
  coord_sf()+
  guides(fill = "none")+
  facet_wrap(~year(sampledate))+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "left",
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank())+

  labs(title = "Orthocladiinae",
       color = element_blank())




tanypplot <- cc %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  # geom_point(aes(easting, northing))+
  geom_point(aes(easting, northing, size = tanyp), color = "#7570b3")+
  coord_sf()+
  facet_wrap(~year(sampledate))+
  guides(fill = "none")+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "left",
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank())+

  scale_color_viridis_c()+
  labs(title = "Tanypodinae",
       color = element_blank())



tubifexplot <- cc %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  # geom_point(aes(easting, northing))+
  geom_point(aes(easting, northing, size = tubifex), color = "gray20")+
  coord_sf()+
  guides(fill = "none")+
  facet_wrap(~year(sampledate))+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "left",
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank())+

  scale_color_viridis_c()+
  labs(title = "Tubifex",
       color = element_blank())


plot_grid(tanytplot, chiroplot, orthoplot, tanypplot)

ggpreview(plot = last_plot(), width = 6, height = 4, units = "in", dpi = 650)

plot_grid(tanytplot, chiroplot, orthoplot, tanypplot, tubifexplot)


midgecomm <-  cc %>% 
  select(sampledate, easting, northing, tanyt, chiro, ortho, tanyp) %>% 
  mutate(total = tanyt + chiro + ortho+ tanyp) %>% 
  rename(Chironomini = chiro,
         Tanytarsini = tanyt,
         Orthocladiinae = ortho,
         Tanypodinae = tanyp)


cc %>% 
  select(spot2021, year, tanyt, chiro, ortho, tanyp) %>% 
  rename(Chironomini = chiro,
         Tanytarsini = tanyt,
         Orthocladiinae = ortho,
         Tanypodinae = tanyp) %>% 
  gather(taxon, count, -spot2021, -year) %>% 
  spread(year, count) %>% 
  group_by(taxon, `2021`, `2022`) %>% 
  count() %>% 
  ggplot(aes(x = `2021`, y = `2022`))+
  facet_wrap(~taxon, scales = "free")+
  geom_point(aes(color = log10(n)), size = 3)+
  scale_color_viridis_c()


commplot <- ggplot()+
  geom_polygon(aes(long, lat, group = piece), color = "black", alpha = 0.3, data = myvatnshp)+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())+
  geom_scatterpie(aes(x = easting, y = northing, r = total*100+1), legend_name = "taxa", cols = c("Tanytarsini", "Chironomini", "Orthocladiinae", "Tanypodinae"),  data = midgecomm, alpha = 0.8)+
  geom_scatterpie_legend(midgecomm$total*100+1, x = 404500, 7272500, n = 3, labeller = function(x){round(x/100)})+
  theme_bw()+
  facet_wrap(~year(sampledate))+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank())+
  scale_fill_manual(values = c("Tanytarsini" = "#e7298a",
                               "Chironomini" = "#1b9e77",
                               "Orthocladiinae" = "#d95f02",
                               "Tanypodinae" = "#7570b3"))+
  coord_sf()

ggpreview(plot = commplot, width = 6, height = 4, units = "in", dpi = 650)


ggplot()+
  geom_polygon(aes(long, lat, group = piece), color = "black", fill = "dodgerblue", alpha = 0.3, data = myvatnshp)+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())+
  geom_scatterpie(aes(x = easting, y = northing, r = total*100+1), legend_name = "taxa", cols = c("Tanytarsini", "Chironomini", "Orthocladiinae", "Tanypodinae"),  data = midgecomm, alpha = 0.8)+
  geom_scatterpie_legend(midgecomm$total*100+1, x = 404500, 7272500, n = 3, labeller = function(x){round(x/100)})+
  theme(legend.title = element_blank())+
  scale_fill_brewer(palette = "Dark2")+
  coord_sf()


tubifexplot
p1 <- plot_grid(tanytplot, chiroplot, orthoplot, tanypplot)

plot_grid(p1, commplot, rel_widths = c(1.25,1))
ggpreview(last_plot(), dpi = 650, width = 12, height = 6, units = "in")


cm %>% 
  group_by(spot, species_name) %>% 
  count() 

cc %>% 
  select(spot, tanyt:ortho) %>% 
  gather(species, n_id, -spot) %>% 
  mutate(species_name = ifelse(species == "tanyt", "tt",
                               ifelse(species == "chiro","ch",
                                      ifelse(species == "tanyp", "tp",
                                             ifelse(species == "ortho", "or", NA))))) %>% 
  full_join(cm %>% 
              group_by(spot, species_name) %>% 
              count()) %>% 
  mutate(n = ifelse(is.na(n), 0, n)) %>% 
  filter(n_id!=n)

#tanytarsini are all second instar except one found in spot 6, which had a head size of 9 (3rd instar)
cm %>%
  filter(species_name == "tt") %>% 
  ggplot()+
  geom_histogram(aes(x = head_size), center = TRUE)+
  geom_vline(xintercept = c(4.3, 7.3, 12, 19.5))

cm %>%
  filter(species_name == "ch") %>% 
  ggplot(aes(x = head_size))+
  geom_histogram(center = TRUE)+
  geom_vline(xintercept = c(7, 12, 25.5, 44))
#im surprised there are any second instar chironomini (right on the edge of what we call first instar.)
cm %>%
  filter(species_name == "ch") %>% 
  select(spot, species_name, head_size) %>% 
  arrange(head_size)

mg_instars = c(0,4.3,7.3,12, 19.5)/55*1000 #divide by ocular micrometer units to get mm then 1000 to get micrometers
ch_instars = c(3,7,12,25.5,44)/55*1000

instarred <- cm %>% 
  mutate(conversion = head_size_units_to_mm, 
         head_size_micrometer = head_size,
         #convert ocular micrometer units to mm then micrometers
         head_size = head_size/conversion*1000, 
         #assign instar based on head capsule widths
         instar = ifelse(species_name == "tt",
                         cut(head_size, 
                             breaks = mg_instars, 
                             include.lowest =  FALSE),
                         ifelse(species_name == "ch",
                                cut(head_size,
                                    breaks = ch_instars,
                                    include.lowest = FALSE),
                                NA)),
         instar = paste0("s", instar)) 

tts <- instarred %>% 
  filter(species_name == "tt") %>% 
  select(spot, instar) %>% 
  group_by(spot, instar) %>% 
  count() %>% 
  spread(instar, n) %>% 
  full_join(grid_sites) %>% 
  select(sampledate, spot, northing, easting, s2, s3, s4) %>% 
  rowwise() %>% 
  mutate(total = sum(s2, s3, s4, na.rm = TRUE),
         s2 = ifelse(is.na(s2), 0, s2),
         s3 = ifelse(is.na(s3), 0, s3), 
         s4 = ifelse(is.na(s4), 0, s4)) %>% 
  add_row(spot = NA, northing = NA, easting = NA, s2 = 0, s3 = 0, s4 =0, total = 1)

tins_plot <- ggplot()+
  geom_polygon(aes(long, lat, group = piece), color = "black", fill = "dodgerblue", alpha = 0.3, data = myvatnshp)+
  facet_wrap(~year(sampledate))+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())+
  geom_scatterpie(aes(x = easting, y = northing, r = total*150+1), legend_name = "instar", cols = c("s2", "s3", "s4"),  data = tts, alpha = 0.8)+
  # geom_scatterpie_legend(tts$total*100+1, x = 404500, 7272500, n = 3, labeller = function(x){round(x/100)})+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))+
  labs(title = "Tanytarsini")+
  scale_fill_manual(values = c("white", "gray50", "black"))+
  coord_sf()


chs <- instarred %>% 
  filter(species_name == "ch") %>% 
  select(spot, instar) %>% 
  group_by(spot, instar) %>% 
  count() %>% 
  spread(instar, n) %>% 
  full_join(grid_sites) %>% 
  select(sampledate, spot, northing, easting,s1, s2, s3, s4) %>% 
  rowwise() %>% 
  mutate(total = sum(s1, s2, s3, s4, na.rm = TRUE),
         s1 = ifelse(is.na(s2), 0, s1),
         s2 = ifelse(is.na(s2), 0, s2),
         s3 = ifelse(is.na(s3), 0, s3),
         s4 = ifelse(is.na(s4),0, s4))

cins_plot <- ggplot()+
  geom_polygon(aes(long, lat, group = piece), color = "black", fill = "dodgerblue", alpha = 0.3, data = myvatnshp)+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())+
  facet_wrap(~year(sampledate))+
  geom_scatterpie(aes(x = easting, y = northing, r = total*150+1), legend_name = "instar", cols = c("s2", "s3", "s4"),  data = chs, alpha = 0.8)+
  geom_scatterpie_legend(chs$total*150+1, x = 404500, 7272500, n = 3, labeller = function(x){round(x/150)})+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(title = "Chironomini")+
  scale_fill_manual(values = c("white", "gray50", "black"))+
  coord_sf()


ggpubr::ggarrange(cins_plot, tins_plot, common.legend = TRUE, legend = "bottom")
ggpreview(last_plot(), dpi = 650, width = 8, height = 6, units = "in")


#====Plot Turner Probe Data ====
turner <- full_join(turner, grid_sites) 

turner.avg <- turner %>% 
  group_by(spot, spot2021, year, easting, northing) %>% 
  summarise(chl = mean(chl),
            phyc = mean(phyc))

turner_chl <- turner.avg %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  # geom_point(aes(easting, northing))+
  geom_point(aes(easting, northing, color = chl))+
  facet_wrap(~year, nrow = 2)+
  coord_sf()+
  guides(fill = "none")+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "left",
        plot.title = element_text(hjust = 0.5))+

  scale_color_viridis_c()+
  NULL

turner_phyc <- turner.avg %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  # geom_point(aes(easting, northing))+
  geom_point(aes(easting, northing, color = phyc))+
  coord_sf()+
  guides(fill = "none")+
  facet_wrap(~year, nrow = 2)+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "left",
        plot.title = element_text(hjust = 0.5))+

  scale_color_viridis_c()+
  NULL


plot_grid(turner_chl, turner_phyc)
ggpreview(last_plot(), dpi = 650, width = 4, height = 4, units = "in")

turner.avg %>% 
  mutate(phycfrac = phyc/chl) %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  # geom_point(aes(easting, northing))+
  geom_point(aes(easting, northing, color = phycfrac, size = log10(chl+phyc+0.1)))+
  coord_sf()+
  guides(fill = "none")+
  facet_wrap(~year(sampledate))+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "left",
        plot.title = element_text(hjust = 0.5))+

  scale_color_viridis_c(trans = "log10")+
  NULL

turner.avg %>% 
  ungroup %>% 
  summarise(minc = min(chl),
            maxc = max(chl),
            medianc = median(chl),
            minp = min(phyc),
            maxp = max(phyc),
            medianp = median(phyc))

turner.avg %>% 
  gather(pig, read, chl, phyc) %>% 
  ggplot(aes(read))+
  facet_wrap(~pig, scales = "free")+
  geom_histogram()+
  theme_classic()

turner.avg %>% 
  ggplot(aes(x= phyc, y = chl))+
  geom_point()

ggpreview(last_plot(), dpi = 650, width = 2, height = 2, units = "in")


lm(turner.avg$chl~ turner.avg$phyc) %>% summary()
#====Plot BenthoTorch Data====
bt <- full_join(bt, grid_sites) 


pc_plot <- bt %>% 
  filter(year == 2021) %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  # geom_point(aes(easting, northing))+
  geom_point(aes(easting, northing, color = pel_cyano))+
  coord_sf()+
  guides(fill = "none")+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "left",
        plot.title = element_text(hjust = 0.5))+

  scale_color_viridis_c(breaks = c(0, 0.05, 0.1, 0.15))+
  theme(legend.position = "bottom")



pg_plot <- bt %>% 
  filter(year == 2021) %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  # geom_point(aes(easting, northing))+
  geom_point(aes(easting, northing, color = pel_greens))+
  coord_sf()+
  guides(fill = "none")+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "left",
        plot.title = element_text(hjust = 0.5))+

  scale_color_viridis_c(breaks = c(0, 0.05, 0.1, 0.15))+
  theme(legend.position = "bottom")
  

pd_plot <- bt %>% 
  filter(year == 2021) %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  # geom_point(aes(easting, northing))+
  geom_point(aes(easting, northing, color = pel_diatoms))+
  coord_sf()+
  guides(fill = "none")+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "left",
        plot.title = element_text(hjust = 0.5))+

  scale_color_viridis_c(breaks = c(0, 0.05, 0.1, 0.15))+
  theme(legend.position = "bottom")
  
#plot pelagic benthotorch measures (300 mL filtered)
plot_grid(pc_plot, pg_plot, pd_plot, nrow = 1)


pelsp <- bt %>%
  filter(year == 2021) %>% 
  select(spot2021, year, easting, northing, contains("pel")) %>% 
  mutate(chlorophyll = pel_diatoms + pel_cyano + pel_greens) %>% 
  rename(diatoms = pel_diatoms, cyanobacteria = pel_cyano, greens = pel_greens)


ggplot()+
  geom_polygon(aes(long, lat, group = piece), color = "black", alpha = 0.3, data = myvatnshp)+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())+
  geom_scatterpie(aes(x = easting, y = northing, r = chlorophyll*2000), legend_name = "taxa", cols = c("cyanobacteria", "diatoms", "greens"),  data = pelsp)+
  geom_scatterpie_legend(pelsp$chlorophyll*2000, x = 404500, 7272500, n = 2, labeller = function(x) x/2000)+
  scale_fill_manual(values = c("darkorchid", "yellow", "darkolivegreen4"))+
  coord_sf()

ggpreview(last_plot(), dpi = 650, width = 5, height = 4, units = "in")


bt %>% 
  filter(year == 2021) %>% 
  select(contains("pel")) %>% 
  gather(pig, read) %>% 
  ggplot(aes(read))+
  facet_wrap(~pig, scales = "free")+
  geom_histogram()+
  theme_classic()
ggpreview(last_plot(), dpi = 650, width = 6, height = 2, units = "in")


bt %>% 
  filter(year == 2021) %>% 
  select(contains("pel")) %>% 
  gather(pig, read) %>% 
  group_by(pig) %>% 
  summarise(median = median(read),
            min = min(read),
            max = max(read))

bt %>% 
  ggplot(aes(x = pel_cyano, y = pel_greens))+
  geom_point()

bt %>% 
  lm(pel_cyano~pel_greens, data = .) %>% summary()

bc_plot <- bt %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  # geom_point(aes(easting, northing))+
  geom_point(aes(easting, northing, color = sed_cyano))+
  facet_wrap(~year)+
  coord_sf()+
  guides(fill = "none")+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "left",
        plot.title = element_text(hjust = 0.5))+

  scale_color_viridis_c()+
  theme(legend.position = "bottom")

max(bt$sed_greens)


bd_plot <- bt %>% 
  ggplot()+
  facet_wrap(~year)+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  # geom_point(aes(easting, northing))+
  geom_point(aes(easting, northing, color = sed_diatoms))+
  coord_sf()+
  guides(fill = "none")+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "left",
        plot.title = element_text(hjust = 0.5))+

  scale_color_viridis_c()+
  theme(legend.position = "bottom")

#units are micrgrams of chlorophylla/ cm2
br_plot <- bt %>% 
  ggplot()+
  facet_wrap(~year)+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  # geom_point(aes(easting, northing))+
  geom_point(aes(easting, northing, color = sed_cyano/sed_diatoms))+
  coord_sf()+
  guides(fill = "none")+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "left",
        plot.title = element_text(hjust = 0.5))+

  scale_color_viridis_c()+
  theme(legend.position = "bottom")


plot_grid(bc_plot, bd_plot, br_plot, nrow = 1)

bt %>% 
  mutate(cdr)

bensp <-  bt %>% 
  select(year, spot2021, easting, northing, contains("sed")) %>% 
  mutate(chlorophyll = sed_diatoms + sed_cyano + sed_greens) %>% 
  rename(diatoms = sed_diatoms, cyanobacteria = sed_cyano, greens = sed_greens)


ggplot()+
  geom_polygon(aes(long, lat, group = piece), color = "black", alpha = 0.3, data = myvatnshp)+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())+
  facet_wrap(~year)+
  geom_scatterpie(aes(x = easting, y = northing, r = chlorophyll*75), legend_name = "taxa", cols = c("cyanobacteria", "diatoms", "greens"),  data = bensp)+
  geom_scatterpie_legend(bensp$chlorophyll*75, x = 404500, 7272500, n = 3, labeller = function(x) round(x/75, digits = 1))+
  scale_fill_manual(values = c("darkorchid", "yellow", "darkolivegreen4"))+
  theme(legend.title = element_blank())+
  coord_sf()
ggpreview(last_plot(), dpi = 650, width = 5, height = 3, units = "in")

bensp %>% 
  mutate(cdr = cyanobacteria/diatoms) %>% 
  select(spot2021, year, cdr) %>% 
  spread(year, cdr) %>% 
  ggplot(aes(x = `2021`, y = `2022`))+
  geom_point()

bensp %>% 
  mutate(cdr = cyanobacteria/diatoms) %>% 
  select(year, spot2021, cdr) %>% 
  left_join(pelsp %>% filter(year == 2021) %>% select(-year)) %>% 
  ggplot(aes(x = cyanobacteria, y = cdr))+
  facet_wrap(~year)+
  geom_point(alpha = 0.5)+
  labs(x = "Pelagic Cyanobacteria (Aug 2021)",
       y = "Benthic Cyanobacteria:Diatoms")+
  theme_bw()

bsp <- bensp %>% 
  left_join(pelsp %>% filter(year == 2021) %>% select(spot2021, pel_cyano = cyanobacteria)) %>% 
  mutate(cdr = cyanobacteria/diatoms,
         eaststd = (easting-min(easting))/1000,
         northstd = (northing-min(northing))/1000,
         sed_type2 = ifelse(sed_type == "cladophorales", "cladophora", "sediment"),
         sed_type2 = fct_relevel(sed_type2, c("sediment", "cladophora")),
         year = as.factor(year))


gls(cdr~pel_cyano:year,
    data = bsp,
    correlation = corExp(form = ~ easting+northing, nugget = TRUE)) %>% 
  summary()
  
bsp %>% 
  ggplot(aes(x = sed_type2, y = cdr, color = sed_type, shape = year))+
  geom_jitter(width = 0.2)

gls(cdr~sed_type2 + eaststd + northstd+ eaststd:year + northstd:year,
    data = bsp,
    correlation = corExp(form = ~ easting+northing, nugget = TRUE)) %>% 
  summary()

#====Zooplankton====

zoops.long <- zoops %>% 
  select(-rotifers, -alona_total, -comments) %>% 
  gather(species, count, cycl:aspl) %>% 
  mutate(count = count/fract_count) %>% 
  group_by(spot, species) %>% 
  summarise(count = sum(count, na.rm = TRUE))

rotifers <- c("kera", "brac", "poly", "fili", "aspl") 

zooplump <- zoops.long %>%
  group_by(species) %>% 
  mutate(species2 = ifelse(sum(count)>20, species, "Other"),
         species2 = ifelse(species %in% rotifers, "rotfiers", species2)) %>% 
  group_by(spot, species2) %>% 
  summarise(count = sum(count)) %>% 
  left_join(grid_sites) %>% 
  select(easting, northing, species2, count)


zl <-zooplump %>% group_by(spot) %>% mutate(total = sum(count)) %>% spread(species2, count) %>% ungroup %>%  select(-spot)


totplot <- ggplot()+
  geom_polygon(aes(long, lat, group = piece), color = "black", alpha = 0.3, data = myvatnshp)+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())+
  geom_scatterpie(aes(x = easting, y = northing, r = total), alpha = 0.8, legend_name = "taxa", cols = colnames(zl %>% select(-easting, -northing, -total)),  data = zl)+
  geom_scatterpie_legend(zl$total, x = 404500, 7272500, n = 3, labeller = function(x) {round(x/5)*5})+
  theme(legend.title = element_blank())+
  scale_fill_viridis_d()+
  coord_sf(xlim = c(404000, 408000),
           ylim = c(7272000, 7280000))  

ggpreview(last_plot(), dpi = 650, width = 5, height = 4, units = "in")




zooplump2 <- zoops.long %>%
  ungroup %>% 
  group_by(species) %>% 
  filter(!species %in% rotifers,
         sum(count) >0) %>% 
  group_by(species) %>% 
  mutate(species2 = ifelse(species %in% c("alre", "alsp"), "alona", species),
        species2 = ifelse(sum(count)<10, "others", species2)) %>% 
  group_by(spot, species2) %>% 
  summarise(count = sum(count)) %>% 
  left_join(grid_sites) %>% 
  select(easting, northing, species2, count)


zl2 <-zooplump2 %>% group_by(spot) %>% mutate(total = sum(count)) %>% spread(species2, count) %>% ungroup %>%  select(-spot)


rlessplot <- ggplot()+
  geom_polygon(aes(long, lat, group = piece), color = "black", alpha = 0.3, data = myvatnshp)+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())+
  geom_scatterpie(aes(x = easting, y = northing), pie_scale = 5, alpha = 0.8, legend_name = "taxa", cols = colnames(zl2 %>% select(-easting, -northing, -total)),  data = zl2)+
  theme(legend.title = element_blank())+
  scale_fill_viridis_d()+
  coord_sf(xlim = c(404000, 408000),
           ylim = c(7272000, 7280000))  

ggpreview(last_plot(), dpi = 650, width = 5, height = 4, units = "in")



zoops.wide <- zoops.long %>% 
  spread(species, count)


#plot Daphnia

zoops.wide %>% 
  left_join(grid_sites) %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  # geom_point(aes(easting, northing))+
  geom_point(aes(easting, northing, color = log10(dalo)))+
  coord_sf()+
  guides(fill = "none")+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "left",
        plot.title = element_text(hjust = 0.5))+

  scale_color_viridis_c()+
  theme(legend.position = "bottom")

zoops.wide %>% 
  left_join(grid_sites) %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  # geom_point(aes(easting, northing))+
  geom_point(aes(easting, northing, color = (cycl)))+
  coord_sf()+
  guides(fill = "none")+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "left",
        plot.title = element_text(hjust = 0.5))+

  scale_color_viridis_c()+
  theme(legend.position = "bottom")

zoops.wide %>% 
  mutate(rotifers = kera + brac + poly + fili + aspl) %>% 
  left_join(grid_sites) %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  # geom_point(aes(easting, northing))+
  geom_point(aes(easting, northing, color = rotifers))+
  coord_sf()+
  guides(fill = "none")+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "left",
        plot.title = element_text(hjust = 0.5))+

  scale_color_viridis_c()+
  theme(legend.position = "bottom")

zoops.wide %>% 
  left_join(bt) %>% 
  ggplot(aes(x = log10(dalo), y = pel_cyano))+
  geom_point()

zoops.wide %>% 
  mutate(rotifers = kera + brac + poly + fili + aspl) %>% 
  left_join(bt) %>% 
  ggplot(aes(x = rotifers, y = pel_cyano))+
  geom_point()

rotifers <- zoops.wide %>% 
  select(spot, kera, brac, poly, fili, aspl) %>% 
  left_join(grid_sites) %>%
  mutate(total_rotifers = kera+brac + poly + fili + aspl)

rplot <- ggplot()+
  geom_polygon(aes(long, lat, group = piece), color = "black", alpha = 0.3, data = myvatnshp)+
  geom_scatterpie(aes(x = easting, y = northing, r = total_rotifers), legend_name = "rotifer", cols = c("kera", "brac", "poly", "fili", "aspl"),  data = rotifers)+
  geom_scatterpie_legend(rotifers$total_rotifers, x = 404500, 7272500, n = 3, labeller = function(x){round(x/5)*5})+
  coord_sf(xlim = c(404000, 408000),
           ylim = c(7272000, 7280000))+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())
  

plot_grid(totplot, rlessplot, rplot, nrow = 1)
ggpreview(last_plot(), dpi = 650, width = 10, height = 4, units = "in")


zoopmat <- as.matrix(zoops.wide[,2:19])

row.names(zoopmat) <- zoops.wide[,1]

plot_grid()

#calculate bray curtis dissimilarity
zoopdissim <- vegan::vegdist(zoopmat, method = "bray") %>% as.matrix() 



zdissim_long <- zoopdissim %>% 
  as.data.frame() %>% 
  mutate(spot1 = row.names(.)) %>% 
  gather(spot2, dissim, -spot1) %>% 
  mutate(spot2 = as.numeric(spot2),
         spot1 = as.numeric(spot1)) %>% 
  filter(spot1 != spot2)

zdissim_long %>% 
  ggplot(aes(x = spot1, y = spot2, fill = 1-dissim))+
  geom_tile()+
  scale_y_reverse()+
  scale_fill_viridis_c()


left_join(zdissim_long, site_dist_long) %>% 
  filter(spot1>=spot2) %>% 
  ggplot(aes(x = distance/1000, y = 1-dissim))+
  geom_point()+
  labs(x = "Euclidian Distance (km)",
       y = "Community Similarity\n(1-Bray Curtis Dissimilarity)")+
  theme_classic()

ggpreview(last_plot(), dpi = 650, width = 3, height = 3, units = "in")



