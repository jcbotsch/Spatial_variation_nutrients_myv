#====Load Packages====
library(tidyverse)
library(ggspatial)
library(tidyterra)
library(terra)
library(cowplot)
library(scatterpie)
library(lubridate)
library(mgcv)
library(ape)

#====Set aesthetics====
theme_set(  theme_bw()+
              theme(panel.grid = element_line(color = "white"),
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

#=====Load data====
file <- "Lattice_MASTER_19Jan24.xlsx"

site_orig <- read_csv("1km_grid.csv")

grid_sites <- readxl::read_xlsx(file, sheet ="Coordinates", na = "NA") %>%  mutate(year = as.character(year(sampledate)),
                                                                                   spot2021 = ifelse(year == "2021", spot,
                                                                                                     ifelse(spot>=27+31, spot-30, spot-31)))
profiles <- readxl::read_xlsx(file, sheet = "profiles", na = "NA")
turner <- readxl::read_xlsx(file, sheet = "Turner", na = "NA")
cc <- readxl::read_xlsx(file, sheet = "chiro_counts", na = "NA")
nutrients <- readxl::read_xlsx(file, sheet = "nutrients", na = "NA")
detectionlimits <- readxl::read_xlsx(file, sheet = "detection_limits", na = "NA")
bt <- readxl::read_xlsx(file, sheet = "benthotorch", na = "NA")

#=====combine data=====
# pull midges
midges <- cc %>% 
  mutate(chironominae = (chiro + tanyt)/fract_count) %>% 
  dplyr::select(spot, chironominae)

# pull fluorometer data
pel_algae <- turner %>% 
  group_by(spot) %>% 
  summarise(pel_chl = mean(chl),
            pel_phyc = mean(phyc))

# pull site depths
site_depth <- profiles %>% 
  dplyr::select(spot, depth) %>% 
  unique()

# get summary statistics
summary(site_depth)

# prep nutrients
bothyears <- c("nh4", "no3", "po4") # dplyr::select nutrients in both years

nut1 <- nutrients %>% 
  mutate(analyte2 = ifelse(analyte == "nh4_nh3", "nh4", analyte),
         layer_analyte = paste(layer, analyte2, sep = "_")) %>% 
  filter(analyte2 %in% bothyears) %>%                  
  dplyr::select(spot, layer_analyte, mgl) %>%
  spread(layer_analyte, mgl)

# prep site information
sites <- grid_sites %>% 
  mutate(east.std = easting/1000 - min(easting/1000), #km distance from most western site
         north.std = northing/1000 - min(northing/1000), # km distance from most southern site
         basin = ifelse(spot2021 %in% c(21:23), "north", "main")) %>% 
  dplyr::select(site = spot2021, spot, year, sed_type, east.std, north.std, easting, northing, basin)

# join all data together
nut <- sites %>% 
  full_join(site_depth) %>% 
  full_join(nut1) %>% 
  full_join(midges) %>% 
  full_join(pel_algae) %>% 
  full_join(bt %>% dplyr::select(spot, contains("sed"))) %>% 
  mutate(site = as.character(site),
         spot = as.character(spot)) 

#====Get information on conditions during sampling events====
# Limnological data
limno <- profiles %>% 
  filter(sampledepth == 0.5) %>% # take only samples from half a meter below surface
  mutate(secchi = (as.numeric(secchi_dn) + as.numeric(secchi_up))/2) %>% 
  dplyr::select(spot, depth, secchi, wtemp, o2sat) %>% 
  full_join(grid_sites) 

limno %>% 
  dplyr::select(year, spot, depth, wtemp, o2sat, secchi) %>% 
  gather(var, val, depth:secchi) %>% 
  group_by(var, year) %>% 
  reframe(quantiles = quantile(val, c(0, 0.25, 0.5, 0.75, 1)),
            mean = mean(val),
            se = sd(val)/sqrt(n())) %>% 
  print(n = 100)

# midges
midges %>% 
  left_join(grid_sites) %>% 
  # group_by(year) %>% 
  filter(!is.na(chironominae)) %>% 
  summarise(mean = mean(chironominae/ (pi*2.5^2) * 1000),
            se = sd(chironominae/ (pi*2.5^2) *1000)/sqrt(n()))

# pelagic algae
turner.avg <- pel_algae %>% 
  left_join(grid_sites)

turner.avg %>% 
  dplyr::select(year, pel_chl, pel_phyc) %>% 
  gather(var, val, pel_chl, pel_phyc) %>% 
  group_by(var, year) %>% 
  reframe(quantiles = quantile(val, c(0, 0.25, 0.5, 0.75, 1)),
          mean = mean(val),
          se = sd(val)/sqrt(n())) %>% 
  print(n = 100)

# sediment type
sites %>% 
  group_by(year) %>% 
  count(sed_type) %>% 
  spread(year, n)

#=====Global Moran's I====
# create function to calculate global Moran's I 
getMoran.I <- function(sitedata, variable, returnmat =F ){
  distmat <- as.matrix(dist(cbind(sitedata$easting, sitedata$northing)))
  distmat.inv <- 1/distmat
  diag(distmat.inv) <- 0
  
  # print to check
  if(returnmat){
    return(list(distmat.inv <- distmat.inv, Moran = Moran.I(variable, distmat.inv)))
    
  }
  return(Moran = Moran.I(variable, distmat.inv))
  
}

# create function to apply Moran's I to each analyte
yearmorans <- function(data){
  out <- list()
  for(i in 2:ncol(nut1)){
    response <- data %>% 
      pull(colnames(nut1[i]))
    if(length(unique(response))==1){
      NA
    }else{
      out[[i-1]] <- getMoran.I(data, response)
      names(out)[[i-1]] <- colnames(nut1[i])
    }
  }
  return(out)
}


# calculate moran's I to each analyte in each year
sites21 <- sites %>% filter(year == 2021)
sites22 <- sites %>% filter(year == 2022)
nut21 <- sites21 %>% left_join(nut1)
nut22 <- sites22 %>% left_join(nut1)

# print output
y21 <- yearmorans(nut21) 
y22 <- yearmorans(nut22)

y21
y22

names(y21) <- paste0("y21_", names(y21))
names(y22) <- paste0("y22_", names(y22))
# correct pvalues
p.adjust(lapply(c(y21, y22[c(1:3, 5:6)]), FUN = function(x){x$p.value}), method = "fdr") %>% 
  data.frame %>% 
  rename("p.value" = 1) %>% 
  rownames_to_column("details") %>% 
  mutate(moran.i = lapply(c(y21, y22[c(1:3, 5:6)]), FUN = function(x){x$observed}),
         moran.i2 = signif(as.numeric(moran.i), digits = 2),
         p.value2 = signif(p.value, digits = 2)) %>% 
  separate(details, into = c("year", "layer", "analyte")) %>% 
  arrange(analyte, year, layer) 

#====Rank correlations between years and layers====
nut %>%
  gather(analyte, conc, contains("nh4"), contains("no3"), contains("po4")) %>% 
  separate(analyte, into = c("layer", "analyte"), sep = "_") %>% 
  dplyr::select(year, analyte, conc, site, layer) %>% 
  mutate(year = paste0("yr", year)) %>% 
  spread(year, conc) %>% 
  group_by(analyte, layer) %>%
  summarise(cor = cor(yr2022, yr2021, method = "spearman")) %>% 
  ungroup %>% 
  spread(layer, cor)


nut %>% 
  gather(analyte, conc, contains("nh4"), contains("no3"), contains("po4")) %>% 
  separate(analyte, into = c("layer", "analyte"), sep = "_") %>% 
  dplyr::select(year, analyte, conc, site, layer) %>% 
  spread(layer, conc) %>% 
  group_by(analyte, year) %>%
  summarise(cor = cor(ben, pel, method = "spearman")) %>% 
  spread(year, cor)

#====Factors influencing nutrient concentrations=====
nut_nona <- nut %>% 
  filter(!is.na(chironominae)) %>%  # remove NAs
  mutate(sed_type= factor(sed_type, levels = c("sediment", "macrophyte", "cladophorales")))

# benthic ammonium
bnh4 <- gls(log10(ben_nh4 + 5e-4) ~ year + 
                  depth + 
                  east.std + basin + 
                  sed_cyano + sed_diatoms + 
              sed_type + 
                   chironominae, 
                correlation = corExp(form = ~easting + northing| year, nugget = TRUE),
                data = nut_nona,
            method = "ML")


summary(bnh4)
plot(bnh4)
qqnorm(residuals(bnh4, type = "pearson"))
qqline(residuals(bnh4, type = "pearson"))

#====Benthic NH4 forward selection====
bnh4.0 <- gls(log10(ben_nh4 + 5e-4) ~ 1, 
              correlation = corExp(form = ~easting + northing| year, nugget = TRUE),
              data = nut_nona,
              method = "ML")


AIC(bnh4.0, 
      update(bnh4.0, ~ + year),
      update(bnh4.0, ~ + pel_phyc),
      update(bnh4.0, ~ + pel_chl),
      update(bnh4.0, ~ + depth),
      update(bnh4.0, ~ + east.std),
      update(bnh4.0, ~ + basin),
      update(bnh4.0, ~ + sed_cyano),
      update(bnh4.0, ~ + sed_diatoms),
      update(bnh4.0, ~ + sed_type),
      update(bnh4.0, ~ + chironominae))# model with depth is best (AIC = 121)

bnh4.1 <- update(bnh4.0, ~ + depth)

AIC(bnh4.1, 
    update(bnh4.1, ~ . + year),
    update(bnh4.1, ~ . + pel_phyc),
    update(bnh4.1, ~ . + pel_chl),
    update(bnh4.1, ~ . + east.std),
    update(bnh4.1, ~ . + basin),
    update(bnh4.1, ~ . + sed_cyano),
    update(bnh4.1, ~ . + sed_diatoms),
    update(bnh4.1, ~ . + sed_type),
    update(bnh4.1, ~ . + chironominae)) # model with depth and year best (AIC = 117)

bnh4.2 <- update(bnh4.1, ~ . + year)

AIC(bnh4.2, 
    update(bnh4.2, ~ . + pel_phyc),
    update(bnh4.2, ~ . + pel_chl),
    update(bnh4.2, ~ . + year:depth),
    update(bnh4.2, ~ . + east.std),
    update(bnh4.2, ~ . + basin),
    update(bnh4.2, ~ . + sed_cyano),
    update(bnh4.2, ~ . + sed_diatoms),
    update(bnh4.2, ~ . + sed_type),
    update(bnh4.2, ~ . + chironominae)) # model with depth, year, and diatoms best (AIC = 113)


bnh4.3 <- update(bnh4.2, ~ . + sed_diatoms)

AIC(bnh4.3, 
    update(bnh4.3, ~ . + year:depth),
    update(bnh4.3, ~ . + sed_diatoms:depth),
    update(bnh4.3, ~ . + year:sed_diatoms),
    update(bnh4.3, ~ . + pel_phyc),
    update(bnh4.3, ~ . + pel_chl),
    update(bnh4.3, ~ . + east.std),
    update(bnh4.3, ~ . + basin),
    update(bnh4.3, ~ . + sed_cyano),
    update(bnh4.3, ~ . + sed_type),
    update(bnh4.3, ~ . + chironominae)) # model with year diatom interaction best (AIC = 109, model with midges approximately similar)

bnh4.4 <- update(bnh4.3, ~ . + year:sed_diatoms)

AIC(bnh4.4, 
    update(bnh4.4, ~ . + year:depth),
    update(bnh4.4, ~ . + sed_diatoms:depth),
    update(bnh4.4, ~ . + pel_phyc),
    update(bnh4.4, ~ . + pel_chl),
    update(bnh4.4, ~ . + east.std),
    update(bnh4.4, ~ . + basin),
    update(bnh4.4, ~ . + sed_cyano),
    update(bnh4.4, ~ . + sed_type),
    update(bnh4.4, ~ . + chironominae)) # year depth interaction best (AIC = 106)

bnh4.5 <- update(bnh4.4, ~. + year:depth)

AIC(bnh4.5, 
    update(bnh4.4, ~ . + year:depth:sed_diatoms),
    update(bnh4.4, ~ . + sed_diatoms:depth),
    update(bnh4.4, ~ . + pel_phyc),
    update(bnh4.4, ~ . + pel_chl),
    update(bnh4.4, ~ . + east.std),
    update(bnh4.4, ~ . + basin),
    update(bnh4.4, ~ . + sed_cyano),
    update(bnh4.4, ~ . + sed_type),
    update(bnh4.4, ~ . + chironominae)) #base model best

AIC(bnh4.5,
    update(bnh4.5, ~. - year),
    update(bnh4.5, ~. -depth),
    update(bnh4.5, ~. -sed_diatoms),
    update(bnh4.5, .~ - sed_diatoms - depth))


bnh4.best <- update(bnh4.5, method = "REML")
summary(bnh4.best)
plot(bnh4.best)
qqnorm(residuals(bnh4.best, type = "pearson"))
qqline(residuals(bnh4.best, type = "pearson"))

# benthic phosphate
bpo4 <- gls(log10(ben_po4 + 5e-4) ~ year + 
              depth + 
              east.std + basin + 
              sed_cyano + sed_diatoms + 
              chironominae,
                correlation = corExp(form = ~easting + northing| year, nugget = TRUE),
                data = nut_nona)
summary(bpo4)
plot(bpo4)
qqnorm(residuals(bpo4, type = "pearson"))
qqline(residuals(bpo4, type = "pearson"))

#====Benthic phosphate forward selection====
bpo4.0 <- gls(log10(ben_po4 + 5e-4) ~ 1, 
              correlation = corExp(form = ~easting + northing| year, nugget = TRUE),
              data = nut_nona,
              method = "ML")


AIC(bpo4.0, 
    update(bpo4.0, ~ + year),
    update(bpo4.0, ~ + pel_phyc),
    update(bpo4.0, ~ + pel_chl),
    update(bpo4.0, ~ + depth),
    update(bpo4.0, ~ + east.std),
    update(bpo4.0, ~ + basin),
    update(bpo4.0, ~ + sed_cyano),
    update(bpo4.0, ~ + sed_diatoms),
    update(bpo4.0, ~ + sed_type),
    update(bpo4.0, ~ + chironominae)) # model with depth is best (AIC = 132), but effect is weak dAIC = 1.3)

bpo4.1 <- update(bpo4.0, ~ + depth)

AIC(bpo4.1, 
    update(bpo4.1, ~ . + year),
    update(bpo4.1, ~ . + pel_phyc),
    update(bpo4.1, ~ . + pel_chl),
    update(bpo4.1, ~ . + east.std),
    update(bpo4.1, ~ . + basin),
    update(bpo4.1, ~ . + sed_cyano),
    update(bpo4.1, ~ . + sed_diatoms),
    update(bpo4.1, ~ . + sed_type),
    update(bpo4.1, ~ . + chironominae)) # no model improves AIC

bpo4.best <- update(bpo4.1, method = "REML")
summary(bpo4.best)
qqnorm(residuals(bpo4.best, type = "pearson"))
qqline(residuals(bpo4.best, type = "pearson"))

#====Pelagic phosphate forward selection====
ppo4.0 <- gls(log10(pel_po4 + 5e-4) ~  1,
                 correlation = corExp(form = ~easting + northing| year, nugget = TRUE),
                 data = nut_nona,
              method = "ML")

AIC(ppo4.0, 
    update(ppo4.0, ~ + year),
    update(ppo4.0, ~ + pel_phyc),
    update(ppo4.0, ~ + pel_chl),
    update(ppo4.0, ~ + depth),
    update(ppo4.0, ~ + east.std),
    update(ppo4.0, ~ + basin),
    update(ppo4.0, ~ + sed_cyano),
    update(ppo4.0, ~ + sed_diatoms),
    update(ppo4.0, ~ + sed_type),
    update(ppo4.0, ~ + chironominae)) # sed type best, although effect weak dAIC = 1.3

ppo4.1 <- update(ppo4.0, ~ + sed_type)

AIC(ppo4.1, 
    update(ppo4.1, ~ . + year),
    update(ppo4.1, ~ . + pel_phyc),
    update(ppo4.1, ~ . + pel_chl),
    update(ppo4.1, ~ . + depth),
    update(ppo4.1, ~ . + east.std),
    update(ppo4.1, ~ . + basin),
    update(ppo4.1, ~ . + sed_cyano),
    update(ppo4.1, ~ . + sed_diatoms),
    update(ppo4.1, ~ . + chironominae)) # sed diatoms best, but effect weak dAIC = 1.8

ppo4.2 <- update(ppo4.1, ~ . + sed_diatoms)

AIC(ppo4.2, 
    update(ppo4.2, ~ . + sed_type:sed_diatoms),
    update(ppo4.2, ~ . + year),
    update(ppo4.2, ~ . + pel_phyc),
    update(ppo4.2, ~ . + pel_chl),
    update(ppo4.2, ~ . + depth),
    update(ppo4.2, ~ . + east.std),
    update(ppo4.2, ~ . + basin),
    update(ppo4.2, ~ . + sed_cyano),
    update(ppo4.2, ~ . + chironominae)) # pel chl best, but effect weak dAIC = 0.8

ppo4.3 <- update(ppo4.2, ~ . + pel_chl)

AIC(ppo4.3, 
    update(ppo4.3, ~ . + sed_type:sed_diatoms),
    update(ppo4.3, ~ . + sed_type:pel_chl),
    update(ppo4.3, ~ . + pel_chl:sed_diatoms),
    update(ppo4.3, ~ . + year),
    update(ppo4.3, ~ . + pel_phyc),
    update(ppo4.3, ~ . + depth),
    update(ppo4.3, ~ . + east.std),
    update(ppo4.3, ~ . + basin),
    update(ppo4.3, ~ . + sed_cyano),
    update(ppo4.3, ~ . + chironominae)) # depth best, but effect weak dAIC = 0.1

ppo4.4 <- update(ppo4.3, ~ . + depth)

AIC(ppo4.4, 
    update(ppo4.4, ~ . + sed_type:sed_diatoms),
    update(ppo4.4, ~ . + sed_type:pel_chl),
    update(ppo4.4, ~ . + pel_chl:sed_diatoms),
    update(ppo4.4, ~ . + year),
    update(ppo4.4, ~ . + pel_phyc),
    update(ppo4.4, ~ . + east.std),
    update(ppo4.4, ~ . + basin),
    update(ppo4.4, ~ . + sed_cyano),
    update(ppo4.4, ~ . + chironominae)) # no model improves 
    
ppo4.best <- update(ppo4.4, method = "REML")

summary(ppo4.best)
qqnorm(residuals(ppo4.best, type = "pearson"))
qqline(residuals(ppo4.best, type = "pearson"))


#====plots====

nutrient_labels <- data.frame(analyte = c("nh4", "no3", "po4"),
                              analytelabel = c("NH[4]", "NO[3]", "PO[4]"))

layer_labels <- data.frame(layer = c("ben", "pel"),
                           layerlab = c("Benthic", "Pelagic"))

year_labels <- data.frame(year = c("2021", "2022"),
                          yearlab = c("August 2021", "June 2022"))


 # Fig 1
n21 <- nut %>% 
  mutate(layer = ifelse(layer == "ben", "Benthic", "Pelagic")) %>% 
  filter(analyte2 %in% bothyears) %>% 
  mutate(analyte2 = toupper(analyte2),
         analyte2 = paste0(substr(analyte2, 1,2), "[", substr(analyte2,3,3), "]")) %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  geom_point(aes(easting, northing, color = log10(mgl)),  data = . %>% filter(year == 2021))+
  coord_sf()+
  facet_grid(fct_rev(layer)~analyte2, switch = "y", labeller = label_parsed)+
  guides(fill = "none")+
  labs(color = expression(~log[10](mg/L)),
       title = "August 2021")+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "top",
        legend.key.height = unit(0.5, "lines"),
        legend.title = element_text(vjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  scale_color_viridis_c(limits = c(min(log10({nut %>% filter(mgl>0)}$mgl), na.rm = T), max(log10(nut$mgl), na.rm = T)))


n22.2 <- nut %>% 
  mutate(layer = ifelse(layer == "ben", "Benthic", "Pelagic")) %>% 
  filter(analyte2 %in% bothyears) %>% 
  mutate(analyte2 = toupper(analyte2),
         analyte2 = paste0(substr(analyte2, 1,2), "[", substr(analyte2,3,3), "]")) %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  geom_point(aes(easting, northing, color = log10(mgl)),  data = . %>% filter(year == 2022))+
  coord_sf()+
  facet_grid(fct_rev(layer)~analyte2, switch = "y", labeller = label_parsed)+
  guides(fill = "none")+
  labs(color = expression(~log[10](mg/L)),
       title = "June 2022")+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom",
        legend.key.height = unit(0.5, "lines"),
        legend.title = element_text(vjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  scale_color_viridis_c(limits = c(min(log10({nut %>% filter(mgl>0)}$mgl), na.rm = T), max(log10(nut$mgl), na.rm = T)))

plot_grid(nz21, nz22)

long <- plot_grid(n21, n22.2, ncol = 1)

ggpreview(plot = long, width = 84, height= 200, units = "mm", dpi = 650)


ggsave(plot = long, width = 84, height= 200, units = "mm", dpi = 650, filename = "spat_nuts_fig2.jpeg")


# Fig 2
ncor1 <- nut %>% 
  left_join(year_labels) %>% 
  gather(analyte, conc, contains("ben"), contains("pel")) %>% 
  separate(analyte, into = c("layer", "analyte"), sep = "_") %>% 
  left_join(layer_labels) %>% 
  dplyr::select(yearlab, analyte, conc, site, layerlab) %>% 
  spread(layerlab, conc) %>%
  left_join(nutrient_labels) %>% 
  mutate(analytelabel = paste(yearlab, analytelabel, sep = "~"),
         analytelabel = str_replace(analytelabel, " ", "~")) %>% 
  ggplot(aes(x = Benthic, y = Pelagic, fill = yearlab))+
  facet_wrap(~analytelabel, scales = "free", labeller = label_parsed) +
  geom_point(shape = 21) +
  scale_fill_manual(values = c("black", "gray"))+
  labs(fill = element_blank()) + 
  theme(legend.position = "none")


ncor2 <- nut %>% 
  gather(analyte, conc, contains("ben"), contains("pel")) %>% 
  separate(analyte, into = c("layer", "analyte"), sep = "_") %>% 
  dplyr::select(year, analyte, conc, site, layer) %>% 
  left_join(nutrient_labels) %>% 
  left_join(year_labels) %>% 
  dplyr::select(yearlab, analytelabel, conc, site, layer) %>% 
  spread(yearlab, conc) %>%  
  left_join(layer_labels) %>% 
  mutate(analytelabel = paste(layerlab, analytelabel, sep = "~"),
         analytelabel = fct_relevel(analytelabel, c("Pelagic~NH[4]", "Pelagic~NO[3]", "Pelagic~PO[4]", "Benthic~NH[4]", "Benthic~NO[3]", "Benthic~PO[4]"))) %>%
  ggplot(aes(x = `August 2021`, y = `June 2022`, fill = layerlab))+
  facet_wrap(~analytelabel, scales = "free", labeller = label_parsed) +
  geom_point(shape = 21) + 
  scale_fill_manual(values = c("darkorange", "dodgerblue"))+
  theme(legend.position = "none")


ncor <- plot_grid(ncor1, ncor2, align = "v", ncol = 1, labels = c("a", "b"))

ggpreview(plot = ncor, width = 174, height= 174, units = "mm", dpi = 650)


#############################################
############## SUPPLEMENT ###################
#############################################

# N:P
totnut <- nutrients %>% 
  filter(analyte %in% c("tn", "tp")) %>%                  
  dplyr::select(spot, layer, analyte, mgl) %>% 
  spread(analyte, mgl) %>% 
  mutate(spot = as.character(spot),
         np = tn/tp,
         np_molar = (tn/14.01)/(tp/30.97)) %>% 
  left_join(nut)

maptots <- nutrients %>% 
  filter(analyte %in% c("tn", "tp")) %>%                  
  dplyr::select(spot, layer, analyte, mgl) %>% 
  left_join(sites) %>% 
  filter(year == "2021") %>% 
  mutate(layer = ifelse(layer == "ben", "Benthic", "Pelagic")) %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  geom_point(aes(easting, northing, color = mgl))+
  coord_sf()+
  facet_grid(toupper(analyte)~fct_rev(layer), switch = "y", labeller = label_parsed)+
  guides(fill = "none")+
  scale_color_viridis_c(trans = "log10", breaks = c(0.01, 0.1, 1,10))+
  labs(color = "Concentration (mg/L)",
       title = "August 2021")+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "top",
        legend.key.height = unit(0.5, "lines"),
        legend.title = element_text(vjust = 1),
        plot.title = element_text(hjust = 0.5))

histtots <- nutrients %>% 
  filter(analyte %in% c("tn", "tp")) %>%                  
  dplyr::select(spot, layer, analyte, mgl) %>% 
  left_join(sites) %>% 
  filter(year == "2021") %>% 
  mutate(layer = ifelse(layer == "ben", "Benthic", "Pelagic")) %>% 
  ggplot(aes(mgl)) + 
  facet_grid(toupper(analyte)~fct_rev(layer), scales = "free_x") + 
  geom_histogram() + 
  labs(x = "Concentration (mg/L)")

nplim1 <- totnut %>% 
  mutate(layer = ifelse(layer == "ben", "Benthic", "Pelagic")) %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  geom_point(aes(easting, northing, color = np_molar))+
  coord_sf()+
  facet_grid(~fct_rev(layer), switch = "y", labeller = label_parsed)+
  guides(fill = "none")+
  labs(color = "N:P",
       title = "August 2021")+
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "top",
        legend.key.height = unit(0.5, "lines"),
        legend.title = element_text(vjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  scale_color_gradient2(low = "red", high = "darkblue", midpoint = log10(16), mid = "gray", trans = "log10")


nplim2 <- totnut %>% 
  mutate(layer = ifelse(layer == "ben", "Benthic", "Pelagic")) %>% 
  ggplot(aes(np_molar)) + 
  facet_wrap(~fct_rev(layer)) + 
  geom_histogram() + 
  geom_vline(xintercept = 16) +
  labs(x = "TN : TP (Molar)")

plot_grid(nplim1, maptots, nplim2, histtots, ncol = 2, align = "v", axis = "lr", rel_heights = c(1, 0.5))

ggpreview(plot = last_plot(), width = 160, height= 220, units = "mm", dpi = 650)

totnut %>% 
  group_by(layer) %>% 
  count(plim = np_molar>16)

totnut %>% 
  group_by(layer) %>% 
  reframe(mean = mean(np_molar),
          se = sd(np_molar)/ sqrt(n()),
          quantiles = quantile(np_molar, c(0.5, 0.25, 0.5, 0.75, 0.95)))

# depth
limno %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  # geom_point(aes(easting, northing))+
  geom_point(aes(easting, northing, color = depth), size = 3)+
  facet_wrap(~year) +
  coord_sf()+
  guides(fill = "none")+
  labs(color = "Depth (m)") + 
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "left",
        plot.title = element_text(hjust = 0.5))+
  
  scale_color_gradient(low = "white", high = "black")


# Oxygen
limno %>% 
  ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  # geom_point(aes(easting, northing))+
  geom_point(aes(easting, northing, color = o2sat/100), size = 3)+
  facet_wrap(~year) +
  coord_sf()+
  guides(fill = "none")+
  labs(color = "Dissolved oxygen\n% saturation at 50cm") + 
  theme_bw()+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "left",
        plot.title = element_text(hjust = 0.5))+
  
  scale_color_gradient2(low = "red", high = "blue", midpoint= 1, labels= scales::percent)


# midges
chircomm <-  cc %>% 
  left_join(sites) %>% 
  dplyr::select(year, easting, northing, tanyt, chiro, ortho, tanyp) %>% 
  mutate(total = tanyt + chiro) %>% 
  rename(Chironomini = chiro,
         Tanytarsini = tanyt,
         Orthocladiinae = ortho,
         Tanypodinae = tanyp) %>% 
  filter(!is.na(year))

commplot <- ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())+
  geom_point(aes(easting, northing), shape = 21, size = 0.1, color = "gray20", ,  data = chircomm %>% filter(total==0, year!="NA"))+
  geom_scatterpie(aes(x = easting, y = northing, r = total*100+10), legend_name = "taxa", cols = c("Tanytarsini", "Chironomini"),  data = chircomm, alpha = 0.8)+
  geom_scatterpie_legend(chircomm$total*100+1, x = 404500, 7272500, n = 3, labeller = function(x){round(x/100)})+
  theme_bw()+
  facet_wrap(~year)+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "bottom",
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank())+
  scale_fill_manual(values = c("Tanytarsini" = "#e7298a",
                               "Chironomini" = "#1b9e77"))+
  coord_sf()

ggpreview(plot = commplot, width = 6, height = 4, units = "in", dpi = 650)

# benthotorch
bt_com <-  bt %>% 
  left_join(sites) %>% 
  dplyr::select(year, easting, northing, sed_diatoms, sed_cyano, sed_greens) %>% 
  mutate(total = sed_diatoms, sed_cyano, sed_greens) %>% 
  rename(Diatoms = sed_diatoms,
         `Green algae` = sed_greens,
         Cyanobacteria = sed_cyano) %>% 
  filter(!is.na(year))

benalgaeplot <- ggplot()+
  geom_spatvector(fill = "dodgerblue", alpha = 0.2, data = myvatnshp)+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())+
  geom_point(aes(easting, northing), shape = 21, size = 0.1, color = "gray20", ,  data = bt_com)+
  geom_scatterpie(aes(x = easting, y = northing, r = total*100+1), legend_name = "taxa", cols = c("Diatoms", "Green algae", "Cyanobacteria"),  data = bt_com, alpha = 0.8)+
  geom_scatterpie_legend(bt_com$total*100+1, x = 404500, 7272500, n = 3, labeller = function(x){round(x/100)})+
  theme_bw()+
  facet_wrap(~year)+
  theme(panel.grid = element_line(color = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "bottom",
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank())+
  scale_fill_manual(values = c("Cyanobacteria" = "darkorchid",
                               "Diatoms" = "goldenrod",
                               "Green algae" = "green3"))+
  coord_sf()

ggpreview(plot = benalgaeplot, width = 6, height = 4, units = "in", dpi = 650)


