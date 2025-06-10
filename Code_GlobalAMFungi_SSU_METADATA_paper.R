# =====================================================
# LOAD REQUIRED LIBRARIES
# =====================================================
library(raster)
library(sp)
library(sf)
library(dplyr)
library(ggplot2)
library(terra)
library(viridis)
library(stringr)  # For text wrapping
library(CoordinateCleaner)

# =====================================================
# SET WORKING DIRECTORY
# =====================================================
setwd("")

# =====================================================
# LOAD DATA: Ecoregion Shapefile, Raster, Sample Coordinates
# =====================================================
ecoregions <- st_read("Ecoregions2017.shp")
ecoregions_raster <- raster("Ecoregion_GEE.tif")

# Metadata table
EcoregionMetadata <- data.frame(
  EcoID = ecoregions$ECO_ID,
  BIOME_NAME = ecoregions$BIOME_NAME,
  ECO_NAME = ecoregions$ECO_NAME,
  COLOR = ecoregions$COLOR
)

# Load and filter samples
GF_Samples <- read.csv("GlobalAMFungi_SSU_METADATA_paper.csv", header = TRUE) %>%
  filter(!X...sample_type..TEXT.NOT.NULL. %in% c("root", "sediment"))

Samples_latlong <- na.omit(GF_Samples %>%
                             dplyr::select(longitude = longitude, latitude = latitude) %>%
                             mutate(across(c(longitude, latitude), as.numeric)))

# =====================================================
# ASSIGN SAMPLES TO ECOREGIONS
# =====================================================
extract_values <- terra::extract(rast(ecoregions_raster), Samples_latlong)
Samples_ecoregions <- cbind.data.frame(Samples_latlong, extract_values)
Samples_ecoregions$EcoID <- Samples_ecoregions$Ecoregion_GEE

# Merge ecoregion metadata with sample assignments
ecoregion.full.list <- full_join(EcoregionMetadata, Samples_ecoregions, by = "EcoID")

# =====================================================
# SUMMARIZE SAMPLING PER ECOREGION
# =====================================================
ecoregion.sampled.list <- Samples_ecoregions %>%
  group_by(EcoID) %>%
  summarise(total_count = n())

ecoregion.full.list <- full_join(EcoregionMetadata, ecoregion.sampled.list, by = "EcoID")
ecoregion.full.list$total_count[is.na(ecoregion.full.list$total_count)] <- 0  # Set NA sample counts to 0
ecoregion.full.list$ECO_NAME[is.na(ecoregion.full.list$ECO_NAME)] <- 'Rock&Ice'

ecoregion.full.list.narm <- ecoregion.full.list %>% filter(!is.na(EcoID))

# =====================================================
# COMPUTE SAMPLING STATISTICS
# =====================================================
total_ecoregions <- nrow(ecoregion.full.list.narm)
num_sampled_ecoregions <- sum(ecoregion.full.list.narm$total_count > 0)
mean_total_count <- mean(ecoregion.full.list.narm$total_count)
range_total_count <- range(ecoregion.full.list.narm$total_count)

percent_gt_5 <- sum(ecoregion.full.list.narm$total_count > 5) / total_ecoregions * 100
percent_gt_10 <- sum(ecoregion.full.list.narm$total_count > 10) / total_ecoregions * 100
percent_gt_50 <- sum(ecoregion.full.list.narm$total_count > 50) / total_ecoregions * 100

# Display results
cat("Total ecoregions:", total_ecoregions, "\n")
cat("Number of sampled ecoregions:", num_sampled_ecoregions, "\n")
cat("Mean sample count:", mean_total_count, "\n")
cat("Range of sample count:", range_total_count, "\n")
cat("Percentage of ecoregions with >5 samples:", round(percent_gt_5, 2), "%\n")
cat("Percentage of ecoregions with >10 samples:", round(percent_gt_10, 2), "%\n")
cat("Percentage of ecoregions with >50 samples:", round(percent_gt_50, 2), "%\n")

# =====================================================
# BAR PLOTS: PERCENTAGE OF ECOREGIONS SAMPLED PER BIOME
# =====================================================
biome_summary <- ecoregion.full.list.narm %>%
  filter(!is.na(BIOME_NAME) & BIOME_NAME != "N/A") %>%
  mutate(sampled = total_count > 0) %>%
  group_by(BIOME_NAME) %>%
  summarise(
    ecoregions_sampled = sum(sampled),
    total_ecoregions = n(),
    percent_sampled = (ecoregions_sampled / total_ecoregions) * 100
  ) %>%
  arrange(desc(percent_sampled))

biome_summary$BIOME_NAME_wrapped <- factor(
  str_wrap(biome_summary$BIOME_NAME, width = 20),
  levels = str_wrap(biome_summary$BIOME_NAME, width = 20)
)

sampled_percents <- ggplot(biome_summary, aes(x = BIOME_NAME_wrapped, y = percent_sampled)) +
  geom_col(fill = "forestgreen") +
  geom_text(aes(label = sprintf("%.1f%%", percent_sampled)), vjust = -0.5, size = 4, color = "black") +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  theme_minimal() +
  labs(x = "Biome", y = "Percentage of Sampled Ecoregions", title = "Percentage of Ecoregions Sampled per Biome") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Manuscript/Figures/Figure3_FEMS.pdf", plot = sampled_percents, device = "pdf", scale = 0.6)

# =====================================================
# BAR PLOT: TOTAL SAMPLES PER BIOME
# =====================================================
biome_summary_counts <- ecoregion.full.list.narm %>%
  filter(!is.na(BIOME_NAME) & BIOME_NAME != "N/A") %>%
  group_by(BIOME_NAME) %>%
  summarise(total_samples = sum(total_count)) %>%
  arrange(desc(total_samples))

sampled_counts <- ggplot(biome_summary_counts, aes(x = BIOME_NAME, y = total_samples)) +
  geom_col(fill = "cornflowerblue") +
  geom_text(aes(label = total_samples), vjust = -0.5, size = 4, color = "black") +
  theme_minimal() +
  labs(x = "Biome", y = "Total Sample Count", title = "Total Samples per Biome") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Manuscript/Figures/Figure4_FEMS.pdf", plot = sampled_counts, device = "pdf", scale = 0.6)

# =====================================================
# IDENTIFY SAMPLES IN URBAN AREAS
# =====================================================
samples_df <- as.data.frame(Samples_latlong)
colnames(samples_df) <- c("decimalLatitude", "decimalLongitude")

samples_df$urban <- clean_coordinates(
  x = samples_df, lon = "decimalLongitude", lat = "decimalLatitude",
  tests = "urban", value = "flagged"
)

num_samples_in_cities <- sum(samples_df$urban == TRUE, na.rm = TRUE)
percent_in_cities <- (num_samples_in_cities / nrow(samples_df)) * 100

cat("Samples found in urban areas:", num_samples_in_cities, "\n")
cat("Percentage of samples in urban areas:", round(percent_in_cities, 2), "%\n")
