# Clean Global Environment
rm(list = ls())

# Load necessary packages using 'pacman'
library(pacman)
pacman::p_load(tidyverse, FactoMineR, factoextra, readxl)
source("https://raw.githubusercontent.com/Cassava2050/PPD/main/utilities_tidy.R")

# Define the output folder and list files
folder_output <- here::here("output//")
list_file <- list.files(folder_output)

# Select files containing "BLUPs" in their names
sel_file <- list_file[str_detect(list_file, "BLUPs")]
sel_file

# Select the first file from the filtered list
sel_file[1]

# Read data from the selected Excel file
blupDF_kp <- read_excel(
  paste(folder_output,
        sel_file[1],
        sep = ""
  ),
  sheet = paste0("BLUPs_", "gxe")
)

# Check column names
colnames(blupDF_kp)

# Define traits of interest
index_traits <- colnames(blupDF_kp)[-1]

# Select relevant columns and drop rows with missing values
index_dat <- blupDF_kp %>%
  select("accession_name", all_of(index_traits)) %>% 
  drop_na()


# Define a function for multi-trait PCA-based index selection
pca_index <- function(data, id, variables = NULL, percentage = 0.60, b) {
  # Convert data to a data frame
  data <- as.data.frame(data)
  
  # Set row names
  rownames(data) <- data[, id]
  
  # If 'variables' is not specified, use all columns except 'id'
  if (is.null(variables)) variables <- names(data)[names(data) != id]
  
  # Select relevant columns
  data <- data[, variables]
  
  # Calculate the selection index
  index <- selIndex(Y = as.matrix(data), b = b, scale = TRUE)
  index <- c(index)
  
  # Add the selection index to the data frame and sort by descending index values
  data$index <- index
  data <- data %>% arrange(desc(index))
  
  # Mark the top 'percentage' of genotypes as selected
  data$selected <- NA
  data$selected[1:(round(percentage * nrow(data)))] <- TRUE
  
  # Set 'selected' to FALSE for genotypes not selected
  data$selected <- ifelse(is.na(data$selected), FALSE, data$selected)
  
  # Perform PCA
  res.pca <- PCA(data, graph = TRUE, scale.unit = TRUE, quali.sup = ncol(data))
  
  # Create a PCA biplot
  final <- fviz_pca_biplot(res.pca,
                           habillage = data$selected,
                           geom = c("point"),
                           addEllipses = TRUE,
                           col.var = "black",
                           ggtheme = theme_minimal()
  )
  
  # Filter selected genotypes
  selection <- data %>% filter(selected == TRUE)
  
  # Return results as a list
  return(list(res.pca = res.pca, final = final, results = data, selection = selection))
}

# Function to calculate the selection index
selIndex <- function(Y, b, scale = FALSE) {
  if (scale) {
    return(scale(Y) %*% b)
  }
  return(Y %*% b)
}



# Perform PCA-based index selection
res.pca <- pca_index(
  data = index_dat, id = "accession_name",
  variables = index_traits,
  b = c(-2.47, 1.044, 2.523, -.04, 1.412, 1.615), percentage = 0.60
)

# Save the PCA biplot as an image
res.pca_final <- res.pca$final + theme_xiaofei()
ggsave(paste("images/selection", Sys.Date(), ".png"),
       plot = res.pca_final, units = "in", dpi = 300, width = 7, height = 6
)

# Print selected genotypes
res.pca$selection

# Convert 'results' data frame to a format that can be joined with BLUP data
selections <- res.pca$results %>% rownames_to_column(var = "accession_name")

# Join the selection data with BLUP data and select relevant columns
selections <- blupDF_kp %>% 
  left_join(selections %>% 
              select(accession_name, index, selected), by = "accession_name") %>% 
  select(accession_name, index, selected, DM_gravity, yield_ha, plant_type, everything()) %>% 
  arrange(desc(index))

# Save the selections to a CSV file
write.csv(selections, paste0(folder_output, "selections_by_index.csv"), row.names = FALSE)


