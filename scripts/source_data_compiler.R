# compile_source_data.R
if (!require("openxlsx")) install.packages("openxlsx")
library(openxlsx)

# Configuration
input_dir <- "data/source_data"
output_file <- "data/Source_Data.xlsx"

# Initialize Workbook
wb <- createWorkbook()

# --- 1. Helper function to load RDS and write to sheet ---
add_sheet_from_file <- function(wb, filename) {
  file_path <- file.path(input_dir, filename)
  
  # Derive sheet name from filename (e.g., "Fig2a.Rds" -> "Fig2a")
  sheet_name <- gsub("\\.Rds$", "", filename)
  
  data <- readRDS(file_path)
  
  # Special handling for Lists (like Network nodes/edges)
  if (inherits(data, "list") && !is.data.frame(data)) {
    for (sub_name in names(data)) {
      # Create name like "Fig2a_nodes"
      sub_sheet <- paste0(sheet_name, "_", sub_name)
      # Truncate to 31 chars (Excel limit)
      sub_sheet <- substr(sub_sheet, 1, 31)
      
      addWorksheet(wb, sub_sheet)
      writeData(wb, sub_sheet, data[[sub_name]])
      message(paste("  Added sheet:", sub_sheet))
    }
  } else {
    # Standard Dataframe
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, data)
    message(paste("  Added sheet:", sheet_name))
  }
}

# --- 2. Main Execution Loop ---

# Get all RDS files
files <- list.files(input_dir, pattern = "\\.Rds$", full.names = FALSE)

# Sort them just to be safe (Fig2a, Fig2b, Fig2c...)
# Note: If you have Fig1 and Fig10, use stringr::str_sort(files, numeric=TRUE)
files <- sort(files) 

if(length(files) == 0) stop("No .Rds files found in data/source_data/")

message(paste("Found", length(files), "files. Compiling..."))

# Loop through every file and add it
for (f in files) {
  add_sheet_from_file(wb, f)
}

# --- 3. Save ---
saveWorkbook(wb, output_file, overwrite = TRUE)
message(paste("\nDONE! Saved to:", output_file))