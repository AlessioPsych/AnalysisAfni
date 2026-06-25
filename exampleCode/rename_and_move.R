library(tools)

# Rename files so that number is in front

# Set the directory path where your files are located
directory_path <- "I:\\20230620_CCY20"
setwd(directory_path)

# Get a list of all files in the directory
files <- list.files(directory_path)

# Iterate through each file and rename them
for (file in files) {
  # Check if the file has the desired extensions
  if (file_ext(file) %in% c("nii", "json")) {
    # Extract the number from the file name
    number <- sub(".*_(\\d+)_.*", "\\1", file)
    
    # Create the new file name by concatenating the number and the original file name
    new_file_name <- paste0(number, "_", file)
    
    # Construct the full paths for the old and new file names
    old_path <- file.path(directory_path, file)
    new_path <- file.path(directory_path, new_file_name)
    
    # Rename the file
    file.rename(old_path, new_path)
  }
}

# Move all TOPUP Amplitude files into TOPUP folder and all TOPUP Phase files into TOPUP_PHASE folder

# Create the destination folder paths
destination_folder_odd <- file.path(directory_path, "TOPUP")
destination_folder_even <- file.path(directory_path, "TOPUP_PHASE")

# Get a list of all files in the directory
files <- list.files(directory_path)

# Iterate through each file and copy eligible files to the respective destination folders
for (file in files) {
  # Check if the file starts with an odd number and contains '_TOPUP'
  if (grepl("^\\d+.*_TOPUP.*", file)) {
    # Extract the leading digits from the file name
    leading_digits <- as.integer(str_extract(file, "^\\d+"))
    
    # Construct the full paths for the source and destination files
    source_path <- file.path(directory_path, file)
    
    # Determine the destination folder based on whether the leading digits are even or odd
    if (leading_digits %% 2 == 0) {
      destination_folder <- destination_folder_even
    } else {
      destination_folder <- destination_folder_odd
    }
    
    destination_path <- file.path(destination_folder, file)
    
    # Move the file to the destination folder
    file.rename(source_path, destination_path)
  }
}


# Iterate through each file and copy eligible files to the respective destination folders
for (file in files) {
  # Check if the file starts with an odd number and contains '_TOPUP'
  if (grepl("^\\d+.*_TOPUP.*", file)) {
    # Extract the leading digits from the file name
    leading_digits <- as.integer(str_extract(file, "^\\d+"))
    
    # Construct the full paths for the source and destination files
    source_path <- file.path(directory_path, file)
    
    # Determine the destination folder based on whether the leading digits are even or odd
    if (leading_digits %% 2 == 0) {
      destination_folder <- destination_folder_even
    } else {
      destination_folder <- destination_folder_odd
    }
    
    destination_path <- file.path(destination_folder, file)
    
    # Move the file to the destination folder
    file.rename(source_path, destination_path)
  }
}

# Move ANATOMY files

destination_folder <- file.path(directory_path, "ANATOMY")

# Get a list of all files in the directory
files <- list.files(directory_path)

# Iterate through each file and move files containing 'mp2rage' to the destination folder
for (file in files) {
  # Check if the file name contains 'mp2rage'
  if (grepl("mp2rage", file, ignore.case = TRUE)) {
    # Construct the full paths for the source and destination files
    source_path <- file.path(directory_path, file)
    destination_path <- file.path(destination_folder, file)
    
    # Move the file to the destination folder
    file.rename(source_path, destination_path)
  }
}

## Finally move all nifti encoding and retrievals to EPI and EPI_PHASE

destination_folder_odd <- file.path(directory_path, "EPI")
destination_folder_even <- file.path(directory_path, "EPI_PHASE")

# Get a list of all files in the directory
files <- list.files(directory_path)

# Iterate through each file and move eligible files to the respective destination folders
for (file in files) {
  # Check if the file has the extension "*.nii"
  if (endsWith(file, ".nii")) {
    # Extract the leading digits from the file name
    leading_digits <- as.integer(str_extract(file, "^\\d+"))
    
    # Construct the full paths for the source and destination files
    source_path <- file.path(directory_path, file)
    
    # Determine the destination folder based on whether the leading digits are even or odd
    if (leading_digits %% 2 == 0) {
      destination_folder <- destination_folder_even
    } else {
      destination_folder <- destination_folder_odd
    }
    
    destination_path <- file.path(destination_folder, file)
    
    # Move the file to the destination folder
    file.rename(source_path, destination_path)
  }
}


## ... same for jsons encoding and retrievals to EPI and EPI_PHASE

destination_folder_odd <- file.path(directory_path, "EPI_jsons")
destination_folder_even <- file.path(directory_path, "EPI_PHASE_jsons")

# Get a list of all files in the directory
files <- list.files(directory_path)

# Iterate through each file and move eligible files to the respective destination folders
for (file in files) {
  # Check if the file has the extension "*.json"
  if (endsWith(file, ".json")) {
    # Extract the leading digits from the file name
    leading_digits <- as.integer(str_extract(file, "^\\d+"))
    
    # Construct the full paths for the source and destination files
    source_path <- file.path(directory_path, file)
    
    # Determine the destination folder based on whether the leading digits are even or odd
    if (leading_digits %% 2 == 0) {
      destination_folder <- destination_folder_even
    } else {
      destination_folder <- destination_folder_odd
    }
    
    destination_path <- file.path(destination_folder, file)
    
    # Move the file to the destination folder
    file.rename(source_path, destination_path)
  }
}