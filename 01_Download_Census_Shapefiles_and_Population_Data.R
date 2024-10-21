# With distributed computing, several iterations of your script are running
# simultaneously. If the code includes API calls -- for instance, downloading
# a Census shapefile for your area of interest -- these simultaneous requests
# can be blocked by the server, resulting in your code failing. For this tutorial,
# we need census-block shapefiles for every state, as well as a ZCTA shapefile
# for the entire country. Instead of downloading these files within the same code
# as the processing steps, we will instead download all of the necessary files
# locally first and have the next step read in those downloaded files. 
#
# NOTE: These files will consume a substantial amount of storage space on your machine.
#       If you want to run just a sample, set the states variable assignment below
#       to include only "11" (state FIPS code for Washington, D.C.), rather than
#       for all 50 states + DC.
#
library("tigris")     # Used to download census shapefiles; no API needed
library("sf")         # Used for reading/write vector files and performing spatial analysis
library("tidycensus") # Used to download census data; API required (see below)

# %%%%%%%%%%%%%%%%%%%%%%% USER-DEFINED PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#
# Obtain a Census API key here: https://www.census.gov/data/developers.html
# Run the following commented-out lines of code ONCE prior to submitting via Bash script:
#
# my_census_API_token <- "" # enter your API token here -- do NOT save this in your code as plain text
# census_api_key(my_census_API_token, install = TRUE) # saves your API to a .Renviron file
#
# NOTE: The steps above will save your Census API token in a plain-text file on your machine.
#       Anyone with access to the directory where this is saved will be able to see
#       your API token. Alternatively, you can use the keyring() package to store credentials
#       securely in R, but this presents challenges if you ever choose to process code
#       that downloads census data via bash script, as it requires user input of the password to 
#       unlock the keyring. For more information on securely storing credentials
#       in R using keyring(), see: https://www.infoworld.com/article/3320999/r-tip-keep-your-passwords-and-tokens-secure-with-the-keyring-package.html
#
# Set the input data directory
#
input_data_dir <- "In_Dir/" # You may need to change this to wherever you want to save the shapefiles

# Define the year for which you want shapefiles. Major changes to ZCTAs and blocks
# occur during the Decennial Census (2000, 2010, 2020, etc.), but some small corrections
# and changes may occur between them as well. Note, however, that the Census Bureau
# will not necessarily update FIPS codes for block shapefiles between Decennial
# Censuses, even if there are changes to block groups, tracts, and/or counties.
# For example, as of January 2024, changes to Connecticut county-equivalents
# (which went into effect starting with 2022 census data products) are not reflected
# in the GEOID codes for 2022 or 2023 block shapefiles, but they *are* changed in
# the tract- and county-level shapefiles.
#
year <- 2020

# Identify the FIPS codes (GEOIDs) of the 50 states + DC
#
states <- c("01", "02", "04", "05", "06", "08", "09", "10", "11", "12", "13", "15", "16", "17", "18", "19", "20", "21", "22", "23",
            "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40", "41", "42", "44",
            "45", "46", "47", "48", "49", "50", "51", "53", "54", "55", "56")

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# %%%%%%%%%%%%%%%%%%%%% DOWNLOAD SHAPEFILES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

# The ZCTA shapefile cannot be downloaded by state for 2020. Note also that
# ZCTAs occassionally cross state lines. Download the national file once and
# save it to the input_data_directory as a geopackage (.gpkg):
#
zctas <- tigris::zctas(year = year)
st_write(zctas, paste0(input_data_dir, "ZCTA_", year, "_US.gpkg"))

# Download the block shapefiles by state and save them to the input_data_directory
#
for (i in 1:length(states)) {
  
  cat("--------------------------------------------------------------------- \n")
  cat("Beginning download of state FIPS", states[i], "\n")
  blocks <- tigris::blocks(state = states[i], year = year)
  st_write(blocks, paste0(input_data_dir, "Blocks_", year, "_", states[i], ".gpkg"))
  
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# %%%%%%%%%%%%%%%%%%%%% DOWNLOAD POPULATION DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#
# Obtain block-level population data. Note that the Census Bureau only releases
# population counts at the block level in the Decennial Census (2000, 2010, 2020, etc.).
# Note also that block identifiers (i.e., GEOID or FIPS code) can change over time,
# even between Decennial Censuses -- for example, a county may merge with an
# adjacent county and the Census Bureau will have to adjust the GEOIDs for the
# affected Census units within the relevant counties. However, these GEOID changes
# may not necessarily be reflected in the shapefiles downloaded from the Census Bureau.
#
# You can download Census data in a variety of ways. Here, we are using the
# tidycensus package, which uses the Census Bureau's API. You must first obtain 
# an API key. Visit https://www.census.gov/data/developers.html and select the
# large "Request a Key" image on the left side of the screen to request an API key.
# Enter your API key at the top of the script under "USER-DEFINED PARAMETERS";
# it is assigned via the census_api_key() function.
# Block populations are only available for Decennial Censuses.
#
decennial_year <- ifelse(year %in% 2000:2009, 2000,
                         ifelse(year %in% 2010:2019, 2010,
                                ifelse(year %in% 2020:2029, 2020, NA)))

for (i in 1:length(states)) {
  
  cat("----------------------------------------------------------------------\n")
  cat("Processing data for state FIPS", states[i], "\n")
  
  counties <- tigris::counties(state = states[i], year = year, cb = TRUE)
  counties <- unique(counties[[grep("^COUNTYFP", names(counties), ignore.case = TRUE)]])
  
  pop_var_name <- "P1_001N"
  geography <- "block"
  
  blockpop <- get_decennial(geography = geography,
                            variables = pop_var_name, # Use load_variables("pl", year = 2020) to see available vars
                            year = decennial_year,
                            state = states[i],
                            county = counties,
                            sumfile = "pl",
                            key = Sys.getenv("CENSUS_API_KEY"))
  
  saveRDS(blockpop, paste0(input_data_dir, "Census_", pop_var_name, "_", geography, "_", decennial_year, "_", states[i], ".Rds"))
}

