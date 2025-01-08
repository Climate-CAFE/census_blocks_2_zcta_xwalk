library("sf")
library("tigris")
library("tidycensus")
library("doBy")
library("dplyr")

# This script takes all of the state-level, block-to-ZCTA crosswalks created in 
# script "02..." and combines them into a single nationwide crosswalk. It also checks
# the error logs from the distributed computing to ensure that no errors were thrown
#
# This script has been developed for 2000, 2010, or 2020 decennial census 
# geographies. You can update the year below for your project
#
year <- 2020

input_data_dir <- "In_Dir/"
output_data_dir <- "Out_Dir/State_Level_Files/"
final_output_data_dir <- "Out_Dir/"

# Function for summing values that returns NA if all values are NA and NA values are removed
#
sumfun <- function(x) { ifelse(all(is.na(x)), return(NA), return(sum(x, na.rm = TRUE))) }

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# %%%%%%%%%%%%%%%%%%%%%% MERGE THE OUTPUT FILES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

state_files <- paste0(output_data_dir, list.files(output_data_dir, pattern = paste0("Block_to_ZCTA_", year, ".*.Rds")))

# Automated QC check -- ensure 51 files are present
#
if (length(state_files) != 51) { cat("ERROR: there are not 51 files present \n") } else { cat(":) 51 files present \n") }

national_xwalk <- do.call(rbind, lapply(state_files, readRDS))

# Get applicable variables
#
block_geoid <- names(national_xwalk)[grep("^GEOID", names(national_xwalk))]
zcta_geoid <- names(national_xwalk)[grep("^ZCTA5", names(national_xwalk))]
blockpopvar <- names(national_xwalk)[grep("^Pop_20", names(national_xwalk), ignore.case = TRUE)]

# Automated QC check -- confirm that there are no double-counted blocks & that all ZCTAs are present
#
if (length(unique(national_xwalk[[block_geoid]])) != dim(national_xwalk)[1]) {
  cat("ERROR: duplicate block GEOIDs! \n") } else { cat(":) no duplicate block GEOIDs \n") }

zctas <- st_read(paste0(input_data_dir, "ZCTA_", year, "_US.gpkg"))
zcta_id <- names(zctas)[grep("^ZCTA5", names(zctas))]
zcta_lat_var <- names(zctas)[grep("^INTPTLAT", names(zctas))]

# Subset to 50 states + DC
#
zctas <- zctas[which(substr(zctas[[zcta_id]], 1, 2) != "00" &
                       as.numeric(zctas[[zcta_lat_var]]) > 18.5),]

all_zctas <- unique(zctas[[zcta_id]])

if (length(which( !(all_zctas %in% national_xwalk[[zcta_geoid]]) )) > 0) {
  cat("ERROR: some ZCTAs are missing from the national crosswalk! \n") 
} else { cat(":) all ZCTAs are present in the national crosswalk \n") }

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# %%%%%%%%%%%%%%%%%%%%%% CORRECT SPATIAL WEIGHTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

# Merge the total population of the ZCTAs according to the 2020 Decennial Census
# that were downloaded in the first script. We need to ensure that the sum of the 
# block-level populations within each ZCTA match the ZCTA population counts
#
# Total population variable will vary by Decennial Census
#
pop_var_name <- ifelse(year %in% 1990:2009, "P001001",
                       ifelse(year %in% 2010:2019, "P001001",
                              ifelse(year %in% 2020:2029, "P1_001N", NA)))

geography <- "zcta"
zcta_pop_census <- readRDS(paste0(input_data_dir, "Census_", pop_var_name, "_", geography, "_", year, "_US", ".Rds"))
zcta_pop_census <- zcta_pop_census[,c("GEOID", "value")]
names(zcta_pop_census) <- c(zcta_id, "PopZCTA_Census")

# Merge in the Decennial Census population counts by ZCTA to compare
#
national_xwalk <- merge(national_xwalk, zcta_pop_census, by = zcta_id, all.x = TRUE)
if (!(isTRUE(all.equal(national_xwalk$Pop_ZCTA, national_xwalk$PopZCTA_Census)))) {
  cat("ERROR: populations do not match! \n") } else { cat(":) populations match \n") }

# Population from block aggregation and census doesn't always match due to 
# unassigned blocks. This also will arise if  if there are ZCTAs that had block 
# subpopulations summed separately in different states the ZCTA crosses. For 
# example the Census may record 

# Stratify ZCTA with population in multiple states and adjust weights accordingly
# Subset to address population differences across states
#
# Add state marker 
#
national_xwalk$ST <- substr(national_xwalk[[block_geoid]], 1, 2)

# Assess ZCTAs that are in multiple states - first subset to single state-ZCTA
# instead of having block-level results. Our checks for this set of code
# are only based on the ZCTA population, which will be consistent across 
# blocks within a state
#
national_xwalk_cross_state <- national_xwalk %>% 
  group_by(!!sym(zcta_id), ST) %>% 
  slice(1)

# Assess weight across state. Calculate the sum population for a ZCTA regardless
# of state, then compute the weight for that ZCTA at the state level. We 
# will output these state-level ZCTA weights so that when we process exposure
# based on the state-level ZCTA weight, we can then re-weight to the cross-state
# ZCTA total population
#
national_xwalk_state_weight <- national_xwalk_cross_state %>%
  group_by(!!sym(zcta_id)) %>%
  dplyr::mutate(Pop_ZCTA_Total = sum(Pop_ZCTA),
                State_Wt = Pop_ZCTA / Pop_ZCTA_Total) %>%
  select(!!sym(zcta_id), Pop_ZCTA, Pop_ZCTA_Total, State_Wt, ST)

# Save ZCTAs with state-based weights
#
saveRDS(national_xwalk_state_weight, paste0(final_output_data_dir, "ZCTA_StateWt_", year, "_US.Rds"))

# Restrict to those ZCTAs that appear in multiple states
#
national_xwalk_cross_state <- national_xwalk_cross_state %>%
  group_by(!!sym(zcta_id)) %>%
  dplyr::mutate(n = n()) %>% 
  filter(n > 1) 

# Get full crosswalk for those that have different pop (including all blocks)
#
national_xwalk_cross_state_all <- filter(national_xwalk, !!sym(zcta_id) %in% national_xwalk_cross_state[[zcta_id]])

# In these ZCTAs appearing across states, recalculate the spatial weight based
# on the sum population across states
#
eqn <- as.formula(paste0("Pop_ZCTA", " ~ ", zcta_id))
cross_state_sum_ZCTApop <- doBy::summaryBy(eqn, data = national_xwalk_cross_state, FUN = sumfun)

# Join the sum across ZCTAs to the cross state data
#
national_xwalk_cross_state_all <- left_join(national_xwalk_cross_state_all, cross_state_sum_ZCTApop, by = zcta_id)

# Update weights based on full population
#
national_xwalk_cross_state_all$Spatial_Weight_upd <- national_xwalk_cross_state_all[[blockpopvar]] / national_xwalk_cross_state_all$Pop_ZCTA.sumfun

# Replace columns
#
national_xwalk_cross_state_all$Pop_ZCTA <- national_xwalk_cross_state_all$Pop_ZCTA.sumfun
national_xwalk_cross_state_all$Spatial_Weight <- national_xwalk_cross_state_all$Spatial_Weight_upd

# Subset to columns needed
#
national_xwalk_cross_state_all <- national_xwalk_cross_state_all[names(national_xwalk)]

# Get data without cross-state summation issue
#
national_xwalk_single_state <- national_xwalk[!national_xwalk[[zcta_id]] %in% national_xwalk_cross_state[[zcta_id]], ]

# Rejoin data together
#
national_xwalk <- rbind(national_xwalk_cross_state_all, national_xwalk_single_state)

# After addressing ZCTAs across states is there still a population difference?
#
if (!(isTRUE(all.equal(national_xwalk$Pop_ZCTA, national_xwalk$PopZCTA_Census)))) {
  cat("ERROR: populations do not match! \n") } else { cat(":) populations match \n") }

# Save national crosswalk
#
saveRDS(national_xwalk, paste0(final_output_data_dir, "Block_to_ZCTA_", year, "_US.Rds"))

# For some years, a discrepancy may remain betweent the census and block
# aggregate population. Depending on your research question and geographic
# extent, the level error introduced by this discrepancy may require 
# modifying the block -- ZCTA linkage process. The discrepancies can be further
# characterized using the code below.
#
# If the populations match, this code is not needed
#
check_national_xwalk_pop <- national_xwalk[national_xwalk$Pop_ZCTA != national_xwalk$PopZCTA_Census, ]

# Check how pervasive this is
#
length(unique(national_xwalk[[zcta_id]]))           # Total unique ZCTAs
length(unique(check_national_xwalk_pop[[zcta_id]])) # Unique ZCTAs with mismatch in block agg. and census population

# How large are these differences?
#
# Subset to single ZCTA (instead of having all blocks included for each ZCTA)
#
check_national_xwalk_pop_slice <- check_national_xwalk_pop %>%
  group_by(!!sym(zcta_id)) %>%
  slice(1) %>%
  mutate(prop_diff = PopZCTA_Census / Pop_ZCTA,
         tot_diff = PopZCTA_Census - Pop_ZCTA) %>%
  filter(!is.na(PopZCTA_Census))

# Note: This will cause errors if the populations match!
# Plot Census v. Block Aggregated population and assess correlation
#
plot(check_national_xwalk_pop_slice$Pop_ZCTA, check_national_xwalk_pop_slice$PopZCTA_Census, 
     xlab = "ZCTA Block Agg Pop",
     ylab = "ZCTA Census Pop")
cor(check_national_xwalk_pop_slice$Pop_ZCTA, check_national_xwalk_pop_slice$PopZCTA_Census)
