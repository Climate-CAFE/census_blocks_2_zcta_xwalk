# blocks2zcta_xwalk
R scripts to generate a nation wide census block to zcta crosswalk

## Background

Census blocks are the smallest enumeration units defined by the US Census Bureau, and they are the only administrative unit that "nests" within **every** other administrative unit. That is, census blocks can neatly be aggregated up to any census geography (e.g., census tracts, counties, states, congressional districts, etc.). 

ZIP Code Tabulation Areas (ZCTAs) are drawn by the US Census Bureau as areal _approximations_ of ZIP Codes. More information can be found [here](https://github.com/Climate-CAFE/zip_codes_and_zctas). Since ZCTAs can cross counties and even state boundaries, the only sub-unit geography aligning within them is census blocks. Similarly, ZCTAs cannot be _directly_ aggregated up to counties or states without a relationship file to account for the ZCTAs that cross county or state boundaries.

Although census blocks share geographic identifier (GEOID) codes with census block groups, census tracts, counties, and states — the first 12, 11, 5, and 2 characters of the census block GEOID constitute the GEOID for each of these other geographies, respectively — they do *not* share common GEOIDs with ZCTAs. Therefore, to aggregate census-block-level data up to ZCTAs, users must employ a crosswalk/relationship file between them.


### Code Overview

The code in this directory creates a nationwide crosswalk file that shows the relationship between each 2020 census block and the 2020 ZCTA in which it is located. The proportion of the block population relative to the total population of the ZCTA constitutes the "spatial weight" of that block. The spatial weight of these blocks can be used to create population-weighted mean values at the ZCTA level based on block-level data.

### Data Dictionary

The final output contains the following variables:

| **— Variable —** | **— Description —** |
|--------:|:-----------|
| **GEOID20** | Unique GEOID for the 2020 *census block* from the US Census Bureau |
| **ZCTA5CE20** | Unique GEOID for the 2020 *ZCTA* from the US Census Bureau |
| **Pop_Block_2020** | Block population from the 2020 Decennial Census (variable P1_001N) |
| **Pop_ZCTA_2020** | ZCTA population from the 2020 Decennial Census (variable P1_001N) |
| **Weight** | Weighting term used to population-weight from block to ZCTA. This value represents the proportion of the block's population relative to the total ZCTA population in which it is located |

### How to Use this Crosswalk


#### To calculate _population-weighted means_ of block-level data:

1. Merge your block-level data with this crosswalk by the GEOID20 value. **Make sure your input data has retained the leading 0 of the GEOID**.
2. If your input data contains **any** missing data, you should re-calculate the weights based on the availability of data. For example, if a hypothetical ZCTA is comprised of two blocks (Block A with Weight=0.75 and Block B with Weight=0.25), and your input data has an NA value for Block B, then you need to change the weight of Block A to 1.0 and the weight of Block B to 0. This is most easily achieved by summing the Weight values of the _non-missing_ blocks within each ZCTA and dividing each weight by the total of the available ZCTA-level weights. In the example above, the total sum of weights among the non-missing data is 0.75, so you would calculate the new weights as 0.75/0.75 = 1.0 for Block A and 0/0.75 = 0.0 for Block B. See the note on missing data in the Important Considerations section below.
3. Multiply your block-level input data value(s) by the Weight
4. Sum these new values within unique ZCTAs

#### To calculate _total counts_ of **enumerable** data:

- Simply merge with the crosswalk by the unique GEOID and sum within the unique ZCTAs. If you are ultimately calculating proportions at the ZCTA level, **be sure that the denominator reflects only the non-missing blocks in your input data**. This is similar to what is described above when calculating population-weighted means.

### Important Considerations

- **Missing data**: When using this crosswalk, you should consider defining a threshold of missingness in your block-level input data that is tolerable at the ZCTA level. For example, you might decide that you only want to have a ZCTA average or total if at least 50% of the block-level data is non-missing. The user needs to decide what is an acceptable amount of missingness in their particular data that they feel can still be representative of the ZCTA overall. **If the block-level data are missing systematically (as opposed to randomly), then your ZCTA level averages will likely be biased.**
- Aggregating _land-area_ variables: the Weight variable in this crosswalk reflects population proportions rather than land area. Users should consider whether population-weighting their input data is appropriate. If land-area weighting is needed, the code can be modified to create weights based on the land areas of the census blocks rather than the populations.
- 2020 Census data includes differential privacy, which particularly affects small areal units (such as census blocks) with low population counts. You can read more about the Census Bureau's use of differential privacy [here](https://www.census.gov/programs-surveys/decennial-census/decade/2020/planning-management/process/disclosure-avoidance/differential-privacy.html).

