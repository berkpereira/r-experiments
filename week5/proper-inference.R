# setwd("repos/r-experiments/week5")
getwd()  # Displays the current working directory

library(fluEvidenceSynthesis)
library(pander)
library(ggplot2)
library(dplyr)
library(tidyr)
library(epiparameter)
library(qs)

data("demography")
data("polymod_uk")


################################################################################
# SET STRAIN FOR INFERENCE
################################################################################
# Either H3N2, B, or H1N1
strain <- "H3N2"



################################################################################
# PREPARE THE DATA
################################################################################

# De-serialise Edwin-sent data
all_data <- qread("2017-18_EN_viro_ili.qs")
if (strain == "H3N2") {
    infer_data <- all_data$H3N2
} else if (strain == "B") {
    infer_data <- all_data$B
} else if (strain == "H1N1") {
    infer_data <- all_data$H1N1
} else {
    stop("Invalid strain. Choose either H3N2, B, or H1N1.")
}

ili <- infer_data$ili
viro <- infer_data$viro

# Group all the example data into two age groups:
polymod <- polymod_uk[,c(1,2)]
polymod[,3] <- rowSums(polymod_uk[,c(3,4,5,6,7,8)])
polymod[,4] <- polymod_uk[,9]

# Initialize ili_df as a list
ili_df <- list()

# Manipulate ili$ili according to the specifications with renamed columns
ili_df$ili <- ili$ili %>%
    transmute(
        Date = Date,
        `[0,13)` = rowSums(select(., `[0,2)`, `[2,5)`, `[5,13)`), na.rm = TRUE),
        `[13,65)` = rowSums(select(., `[13,17)`, `[17,45)`, `[45,65)`), na.rm = TRUE),
        `[65,+)` = `[65,+)`
    )

# Manipulate ili$mon_pop according to the same specifications with renamed columns
ili_df$mon_pop <- ili$mon_pop %>%
    transmute(
        Date = Date,
        `[0,13)` = rowSums(select(., `[0,2)`, `[2,5)`, `[5,13)`), na.rm = TRUE),
        `[13,65)` = rowSums(select(., `[13,17)`, `[17,45)`, `[45,65)`), na.rm = TRUE),
        `[65,+)` = `[65,+)`
    )


# Assuming viro_df is the desired result structure
viro_df <- list()

# Manipulate viro$positive according to the specifications
viro_df$positive <- viro$positive %>%
    transmute(
        Date = Date,
        `[0,13)` = rowSums(select(., `[0,2)`, `[2,5)`, `[5,13)`), na.rm = TRUE),
        `[13,65)` = rowSums(select(., `[13,17)`, `[17,45)`, `[45,65)`), na.rm = TRUE),
        `[65,+)` = `[65,+)`
    )

# Manipulate viro$total according to the same specifications
viro_df$total <- viro$total %>%
    transmute(
        Date = Date,
        `[0,13)` = rowSums(select(., `[0,2)`, `[2,5)`, `[5,13)`), na.rm = TRUE),
        `[13,65)` = rowSums(select(., `[13,17)`, `[17,45)`, `[45,65)`), na.rm = TRUE),
        `[65,+)` = `[65,+)`
)



################################################################################
# DEMOGRAPHY AND CONTACTS DATA
################################################################################
polymod <- polymod_uk[,c(1,2,3)]
polymod[,3] <- rowSums(polymod_uk[,c(3,4,5)])
polymod[,4] <- rowSums(polymod_uk[,c(6,7,8)])
polymod[,5] <- polymod_uk[,9]
names(polymod)[3:5] <- c("[0,15)", "[15,65)", "[65,+)")

rm(polymod_uk) # Get rid of this to avoid confusion

contacts <- contact_matrix(as.matrix(polymod), demography, c(15,65))

# Separate into 3 age groups
ag <- stratify_by_age(demography, limits=c(15, 65)) # c( 43670500, 8262600 )

# THE BELOW TAKES NUMBERS FROM EDWIN'S VIGNETTES
# Separate in risk groups. In this case we assume one additional (high) risk groups and that respectively 
# 1 percent of non elderly and 40% of elderly (65+) are in this risk group.
population <- stratify_by_risk(ag, matrix(c(0.01, 0.01, 0.4),nrow=1), 
                               labels = c("LowRisk", "HighRisk"))


################################################################################
# VACCINATION CALENDAR
################################################################################

# LAYOUT OF COVERAGE MATRIX:
# low risk young, middle, old
# high risk young, middle, old


# The coverage data below is somewhat loosely based on the real coverage data
# for the 2022-23 season.

dates_vector <- as.Date(c("2022-09-01", "2022-09-30", "2022-10-31",
                          "2022-11-30", "2022-12-31", "2023-01-31", "2023-02-28"))

coverage_matrix = rbind(c(0, 0, 0, 0, 0, 0),
                        c(0.015, 0.013,  0.07, 0.030, 0.067, 0.14),
                        c(0.073, 0.061, 0.218, 0.146, 0.308, 0.436),
                        c(0.133, 0.085, 0.254, 0.266, 0.427, 0.508),
                        c(0.150, 0.093, 0.261, 0.300, 0.463, 0.523),
                        c(0.155, 0.097, 0.265, 0.315, 0.486, 0.529),
                        c(0.166, 0.098, 0.266, 0.333, 0.491, 0.533))

# VACCINE EFFICACY ASSUMPTIONS HELPED BY Baguelin2013.pdf, page 6
vaccine_calendar <- as_vaccination_calendar(
    efficacy = c(0.7, 0.7, 0.46, 0.7, 0.7, 0.46),
    coverage = as.data.frame(coverage_matrix),
    dates = dates_vector,
    no_age_groups = 3,
    no_risk_groups = 2
)
