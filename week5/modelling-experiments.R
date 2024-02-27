# setwd("repos/r-experiments/week5")


library(fluEvidenceSynthesis)
library(pander)
library(ggplot2)
library(dplyr)
library(tidyr)
library(epiparameter)
library(qs)

data(polymod_uk)
data(demography)

load("inference-data-results.RData")


run_model <- function(population, initial_infected, vaccine_calendar, contacts, 
                               susceptibility, transmissibility, infection_delays, interval) {
    
    # Run model and store results
    odes <- infectionODEs( population, initial_infected, vaccine_calendar, contacts, 
                           susceptibility, transmissibility, infection_delays, interval)
    
    # Compute fraction of each subpopulation that is infective at any point in time
    fraction.infected <- odes %>%
        gather(Group, Incidence, -Time) %>%
        mutate(fraction = Incidence/population[Group])
    
    return(odes)
}


####################################################################################
# INTEGRATION dt
####################################################################################

interval <- 7



####################################################################################
# INITIAL INFECTED
####################################################################################

ag <- c(10 ^ 1.5, 10 ^ 1.5)
initial_infected <- stratify_by_risk( ag, matrix(c(0.01,0.4),nrow=1)) # c(990, 600, 10, 400)



####################################################################################
# VACCINE CALENDAR
####################################################################################

PLOT_COVERAGE <- T

if (PLOT_COVERAGE) {
    plot_coverage_time_series(baseline_dates_vector, baseline_coverage_matrix)
}


modify_coverage_data <- function(baseline_dates, baseline_coverage, start_date_shift = 0, coverage_scaling = 1) {
    # Check if coverage_scaling is within the valid range
    max_coverage_value <- max(baseline_coverage)
    if(coverage_scaling <= 0 || coverage_scaling > 1 / max_coverage_value) {
        stop("coverage_scaling must be between 0 and the reciprocal of the maximum entry in baseline_coverage_matrix!")
    }
    if (start_date_shift < 0) {
        stop("Vaccine calendar shift must be positive (delay)!")
    }
    
    # Convert baseline_dates to Date class if not already
    baseline_dates <- as.Date(baseline_dates)
    
    # Shift the baseline dates
    new_dates <- baseline_dates + start_date_shift
    
    # If start_date_shift is not zero, add an extra date and a row of zeros
    if (start_date_shift != 0) {
        extra_date <- as.Date("2022-09-01") # or any other specific date you need to add
        new_dates <- c(extra_date, new_dates)
        
        # Add a row of zeros at the start of the coverage matrix
        zero_row <- matrix(rep(0, ncol(baseline_coverage)), nrow = 1)
        new_coverage <- rbind(zero_row, baseline_coverage)
    } else {
        new_coverage <- baseline_coverage
    }
    
    # Scale the baseline (or modified) coverage matrix
    new_coverage <- new_coverage * coverage_scaling
    
    # Return a list containing the modified date vector and coverage matrix
    return(list(new_dates, new_coverage))
}




coverage_scaling = 1.1


# Assuming baseline_dates_vector and baseline_coverage_matrix are already defined

temp <- modify_coverage_data(baseline_dates_vector,
                             baseline_coverage_matrix,
                             coverage_scaling = coverage_scaling)
    
new_dates_vector <- temp[[1]]
new_coverage_matrix <- temp[[2]]


new_vaccine_calendar <- as_vaccination_calendar(efficacy = baseline_vaccine_efficacy,
                                                dates = new_dates_vector,
                                                coverage = new_coverage_matrix,
                                                no_age_groups = 3,
                                                no_risk_groups = 2)

if (PLOT_COVERAGE) {
    plot_coverage_time_series(new_dates_vector, new_coverage_matrix)
}



####################################################################################
# EPIDEMIOLOGICAL PARAMETERS
####################################################################################

susceptibility <- c(0.7, 0.3) # Different for different ages
transmissibility <- 0.14 # Same for all ages


# LOOK INTO EPIPARAMETER OR OTHER LITERATURE TO OBTAIN THESE VALUES HERE
# values from paper: https://doi.org/10.1038/nature04017, as also used in
# baguelin2013.pdf (see page 6)
infection_delays <- c( 0.8, 1.8 ) 



####################################################################################
# ITERATE OVER PARAMETER VALUES
####################################################################################

# Define a range of transmissibility values to iterate over
transmissibility_values <- seq(0.12, 0.16, by = 0.01) 






# Save the plot to a PDF file
#ggsave("output.pdf", plot = plot_object, width = 7, height = 5)


