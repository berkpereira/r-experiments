# setwd("repos/r-experiments/week5")


library(fluEvidenceSynthesis)
library(pander)
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(epiparameter)
library(qs)
library(ggplot2)
library(tidyverse)
library(beepr)

data(polymod_uk)
data(demography)





####################################################################################
# LOAD DATA
####################################################################################

# load("inference-data-results.RData")

load("GOOD-INFERENCE.RData")
batch_tibble <- as_tibble(inference_results$batch) # just casting to a more conveniently structured object
mean_inferred_params <- sapply(batch_tibble, mean)


# Values below taken from Baguelin 2013, page 6
# In DAYS
T_latent <- 2 / 2.5
T_infectious <- 2 / 1.1




####################################################################################
# REDEFINE THIS TO HAVE PROPER SERIES LABELS
####################################################################################


plot_coverage_time_series <- function(dates, coverage) {
    # Convert the coverage matrix to a long format data frame suitable for ggplot2
    coverage_df <- as.data.frame(coverage)
    names(coverage_df) <- c("Low Risk [0,15)",
                            "Low Risk [15,65)",
                            "Low Risk [65,+)",
                            "High Risk [0,15)",
                            "High Risk [15,65)",
                            "High Risk [65,+)")
    
    coverage_df$Date <- dates
    long_coverage_df <- reshape2::melt(coverage_df, id.vars = "Date", variable.name = "Series", value.name = "Coverage")
    
    # Plotting
    p <- ggplot(data = long_coverage_df, aes(x = Date, y = Coverage, color = Series)) +
        geom_line() +
        geom_point() +
        theme_minimal() +
        labs(title = "Coverage Over Time", x = "Date", y = "Coverage (%)", color = "Series") +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1)) # Display y-axis labels as percentages
    
    # Customize the date breaks and labels on the x-axis
    p <- p + scale_x_date(date_breaks = "1 month", date_labels = "%b %d")
    # Adjust the date_breaks and date_labels arguments according to your needs
    # "%b %d" will format dates as 'Month day', e.g., 'Jan 01'
    
    print(p)
}



####################################################################################
# INTEGRATION dt
####################################################################################

interval <- 7



####################################################################################
# INITIAL INFECTED
####################################################################################

# same in all age groups
ag <- rep(10^3, 3)

# for each age group, distribute initial cases uniformly across risk groups.
# i.e., risk groups play no role in distributing initial infections.
initial_infected <- stratify_by_risk(ag, matrix(risk_ratios,nrow=1))



####################################################################################
# VACCINE CALENDAR
####################################################################################

PLOT_COVERAGE <- FALSE

if (PLOT_COVERAGE) {
    plot_coverage_time_series(baseline_dates_vector, baseline_coverage_matrix)
}


modify_coverage_data <- function(baseline_dates, baseline_coverage,
                                 start_date_shift = 0, coverage_scaling = 1,
                                 uptake_speedup = 1) {
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
    
    # Adjust dates for uptake speedup before shifting
    if (uptake_speedup != 1 && length(baseline_dates) > 1) {
        # Calculate the differences between consecutive dates
        date_diffs <- diff(baseline_dates)
        # Convert 'difftime' object to numeric days
        date_diffs_days <- as.numeric(date_diffs, units = "days")
        # Adjust the gaps by the uptake_speedup factor
        adjusted_diffs <- date_diffs_days / uptake_speedup
        # Reconstruct the dates vector starting from the first date
        adjusted_dates <- c(baseline_dates[1], baseline_dates[1] + cumsum(adjusted_diffs))
        baseline_dates <- adjusted_dates
    }
    
    # Shift the adjusted dates
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
# DEFINE PLOTTING FUNCTIONS
####################################################################################
normalise_odes <- function(odes, population_vector) {
    normalised_odes <- odes
    for(col_name in names(population_vector)) {
        normalised_odes[[col_name]] <- odes[[col_name]] / population_vector[col_name]
    }
    return(normalised_odes)
}



plot_odes <- function(odes, normalised=FALSE) {
    # Melt the data frames for ggplot
    odes_long <- melt(odes, id.vars = "Time", variable.name = "Group", value.name = "Cases")
    
    # Plot raw numbers
    if (!normalised) {
        ggplot(odes_long, aes(x = Time, y = Cases, color = Group)) +
            geom_line() +
            labs(title = "Weekly New Infection Cases by Group", x = "Week", y = "Number of Cases") +
            theme_minimal()
    } else {
        ggplot(odes_long, aes(x = Time, y = Cases, color = Group)) +
            geom_line() +
            labs(title = "Weekly New Infection Cases by Group, Normalised", x = "Week", y = "Newly Infected Fraction of Group") +
            theme_minimal()
    }
    
}





####################################################################################
# ITERATE OVER PARAMETER VALUES
####################################################################################


PLOT_CASE_SERIES <- TRUE



# Define a range of vaccine calendar changes to iterate over
vaccine_delays   <- seq(from = 0, to = 30, length.out = 7)
vaccine_scalings <- seq(from = 1, to = 1, length.out = 7)
vaccine_speedups <- seq(from = 1, to = 1, length.out = 7)

# Iterate over the index range of vaccine_delays
for(i in 1:length(vaccine_delays)) {
    # Access each element by its index
    
    cov_scaling      <- vaccine_scalings[i]
    delay            <- vaccine_delays[i]
    calendar_speedup <- vaccine_speedups[i]
    
    temp <- modify_coverage_data(baseline_dates_vector,
                                 baseline_coverage_matrix,
                                 coverage_scaling = cov_scaling,
                                 start_date_shift = delay,
                                 uptake_speedup = calendar_speedup)
    
    new_dates_vector <- temp[[1]]
    new_coverage_matrix <- temp[[2]]
    
    new_calendar <- as_vaccination_calendar(efficacy = baseline_vaccine_efficacy,
                                            dates = new_dates_vector,
                                            coverage = new_coverage_matrix,
                                            no_risk_groups = 2,
                                            no_age_groups = 3)
    
    plot_coverage_time_series(new_dates_vector, new_coverage_matrix)
    
    
    # RUN MODEL
    odes <- infectionODEs(population = population,
                          initial_infected = initial_infected,
                          vaccine_calendar = new_calendar,
                          contact_matrix = contacts,
                          susceptibility = c(mean_inferred_params['susceptibility_1'], mean_inferred_params['susceptibility_2'], mean_inferred_params['susceptibility_3']),
                          transmissibility = mean_inferred_params['transmissibility'],
                          infection_delays = c(T_latent, T_infectious),
                          interval = 7)
    
    
    normalised_odes <- normalise_odes(odes, population)
    
    print(plot_odes(normalised_odes, normalised=TRUE))
}








####################################################################################
# SAVE PLOT TO PDF
####################################################################################

SAVE_PLOT = FALSE

if (SAVE_PLOT) {
    ggsave("output.pdf", plot = plot_object, width = 7, height = 5)
}



