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




PLOT_WIDTH <- 4
PLOT_HEIGHT <- 4.5


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


plot_coverage_time_series <- function(dates, coverage, delay=NULL, speedup=NULL) {
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
    
    # Determine additional title text based on provided delay and speedup
    additional_title <- ""
    if (!is.null(delay)) {
        additional_title <- paste(additional_title, "Delay:", delay, "days")
    }
    if (!is.null(speedup)) {
        additional_title <- paste(additional_title, "Speedup:", speedup, "x")
    }
    full_title <- paste("Coverage Over Time", additional_title)
    
    # Plotting
    p <- ggplot(data = long_coverage_df, aes(x = Date, y = Coverage, color = Series)) +
        geom_line() +
        geom_point() +
        theme_minimal() +
        labs(title = full_title, x = "Date", y = "Coverage (%)", color = "Series") +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
        scale_x_date(date_breaks = "2 months", date_labels = "%b") +  # Only show month name
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))  # Rotate and adjust size of tick labels
    
    print(p)
}




plot_param_hists <- function(inference_results, all_params = FALSE) {
    # Convert the batch results to a tibble for better handling with tidyverse functions
    batch_tibble <- as_tibble(inference_results$batch)
    
    # Check if all parameters should be included
    if (!all_params) {
        # If not all parameters, select only the specified parameters by their indices
        batch_tibble <- batch_tibble %>%
            select(5, 6, 7, 8)
    }
    
    # Pivot the data to a long format
    batch_long <- batch_tibble %>%
        pivot_longer(cols = everything(), names_to = "Parameter", values_to = "Value")
    
    # Plot histograms for each parameter with adjusted layout if not all parameters are included
    number_of_columns <- if (all_params) 3 else 2  # Adjust the number of columns based on the number of plots
    
    ggplot(batch_long, aes(x = Value)) +
        geom_histogram(bins = 25, fill = "blue", color = "black") +
        facet_wrap(~ Parameter, ncol = number_of_columns, scales = "free") +
        theme_minimal() +
        labs(x = "Parameter Value", y = "Frequency", title = "Histograms of Infererred Parameters")
}


####################################################################################
# INTEGRATION dt
####################################################################################

interval <- 7



####################################################################################
# INITIAL INFECTED
####################################################################################

# same proportion in all age groups.
# E.g., this amounts to 10110 initial infected in the [15, 65) age group.
ag <- 3e-4 * c(population[1] + population[4],
               population[2] + population[5],
               population[3] + population[6])

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



plot_odes <- function(odes, normalised=FALSE, delay=NULL,
                      calendar_speedup=NULL, cutoff_date=NULL, y_max=NULL, labs=TRUE) {
    # Melt the data frames for ggplot
    odes_long <- melt(odes, id.vars = "Time", variable.name = "Group", value.name = "Cases")
    
    # Filter the data if a cutoff_date is provided
    if (!is.null(cutoff_date)) {
        cutoff_date <- as.Date(cutoff_date)  # Ensure cutoff_date is in Date format
        odes_long <- odes_long[odes_long$Time <= cutoff_date, ]
    }
    
    # Dynamically generate the title based on delay and calendar_speedup
    # Insert a newline character (\n) to split the title across two lines
    title_text <- paste("Weekly New Infections",
                        if(normalised) {", Normalised"} else {""}, "\n",
                        if(!is.null(delay) && !is.null(calendar_speedup)) {
                            paste("Delay:", delay, "days, Speedup:", calendar_speedup)
                        } else {""})
    
    # Initialize ggplot object
    p <- ggplot(odes_long, aes(x = Time, y = Cases, color = Group)) +
        geom_line() +  # Draw lines connecting the points
        geom_point() +  # Add points at each data point
        labs(title = title_text, x = "Date", y = if(normalised) {"Newly Infected Fraction of Group"} else {"Number of Cases"}) +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5))  # Center-align the title
    
    # Conditionally set the y-axis limit if y_max is provided
    if (!is.null(y_max)) {
        p <- p + scale_y_continuous(limits = c(NA, y_max))
    }
    
    # Conditionally remove the legend if labs is FALSE
    if (!labs) {
        p <- p + theme(legend.position = "none") +
            theme(plot.margin = unit(c(1,1,1,1), "lines"))  # Adjust plot margins
    }
    
    print(p)
}











####################################################################################
# ITERATE OVER PARAMETER VALUES
####################################################################################


PLOT_CASE_SERIES <- TRUE



# Define a range of vaccine calendar changes to iterate over
vaccine_scalings <- c(1, 1, 1, 1)

vaccine_delays   <- c(0, 30, 0, 30)
vaccine_speedups <- c(1, 1, 1.3, 1.3)

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
    
    cov_plot <- plot_coverage_time_series(new_dates_vector, new_coverage_matrix)
    print(cov_plot)
    
    cov_filename <- paste("coverage-delay", delay,
                          "-speedup", calendar_speedup, ".pdf", sep = "")
    ggsave(cov_filename, plot = cov_plot, width = PLOT_WIDTH, height = PLOT_HEIGHT)
    
    
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
    
    truncation_date = "2023-04-01"
    
    cases_plot <- plot_odes(normalised_odes, normalised=TRUE, delay=delay,
                            calendar_speedup=calendar_speedup, cutoff_date = truncation_date,
                            y_max=0.009, labs = FALSE)
    
    print(cases_plot)
    
    cases_filename <- paste("cases-delay", delay,
                            "-speedup", calendar_speedup, ".pdf", sep = "")
    ggsave(cases_filename, plot = cases_plot, width = PLOT_WIDTH, height = PLOT_HEIGHT)
    
    # cases <- rowSums(vaccination_scenario(demography = demography,
    #                                       vaccine_calendar = new_calendar,
    #                                       polymod_data = polymod,
    #                                       contact_ids = inference_results$contact.ids,
    #                                       parameters = inference_results$batch,
    #                                       verbose = F))
    # cases_df <- data.frame(value = cases, scenario = "Original")
    # ggplot(data = cases_df) + geom_histogram(aes(x = value), bins = 25)
}








####################################################################################
# SAVE PLOT TO PDF
####################################################################################

SAVE_PLOT = FALSE

if (SAVE_PLOT) {
    ggsave("output.pdf", plot = plot_object, width = 7, height = 5)
}



