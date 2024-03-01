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


strain <- "H3N2"

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




# Values below taken from Baguelin 2013, page 6
# In DAYS
T_latent <- 2 / 2.5
T_infectious <- 2 / 1.1




####################################################################################
# REDEFINE THIS TO HAVE PROPER SERIES LABELS
####################################################################################


age_group_colours <- c("[0,15)" = "#377eb8", "[15,65)" = "#ff7f00", "[65,+)" = "#4daf4a")
risk_group_shapes <- c("Low Risk" = 16, "High Risk" = 17)  # 16 = circle, 17 = triangle
marker_size <- 2



plot_coverage_time_series <- function(dates, coverage, delay=NULL, speedup=NULL, cutoff_date=NULL) {
    coverage_df <- as.data.frame(coverage)
    colnames(coverage_df) <- c("Low Risk [0,15)", "Low Risk [15,65)", "Low Risk [65,+)",
                               "High Risk [0,15)", "High Risk [15,65)", "High Risk [65,+)")
    coverage_df$Date <- dates
    
    # Check and process cutoff_date if provided
    if (!is.null(cutoff_date)) {
        cutoff_date <- as.Date(cutoff_date)
        if(cutoff_date < min(dates)) {
            stop("cutoff_date cannot be earlier than any date in the dates vector.")
        }
        
        # Append a new row with the cutoff_date and the last values if cutoff_date is after the last date
        if(cutoff_date > max(dates)) {
            last_values <- tail(coverage_df, 1)[,-ncol(coverage_df)]
            new_row <- cbind(Date = cutoff_date, last_values)
            coverage_df <- rbind(coverage_df, new_row)
        }
    }
    
    long_coverage_df <- reshape2::melt(coverage_df, id.vars = "Date", variable.name = "Series", value.name = "Coverage")
    long_coverage_df$AgeGroup <- sub(".*(\\[.*\\)).*", "\\1", long_coverage_df$Series)
    long_coverage_df$RiskGroup <- ifelse(grepl("Low Risk", long_coverage_df$Series), "Low Risk", "High Risk")
    
    ggplot(long_coverage_df, aes(x = Date, y = Coverage, color = AgeGroup, shape = RiskGroup)) +
        geom_line(aes(group = Series)) + 
        geom_point(size = 3) + 
        scale_color_manual(values = age_group_colours) +
        scale_shape_manual(values = risk_group_shapes) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
        scale_x_date(limits = c(min(dates), ifelse(is.null(cutoff_date), max(dates), cutoff_date)), 
                     date_breaks = "1 month", date_labels = "%b") +
        labs(title = paste("Vaccine Coverage Over Time", ifelse(!is.null(delay) && !is.null(speedup), 
                                                                paste("\nDelay:", delay, "days, Speedup:", speedup), 
                                                                "")), 
             x = "Date", y = "Coverage (%)") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5),  
              axis.text.x = element_text(angle = 45, hjust = 1))
}

# Adjust plot_odes function
plot_odes <- function(odes, normalised=FALSE, delay=NULL, calendar_speedup=NULL, cutoff_date=NULL, y_max=NULL, labs=TRUE) {
    odes_long <- melt(odes, id.vars = "Time", variable.name = "Group", value.name = "Cases")
    
    # Extract age group and risk group from the 'Group' column


    odes_long$AgeGroup <- sub(".*(\\[.*\\)).*", "\\1", odes_long$Group)
    odes_long$RiskGroup <- ifelse(grepl("LowRisk", odes_long$Group), "Low Risk", "High Risk")
    
    p <- ggplot(odes_long, aes(x = Time, y = Cases, color = AgeGroup, shape = RiskGroup)) +
        geom_line(aes(group = Group)) +
        geom_point(size = marker_size) +
        scale_color_manual(values = age_group_colours) +
        scale_shape_manual(values = risk_group_shapes) +
        labs(title = paste("Weekly New Infections,", if(normalised) {"Normalised"} else {""}, "\n", 
                           if(!is.null(delay) && !is.null(calendar_speedup)) {
                               paste("Delay:", delay, "days, Speedup:", calendar_speedup)
                           } else {""}),
             x = "Date", y = if(normalised) {"Newly Infected Fraction of Group"} else {"Number of Cases"}) +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5),  # Center-align the title
              legend.position = if(labs) "right" else "none",  # Conditionally remove the legend if labs is FALSE
              axis.text.x = element_text(angle = 45, hjust = 1)) +  # Apply 45-degree inclination to x-axis labels
        scale_y_continuous(limits = c(NA, y_max)) +
        scale_x_date(date_breaks = "1 month", date_labels = "%b")  # Generate ticks for each month
    
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






####################################################################################
# ITERATE OVER PARAMETER VALUES
####################################################################################


PLOT_CASE_SERIES <- TRUE

SAVE_PLOT = FALSE


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
    
    
    
    truncation_date = "2023-04-01"
    
    
    
    cov_plot <- plot_coverage_time_series(new_dates_vector, new_coverage_matrix,
                                          cutoff_date = truncation_date,
                                          delay = delay,
                                          speedup = calendar_speedup)
    print(cov_plot)
    
    if (SAVE_PLOT) {
        cov_filename <- paste("coverage-delay", delay,
                              "-speedup", calendar_speedup, ".pdf", sep = "")
        ggsave(cov_filename, plot = cov_plot, width = PLOT_WIDTH, height = PLOT_HEIGHT)
    }

    
    
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
    
    
    
    cases_plot <- plot_odes(normalised_odes, normalised=TRUE, delay=delay,
                            calendar_speedup=calendar_speedup, cutoff_date = truncation_date,
                            y_max=0.009, labs = FALSE)
    
    print(cases_plot)
    
    if (SAVE_PLOT) {
        cases_filename <- paste("cases-delay", delay,
                                "-speedup", calendar_speedup, ".pdf", sep = "")
        ggsave(cases_filename, plot = cases_plot, width = PLOT_WIDTH, height = PLOT_HEIGHT)
    }

    
    # cases <- rowSums(vaccination_scenario(demography = demography,
    #                                       vaccine_calendar = new_calendar,
    #                                       polymod_data = polymod,
    #                                       contact_ids = inference_results$contact.ids,
    #                                       parameters = inference_results$batch,
    #                                       verbose = F))
    # cases_df <- data.frame(value = cases, scenario = "Original")
    # ggplot(data = cases_df) + geom_histogram(aes(x = value), bins = 25)
}



