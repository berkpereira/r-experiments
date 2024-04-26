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
library(plotly)
library(latex2exp)
library(extrafont)

data(polymod_uk)
data(demography)

PLOT_WIDTH <- 5.5 # inches
PLOT_HEIGHT <- 3 # inches
FONT_SIZE <- 12 # points

# Set a global theme for all ggplot2 plots
theme_set(theme_minimal(base_size = 12) + 
              theme(
                  plot.title = element_text(size = rel(1.2)),  # 120% of the base size
                  axis.text = element_text(size = rel(0.8)),   # 80% of the base size
                  axis.title = element_text(size = rel(1))     # 100% of the base size
              ))


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


####################################################################################
# MISC FUNCTIONS DEFINITIONS
####################################################################################


age_group_colours <- c("[0,15)" = "#377eb8", "[15,65)" = "#ff7f00", "[65,+)" = "#4daf4a")
risk_group_shapes <- c("Low Risk" = 16, "High Risk" = 17)  # 16 = circle, 17 = triangle
marker_size <- 2



plot_coverage_time_series <- function(dates, coverage, delay=NULL, speedup=NULL,
                                      cutoff_date=NULL, title=FALSE) {
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
        labs(title = if (title) {paste("Vaccine Coverage Over Time", ifelse(!is.null(delay) && !is.null(speedup), 
                                                                paste("\nDelay:", delay, "days, Speedup:", speedup), 
                                                                "")) } else NULL, 
             x = "Date", y = "Coverage (%)") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5),  
              axis.text.x = element_text(angle = 45, hjust = 1))
}

# Adjust plot_odes function
plot_odes <- function(odes, normalised=FALSE, delay=NULL, calendar_speedup=NULL,
                      cutoff_date=NULL, y_max=NULL, labs=TRUE, title=FALSE) {
    odes_long <- melt(odes, id.vars = "Time", variable.name = "Group", value.name = "Cases")
    
    # Extract age group and risk group from the 'Group' column
    odes_long$AgeGroup <- sub(".*(\\[.*\\)).*", "\\1", odes_long$Group)
    odes_long$RiskGroup <- ifelse(grepl("LowRisk", odes_long$Group), "Low Risk", "High Risk")
    
    # Ensure cutoff_date is a Date object and filter the data if provided
    if (!is.null(cutoff_date)) {
        cutoff_date <- as.Date(cutoff_date)
        odes_long <- odes_long[odes_long$Time <= cutoff_date, ]
    }
    
    p <- ggplot(odes_long, aes(x = Time, y = Cases, color = AgeGroup, shape = RiskGroup)) +
        geom_line(aes(group = Group)) +
        geom_point(size = 3) +  # Adjust size as needed
        scale_color_manual(values = age_group_colours) +
        scale_shape_manual(values = risk_group_shapes) +
        labs(title = if (title) {paste("Weekly New Infections,", if(normalised) {"Normalised"} else {""}, "\n", 
                           if(!is.null(delay) && !is.null(calendar_speedup)) {
                               paste("Delay:", delay, "days, Speedup:", calendar_speedup)
                           } else {""})} else NULL,
             x = "Date", y = if(normalised) {"Newly Infected Fraction of Group"} else {"Number of Cases"}) +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5),  # Center-align the title
              legend.position = if(labs) "right" else "none",  # Conditionally remove the legend if labs is FALSE
              axis.text.x = element_text(angle = 45, hjust = 1)) +  # Apply 45-degree inclination to x-axis labels
        scale_y_continuous(limits = c(NA, y_max)) +
        scale_x_date(limits = c(min(odes_long$Time), cutoff_date), date_breaks = "1 month", date_labels = "%b")  # Adjusted for cutoff_date
    
    print(p)
}

total_cases <- function(odes) {
    # Ensure 'odes' is a data frame and remove the 'Time' column
    odes <- as.data.frame(odes)
    odes$Time <- NULL
    
    # Sum all values for each group/column and round up to the nearest integer
    total_cases_sum <- ceiling(colSums(odes, na.rm = TRUE))
    
    # Create a new data frame with the rounded summed values
    total_cases_output <- as.data.frame(t(total_cases_sum))
    colnames(total_cases_output) <- names(total_cases_sum)
    
    return(total_cases_output)
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

peak_dates <- function(odes) {
    # Exclude the 'Time' column for calculation
    data_without_time <- odes[, -1]
    
    # Initialize an empty vector to store peak dates
    peak_dates_vec <- vector("list", length = ncol(data_without_time))
    
    # Iterate over each column to find the index of the max value and map it to its date
    for (i in seq_along(data_without_time)) {
        max_index <- which.max(data_without_time[[i]])
        peak_dates_vec[[i]] <- odes$Time[max_index]
    }
    
    # Return the list of peak dates
    names(peak_dates_vec) <- colnames(data_without_time)
    return(peak_dates_vec)
}

peak_cases <- function(odes) {
    # Exclude the 'Time' column for calculation
    data_without_time <- odes[, -1]
    
    # Apply max function column-wise and return the result
    peak_cases_vec <- apply(data_without_time, 2, max, na.rm = TRUE)
    return(peak_cases_vec)
}

sensitivity_contours <- function(matrix, delay_vector, speedup_vector, title = "Contour Plot") {
    
    # CAREFUL TO FLIP X Y RELATIVE TO ROW, COL IN THE DATA MATRIX INDICES!!! 
    fig <- plot_ly(
        y = delay_vector, 
        x = speedup_vector, 
        z = matrix, 
        type = "contour",
        colorscale = 'Hot',
        contours = list(showlabels = TRUE)
    ) %>% layout(title = title,
                 xaxis = list(title = "Vaccination Calendar Speedup"),
                 yaxis = list(title = "Vaccination Calendar Delay (days)"))
    
    return(fig)
}

# Function to generate contour plots with ggplot2 and aesthetic adjustments
# Function to generate contour plots with ggplot2
sensitivity_contours_ggplot <- function(matrix, delay_vector, speedup_vector, title = "Contour Plot") {
    data_melted <- melt(matrix)
    
    names(data_melted) <- c("X1", "Y1", "value")
    
    data_melted$X <- delay_vector[data_melted$X1]
    data_melted$Y <- speedup_vector[data_melted$Y1]
    
    # Generate the plot using the corrected data frame
    p <- ggplot(data_melted, aes(x = Y, y = X, z = value)) +
        + stat_contour(aes(fill=..level..))
    
    print(p)
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

# FUNCTION TO GENERATE MODIFIED VACCINATION CALENDARS
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


####################################################################################
# EPIDEMIOLOGICAL PARAMETERS
####################################################################################

# susceptibility <- c(0.7, 0.3) # (0.7, 0.3) from MODELLING VIGNETTE 
susceptibility <- c(0.57, 0.72, 0.56) # Different for different ages
transmissibility <- 0.14 # 0.14 taken from leeuwen2023 supplementary material


# LOOK INTO EPIPARAMETER OR OTHER LITERATURE TO OBTAIN THESE VALUES HERE
# values from paper: https://doi.org/10.1038/nature04017, as also used in
# baguelin2013.pdf (see page 6)
infection_delays <- c( 2/2.5, 2/1.1 ) 



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
# CASE SUMMARY DATA SENSITIVITY TO VACCINATION CALENDAR DELAY AND SPEEDUP
####################################################################################

# THIS FUNCTION runs a few vaccination scenarios for exlpicit plotting of each calendar
# and resulting epidemic wave
# NOTE: vaccine_delays in days; speedups in factor stretching (applied before the delay, of course);
# vaccine_scalings is the vertical scaling factor (vertical when looking at a coverage time series plot)
run_calendar_scenarios <- function(vaccine_delays, vaccine_speedups, vaccine_scalings, infection_delays,
                                   susceptibility, transmissibility, initial_infected, savefig) {
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
        
        # Show coverage plot for all groups and save plot if requested
        if (savefig) {
            cov_plot <- plot_coverage_time_series(new_dates_vector, new_coverage_matrix,
                                                  cutoff_date = truncation_date,
                                                  delay = delay,
                                                  speedup = calendar_speedup,
                                                  title=FALSE)
            cov_filename <- paste("coverage-delay", delay,
                                  "-speedup", calendar_speedup, ".pdf", sep = "")
            ggsave(cov_filename, plot = cov_plot, width = PLOT_WIDTH, height = PLOT_HEIGHT)
        } 
        else {
            cov_plot <- plot_coverage_time_series(new_dates_vector, new_coverage_matrix,
                                                  cutoff_date = truncation_date,
                                                  delay = delay,
                                                  speedup = calendar_speedup,
                                                  title=TRUE)
        }
        
        print(cov_plot)
        
        
        # Run model and normalise results
        odes <- infectionODEs(population = population,
                              initial_infected = initial_infected,
                              vaccine_calendar = new_calendar,
                              contact_matrix = contacts,
                              susceptibility = c(mean_inferred_params['susceptibility_1'], mean_inferred_params['susceptibility_2'], mean_inferred_params['susceptibility_3']),
                              transmissibility = mean_inferred_params['transmissibility'],
                              infection_delays = infection_delays,
                              interval = 7)
        
        normalised_odes <- normalise_odes(odes, population)
        
        if (savefig) {
            cases_plot <- plot_odes(normalised_odes, normalised=TRUE, delay=delay,
                                    calendar_speedup=calendar_speedup, cutoff_date = truncation_date,
                                    y_max=0.009, labs = TRUE, title=FALSE)
            cases_filename <- paste("cases-delay", delay,
                                    "-speedup", calendar_speedup, ".pdf", sep = "")
            ggsave(cases_filename, plot = cases_plot, width = PLOT_WIDTH, height = PLOT_HEIGHT)
        }
        else {
            cases_plot <- plot_odes(normalised_odes, normalised=TRUE, delay=delay,
                                    calendar_speedup=calendar_speedup, cutoff_date = truncation_date,
                                    y_max=0.009, labs = TRUE, title=TRUE)
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
    
}


# THIS FUNCTION is used in the big sensitivity analysis which we then
# went on to plot in MATLAB
calendar_sensitivity <- function(delay_vector, speedup_vector, infection_delays) {
    if (length(delay_vector) != length(speedup_vector)) {
        stop("Delay and speedup vectors must have the same length!")
    }
    
    # Initialize matrices and lists for storing the results
    total_ili_numbers_list <- list()
    peak_case_dates_list   <- list()
    peak_case_numbers_list <- list()
    
    # Initialize lists for group names
    group_names <- names(population)
    
    for (group in group_names) {
        total_ili_numbers_list[[group]] <- matrix(nrow = length(delay_vector), ncol = length(speedup_vector))
        peak_case_dates_list[[group]]   <- matrix(nrow = length(delay_vector), ncol = length(speedup_vector))
        peak_case_numbers_list[[group]] <- matrix(nrow = length(delay_vector), ncol = length(speedup_vector))
    }
    
    for(i in seq_along(delay_vector)) {
        for (j in seq_along(speedup_vector)) {
            # Create the vaccination calendar
            delay            <- delay_vector[i]
            calendar_speedup <- speedup_vector[j]
            
            temp <- modify_coverage_data(baseline_dates = baseline_dates_vector,
                                         baseline_coverage = baseline_coverage_matrix,
                                         start_date_shift = delay,
                                         uptake_speedup = calendar_speedup)
            
            new_dates_vector <- temp[[1]]
            new_coverage_matrix <- temp[[2]]
            
            
            new_calendar <- as_vaccination_calendar(efficacy = baseline_vaccine_efficacy,
                                                    dates = new_dates_vector,
                                                    coverage = new_coverage_matrix,
                                                    no_risk_groups = 2,
                                                    no_age_groups = 3)
            
            # RUN MODEL
            odes <- infectionODEs(population = population,
                                  initial_infected = initial_infected,
                                  vaccine_calendar = new_calendar,
                                  contact_matrix = contacts,
                                  susceptibility = c(mean_inferred_params['susceptibility_1'], mean_inferred_params['susceptibility_2'], mean_inferred_params['susceptibility_3']),
                                  transmissibility = mean_inferred_params['transmissibility'],
                                  infection_delays = infection_delays,
                                  interval = 7)
            
            # Calculate the results for each combination of delay and speedup
            # Placeholder for actual function calls
            total_ili_numbers <- total_cases(odes)  # Placeholder
            peak_case_dates   <- peak_dates(odes)   # Placeholder
            peak_case_numbers <- peak_cases(odes)   # Placeholder
            
            # Store the results
            for (k in 1:length(group_names)) {
                group_name <- group_names[k]
                
                total_ili_numbers_list[[group_name]][i, j] <- total_ili_numbers[[group_name]]
                peak_case_dates_list[[group_name]][i, j]   <- peak_case_dates[[group_name]]
                peak_case_numbers_list[[group_name]][i, j] <- peak_case_numbers[[group_name]]
            }
        }
    }
    
    # Consolidate results into a list
    results <- list(
        total_ili_numbers = total_ili_numbers_list,
        peak_case_dates   = peak_case_dates_list,
        peak_case_numbers = peak_case_numbers_list
    )
    
    names(results) <- c("TotalIliNumbers", "PeakCaseDates", "PeakCaseNumbers")
    
    return(results)
}




######################## RUN STUFF PROPER


vaccine_delays   <- c(0, 30, 0, 30)
vaccine_speedups <- c(1, 1, 1.3, 1.3)
vaccine_scalings <- c(1, 1, 1, 1) # keep it to ones, not very realistic to upscale as per Jasmina

PLOT_BASELINE_COVERAGE <- TRUE

SAVE_PLOT <- TRUE

SENSITIVITY_ANALYSIS <- FALSE

RUN_CALENDAR_SCENARIOS <- TRUE

# SHOULD IMPROVE FILE NAMING FOR MATLAB BEFORE SAVING ANY FURTHER STUFF!
SAVE_FOR_MATLAB <- FALSE

# Plots baseline vaccine calendar coverage
if (PLOT_BASELINE_COVERAGE) {
    plot_coverage_time_series(baseline_dates_vector, baseline_coverage_matrix)
}


# THIS RUN runs through the few select cases for explicit plotting of vaccination
# calendars and consequent epidemic curves.
if (RUN_CALENDAR_SCENARIOS) {
    run_calendar_scenarios(vaccine_delays=vaccine_delays,
                           vaccine_speedups=vaccine_speedups,
                           vaccine_scalings=vaccine_scalings,
                           infection_delays=infection_delays,
                           susceptibility=susceptibility,
                           transmissibility=transmissibility,
                           initial_infected=initial_infected,
                           savefig=SAVE_PLOT)
}

# THIS RUN makes a large run of many scenarios and stores only summary data
# e.g., peak cases or total cases in a group.
if (SENSITIVITY_ANALYSIS) {
    delay_vector   <- seq(0, 60, length.out = 50)
    speedup_vector <- seq(0.8, 2, length.out = 50)
    
    calendar_sensitivity_results <- calendar_sensitivity(delay_vector=delay_vector,
                                                         speedup_vector=speedup_vector,
                                                         infection_delays=infection_delays,
                                                         )
    
    # Plot for elderly or the high risk
    contour_group <- "elderly"
    quantity <- "PeakCaseNumbers"
    
    if (contour_group == 'elderly') {
        if (quantity == 'TotalIliNumbers') {
            matrix_to_plot <- calendar_sensitivity_results$TotalIliNumbers[["LowRisk [65,+)"]] +
                calendar_sensitivity_results$TotalIliNumbers[["HighRisk [65,+)"]]
        } else if (quantity == 'PeakCaseNumbers') {
            matrix_to_plot <- calendar_sensitivity_results$PeakCaseNumbers[["LowRisk [65,+)"]] +
                calendar_sensitivity_results$PeakCaseNumbers[["HighRisk [65,+)"]]
        }
    } else if (contour_group == 'high-risk') {
        if (quantity == 'TotalIliNumbers') {
            matrix_to_plot <- calendar_sensitivity_results$TotalIliNumbers[["HighRisk [0,15)"]] +
                calendar_sensitivity_results$TotalIliNumbers[["HighRisk [15,65)"]] +
                calendar_sensitivity_results$TotalIliNumbers[["HighRisk [65,+)"]]
        }
    }
    

    if (SAVE_FOR_MATLAB) {
        
        datamatrix_filename <- paste("~/Documents/MATLAB/oxford/case-study-modelling/data_matrix-",
                                     quantity, "-", contour_group, ".csv", sep = "")
        delay_vec_filename <- paste("~/Documents/MATLAB/oxford/case-study-modelling/delay_vector-",
                                    quantity, "-", contour_group, ".csv", sep = "")
        speedup_vec_filename <- paste("~/Documents/MATLAB/oxford/case-study-modelling/speedup_vector-",
                                      quantity, "-", contour_group, ".csv", sep = "")
        
        # Save the matrix
        write.csv(matrix_to_plot, file = datamatrix_filename, sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)
        
        # Save the delay vector
        write.csv(delay_vector, file = delay_vec_filename, row.names = FALSE, col.names = FALSE)
        
        # Save the speedup vector
        write.csv(speedup_vector, file = speedup_vec_filename, row.names = FALSE, col.names = FALSE)
    }
    
    # CONTOUR STUFF BELOW IS FAIRLY BROKEN, SO I DID IT IN MATLAB INSTEAD
    # contour_fig <- sensitivity_contours_ggplot(matrix_to_plot, delay_vector, speedup_vector,
    #                     title = paste("Total Season Cases in ", contour_group, " Population", sep = ""))
    
    # print(contour_fig)
    
    if (SAVE_CONTOUR) {
        contour_filename <- paste("contour-group", contour_group,
                                  "-quantity-", quantity, ".pdf", sep = "")
        save_image(contour_fig, file=contour_filename, width = PLOT_WIDTH, height = PLOT_HEIGHT)
    }
}
