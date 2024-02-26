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


# Implement the vignette modelling runs in an easy function


run_model_vignette <- function(population, initial.infected, vaccine_calendar, contacts, 
                               susceptibility, transmissibility, infection_delays, interval) {
    
    # Run model and store results
    odes <- infectionODEs( population, initial.infected, vaccine_calendar, contacts, 
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
initial.infected <- stratify_by_risk( ag, matrix(c(0.01,0.4),nrow=1)) # c(990, 600, 10, 400)



####################################################################################
# VACCINE CALENDAR
####################################################################################

PLOT_COVERAGE <- T

if (PLOT_COVERAGE) {
    plot_coverage_time_series(dates_vector, coverage_matrix)
}


####################################################################################
# EPIDEMIOLOGICAL PARAMETERS
####################################################################################

susceptibility <- c(0.7, 0.3) # Different for different ages
transmissibility <- 0.14 # Same for all ages
infection_delays <- c( 0.8, 2.5 ) # 0.8 and 1.8 day.



####################################################################################
# ITERATE OVER PARAMETER VALUES
####################################################################################

# Define a range of transmissibility values to iterate over
transmissibility_values <- seq(0.12, 0.16, by = 0.01)






# Save the plot to a PDF file
#ggsave("output.pdf", plot = plot_object, width = 7, height = 5)


