library(fluEvidenceSynthesis)
library(pander)
library(ggplot2)
library(dplyr)
library(tidyr)
library(epiparameter)

data(polymod_uk)
data(demography)


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
    
    # return(fraction.infected)
    return(odes)
}



####################################################################################
# POPULATION DATA
####################################################################################

population <- c(sum(demography))


####################################################################################
# INTEGRATION dt
####################################################################################

interval <- 7



####################################################################################
# INITIAL INFECTED
####################################################################################

initial.infected <- c(10 ^ 1.5)



####################################################################################
# NO VACCINATION IN SIMPLEST MODEL INSTANCE
####################################################################################


dates_vector <- as.Date(c("2022-09-01", "2023-02-28"))

coverage_matrix = rbind(c(0),
                        c(0.0))

vaccine_calendar <- as_vaccination_calendar(
    efficacy = c(0.0),
    coverage = as.data.frame(coverage_matrix),
    dates = dates_vector,
    no_age_groups = 1,
    no_risk_groups = 1
)


####################################################################################
# CONTACT RATES BETWEEN AGE GROUPS
####################################################################################

# Polymod data is subdivided in seven age groups
# poly <- polymod_uk[,c(1,2,3,9)]
# poly[,3] <- rowSums(polymod_uk[,3:8])
# 
# contacts <- contact_matrix(as.matrix(poly), demography, c(65))

# UNCLEAR HOW TO MAKE THE USUAL WORKFLOW WORK FOR A SINGLE AGE GROUP.
# HELP FROM JASMINA/EDWIN? IN THE MEANTIME SINCE THIS IS FOR A VERY SIMPLE TEST CASE
# ANYWAY, CAN JUST USE SOMETHING WITH THE CORRECT ORDER OF MAGNITUDE.
contacts <- as.matrix(2e-7)


####################################################################################
# EPIDEMIOLOGICAL PARAMETERS
####################################################################################

susceptibility <- c(0.4) # 
transmissibility <- 0.14 # Same for all ages
infection_delays <- c( 0.8, 2.5 ) # 0.8 and 1.8 day.



####################################################################################
# ITERATE OVER TRANSMISSIBILITY VALUES
####################################################################################

# Define a range of transmissibility values to iterate over
transmissibility_values <- seq(0.12, 0.16, by = 0.01)



# ISSUE:
# IT SEEMS THE INTERVAL ARGUMENT HAS A HUGE EFFECT ON THE WAVE THAT IS PREDICTED!
# WHY DO WE USE 7? MY SUSPICION IS THAT THIS IS BECAUSE POLYMOD DATA ACCOUNTS FOR CONTACTS
# OVER A WEEK-LONG PERIOD (EACH SURVEY ENTRY CORRESPONDS TO A WEEK OF CONTACTS), BUT
# MAY BE WORTH CHECKING


# Now loop over the remaining transmissibility values starting from the second one
odes <- run_model_vignette(population, initial.infected, vaccine_calendar, contacts,
                                 susceptibility, transmissibility, infection_delays, interval)

fraction.infected <- odes %>%
    gather(Group, Incidence, -Time) %>%
    mutate(fraction = Incidence/population)

plot <- ggplot(fraction.infected) + geom_line( aes(x=Time, y=fraction, colour = Group) ) + 
    ylab( "Fraction infected") + ggtitle(paste('Transmissibility = ', transmissibility))
print(plot)




# Save the plot to a PDF file
#ggsave("output.pdf", plot = plot_object, width = 7, height = 5)


