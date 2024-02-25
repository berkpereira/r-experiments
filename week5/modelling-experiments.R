setwd("repos/r-experiments/week5")


library(fluEvidenceSynthesis)
library(pander)
library(ggplot2)
library(dplyr)
library(tidyr)
library(epiparameter)
library(qs)

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
    
    return(fraction.infected)
}



####################################################################################
# POPULATION DATA
####################################################################################

# Separate into non elderly (<65) and elderly population
ag <- stratify_by_age(demography, limits=c(65)) # c( 43670500, 8262600 )

# Separate in risk groups. In this case we assume one additional (high) risk groups and that respectively 
# 1 percent of non elderly and 40% of elderly are in this risk group. It is possible to
# define an additional risk group (e.g. pregnant women) if desired
population <- stratify_by_risk(ag, matrix(c(0.01,0.4),nrow=1), 
                               labels = c("LowRisk", "HighRisk")) # c( 43233795, 4957560, 436705, 3305040 )



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

dates_vector <- as.Date(c("2023-10-01", "2023-11-01", "2023-12-01", "2024-01-01", "2024-02-01", "2024-03-01"))

# note: rbind appends row-by-row
# coverage_matrix = rbind(c(0, 0, 0, 0),
#                         c(0, 0.861, 0.123, 0.861),
#                         c(0, 0.861, 0.123, 0.861),
#                         c(0, 0.861, 0.123, 0.861),
#                         c(0, 0.861, 0.123, 0.861),
#                         c(0, 0.861, 0.123, 0.861))
# 
# vaccine_calendar <- as_vaccination_calendar(
#     efficacy = c(0.7,0.4),
#     coverage = as.data.frame(coverage_matrix),
#     dates = dates_vector,
#     no_age_groups = 2,
#     no_risk_groups = 2
# )





########
# THE LAYOUT IN GENERAL is expected to be all age groups for one risk group,
# followed by the age groups for the second risk group etc



####################################################################################
# VACCINE CALENDAR BASED ON ACTUAL 2022-23 UPTAKE DATA
####################################################################################

# We did the following:
# Took 65+ coverage data, put 1/3 into low-risk and 2/3 of it into the high-risk group in this model
# Under 65 and high-risk is given explicitly in the data

# The picture is significantly more complicated for low-risk individuals under 65 (split into school-age children
# and then 50 to 64-year olds), so we're modelling this at an "eyeballed" low rate

dates_vector <- as.Date(c("2022-09-01", "2022-09-30", "2022-10-31",
                          "2022-11-30", "2022-12-31", "2023-01-31", "2023-02-28"))

coverage_matrix = rbind(c(0, 0, 0, 0),
                        c(0.02, 0.07, 0.067, 0.14),
                        c(0.07, 0.218, 0.308, 0.436),
                        c(0.08, 0.254, 0.427, 0.508),
                        c(0.09, 0.261, 0.463, 0.523),
                        c(0.10, 0.265, 0.486, 0.529),
                        c(0.10, 0.266, 0.491, 0.533))

vaccine_calendar <- as_vaccination_calendar(
    efficacy = c(0.7,0.4),
    coverage = as.data.frame(coverage_matrix),
    dates = dates_vector,
    no_age_groups = 2,
    no_risk_groups = 2
)


####################################################################################
# CONTACT RATES BETWEEN AGE GROUPS
####################################################################################

# Polymod data is subdivided in seven age groups
poly <- polymod_uk[,c(1,2,3,9)]
poly[,3] <- rowSums(polymod_uk[,3:8])

contacts <- contact_matrix(as.matrix(poly), demography, c(65))
pander(contacts)


####################################################################################
# EPIDEMIOLOGICAL PARAMETERS
####################################################################################

susceptibility <- c(0.7, 0.3) # Different for different ages
transmissibility <- 0.14 # Same for all ages
infection_delays <- c( 0.8, 2.5 ) # 0.8 and 1.8 day.



####################################################################################
# ITERATE OVER TRANSMISSIBILITY VALUES
####################################################################################

# Define a range of transmissibility values to iterate over
transmissibility_values <- seq(0.12, 0.16, by = 0.01)


# Now loop over the remaining transmissibility values starting from the second one
for(t in transmissibility_values) {  # Exclude the first element since it's already processed
    model_data <- run_model_vignette(population, initial.infected, vaccine_calendar, contacts,
                                     susceptibility, t, infection_delays, interval)
    plot <- ggplot(model_data) + geom_line( aes(x=Time, y=fraction, colour = Group) ) + 
        ylab( "Fraction infected") + ggtitle(paste('Transmissibility = ', t))
    print(plot)
}






# Save the plot to a PDF file
#ggsave("output.pdf", plot = plot_object, width = 7, height = 5)


