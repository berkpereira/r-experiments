# load epidemics
library(epidemics)
library(ggplot2)
library(dplyr)
library(socialmixr)

# load contact and population data from socialmixr::polymod
polymod <- socialmixr::polymod
contact_data <- socialmixr::contact_matrix(
  polymod,
  countries = "United Kingdom",
  age.limits = c(0, 20, 40),
  symmetric = TRUE
)

# prepare contact matrix
contact_matrix <- t(contact_data[["matrix"]])

# prepare the demography vector
demography_vector <- contact_data[["demography"]][["population"]]
names(demography_vector) <- rownames(contact_matrix)

# initial conditions: one in every 1 million is infected
initial_i <- 1e-6
initial_conditions <- c(
  S = 1 - initial_i, E = 0, I = initial_i, R = 0, V = 0
)

# build for all age groups
initial_conditions <- rbind(
  initial_conditions,
  initial_conditions,
  initial_conditions
)
rownames(initial_conditions) <- rownames(contact_matrix)


# prepare the population to model as affected by the epidemic
uk_population <- population(
  name = "UK",
  contact_matrix = contact_matrix,
  demography_vector = demography_vector,
  initial_conditions = initial_conditions
)


# an intervention to close schools
close_schools <- intervention(
  type = "contacts",
  time_begin = 200,
  time_end = 260,
  reduction = matrix(c(0.5, 0.01, 0.01))
)

# view the intervention
close_schools
#>
#>  Intervention name: 
#>  Begins at: 
#> [1] 200
#> 
#>  Ends at: 
#> [1] 260
#> 
#>  Reduction: 
#>              Interv. 1
#> Demo. grp. 1      0.50
#> Demo. grp. 2      0.01
#> Demo. grp. 3      0.01


# run an epidemic model using `epidemic()`
output <- model_default_cpp(
  population = uk_population,
  transmissibility = 1.5 / 7.0,
  infectiousness_rate = 1.0 / 3.0,
  recovery_rate = 1.0 / 7.0,
  intervention = list(contacts = close_schools),
  time_end = 600, increment = 1.0
)


# plot figure of epidemic curve
filter(output, compartment %in% c("exposed", "infectious")) %>%
  ggplot(
    aes(
      x = time,
      y = value,
      col = demography_group,
      linetype = compartment
    )
  ) +
  geom_line() +
  annotate(
    geom = "rect",
    xmin = close_schools$time_begin,
    xmax = close_schools$time_end,
    ymin = 0, ymax = 500e3,
    fill = alpha("red", alpha = 0.2),
    lty = "dashed"
  ) +
  annotate(
    geom = "text",
    x = mean(c(close_schools$time_begin, close_schools$time_end)),
    y = 400e3,
    angle = 90,
    label = "School closure"
  ) +
  scale_y_continuous(
    labels = scales::comma
  ) +
  scale_colour_brewer(
    palette = "Dark2",
    name = "Age group"
  ) +
  expand_limits(
    y = c(0, 500e3)
  ) +
  coord_cartesian(
    expand = FALSE
  ) +
  theme_bw() +
  theme(
    legend.position = "top"
  ) +
  labs(
    x = "Simulation time (days)",
    linetype = "Compartment",
    y = "Individuals"
  )