import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

class BaseModel:
    def __init__(self, initial_conditions, parameters):
        """
        Initialize the base model with initial conditions and parameters.
        
        Parameters:
        - initial_conditions: List of initial conditions for the model compartments.
        - parameters: Dictionary of parameters specific to the model.
        """
        self.initial_conditions = initial_conditions
        self.parameters = parameters

    def model_equations(self, t, y):
        """
        Defines the model's differential equations. Should be overridden by subclasses.
        
        Parameters:
        - t: Time variable
        - y: A list of compartment values at time t
        """
        raise NotImplementedError("Subclass must implement abstract method")

    def solve(self, t_span, t_eval):
        """
        Solves the model over a given time span.
        
        Parameters:
        - t_span: 2-tuple of (start time, end time)
        - t_eval: Array of time points at which to store the computed solutions
        
        Returns:
        - Object returned by scipy.integrate.solve_ivp containing the solution.
        """
        solution = solve_ivp(self.model_equations, t_span, self.initial_conditions, args=(self.parameters,), t_eval=t_eval)
        return solution

    def plot_results(self, solution, compartment_labels):
        """
        Plots the results of the model simulation.
        
        Parameters:
        - solution: Object returned by scipy.integrate.solve_ivp containing the solution.
        - compartment_labels: Labels for the compartments to be plotted.
        """
        plt.plot(solution.t, solution.y.T)
        plt.legend(compartment_labels)
        plt.xlabel('Time / days')
        plt.ylabel('Fraction of population')
        plt.title('Model Dynamics')
        plt.grid(True)
        plt.show()

class SEIRModel(BaseModel):
    def __init__(self, beta, sigma, gamma, S0, E0, I0, R0):
        parameters = {'beta': beta, 'sigma': sigma, 'gamma': gamma}
        initial_conditions = [S0, E0, I0, R0]
        super().__init__(initial_conditions, parameters)

    def model_equations(self, t, y, parameters):
        S, E, I, R = y
        beta, sigma, gamma = parameters.values()
        dSdt = -beta * S * I
        dEdt = beta * S * I - sigma * E
        dIdt = sigma * E - gamma * I
        dRdt = gamma * I
        return [dSdt, dEdt, dIdt, dRdt]

class SEEIIRModel(BaseModel):
    def __init__(self, beta, sigma, gamma, mu, S0, E1_0, E2_0, I1_0, I2_0, R1_0, R2_0):
        parameters = {'beta': beta, 'sigma': sigma, 'gamma': gamma, 'mu': mu}
        initial_conditions = [S0, E1_0, E2_0, I1_0, I2_0, R1_0, R2_0]
        super().__init__(initial_conditions, parameters)

    def model_equations(self, t, y, parameters):
        S, E1, E2, I1, I2, R1, R2 = y
        beta, sigma, gamma, mu = parameters.values()
        dSdt = -beta * S * (I1 + I2)
        dE1dt = beta * S * (I1 + I2) - sigma * E1
        dE2dt = sigma * E1 - sigma * E2
        dI1dt = sigma * E2 - gamma * I1
        dI2dt = mu * I1 - gamma * I2
        dR1dt = (1 - mu) * I1 + gamma * I2
        dR2dt = gamma * I2
        return [dSdt, dE1dt, dE2dt, dI1dt, dI2dt, dR1dt, dR2dt]

if __name__ == '__main__':
    # Example usage for SEIRModel
    seir_model = SEIRModel(beta=0.5, sigma=1/6, gamma=1/10, S0=0.99, E0=0.01, I0=0.0, R0=0.0)
    seir_solution = seir_model.solve(t_span=[0, 160], t_eval=np.linspace(0, 160, 400))
    seir_model.plot_results(seir_solution, ['Susceptible', 'Exposed', 'Infectious', 'Recovered'])

    # Example usage for SEEIIRModel
    seeiir_model = SEEIIRModel(beta=0.5, sigma=1/6, gamma=1/10, mu=0.1, S0=0.99, E1_0=0.005, E2_0=0.005, I1_0=0.0, I2_0=0.0, R1_0=0.0, R2_0=0.0)
    seeiir_solution = seeiir_model.solve(t_span=[0, 160], t_eval=np.linspace(0, 160, 400))
    seeiir_model.plot_results(seeiir_solution, ['Susceptible', 'Exposed 1', 'Exposed 2', 'Infectious 1', 'Infectious 2', 'Recovered 1', 'Recovered 2'])
