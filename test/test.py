import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Constants
E = 1.0  # Electric field strength
mu = 0.1  # Dipole moment
eta = 0.01  # Viscosity of the fluid
frequency = 100  # Frequency of AC field

# Time parameters
t_start = 0
t_end = 10
t_points = 1000
time = np.linspace(t_start, t_end, t_points)

# Initial conditions [theta1, omega1, theta2, omega2, r, dr]
# theta: angle, omega: angular velocity, r: distance, dr: rate of change of distance
initial_conditions = [0.1, 0, -0.1, 0, 5, 0]


def model(t, y):
    theta1, omega1, theta2, omega2, r, dr = y

    # DEP torques and forces (simplified model)
    torque1 = mu * E * np.sin(theta1)
    torque2 = mu * E * np.sin(theta2)

    # Equations of motion for rotational dynamics
    dtheta1_dt = omega1
    domega1_dt = -eta * omega1 + torque1
    dtheta2_dt = omega2
    domega2_dt = -eta * omega2 + torque2

    # Coulomb force approximation for translational motion
    force = mu * E / (r ** 2) if r > 0.5 else 0  # Prevent singularity

    # Equations of motion for translational dynamics
    dr_dt = dr
    ddr_dt = force - eta * dr

    return [dtheta1_dt, domega1_dt, dtheta2_dt, domega2_dt, dr_dt, ddr_dt]


# Solve the differential equations
solution = solve_ivp(model, [t_start, t_end], initial_conditions, t_eval=time)

# Plotting the results
plt.figure(figsize=(12, 6))
plt.subplot(2, 1, 1)
plt.plot(time, solution.y[0], label='Theta 1 (Rotation CNT 1)')
plt.plot(time, solution.y[2], label='Theta 2 (Rotation CNT 2)')
plt.xlabel('Time')
plt.ylabel('Rotation angle (rad)')
plt.legend()

plt.subplot(2, 1, 2)
plt.plot(time, solution.y[4], label='Distance between CNTs')
plt.xlabel('Time')
plt.ylabel('Distance')
plt.legend()

plt.tight_layout()
plt.show()
