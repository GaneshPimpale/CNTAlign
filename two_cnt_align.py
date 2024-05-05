import math
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm


class TwoCNT:
    def __init__(self):
        # Simulation settings
        self.timesteps = 30000
        self.timesteps_per_second = 100
        self.time = self.timesteps_per_second**(-1)

        # CNT i initial settings
        self.a_i = 1.5*10**(-6)  #static
        self.b_i = 0.5*10**(-9)  #static
        self.theta_i = np.zeros((self.timesteps + 2, 1))
        self.x_i = np.zeros((self.timesteps + 2, 1))
        self.y_i = np.zeros((self.timesteps + 2, 1))
        self.theta_i[0:2] = np.pi/6
        self.x_i[0:2] = 0.25*10**(-4)
        self.y_i[0:2] = 0.5*10**(-4)
        
        # CNT j initial settings
        self.a_j = 1.5*10**(-6)  #static
        self.b_j = 0.5*10**(-9)  #static
        self.theta_j = np.zeros((self.timesteps + 2, 1))
        self.x_j = np.zeros((self.timesteps + 2, 1))
        self.y_j = np.zeros((self.timesteps + 2, 1))
        self.theta_j[0:2] = np.pi
        self.x_j[0:2] = 0.75*10**(-4)
        self.y_j[0:2] = 0.5*10**(-4)
        
        # CNT properties
        self.cnt_permittivity = 1*10**(5)
        self.cnt_conductivity = 1*10**(4)  #S/m

        # Fluid and interphase properties
        self.viscosity = 25  #mPa*s
        self.fluid_permittivity = 11
        self.fluid_conductivity = 3*10**(-5)  #S/m
        self.interphase_permittivity = 1*10**(4)
        self.interphase_conductivity = 1*10**(-4)  #S/m
        self.interphase_thickness = 10*10**(-9)  #m

        # Electric field properties
        self.E_mag = 6.6  #V/mm

        # Derived values
        self.eq_permittivity_i = self.interphase_permittivity*((self.cnt_permittivity+(self.interphase_thickness /
                            (2*self.a_i))*(self.cnt_permittivity-self.interphase_permittivity)) /
                            (self.interphase_permittivity+(self.interphase_thickness/(2*self.a_i)) *
                            (self.cnt_permittivity-self.interphase_permittivity)))
        self.eq_permittivity_j = self.interphase_permittivity*((self.cnt_permittivity+(self.interphase_thickness /
                            (2*self.a_j))*(self.cnt_permittivity-self.interphase_permittivity)) /
                            (self.interphase_permittivity+(self.interphase_thickness/(2*self.a_j)) *
                            (self.cnt_permittivity-self.interphase_permittivity)))
        self.V_i = (4/3) * np.pi * self.a_i * self.b_i**2
        self.V_j = (4/3) * np.pi * self.a_j * self.b_j**2
        self.m_i = 2 * self.b_i * self.V_i
        self.m_j = 2 * self.b_j * self.V_j
        self.I_i = self.m_i * (self.a_i ** 2 + (2 * self.b_i) ** 2) / 5
        self.I_j = self.m_j * (self.a_j ** 2 + (2 * self.b_j) ** 2) / 5

        # Post sim results, #TODO: remove all of this
        self.t_dep_i = []
        self.t_fr_i = []
        
    def show_cnts(self, timestep):
        # CNT i
        start_x_i = self.x_i[timestep] + self.a_i * np.cos(self.theta_i[timestep])
        start_y_i = self.y_i[timestep] + self.a_i * np.sin(self.theta_i[timestep])
        end_x_i = self.x_i[timestep] - self.a_i * np.cos(self.theta_i[timestep])
        end_y_i = self.y_i[timestep] - self.a_i * np.sin(self.theta_i[timestep])
        plt.plot([start_x_i, end_x_i], [start_y_i, end_y_i])

        # CNT j
        start_x_j = self.x_j[timestep] + self.a_j * np.cos(self.theta_j[timestep])
        start_y_j = self.y_j[timestep] + self.a_j * np.sin(self.theta_j[timestep])
        end_x_j = self.x_j[timestep] - self.a_j * np.cos(self.theta_j[timestep])
        end_y_j = self.y_j[timestep] - self.a_j * np.sin(self.theta_j[timestep])
        plt.plot([start_x_j, end_x_j], [start_y_j, end_y_j])
        plt.title('Show CNT at timestep ' + str(timestep))
        print('At timestep ' + str(timestep) + ' CNT i theta: ' + str(self.theta_i[timestep]) +
              ' CNT j theta: ' + str(self.theta_j[timestep]))

        plt.show()

    def show_theta(self):
        steps = range(self.timesteps + 2)
        plt.plot(steps, self.theta_i, label='CNT i theta')
        plt.plot(steps, self.theta_j, label='CNT j theta')
        plt.legend(loc='best')
        plt.xlabel('Timestep')
        plt.ylabel('Angle [radians]')
        plt.title('CNT angle vs time')
        plt.show()
        
    def show_x(self):
        steps = range(self.timesteps + 2)
        plt.plot(steps, self.x_i, label='CNT i x')
        plt.plot(steps, self.x_j, label='CNT j x')
        plt.legend(loc='best')
        plt.xlabel('Timestep')
        plt.ylabel('x [meters]')
        plt.title('CNT x pos vs time')
        plt.show()
        
    def show_y(self):
        steps = range(self.timesteps + 2)
        plt.plot(steps, self.y_i, label='CNT i y')
        plt.plot(steps, self.y_j, label='CNT j y')
        plt.legend(loc='best')
        plt.xlabel('Timestep')
        plt.ylabel('y [meters]')
        plt.title('CNT y pos vs time')
        plt.show()

    def calc_viscosity(self, timestep):
        return self.viscosity * np.exp((1 / (timestep/self.timesteps_per_second - 0.5) ** 2) * np.exp(timestep/self.timesteps_per_second - 3))

    def calculate_r(self, timestep):
        r_i = np.sqrt(self.x_i[timestep] ** 2 + self.y_i[timestep] ** 2)
        r_j = np.sqrt(self.x_j[timestep] ** 2 + self.y_j[timestep] ** 2)
        return r_i, r_j

    def calc_dep_t(self, timestep):
        L_i = (np.log(2 * self.a_i / self.b_i) - 1) / ((self.a_i / self.b_i) ** 2)
        L_j = (np.log(2 * self.a_j / self.b_j) - 1) / ((self.a_j / self.b_j) ** 2)
        alpha_star_i = (((self.eq_permittivity_i-self.fluid_permittivity)**2) /
                      ((self.fluid_permittivity+(self.eq_permittivity_i-self.fluid_permittivity)*L_i) *
                       (self.eq_permittivity_i+self.fluid_permittivity)))
        alpha_star_j = (((self.eq_permittivity_j-self.fluid_permittivity)**2) /
                      ((self.fluid_permittivity+(self.eq_permittivity_j-self.fluid_permittivity)*L_j) *
                       (self.eq_permittivity_j+self.fluid_permittivity)))
        t_dep_i = (1/4)*self.V_i*self.fluid_permittivity*(self.E_mag**2)*alpha_star_i*np.sin(2*self.theta_i[timestep-1])
        t_dep_j = (1/4)*self.V_j*self.fluid_permittivity*(self.E_mag**2)*alpha_star_j*np.sin(2*self.theta_j[timestep-1])
        if math.isnan(t_dep_i) or math.isinf(t_dep_i):
            t_dep_i = 0
        if math.isnan(t_dep_j) or math.isinf(t_dep_j):
            t_dep_j = 0
        return t_dep_i, t_dep_j

    def calc_fr_t(self, timestep):
        theta_dot_i = (self.theta_i[timestep - 1][0] - self.theta_i[timestep - 2][0]) / self.time
        theta_dot_j = (self.theta_j[timestep - 1][0] - self.theta_j[timestep - 2][0]) / self.time
        r_e_i = self.a_i ** (1 / 3) * self.b_i ** (2 / 3)
        r_e_j = self.a_j ** (1 / 3) * self.b_j ** (2 / 3)
        p_i = self.a_i/self.b_i
        p_j = self.a_j/self.b_j
        K_t_i = (np.sqrt(1-p_i**(-2)))/(p_i**(-2/3)*np.log(p_i*(1+np.sqrt(1-p_i**(-2)))))
        K_t_j = (np.sqrt(1-p_j**(-2)))/(p_j**(-2/3)*np.log(p_j*(1+np.sqrt(1-p_j**(-2)))))
        K_r_i = ((4*p_i**2)*(1-p_i**2))/(3*((((2*p_i**(2/3))*(2-p_i**(-2)))/K_t_i)-2))
        K_r_j = ((4*p_j**2)*(1-p_j**2))/(3*((((2*p_j**(2/3))*(2-p_j**(-2)))/K_t_j)-2))
        t_fr_i = 8*np.pi*self.calc_viscosity(timestep-1)*(r_e_i**3)*K_r_i*theta_dot_i
        t_fr_j = 8*np.pi*self.calc_viscosity(timestep-1)*(r_e_j**3)*K_r_j*theta_dot_j
        if math.isnan(t_fr_i) or math.isinf(t_fr_i):
            t_fr_i = 0
        if math.isnan(t_fr_j) or math.isinf(t_fr_j):
            t_fr_j = 0
        return t_fr_i, t_fr_j

    def calc_pm_r(self, timestep):
        if self.theta_i[timestep-1] > np.pi/2:
            loc_plus_i = (self.x_i[timestep-1] - self.a_i*np.cos(self.theta_i[timestep-1]), 
                          self.y_i[timestep-1] + self.a_i*np.sin(self.theta_i[timestep-1]))
            loc_minus_i = (self.x_i[timestep-1] + self.a_i*np.cos(self.theta_i[timestep-1]), 
                          self.y_i[timestep-1] - self.a_i*np.sin(self.theta_i[timestep-1]))
            r_i_plus = np.sqrt((loc_plus_i[0])**2 + (loc_plus_i[1])**2)
            r_i_minus = np.sqrt((loc_minus_i[0])**2 + (loc_minus_i[1])**2)
        else:
            loc_plus_i = (self.x_i[timestep-1] + self.a_i*np.cos(self.theta_i[timestep-1]), 
                          self.y_i[timestep-1] + self.a_i*np.sin(self.theta_i[timestep-1]))
            loc_minus_i = (self.x_i[timestep-1] + self.a_i*np.cos(self.theta_i[timestep-1]), 
                          self.y_i[timestep-1] + self.a_i*np.sin(self.theta_i[timestep-1]))
            r_i_plus = np.sqrt((loc_plus_i[0])**2 + (loc_plus_i[1])**2)
            r_i_minus = np.sqrt((loc_minus_i[0])**2 + (loc_minus_i[1])**2)

        if self.theta_j[timestep-1] > np.pi/2:
            loc_plus_j = (self.x_j[timestep-1] - self.a_j*np.cos(self.theta_j[timestep-1]), 
                          self.y_j[timestep-1] + self.a_j*np.sin(self.theta_j[timestep-1]))
            loc_minus_j = (self.x_j[timestep-1] + self.a_j*np.cos(self.theta_j[timestep-1]), 
                          self.y_j[timestep-1] - self.a_j*np.sin(self.theta_j[timestep-1]))
            r_j_plus = np.sqrt((loc_plus_j[0])**2 + (loc_plus_j[1])**2)
            r_j_minus = np.sqrt((loc_minus_j[0])**2 + (loc_minus_j[1])**2)
        else:
            loc_plus_j = (self.x_j[timestep-1] + self.a_j*np.cos(self.theta_j[timestep-1]),
                          self.y_j[timestep-1] + self.a_j*np.sin(self.theta_j[timestep-1]))
            loc_minus_j = (self.x_j[timestep-1] + self.a_j*np.cos(self.theta_j[timestep-1]), 
                          self.y_j[timestep-1] + self.a_j*np.sin(self.theta_j[timestep-1]))
            r_j_plus = np.sqrt((loc_plus_j[0])**2 + (loc_plus_j[1])**2)
            r_j_minus = np.sqrt((loc_minus_j[0])**2 + (loc_minus_j[1])**2)
        return (r_i_plus, r_i_minus, r_j_plus, r_j_minus), (loc_plus_i, loc_minus_i, loc_plus_j, loc_minus_j)


    def calc_coup(self, timestep):
        r_i, r_j = self.calculate_r(timestep=timestep-1)
        # NOTE: same timestep, self.calc_pm_r calculates prev step internally
        r_pm_ij, r_locs = self.calc_pm_r(timestep=timestep)
        r_i_plus, r_i_minus, r_j_plus, r_j_minus = r_pm_ij
        loc_plus_i, loc_minus_i, loc_plus_j, loc_minus_j = r_locs

        q_i = ((self.fluid_permittivity * self.V_i * self.E_mag) / (2 * self.a_i) * 
               (self.eq_permittivity_i - self.fluid_permittivity) / (self.fluid_permittivity + 
                (self.eq_permittivity_i - self.fluid_permittivity) * 2 * self.a_i))
        q_j = ((self.fluid_permittivity * self.V_j * self.E_mag) / (2 * self.a_j) *
               (self.eq_permittivity_j - self.fluid_permittivity) / (self.fluid_permittivity +
                (self.eq_permittivity_j - self.fluid_permittivity) * 2 * self.a_j))

        slope_i = np.tan(self.theta_i[timestep - 1])
        slope_j = np.tan(self.theta_j[timestep - 1])

        # CNT i
        f_plus_plus_i = (1 / 4 * np.pi * self.fluid_permittivity) * (abs(q_i) * abs(q_j)) / (
                abs(r_i_plus - r_j_plus) ** 2)
        slope_plus_plus = (loc_plus_i[1] - loc_plus_j[1]) / (loc_plus_i[0] - loc_plus_j[0])
        gamma_plus_plus_i = np.arctan((slope_i - slope_plus_plus) / (1 + slope_i * slope_plus_plus))

        f_plus_minus_i = (1 / 4 * np.pi * self.fluid_permittivity) * (abs(q_i) * abs(q_j)) / (
                abs(r_i_plus - r_j_minus) ** 2)
        slope_plus_minus = (loc_plus_i[1] - loc_minus_j[1]) / (loc_plus_i[0] - loc_minus_j[0])
        gamma_plus_minus_i = np.arctan((slope_i - slope_plus_minus) / (1 + slope_i * slope_plus_minus))

        f_minus_minus_i = (1 / 4 * np.pi * self.fluid_permittivity) * (abs(q_i) * abs(q_j)) / (
                abs(r_i_minus - r_j_minus) ** 2)
        slope_minus_minus = (loc_minus_i[1] - loc_minus_j[1]) / (loc_minus_i[0] - loc_minus_j[0])
        gamma_minus_minus_i = np.arctan((slope_i - slope_minus_minus) / (1 + slope_i * slope_minus_minus))

        f_minus_plus_i = (1 / 4 * np.pi * self.fluid_permittivity) * (abs(q_i) * abs(q_j)) / (
                abs(r_i_minus - r_j_plus) ** 2)
        slope_minus_plus = (loc_minus_i[1] - loc_plus_j[1]) / (loc_minus_i[0] - loc_plus_j[0])
        gamma_minus_plus_i = np.arctan((slope_i - slope_minus_plus) / (1 + slope_i * slope_minus_plus))

        # CNT j
        f_plus_plus_j = (1 / 4 * np.pi * self.fluid_permittivity) * (abs(q_j) * abs(q_i)) / (
                abs(r_j_plus - r_i_plus) ** 2)
        slope_plus_plus = (loc_plus_j[1] - loc_plus_i[1]) / (loc_plus_j[0] - loc_plus_i[0])
        gamma_plus_plus_j = np.arctan((slope_j - slope_plus_plus) / (1 + slope_j * slope_plus_plus))

        f_plus_minus_j = (1 / 4 * np.pi * self.fluid_permittivity) * (abs(q_j) * abs(q_i)) / (
                abs(r_j_plus - r_i_minus) ** 2)
        slope_plus_minus = (loc_plus_j[1] - loc_minus_i[1]) / (loc_plus_j[0] - loc_minus_i[0])
        gamma_plus_minus_j = np.arctan((slope_j - slope_plus_minus) / (1 + slope_j * slope_plus_minus))

        f_minus_minus_j = (1 / 4 * np.pi * self.fluid_permittivity) * (abs(q_j) * abs(q_i)) / (
                abs(r_j_minus - r_i_minus) ** 2)
        slope_minus_minus = (loc_minus_j[1] - loc_minus_i[1]) / (loc_minus_j[0] - loc_minus_i[0])
        gamma_minus_minus_j = np.arctan((slope_j - slope_minus_minus) / (1 + slope_j * slope_minus_minus))

        f_minus_plus_j = (1 / 4 * np.pi * self.fluid_permittivity) * (abs(q_j) * abs(q_i)) / (
                abs(r_j_minus - r_i_plus) ** 2)
        slope_minus_plus = (loc_minus_j[1] - loc_plus_i[1]) / (loc_minus_j[0] - loc_plus_i[0])
        gamma_minus_plus_j = np.arctan((slope_j - slope_minus_plus) / (1 + slope_j * slope_minus_plus))

        return (
            (f_plus_plus_i[0], gamma_plus_plus_i[0]),
            (f_plus_minus_i[0], gamma_plus_minus_i[0]),
            (f_minus_minus_i[0], gamma_minus_minus_i[0]),
            (f_minus_plus_i[0], gamma_minus_plus_i[0]),
            (f_plus_plus_j[0], gamma_plus_plus_j[0]),
            (f_plus_minus_j[0], gamma_plus_minus_j[0]),
            (f_minus_minus_j[0], gamma_minus_minus_j[0]),
            (f_minus_plus_j[0], gamma_minus_plus_j[0])
        )

    def calc_coup_t(self, timestep):
        # Note: same timestep, self.calc_coup does timestep-1 internally
        forces = self.calc_coup(timestep=timestep)
        t_coup_i = self.a_i*(forces[0][0]*np.sin(forces[0][1]) +
                             forces[1][0]*np.sin(forces[1][1]) +
                             forces[2][0]*np.sin(forces[2][1]) +
                             forces[3][0]*np.sin(forces[3][1])
        )
        t_coup_j = self.a_i*(forces[4][0]*np.sin(forces[4][1]) +
                             forces[5][0]*np.sin(forces[5][1]) +
                             forces[6][0]*np.sin(forces[6][1]) +
                             forces[7][0]*np.sin(forces[7][1])
        )
        if math.isnan(t_coup_i) or math.isinf(t_coup_i):
            t_coup_i = 0
        if math.isnan(t_coup_j) or math.isinf(t_coup_j):
            t_coup_j = 0
        return t_coup_i, t_coup_j

    def calc_coup_f(self, timestep):
        # Note: same timestep, self.calc_coup does timestep-1 internally
        forces = self.calc_coup(timestep=timestep)
        f_coup_x_i = (forces[0][0]*np.cos(self.theta_i[timestep-1]+forces[0][1]) +
                      forces[1][0]*np.cos(self.theta_i[timestep-1]+forces[1][1]) +
                      forces[2][0]*np.cos(self.theta_i[timestep-1]+forces[2][1]) +
                      forces[3][0]*np.cos(self.theta_i[timestep-1]+forces[3][1])
        )
        f_coup_y_i = (forces[0][0]*np.sin(self.theta_i[timestep-1]+forces[0][1]) +
                      forces[1][0]*np.sin(self.theta_i[timestep-1]+forces[1][1]) +
                      forces[2][0]*np.sin(self.theta_i[timestep-1]+forces[2][1]) +
                      forces[3][0]*np.sin(self.theta_i[timestep-1]+forces[3][1])
        )
        f_coup_x_j = (forces[4][0]*np.cos(self.theta_j[timestep-1]+forces[4][1]) +
                      forces[5][0]*np.cos(self.theta_j[timestep-1]+forces[5][1]) +
                      forces[6][0]*np.cos(self.theta_j[timestep-1]+forces[6][1]) +
                      forces[7][0]*np.cos(self.theta_j[timestep-1]+forces[7][1])
        )
        f_coup_y_j = (forces[4][0]*np.sin(self.theta_j[timestep-1]+forces[4][1]) +
                      forces[5][0]*np.sin(self.theta_j[timestep-1]+forces[5][1]) +
                      forces[6][0]*np.sin(self.theta_j[timestep-1]+forces[6][1]) +
                      forces[7][0]*np.sin(self.theta_j[timestep-1]+forces[7][1])
        )
        if math.isnan(f_coup_x_i) or math.isinf(f_coup_x_i):
            f_coup_x_i = 0
        if math.isnan(f_coup_y_i) or math.isinf(f_coup_y_i):
            f_coup_y_i = 0
        if math.isnan(f_coup_x_j) or math.isinf(f_coup_x_j):
            f_coup_x_j = 0
        if math.isnan(f_coup_y_i) or math.isinf(f_coup_x_j):
            f_coup_y_j = 0
        return f_coup_x_i, f_coup_y_i, f_coup_x_j, f_coup_y_j

    def calc_fr_f(self, timestep):
        x_dot_i = (self.x_i[timestep - 1][0] - self.x_i[timestep - 2][0]) / self.time
        y_dot_i = (self.y_i[timestep - 1][0] - self.y_i[timestep - 2][0]) / self.time
        x_dot_j = (self.x_j[timestep - 1][0] - self.x_j[timestep - 2][0]) / self.time
        y_dot_j = (self.y_j[timestep - 1][0] - self.y_j[timestep - 2][0]) / self.time
        r_e_i = self.a_i ** (1 / 3) * self.b_i ** (2 / 3)
        r_e_j = self.a_j ** (1 / 3) * self.b_j ** (2 / 3)
        p_i = self.a_i/self.b_i
        p_j = self.a_j/self.b_j
        K_t_i = (np.sqrt(1-p_i**(-2)))/(p_i**(-2/3)*np.log(p_i*(1+np.sqrt(1-p_i**(-2)))))
        K_t_j = (np.sqrt(1-p_j**(-2)))/(p_j**(-2/3)*np.log(p_j*(1+np.sqrt(1-p_j**(-2)))))
        f_fr_x_i = 8 * np.pi * self.calc_viscosity(timestep - 1) * r_e_i * K_t_i * x_dot_i
        f_fr_y_i = 8 * np.pi * self.calc_viscosity(timestep - 1) * r_e_i * K_t_i * y_dot_i
        f_fr_x_j = 8 * np.pi * self.calc_viscosity(timestep - 1) * r_e_j * K_t_j * x_dot_j
        f_fr_y_j = 8 * np.pi * self.calc_viscosity(timestep - 1) * r_e_j * K_t_j * y_dot_j
        if math.isnan(f_fr_x_i) or math.isinf(f_fr_x_i):
            f_fr_x_i = 0
        if math.isnan(f_fr_y_i) or math.isinf(f_fr_y_i):
            f_fr_y_i = 0
        if math.isnan(f_fr_x_j) or math.isinf(f_fr_x_j):
            f_fr_x_j = 0
        if math.isnan(f_fr_y_i) or math.isinf(f_fr_x_j):
            f_fr_y_j = 0
        return f_fr_x_i, f_fr_y_i, f_fr_x_j, f_fr_y_j

    def calc_rep_f(self, timestep, f_coups):
        r_i, r_j = self.calculate_r(timestep=timestep-1)
        r_pm_ij, r_locs = self.calc_pm_r(timestep=timestep)
        r_i_plus, r_i_minus, r_j_plus, r_j_minus = r_pm_ij
        Delta = 1*10**(-9)
        f_rep_x_i = 0
        f_rep_y_i = 0
        f_rep_x_j = 0
        f_rep_y_j = 0
        if np.sqrt((self.x_i[timestep - 1] - self.x_j[timestep - 1])**2 + (self.y_i[timestep-1] - self.y_j[timestep-1])**2) < Delta:
            f_rep_x_i = f_coups[0] * (np.exp(-100 * (((r_j_minus - r_i_plus) / Delta) - 1)) + np.exp(-100 * (((r_j_plus - r_i_minus) / Delta) - 1)))
            f_rep_y_i = f_coups[1] * (np.exp(-100 * (((r_j_minus - r_i_plus) / Delta) - 1)) + np.exp(-100 * (((r_j_plus - r_i_minus) / Delta) - 1)))
            f_rep_x_j = f_coups[2] * (np.exp(-100 * (((r_i_minus - r_j_plus) / Delta) - 1)) + np.exp(-100 * (((r_i_plus - r_j_minus) / Delta) - 1)))
            f_rep_y_j = f_coups[3] * (np.exp(-100 * (((r_i_minus - r_j_plus) / Delta) - 1)) + np.exp(-100 * (((r_i_plus - r_j_minus) / Delta) - 1)))
        return f_rep_x_i, f_rep_y_i, f_rep_x_j, f_rep_y_j

    def run(self, show=False):
        for timestep in tqdm(range(2, self.timesteps+2)):
            # Calculate new theta:
            t_dep_i, t_dep_j = self.calc_dep_t(timestep=timestep)
            t_coup_i, t_coup_j = (0, 0) #self.calc_coup_t(timestep=timestep)
            t_fr_i, t_fr_j = (0, 0) #self.calc_fr_t(timestep=timestep)
            theta_dot_dot_i = (t_dep_i + t_coup_i + t_fr_i) / self.I_i
            theta_dot_dot_j = (t_dep_j + t_coup_j + t_fr_j) / self.I_j
            theta_dot_i = (self.theta_i[timestep - 1][0] - self.theta_i[timestep - 2][0]) / self.time
            theta_dot_j = (self.theta_j[timestep - 1][0] - self.theta_j[timestep - 2][0]) / self.time
            self.theta_i[timestep] = theta_dot_i * self.time + 0.5 * theta_dot_dot_i * (self.time ** 2)
            self.theta_j[timestep] = theta_dot_j * self.time + 0.5 * theta_dot_dot_j * (self.time ** 2)

            # Calculate new x and y:
            f_fr_i_x, f_fr_i_y, f_fr_j_x, f_fr_j_y = (0, 0, 0, 0) #self.calc_fr_f(timestep=timestep)
            f_coup_i_x, f_coup_i_y, f_coup_j_x, f_coup_j_y = (0, 0, 0, 0) #self.calc_coup_f(timestep=timestep)
            f_rep_i_x, f_rep_i_y, f_rep_j_x, f_rep_j_y = (0, 0, 0, 0) #self.calc_rep_f(timestep=timestep, f_coups=(f_coup_i_x, f_coup_i_y, f_coup_j_x, f_coup_j_y))
            x_dot_dot_i = (f_fr_i_x + f_coup_i_x + f_rep_i_x) / self.m_i
            x_dot_dot_j = (f_fr_j_x + f_coup_j_x + f_rep_j_x) / self.m_j
            y_dot_dot_i = (f_rep_i_y + f_coup_i_y + f_rep_i_y) / self.m_i
            y_dot_dot_j = (f_rep_j_y + f_coup_j_y + f_rep_j_y) / self.m_j
            x_dot_i = (self.x_i[timestep - 1][0] - self.x_i[timestep - 2][0]) / self.time
            x_dot_j = (self.x_j[timestep - 1][0] - self.x_j[timestep - 2][0]) / self.time
            y_dot_i = (self.y_i[timestep - 1][0] - self.y_i[timestep - 2][0]) / self.time
            y_dot_j = (self.y_j[timestep - 1][0] - self.y_j[timestep - 2][0]) / self.time
            self.x_i[timestep] = x_dot_i * self.time + 0.5 * x_dot_dot_i * (self.time ** 2)
            self.x_j[timestep] = x_dot_j * self.time + 0.5 * x_dot_dot_j * (self.time ** 2)
            self.y_i[timestep] = y_dot_i * self.time + 0.5 * y_dot_dot_i * (self.time ** 2)
            self.y_j[timestep] = y_dot_j * self.time + 0.5 * y_dot_dot_j * (self.time ** 2)

            if math.isnan(self.theta_i[timestep]) or math.isinf(self.theta_i[timestep]):
                self.theta_i[timestep] = self.theta_i[timestep-1]
            if math.isnan(self.theta_j[timestep]) or math.isinf(self.theta_j[timestep]):
                self.theta_j[timestep] = self.theta_j[timestep-1]
            if math.isnan(self.x_i[timestep]) or math.isinf(self.x_i[timestep]):
                self.x_i[timestep] = self.x_i[timestep-1]
            if math.isnan(self.y_i[timestep]) or math.isinf(self.y_i[timestep]):
                self.y_i[timestep] = self.y_i[timestep-1]
            if math.isnan(self.x_j[timestep]) or math.isinf(self.x_j[timestep]):
                self.x_j[timestep] = self.x_j[timestep-1]
            if math.isnan(self.y_j[timestep]) or math.isinf(self.y_j[timestep]):
                self.y_j[timestep] = self.y_j[timestep - 1]

        if show:
            self.show_theta()
            self.show_x()
            self.show_y()

        
    
if __name__ == '__main__':
    sim = TwoCNT()
    sim.run(show=True)
    sim.show_cnts(timestep=0)
    sim.show_cnts(timestep=3)
    sim.show_cnts(timestep=30)
    sim.show_cnts(timestep=300)
    sim.show_cnts(timestep=3000)
    sim.show_cnts(timestep=30000)

