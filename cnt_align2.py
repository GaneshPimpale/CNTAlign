from random import uniform
import statistics
import scipy as sp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm



class CNTsim:
    def __init__(self, space_dim, timesteps, timesteps_per_second, num_cnts, d_range, l_range):
        # Simulation settings
        self.timesteps = timesteps
        self.timesteps_per_second = timesteps_per_second
        self.num_cnts = num_cnts
        self.x_range = (0, space_dim[0])
        self.y_range = (0, space_dim[1])
        self.theta_range = (0, np.pi)
        
        # CNT position, orientation and properties
        self.x = np.zeros((self.timesteps+3, self.num_cnts))
        self.y = np.zeros((self.timesteps+3, self.num_cnts))
        self.theta = np.zeros((self.timesteps+3, self.num_cnts))
        self.d = np.zeros((self.num_cnts, 1))
        self.a = np.zeros((self.num_cnts, 1))

        self.cnt_permittivity = None
        self.cnt_conductivity = None
        
        # Fluid and interphase properties
        self.viscosity = 30  #mPa*s
        self.fluid_permittivity = None
        self.fluid_conductivity = None
        self.interphase_permittivity = None
        self.interphase_conductivity = None
        self.interphase_thickness = None

        # Electric field properties
        self.E_mag = 6.6  #V/mm
        self.E_freq = 50  #Hz

        # Create CNT space
        for cnt_idx in range(self.num_cnts):
            self.x[0:3, cnt_idx] = uniform(self.x_range[0], self.x_range[1])
            self.y[0:3, cnt_idx] = uniform(self.y_range[0], self.y_range[1])
            self.theta[0:3, cnt_idx] = uniform(self.theta_range[0], self.theta_range[1])
            self.d[cnt_idx] = uniform(d_range[0], d_range[1])
            self.a[cnt_idx] = uniform(l_range[0], l_range[1]) / 2

    def set_material_properties(self):
        self.cnt_permittivity = 1*10**(5)
        self.cnt_conductivity = 1*10**(4)  #S/m
        self.fluid_permittivity = 11
        self.fluid_conductivity = 3*10**(-5)  #S/m
        self.interphase_permittivity = 1*10**(4)
        self.interphase_conductivity = 1*10**(-4)  #S/m
        self.interphase_thickness = 10*10**(-9)  #m

    def get_viscosity(self, timestep):
        #TODO: all of this
        return self.viscosity
    
    def dx_dt(self, x):
        dx = np.zeros(x.shape[0] -1)
        for idx in range(1, len(x)):
            dx[idx-1] = (x[idx]-x[idx-1])/(self.timesteps_per_second**(-1))
        return dx

    def calculate_eq_permittivity(self, cnt_idx):
        eq_permittivity = self.interphase_permittivity*((self.cnt_permittivity+(self.interphase_thickness /
                            (2*self.a[cnt_idx]))*(self.cnt_permittivity-self.interphase_permittivity)) /
                            (self.interphase_permittivity+(self.interphase_thickness/(2*self.a[cnt_idx])) *
                            (self.cnt_permittivity-self.interphase_permittivity)))
        return eq_permittivity

    def calculate_coulombic_force(self, timestep, cnt_i, cnt_j):
        r_i = np.sqrt(self.x[timestep, cnt_i] ** 2 + self.y[timestep, cnt_i] ** 2) 
        if self.theta[timestep, cnt_i] > np.pi/2:
            loc_plus_i = (self.x[timestep, cnt_i] - self.a[cnt_i]*np.cos(self.theta[timestep, cnt_i]), 
                          self.y[timestep, cnt_i] + self.a[cnt_i]*np.sin(self.theta[timestep, cnt_i]))
            loc_minus_i = (self.x[timestep, cnt_i] + self.a[cnt_i]*np.cos(self.theta[timestep, cnt_i]), 
                          self.y[timestep, cnt_i] - self.a[cnt_i]*np.sin(self.theta[timestep, cnt_i]))
            r_i_plus = np.sqrt((loc_plus_i[0])**2 + (loc_plus_i[1])**2)
            r_i_minus = np.sqrt((loc_minus_i[0])**2 + (loc_minus_i[1])**2)
        else:
            loc_plus_i = (self.x[timestep, cnt_i] + self.a[cnt_i]*np.cos(self.theta[timestep, cnt_i]), 
                          self.y[timestep, cnt_i] + self.a[cnt_i]*np.sin(self.theta[timestep, cnt_i]))
            loc_minus_i = (self.x[timestep, cnt_i] - self.a[cnt_i]*np.cos(self.theta[timestep, cnt_i]), 
                          self.y[timestep, cnt_i] - self.a[cnt_i]*np.sin(self.theta[timestep, cnt_i]))
            r_i_plus = np.sqrt((loc_plus_i[0])**2 + (loc_plus_i[1])**2)
            r_i_minus = np.sqrt((loc_minus_i[0])**2 + (loc_minus_i[1])**2)
        V_i = (4/3)*np.pi*self.a[cnt_i]*(self.d[cnt_i]/2)**2
        eq_permittivity_i = self.calculate_eq_permittivity(cnt_i)
        q_i = ((self.fluid_permittivity * V_i * self.E_mag) / (2 * self.a[cnt_i]) * 
               (eq_permittivity_i - self.fluid_permittivity) / (self.fluid_permittivity + 
                (eq_permittivity_i - self.fluid_permittivity) * 2 * self.a[cnt_i]))

        r_j = np.sqrt(self.x[timestep, cnt_j] ** 2 + self.y[timestep, cnt_j] ** 2)
        if self.theta[timestep, cnt_j] > np.pi/2:
            loc_plus_j = (self.x[timestep, cnt_j] - self.a[cnt_j] * np.cos(self.theta[timestep, cnt_j]),
                          self.y[timestep, cnt_j] + self.a[cnt_j] * np.sin(self.theta[timestep, cnt_j]))
            loc_minus_j = (self.x[timestep, cnt_j] + self.a[cnt_j] * np.cos(self.theta[timestep, cnt_j]),
                           self.y[timestep, cnt_j] - self.a[cnt_j] * np.sin(self.theta[timestep, cnt_j]))
            r_j_plus = np.sqrt((loc_plus_j[0]) ** 2 + (loc_plus_j[1]) ** 2)
            r_j_minus = np.sqrt((loc_minus_j[0]) ** 2 + (loc_minus_j[1]) ** 2)
        else:
            loc_plus_j = (self.x[timestep, cnt_j] + self.a[cnt_j] * np.cos(self.theta[timestep, cnt_j]),
                          self.y[timestep, cnt_j] + self.a[cnt_j] * np.sin(self.theta[timestep, cnt_j]))
            loc_minus_j = (self.x[timestep, cnt_j] - self.a[cnt_j] * np.cos(self.theta[timestep, cnt_j]),
                           self.y[timestep, cnt_j] - self.a[cnt_j] * np.sin(self.theta[timestep, cnt_j]))
            r_j_plus = np.sqrt((loc_plus_j[0]) ** 2 + (loc_plus_j[1]) ** 2)
            r_j_minus = np.sqrt((loc_minus_j[0]) ** 2 + (loc_minus_j[1]) ** 2)
        V_j = (4/3)*np.pi*self.a[cnt_j]*(self.d[cnt_j]/2)**2
        eq_permittivity_j = self.calculate_eq_permittivity(cnt_j)
        q_j = ((self.fluid_permittivity * V_j * self.E_mag) / (2 * self.a[cnt_j]) * 
               (eq_permittivity_j - self.fluid_permittivity) / (self.fluid_permittivity + 
                (eq_permittivity_j - self.fluid_permittivity) * 2 * self.a[cnt_j]))

        m_i = np.tan(self.theta[timestep, cnt_i])
        F_plus_plus = (1 / 4 * np.pi * self.fluid_permittivity) * (abs(q_i) * abs(q_j)) / (
                    abs(r_i_plus - r_j_plus) ** 2)
        m_plus_plus = (loc_plus_i[1] - loc_plus_j[1])/(loc_plus_i[0] - loc_plus_j[0])
        gamma_plus_plus = np.arctan((m_i - m_plus_plus)/(1 + m_i*m_plus_plus))
        
        F_plus_minus = (1 / 4 * np.pi * self.fluid_permittivity) * (abs(q_i) * abs(q_j)) / (
                    abs(r_i_plus - r_j_minus) ** 2)
        m_plus_minus = (loc_plus_i[1] - loc_minus_j[1])/(loc_plus_i[0] - loc_minus_j[0])
        gamma_plus_minus = np.arctan((m_i - m_plus_minus)/(1 + m_i*m_plus_minus))
        
        F_minus_minus = (1 / 4 * np.pi * self.fluid_permittivity) * (abs(q_i) * abs(q_j)) / (
                    abs(r_i_minus - r_j_minus) ** 2)
        m_minus_minus = (loc_minus_i[1] - loc_minus_j[1])/(loc_minus_i[0] - loc_minus_j[0])
        gamma_minus_minus = np.arctan((m_i - m_minus_minus)/(1 + m_i*m_minus_minus))
        
        F_minus_plus = (1 / 4 * np.pi * self.fluid_permittivity) * (abs(q_i) * abs(q_j)) / (
                    abs(r_i_minus - r_j_plus) ** 2)
        m_minus_plus = (loc_minus_i[1] - loc_plus_j[1])/(loc_minus_i[0] - loc_plus_j[0])
        gamma_minus_plus = np.arctan((m_i - m_minus_plus)/(1 + m_i*m_minus_plus))

        return (
            (F_plus_plus, gamma_plus_plus),
            (F_plus_minus, gamma_plus_minus),
            (F_minus_minus, gamma_minus_minus),
            (F_minus_plus, gamma_minus_plus)
        )

    def calculate_torques(self, cnt_idx, timestep):
        x_pos = self.x[timestep-1, cnt_idx]
        y_pos = self.y[timestep-1, cnt_idx]
        ori = self.theta[timestep-1, cnt_idx]
        d = self.d[cnt_idx]
        b = d / 2
        a = self.a[cnt_idx]

        
        # Calculate T_DEP
        eq_permittivity = self.calculate_eq_permittivity(cnt_idx)
        L = (np.log(2*a/b)-1)/((a/b)**2)
        alpha_star = (((eq_permittivity-self.fluid_permittivity)**2)/
                      ((self.fluid_permittivity+(eq_permittivity-self.fluid_permittivity)*L) *
                       (eq_permittivity+self.fluid_permittivity)))
        V = (4/3)*np.pi*a*b**2
        T_DEP_theta = (1/4)*V*self.fluid_permittivity*(self.E_mag**2)*alpha_star*np.sin(2*ori)

        # Calculate T-fr
        theta_dot = self.dx_dt(self.theta[:, cnt_idx])[-1]
        r_e = a**(1/3)*b**(2/3)
        p = a/b
        K_t = (np.sqrt(1-p**(-2)))/(p**(-2/3)*np.log(p*(1+np.sqrt(1-p**(-2)))))
        K_r = ((4*p**2)*(1-p**2))/(3*((((2*p**(2/3))*(2-p**(-2)))/K_t)-2))
        T_fr_theta = 8*np.pi*self.get_viscosity(timestep)*(r_e**3)*K_r*theta_dot

        # Calculate T-coup
        T_coup_theta = 0
        Torques_on_i = []
        for idx in range(self.num_cnts):
            if idx != cnt_idx:
                forces = (self.calculate_coulombic_force(
                    timestep=timestep-1,
                    cnt_i=cnt_idx,
                    cnt_j=idx
                ))
                torque = a*(forces[0][0]*np.sin(forces[0][1]) +
                            forces[1][0]*np.sin(forces[1][1]) +
                            forces[2][0]*np.sin(forces[2][1]) +
                            forces[3][0]*np.sin(forces[3][1]))
                Torques_on_i.append(torque)
        T_coup_theta = sum(Torques_on_i)

        # Return angular acceleration
        m = d*V
        I = m*(a**2 + d**2)/5
        theta_dot_dot = (-T_DEP_theta-T_fr_theta-T_coup_theta)/I
        return theta_dot_dot
        
    def calculate_forces(self, cnt_idx, timestep):
        x_pos = self.x[timestep-1, cnt_idx]
        y_pos = self.y[timestep-1, cnt_idx]
        ori = self.theta[timestep-1, cnt_idx]
        d = self.d[cnt_idx]
        b = d / 2
        a = self.a[cnt_idx]
        
        # Calculate F-fr
        F_fr_x = 0
        F_fr_y = 0

        # Calculate F-coup
        F_coup_x = 0
        F_coup_y = 0

        # Calculate F-rep
        F_rep_x = 0
        F_rep_y = 0

        # Return x and y acceleration
        V = (4/3)*np.pi*a*b**2
        m = d*V
        x_dot_dot = (-F_fr_x-F_coup_x-F_rep_x)/m
        y_dot_dot = (-F_fr_y-F_coup_y-F_rep_y)/m
        return x_dot_dot, y_dot_dot


    def run(self, show=False):
        for timestep in tqdm(range(3, self.timesteps+3)):
            for cnt_idx in range(self.num_cnts):
                # Get angular and translational accelerations on CNT
                theta_dot_dot = self.calculate_torques(cnt_idx, timestep)
                x_dot_dot, y_dot_dot = self.calculate_forces(cnt_idx, timestep)
                
                time = self.timesteps_per_second**(-1)
                delta_theta = self.dx_dt(self.theta[0:timestep-1, cnt_idx])[-1] * time + 0.5 * theta_dot_dot * (time ** 2)
                delta_x = self.dx_dt(self.x[0:timestep-1, cnt_idx])[-1] * time + 0.5 * x_dot_dot * (time ** 2)
                delta_y = self.dx_dt(self.y[0:timestep-1, cnt_idx])[-1] * time + 0.5 * y_dot_dot * (time ** 2)

                self.theta[timestep, cnt_idx] = self.theta[timestep - 1, cnt_idx] + delta_theta[0]
                self.x[timestep, cnt_idx] = self.x[timestep - 1, cnt_idx] + delta_x[0]
                self.y[timestep, cnt_idx] = self.y[timestep - 1, cnt_idx] + delta_y[0]


    def show_cnts(self, timestep):
        for cnt_idx in range(self.num_cnts):
            pos = (self.x[timestep, cnt_idx], self.y[timestep, cnt_idx])
            ori = self.theta[timestep, cnt_idx]
            d = self.d[cnt_idx]
            b = d/2
            a = self.a[cnt_idx]

            # Create start and 
            start_x = a*np.cos(ori) + pos[0]
            start_y = a*np.sin(ori) + pos[1]
            end_x = pos[0] - a*np.cos(ori)
            end_y = pos[1] - a*np.sin(ori)

            # plot line between points
            plt.plot([start_x, end_x], [start_y, end_y])
        plt.show()

    def show_cnt_stats(self, timestep, show=True):
        theta_mean = statistics.mean(self.theta[timestep]*(180/np.pi))
        sd = statistics.stdev(self.theta[timestep]*(180/np.pi))
        x = np.linspace(theta_mean-3*sd, theta_mean+3*sd, 100)
        ax = plt.plot(x, sp.stats.norm.pdf(x, theta_mean, sd))
        print('Timestep: ' + str(timestep) + ' mean: ' + str(theta_mean) + ' stdev: ' + str(sd))
        if show:
            plt.show()
        return ax

    def compare_cnt_stats(self, timesteps):
        for timestep in timesteps:
            self.show_cnt_stats(timestep, show=False)
        plt.show()

    def show_cnt_trends(self):
        pass

    def show_cnt_pos_ori(self, cnt_idx):
        thetas = self.theta[:, cnt_idx]
        timesteps = range(0, self.timesteps+3)
        plt.plot(timesteps, thetas)
        plt.show()

    def show_cnt_accels(self, cnt_idx):
        pass



if __name__ == '__main__':
    sim = CNTsim(
        space_dim=(1*10**(-4),1*10**(-4)),
        timesteps=300,
        timesteps_per_second=100,
        num_cnts=100,
        d_range=(1 * 10 ** (-9), 1 * 10 ** (-9)),
        l_range=(3 * 10 ** (-6), 3 * 10 ** (-6))
    )
    sim.set_material_properties()
    sim.show_cnts(timestep=0)
    sim.show_cnt_stats(timestep=0)
    sim.run(show=True)
    sim.show_cnts(timestep=300)
    sim.show_cnt_pos_ori(cnt_idx=1)

    

