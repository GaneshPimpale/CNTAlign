import numpy as np
from random import uniform
import matplotlib.pyplot as plt
from tqdm import tqdm

class CNTspace:
    def __init__(self, cnt_space_dim, density, d_cnt_range, L_cnt_range):
        # Simulation settings
        self.timesteps = 300
        self.timesteps_per_second = 100

        # CNT generation settings
        self.cnt_density = density
        self.num_cnts = self.cnt_density
        self.x_range = (0, cnt_space_dim[0])
        self.y_range = (0, cnt_space_dim[1])
        self.angle_range = (0, 180)  #degrees
        self.d_cnt_range = d_cnt_range  #m
        self.L_cnt_range = L_cnt_range  #m

        # CNT position, orientation, and properties
        self.x = np.zeros((self.num_cnts, self.timesteps+1))
        self.y = np.zeros((self.num_cnts, self.timesteps+1))
        self.theta = np.zeros((self.num_cnts, self.timesteps+1))
        self.d = np.zeros((self.num_cnts, 1))
        self.a = np.zeros((self.num_cnts, 1))
        self.cnt_permittivity = None
        self.cnt_conductivity = None

        # Fluid and interphase properties
        self.viscocity = 30  #mPa*s
        self.fluid_permittivity = None
        self.fluid_conductivity = None
        self.interphase_permittivity = None
        self.interphase_conductivity = None
        self.interphase_thickness = None

        # Electric field properties
        self.E_mag = None  #V/mm
        self.E_freq = None  #Hz

    def create_cnt(self):
        x_pos = uniform(self.x_range[0], self.x_range[1])
        y_pos = uniform(self.y_range[0], self.y_range[1])
        theta = uniform(self.angle_range[0], self.angle_range[1])
        d = uniform(self.d_cnt_range[0], self.d_cnt_range[1])
        a = uniform(self.L_cnt_range[0], self.L_cnt_range[1])/2        
        return [(x_pos, y_pos), theta, d, a]

    def create_cnt_space(self):
        for idx in range(self.num_cnts):
            cnt = self.create_cnt()
            self.x[idx][0] = cnt[0][0]
            self.y[idx][0] = cnt[0][1]
            self.theta[idx][0] = cnt[1]
            self.d[idx][0] = cnt[2]
            self.a[idx][0] = cnt[3]

    def set_E_field(self):
        self.E_mag = 6.6  #V/mm
        self.E_freq = 50  #Hz

    def set_material_properties(self):
        self.cnt_permittivity = 1*10**(5)
        self.cnt_conductivity = 1*10**(4)  #S/m
        self.fluid_permittivity = 11
        self.fluid_conductivity = 3*10**(-5)  #S/m
        self.interphase_permittivity = 1*10**(4)
        self.interphase_conductivity = 1*10**(-4)  #S/m
        self.interphase_thickness = 10*10**(-9)  #m

    def show_cnts(self, timestep):
        for cnt_idx in range(self.num_cnts):
            pos = (self.x[cnt_idx][timestep], self.y[cnt_idx][timestep])
            ori = self.theta[cnt_idx][timestep]
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

    def set_viscocity(self):
        pass

    def calculate_numerical_derivative(self, X):
        dx = np.zeros(X.shape[0]-1)
        for idx in range(1, len(X)):
            dx[idx-1] = (X[idx]-X[idx-1])/(self.timesteps_per_second**(-1))
        return dx

    """
    def calculate_torques_on_cnt(self, cnt_idx):
        pos = (self.x[cnt_idx][-1], self.y[cnt_idx][-1])
        ori = self.theta[cnt_idx][-1]
        d = self.d[cnt_idx]
        b = d/2
        a = self.a[cnt_idx]

        # Calculate T_DEP
        eq_permittivity = self.interphase_permittivity*((self.cnt_permittivity+(self.interphase_thickness/(2*a))*
            (self.cnt_permittivity-self.interphase_permittivity))/
            (self.interphase_permittivity+(self.interphase_thickness/(2*a))*(self.cnt_permittivity-self.interphase_permittivity))
        )
        L = (np.log(2*a/b)-1)/((a/b)**2)
        alpha_star = ((eq_permittivity-self.fluid_permittivity)**2)/((self.fluid_permittivity+(eq_permittivity-self.fluid_permittivity)*L)*(eq_permittivity+self.fluid_permittivity))
        V = (4/3)*np.pi*a*b**2
        T_DEP_theta = (1/4)*V*self.fluid_permittivity*(self.E_mag**2)*alpha_star*np.sin(2*ori)

        # Calculate T-fr
        theta_dot = self.calculate_numerical_derivative(self.theta[cnt_idx])[-1]
        r_e = a**(1/3)*b**(2/3)
        p = a/b
        K_t = (np.sqrt(1-p**(-2)))/(p**(-2/3)*np.log(p*(1+np.sqrt(1-p**(-2)))))
        K_r = ((4*p**2)*(1-p**2))/(3*((((2*p**(2/3))*(2-p**(-2)))/K_t)-2))
        T_fr_theta = 8*np.pi*self.viscocity*(r_e**3)*K_r*theta_dot

        # Calculate T-coup
        T_coup_theta = 0


        # Return angular acceleration
        m = d*V
        I = m*(a**2 + d**2)/5
        theta_dot_dot = (-T_DEP_theta-T_fr_theta-T_coup_theta)/I
        return theta_dot_dot
    """

    def calculate_forces_on_cnt(self, cnt_idx):
        pos = (self.x[cnt_idx][-1], self.y[cnt_idx][-1])
        ori = self.theta[cnt_idx][-1]
        d = self.d[cnt_idx]
        b = d/2
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


    def run(self, show=True):

        for timestep in tqdm(range(1, self.timesteps+1)):
            for idx in range(self.num_cnts):
                # Get angular and translational accelerations on CNT
                theta_dot_dot = self.calculate_torques_on_cnt(idx)
                x_dot_dot, y_dot_dot = self.calculate_forces_on_cnt(idx)

                # Calculate deltas in position and orientation
                time = self.timesteps_per_second**(-1)
                delta_theta = self.calculate_numerical_derivative(self.theta[idx])[-1]*time + 0.5*(theta_dot_dot)*(time**2)
                delta_x = self.calculate_numerical_derivative(self.x[idx])[-1]*time + 0.5*(x_dot_dot)*(time**2)
                delta_y = self.calculate_numerical_derivative(self.y[idx])[-1]*time + 0.5*(y_dot_dot)*(time**2)

                # Set new position and orientation
                theta_before = self.theta[idx, timestep-1]
                x_before = self.x[idx, timestep-1]
                y_before = self.y[idx, timestep-1]
                self.theta[idx, timestep] = theta_before + delta_theta[0]
                self.x[idx, timestep] = x_before + delta_x[0]
                self.y[idx, timestep] = y_before + delta_y[0]
            print(self.x[:, timestep])
        
        if show:
            cnt_idx = 0
            plot_time = range(0,self.timesteps+1)
            plt.plot(plot_time, self.theta[cnt_idx])
            plt.show()


if __name__ == '__main__':
    sim = CNTspace(
        cnt_space_dim=(1*10**(-4),1*10**(-4)), 
        density=100, 
        d_cnt_range = (1*10**(-9), 1*10**(-9)), 
        L_cnt_range=(1*10**(-6), 1*10**(-6))
    )
    sim.create_cnt_space()
    sim.show_cnts(timestep=0)
    sim.set_E_field()
    sim.set_material_properties()
    sim.run()
    sim.show_cnts(timestep=200)

