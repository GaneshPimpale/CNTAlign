import numpy as np
import scipy
from matplotlib import pyplot, cm
from mpl_toolkits.mplot3d import Axes3D


A = 2.033*10**4
E_a = 5.5*10**4
n = 1.43
H_r = 6.2*10**5
V_r = 2.7*10**(-4)
R = 8.314
rho_r = 901
rho_c = 1020
C_pc = 2.4
k = 0.22

T_i = 25

n = 1
alpha = 1

n_x = 50
n_y = 50
n_it = 500
x = np.linspace(0, 0.15, n_x)
y = np.linspace(0, 0.15, n_y)
X, Y = np.meshgrid(x, y)

T = np.ones((n_y, n_x))*T_i  #Temperature 
D_c = np.zeros((n_y, n_x))  #Degree of cure
V = np.zeros((n_y, n_x))  #Velocity
P = np.zeros((n_y, n_x))  #Presure


def partial_T(T):
    rhs = k*np.gradient(np.gradient(T)) + rho_r*V_r*H_r*(A*np.exp(-(E_a/(R*T)))*(alpha**(2-n)*(1-alpha)**(n)))
    return rhs/(rho_c*C_pc)

for step in range(n_it):
    pass