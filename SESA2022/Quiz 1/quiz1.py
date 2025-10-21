import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import fsolve

rho_air = 1.25
nu_air =1.5e-5
mu_air = rho_air*nu_air

def shape_boundary_layer(v_profile):

    aa, ab = v_profile.columns # naming columns
    y_mm = v_profile[aa]
    u_ms = v_profile[ab]
    U_inf = u_ms.max() # getting freestream velocity

     # Find boundary layer thickness
    delta_99 = np.interp(0.99, u_ms/U_inf, y_mm) # 99% thickness
    theta = np.trapz(u_ms/U_inf * (1 - u_ms/U_inf), y_mm) # momentum thickness

    # Plot velocity profile in dimensional units
    plt.figure(1)
    plt.plot(u_ms, y_mm, 'o')
    
    print("BL thickness (delta_99): ", np.round(delta_99, 3), " mm")
    print("BL thickness (theta): ", np.round(theta, 3), " mm")

    # Calculate shape factor
    H = delta_99 / theta
    print("Shape factor (H): ", np.round(H, 3))

    # Plot velocity profile in non-dimensional units
    plt.figure(2)
    plt.plot(u_ms/U_inf, y_mm/delta_99, 'o')

    plt.show()
    return delta_99

def drag_per_unit_span(v_profile, c, Re_t):


    aa, ab = v_profile.columns # naming columns
    y_mm = v_profile[aa]
    u_ms = v_profile[ab]
    U_inf = u_ms.max() # getting freestream velocity

    print("Freestream velocity (U_inf): ", U_inf, " m/s")

     # Find boundary layer thickness
    theta = np.trapezoid(u_ms/U_inf * (1 - u_ms/U_inf), y_mm) # momentum thickness

    def virtual_origin(x_0):
        return 0.001*theta - (c-x_0)*0.037*((U_inf*(c-x_0)/nu_air)**(-0.2))
    # Find x_0
    x_t = ( Re_t * nu_air ) / U_inf
    x_0 = fsolve(virtual_origin, 0)[0]

    print("x_t: ", x_t, " m")
    print("x_0: ", x_0, " m")

    D = 2 * rho_air * U_inf**2 * (0.037 * (c - x_0) / (U_inf*(c - x_0)/nu_air)**0.2)

    return D

def total_drag_LBL(b, c, nu, rho, U):
    Re = (U * c) / nu
    D_total = 0.5 * rho * U**2 * b * c * (2.656/np.sqrt(Re))
    return D_total

def drag_coeff(c, U_inf, rho, nu):
    x_T = 0.2 * c

    Re_xt = (U_inf * c) / nu 

    print("Re at transition: ", Re_xt)

    Re_cxt = (U_inf * (c - x_T)) / nu

    x_0 = x_T * (1 - 38.22 * Re_xt**(-3/8)) # Working out the virtual origin
    print("Virtual origin is at: ", np.round(x_0, 3), "m, or ", np.round(x_0/c, 2), "c")

    C_F = 0.074 * (1 - x_0 / c) * Re_cxt**(-1/5)
    C_D = np.round(2*C_F, 4)

    return C_D

def drag_oneside(c, x_0, rho, U_inf, Re_c):
    nu = (U_inf * c) / Re_c
    D_prime = rho * U_inf**2 * (0.037 * (c - x_0) / (U_inf * (c- x_0) / nu)**(1/5))
    return D_prime
# Input data and solve Qs

v_profile = pd.read_excel('/Users/dfps16/Documents/GitHub/diegos-tools/SESA2022/Quiz 1/Q5.xlsx')
c = 0
U_inf = 0
rho = 0
Re_c = 0
x_0 = 0

D_prime = drag_oneside(c, x_0, rho, U_inf, Re_c)
print("D is: ", D_prime)