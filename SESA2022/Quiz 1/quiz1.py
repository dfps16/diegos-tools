import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import fsolve

rho_air = 1.25
nu_air = 1.5e-5
mu_air = rho_air*nu_air


def shape_boundary_layer(v_profile):
    """
    Function to find the 99 and momentum thickness of the BL,
    compute the shape factor, and ascertain the nature of the
    boundary layer. Takes a xlsx file with two columns of data.
    Column 1 is the distance along the airfoil in mm
    and column 2 the velocity along the airfoil in m/s.
    """
    aa, ab = v_profile.columns  # naming columns
    y_mm = v_profile[aa]
    u_ms = v_profile[ab]
    u_inf = u_ms.max()  # getting freestream velocity

    # Find boundary layer thickness
    delta_99 = np.interp(0.99, u_ms/u_inf, y_mm)  # 99% thickness
    theta = np.trapezoid(u_ms/u_inf *
                         (1 - u_ms/u_inf), y_mm)  # momentum thickness

    # Plot velocity profile in dimensional units
    plt.figure(1)
    plt.plot(u_ms, y_mm, 'o')
    
    print("BL thickness (delta_99): ", np.round(delta_99, 3), " mm")
    print("BL thickness (theta): ", np.round(theta, 3), " mm")

    # Calculate shape factor
    H = delta_99 / theta
    print("Shape factor (H): ", np.round(H, 3))
    if H > 2:
        print("BL is Laminar")
    else:
        print("BL is Turbulent")

    # Plot velocity profile in non-dimensional units
    plt.figure(2)
    plt.plot(u_ms/u_inf, y_mm/delta_99, 'o')

    plt.show()
    return delta_99


def drag_per_unit_span(v_profile, c, re_t):
    """
    Function to find the drag per unit span, using a velocity profile,
    a chord value, and the Reynolds number at transition.
    """
    aa, ab = v_profile.columns  # naming columns
    y_mm = v_profile[aa]
    u_ms = v_profile[ab]
    u_inf = u_ms.max()  # getting freestream velocity

    print("Freestream velocity (U_inf): ", u_inf, " m/s")

    # Find boundary layer thickness
    theta = np.trapezoid(u_ms/u_inf *
                         (1 - u_ms/u_inf), y_mm)  # momentum thickness

    def virtual_origin(x_0):
        return 0.001 * theta - (c - x_0) * 0.037 * ((u_inf * (c-x_0)
                                                     / nu_air) ** (- 0.2))
    # Find x_0
    x_t = (re_t * nu_air) / u_inf
    x_0 = fsolve(virtual_origin, 0)[0]

    print("x_t: ", x_t, " m")
    print("x_0: ", x_0, " m")

    # Compute drag and multiply by 2 to account
    # for both sides of the airfoil
    d_prime = 2 * rho_air * u_inf**2 * (0.037 * (c - x_0) /
                                        (u_inf * (c - x_0) / nu_air) ** 0.2)

    return d_prime


def total_drag_lbl(b, c, nu, rho, U):
    """
    Function to find the total drag of an airfoil assuming a fully
    laminar boundary layer.
    """
    re = (U * c) / nu
    d_total = 0.5 * rho * U**2 * b * c * (2.656/np.sqrt(re))
    return d_total


def drag_coeff(c, u_inf, rho, nu):
    x_T = 0.2 * c

    Re_xt = (u_inf * c) / nu

    print("Re at transition: ", Re_xt)

    Re_cxt = (u_inf * (c - x_T)) / nu

    x_0 = x_T * (1 - 38.22 * Re_xt**(-3/8))  # Working out the virtual origin
    print("Virtual origin is at: ", np.round(x_0, 3), "m, or ",
          np.round(x_0/c, 2), "c")

    c_f = 0.074 * (1 - x_0 / c) * Re_cxt**(-1/5)
    c_d = np.round(2*c_f, 4)

    return c_d


def drag_one_side(c, x_0, rho, u_inf, re_c):
    nu = (u_inf * c) / re_c
    d_prime = rho * u_inf**2 * (0.037 * (c - x_0) /
                                (u_inf * (c - x_0) / nu) ** (1 / 5))
    return d_prime


# Input data and solve Qs

data = pd.read_excel(
    '/Users/dfps16/Documents/GitHub/diegos-tools/SESA2022/Quiz 1/Q5.xlsx'
)
