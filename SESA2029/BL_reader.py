import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import fsolve


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
    delta_99 = np.interp(0.99, u_ms / u_inf, y_mm)  # 99% thickness
    theta = np.trapezoid(u_ms / u_inf *
                         (1 - u_ms / u_inf), y_mm)  # momentum thickness
    delta_star = np.trapezoid(1-u_ms/u_inf,y_mm)

    # Plot velocity profile in dimensional units
    plt.figure(1)
    plt.plot(u_ms, y_mm, 'o')

    print("BL thickness (delta_99): ", np.round(delta_99, 3), " mm")
    print("BL thickness (theta): ", np.round(theta, 3), " mm")
    print("BL thickness (delta_star): ", np.round(delta_star, 3), " mm")

    # Calculate shape factor
    H = delta_star / theta
    print("Shape factor (H): ", np.round(H, 3))
    if H > 2:
        print("BL is Laminar")
    else:
        print("BL is Turbulent")
        # Turbulent BL metrics
        L = 1
        c_f = 2 * 0.001 * theta / L
        u_tau = u_inf * np.sqrt(c_f/2)
        print("Skin-Friction Coefficient: ", np.round(c_f, 3))
        print("U_tau is: ", np.round(u_tau, 3), "m/s")

    # Plot velocity profile in non-dimensional units
    plt.figure(2)
    plt.plot(u_ms / u_inf, y_mm / delta_99, 'o')

    plt.show()
    return delta_99, theta, H

data = pd.read_excel(
    '/Users/dfps16/Documents/GitHub/diegos-tools/SESA2029/BLLamAdi.xlsx'
)
delta_99, theta, H = shape_boundary_layer(data)
