import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

rho_air = 1.25
nu_air =1.5e-5
mu_air = rho_air*nu_air

def shape_boundary_layer(v_profile):

    aa, ab = v_profile.columns # naming columns
    y_mm = v_profile[aa]
    u_ms = v_profile[ab]
    U_inf = u_ms.max() # getting freestream velocity

    # Plot velocity profile in dimensional units
    plt.figure(1)
    plt.plot(u_ms, y_mm, 'o')

    # Find boundary layer thickness
    delta_99 = np.interp(0.99, u_ms/U_inf, y_mm) # 99% thickness
    theta = np.trapz(u_ms/U_inf * (1 - u_ms/U_inf), y_mm) # momentum thickness
    
    print("BL thickness (delta_99): ", np.round(delta_99, 2), " mm")
    print("BL thickness (theta): ", np.round(theta, 2), " mm")

    # Calculate shape factor
    H = delta_99 / theta
    print("Shape factor (H): ", np.round(H, 2))

    # Plot velocity profile in non-dimensional units
    plt.figure(2)
    plt.plot(u_ms/U_inf, y_mm/delta_99, 'o')

    plt.show()
    return delta_99
    
# Input velocity profile data
data = pd.read_excel('/Users/dfps16/Documents/GitHub/diegos-tools/SESA2022/Quiz 1/LBLdata2.xlsx')

delta_99 = shape_boundary_layer(data)