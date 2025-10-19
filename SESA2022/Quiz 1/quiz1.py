import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

rho_air = 1.25
nu_air =1.5e-5
mu_air = rho_air*nu_air

def ascertain_boundary_layer(v_profile):



    aa, ab = v_profile.columns # naming columns
    y_mm = v_profile[aa]
    u_ms = v_profile[ab]
    U_inf = u_ms.max() # getting freestream velocity

    # Plot velocity profile in dimensional units
    plt.figure(1)
    plt.plot(u_ms, y_mm, 'o')

    # Find boundary layer thickness
    delta_99 = np.interp(0.99, u_ms/U_inf, y_mm)

    print("BL thickness (delta_99): ", np.round(delta_99, 2), " mm")

    # Plot velocity profile in non-dimensional units
    plt.figure(2)
    plt.plot(u_ms/U_inf, y_mm/delta_99, 'o')

    plt.show()
    return delta_99

# Input velocity profile data
data = pd.read_excel('LBLdata2.xlsx')

delta_99 = ascertain_boundary_layer(data)