import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit

def power_law(Re, a, b):
    return a * Re**b

def shape_boundary_layer(y_mm, u_ms):
    """
    Function to find the 99 and momentum thickness of the BL,
    compute the shape factor, and ascertain the nature of the
    boundary layer. Takes a xlsx file with two columns of data.
    Column 1 is the distance along the airfoil in mm
    and column 2 the velocity along the airfoil in m/s.
    """

    u_inf = u_ms.max()  # getting freestream velocity

    # Find boundary layer thickness
    delta_99 = np.interp(0.99, u_ms / u_inf, y_mm)  # 99% thickness
    theta = np.trapezoid(u_ms / u_inf *
                         (1 - u_ms / u_inf), y_mm)  # momentum thickness

    # Plot velocity profile in dimensional units
    # plt.figure(1)
    # plt.plot(u_ms, y_mm, 'o')

    # print("BL thickness (delta_99): ", np.round(delta_99, 3), " mm")
    #print("BL thickness (theta): ", np.round(theta, 3), " mm")

    # Calculate shape factor
    H = delta_99 / theta
    #print("Shape factor (H): ", np.round(H, 3))
    #if H > 2:
        #print("BL is Laminar")
    #else:
        #print("BL is Turbulent")

    # Plot velocity profile in non-dimensional units
    # plt.figure(2)
    # plt.plot(u_ms / u_inf, y_mm / delta_99, 'o')

    # plt.show()
    return theta, H, delta_99, u_inf


def data_processing(theta, total_manometer,
                    static_manometer, v_profile):

    Ue = np.sqrt(
        2 * 0.001 * (total_manometer - static_manometer) * np.cos(
            theta * np.pi / 180) * rho_water * g / rho)
    ReL = Ue * xloc1 / nu

    aa, ab = v_profile.columns  # naming columns
    y_mm = np.array(v_profile[aa])
    man_read_mm = np.array(v_profile[ab])

    u_ms = np.sqrt(2 * 0.001 * (man_read_mm - static_manometer) * np.cos(
        theta * np.pi / 180) * rho_water * g / rho)


    return ReL, Ue, u_ms, y_mm


# --- Environment Set up ---

Temperature_deg = 22  # Degrees C
Pressure = 100100  # Pa

Temperature_kelvin = Temperature_deg + 273  # K
rho = Pressure / (287 * Temperature_kelvin)  # Kg/m^3

mu = 1.81 * 10**-5  # Pa*s dynamic viscosity
nu = 1.81 * 10**-5 / rho  # m^2/s - kinematic viscosity of air

rho_water = 1000  # kg/m^3
g = 9.81  # m/s^2

theta = 60  # deg - manometer angle

xloc1 = 0.265  # distance in m of the loc of measurement from the leading edge

# --- Data setup ---
Re_array = []
Ue_array = []
u_ms_array = []
y_mm_array = []
CF_array = []
Uinf_array = []
delta_99_array = []
C_f_array = []

# --- Data processing ---

profiles = [
    {'file': 'Bench_1_0.33_throttle.csv', 'total': 191, 'static': 173},
    {'file': 'Bench_1_1_throttle.csv', 'total': 182, 'static': 86},
    {'file': 'Bench_2_0.66_throttle.csv', 'total': 86, 'static': 45},
    {'file': 'Bench_2_0.25_throttle.csv', 'total': 42, 'static': 37},
    {'file': 'Bench_3_0.5_throttle.csv', 'total': 220, 'static': 136},
    {'file': 'Bench_3_0.75_throttle.csv', 'total': 222, 'static': 100},
    {'file': 'Bench_4_0.66_throttle.csv', 'total': 220, 'static': 128},
    {'file': 'Bench_4_0.25_throttle.csv', 'total': 134, 'static': 112},
    {'file': 'Bench_5_0.33_throttle.csv', 'total': 140, 'static': 122},
    {'file': 'Bench_5_1_throttle.csv', 'total': 234, 'static': 130},
    {'file': 'Bench_6_0.5_throttle.csv', 'total': 188, 'static': 146},
    {'file': 'Bench_6_0.75_throttle.csv', 'total': 224, 'static': 154},
]

data_id = 0
for index, profile_data in enumerate(profiles):
    data_id += 1
    bench = index // 2 + 1
    profile = pd.read_csv(
        f'/Users/dfps16/Desktop/LabReportSESA2022/BL/{profile_data["file"]}'
    )
    Re, Ue, u_ms, y_mm = data_processing(
        theta,
        profile_data["total"],
        profile_data["static"],
        profile
    )

    momt, shape_factor, delta99, u_inf = shape_boundary_layer(y_mm, u_ms)
    CF = 2 * momt / xloc1  # Viscous Drag Coefficient

    Re_array.append(Re)
    Ue_array.append(Ue)
    u_ms_array.append(u_ms)
    y_mm_array.append(y_mm)
    CF_array.append(CF)
    Uinf_array.append(u_inf)
    delta_99_array.append(delta99)

    print(
        f"ID: {data_id}, Bench {bench}, w/ Ue: {np.round(Ue, 2)} m/s"
    )

    print(f"Momentum thickness: {np.round(momt, 3)} mm, "
          f"Shape factor: {np.round(shape_factor, 2)}, "
          f"CF: {np.round(CF, 2)}, "
          f"Uinf: {np.round(u_inf, 2)} m/s")

for index in range(len(u_ms_array)):
    u = u_ms_array[index]
    y = y_mm_array[index] / 1000  # Converting mm to m

    n_points = 7  # select number of points to analyse near the wall
    y_near_wall = y[:n_points]
    u_near_wall = u[:n_points]

    # Fiting a 2nd order polynomial (a*y^2 + b*y + c)
    coeffs = np.polyfit(y_near_wall, u_near_wall, 2)

    # At y = 0 du/dy = b
    dudy_wall = coeffs[1]
    C_f = 2 * nu / Uinf_array[index] * dudy_wall
    C_f_array.append(C_f)

# --- Re - CF Process ---
# Fit the curve
params, covariance = curve_fit(power_law, Re_array, CF_array)
a, b = params

# Generate fitted curve for plotting
Re_fit = np.linspace(min(Re_array), max(Re_array), 100)
CF_fit = power_law(Re_fit, a, b)

# Calculate C_f

# --- Plots ---
# Plotting Viscous Drag Coeff
fig1, ax1 = plt.subplots(1, 1, figsize=(10, 6), dpi=300,
                         constrained_layout=True)
ax1.plot(Re_array, CF_array, 'o', label='Experimental data')
ax1.plot(Re_fit, CF_fit, '-', label=f'Fit: $C_F = {a:.2e} \\cdot Re^'
                                       f''f'{{b:.3f}}$')
ax1.set_xlabel('Reynolds Number')
ax1.set_ylabel('Viscous Drag Coefficient $C_F$')
ax1.legend()
print(f"Power-law fit: CF = {a:.6e} * Re^({b:.4f})")
ax1.grid(True)

# Plotting local skin friction coeff
fig2, ax2 = plt.subplots(1, 1, figsize=(10, 6), dpi=300,
                         constrained_layout=True)
ax2.plot(Re_array, C_f_array, 'x')
ax2.set_xlabel('Reynolds Number')
ax2.set_ylabel('Local friction coefficient $C_f$')
ax2.grid(True)

# Plotting velocity profiles
fig3, ax3 = plt.subplots(1, 2, figsize=(24, 10), dpi=300,
                         constrained_layout=True)
# In dimensional coordinates
for index in range(len(u_ms_array)):
    bench = index // 2 + 1
    ax3[0].plot(u_ms_array[index], y_mm_array[index], 'o',
             label=f'Ue {np.round(Ue_array[index],2 )} m/s')

ax3[0].set_xlabel('Velocity (m/s)')
ax3[0].set_ylabel('Distance from wall (mm)')
ax3[0].set_title('Velocity Profiles (Dimensional)')
# plt.legend()
ax3[0].grid(True)

# In non-dimensional coordinates
for index in range(len(u_ms_array)):
    bench = index // 2 + 1
    ax3[1].plot(u_ms_array[index]/Uinf_array[index],
             y_mm_array[index]/delta_99_array[index], 'x',
             label=f'Ue {np.round(Ue_array[index],2 )} m/s')

ax3[1].set_xlabel("$u/U_{\infty}$")
ax3[1].set_ylabel('$y/\delta_{99}$')
ax3[1].set_title('Velocity Profiles (Non-Dimensional)')
ax3[1].grid(True)
ax3[1].grid(True)
ax3[1].legend()
# plt.legend()
ax3[1].grid(True)

plt.show()

