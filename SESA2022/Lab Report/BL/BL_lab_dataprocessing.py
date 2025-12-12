import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.pyplot import savefig
from scipy.optimize import curve_fit
from scipy.integrate import simpson


plt.rcParams["font.family"] = "Helvetica Neue"


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
    ReL = u_inf * 0.265 / nu

    # Find boundary layer thickness
    delta_99 = np.interp(0.99, u_ms / u_inf, y_mm)  # 99% thickness
    theta = simpson(u_ms / u_inf * (1 - u_ms / u_inf), y_mm)  # momentum thickness
    delta_star = simpson(1 - u_ms / u_inf, y_mm)
    # Plot velocity profile in dimensional units
    # plt.figure(1)
    # plt.plot(u_ms, y_mm, 'o')

    # print("BL thickness (delta_99): ", np.round(delta_99, 3), " mm")
    # print("BL thickness (theta): ", np.round(theta, 3), " mm")

    # Calculate shape factor
    H = delta_star / theta
    # print("Shape factor (H): ", np.round(H, 3))
    # if H > 2:
    # print("BL is Laminar")
    # else:
    # print("BL is Turbulent")

    # Plot velocity profile in non-dimensional units
    # plt.figure(2)
    # plt.plot(u_ms / u_inf, y_mm / delta_99, 'o')

    # plt.show()
    return theta, H, delta_99, u_inf, delta_star, ReL


def data_processing(theta, total_manometer, static_manometer, v_profile):
    Ue = np.sqrt(
        2
        * 0.001
        * (total_manometer - static_manometer)
        * np.cos(theta * np.pi / 180)
        * rho_water
        * g
        / rho
    )

    aa, ab = v_profile.columns  # naming columns
    y_mm = np.array(v_profile[aa])
    man_read_mm = np.array(v_profile[ab])

    u_ms = np.sqrt(
        2
        * 0.001
        * (man_read_mm - static_manometer)
        * np.cos(theta * np.pi / 180)
        * rho_water
        * g
        / rho
    )

    return Ue, u_ms, y_mm


# --- Environment Set up ---

Temperature_deg = 23  # Degrees C
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
delta_star_array = []

# --- Data processing ---

profiles = [
    {"file": "Bench_1_1_throttle.csv", "total": 182, "static": 86},
    {"file": "Bench_1_0.33_throttle.csv", "total": 191, "static": 173},
    {"file": "Bench_2_0.66_throttle.csv", "total": 86, "static": 45},
    {"file": "Bench_2_0.25_throttle.csv", "total": 42, "static": 37},
    {"file": "Bench_3_0.5_throttle.csv", "total": 220, "static": 136},
    {"file": "Bench_3_0.75_throttle.csv", "total": 222, "static": 100},
    {"file": "Bench_4_0.66_throttle.csv", "total": 220, "static": 128},
    {"file": "Bench_4_0.25_throttle.csv", "total": 134, "static": 112},
    {"file": "Bench_5_0.33_throttle.csv", "total": 140, "static": 122},
    {"file": "Bench_5_1_throttle.csv", "total": 234, "static": 130},
    {"file": "Bench_6_0.5_throttle.csv", "total": 188, "static": 146},
    {"file": "Bench_6_0.75_throttle.csv", "total": 224, "static": 154},
]

data_id = 0
for index, profile_data in enumerate(profiles):
    data_id += 1
    bench = index // 2 + 1
    profile = pd.read_csv(
        f"/Users/dfps16/Desktop/LabReportSESA2022/BL/{profile_data['file']}"
    )
    Ue, u_ms, y_mm = data_processing(
        theta, profile_data["total"], profile_data["static"], profile
    )
    if data_id == 1:
        print(Ue, u_ms)
    momt, shape_factor, delta99, u_inf, delta_star, Re = shape_boundary_layer(
        y_mm, u_ms
    )
    if data_id == 1:
        print(momt)
    CF = 2 * momt * 0.001 / 0.265  # Viscous Drag Coefficient

    Re_array.append(Re)
    Ue_array.append(Ue)
    u_ms_array.append(u_ms)
    y_mm_array.append(y_mm)
    CF_array.append(CF)
    Uinf_array.append(u_inf)
    delta_99_array.append(delta99)
    delta_star_array.append(delta_star)

    pr = True
    if pr:
        print(
            f"ID: {data_id}, Bench {bench}, w/ Ue: {np.round(Ue, 2)} m/s,"
            f" w/ Re: {np.round(Re, 0)}"
        )

        print(
            f"Momentum thickness: {np.round(momt, 3)} mm, "
            f"Shape factor: {np.round(shape_factor, 2)}, "
            f"CF: {np.round(CF, 4)}, "
            # f"Cf: {np.round(Cf, 4)}, "
            f"Uinf: {np.round(u_inf, 2)} m/s"
        )

    for index in range(len(u_ms_array)):
        u = u_ms_array[index]
        y = y_mm_array[index] / 1000  # Converting mm to m

        n_points = 2  # select number of points to analyse near the wall
        y_near_wall = y[1 : n_points + 1]
        u_near_wall = u[1 : n_points + 1]

        # Fiting a 2nd order polynomial (a*y^2 + b*y + c)
        # coeffs = np.polyfit(y_near_wall, u_near_wall, 2)

        # At y = 0 du/dy = b
        dudy_wall = (u_near_wall[0] - u_near_wall[1]) / (
            y_near_wall[0] - y_near_wall[1]
        )
        C_f = 2 * nu / Uinf_array[index] ** 2 * dudy_wall
        C_f_array.append(C_f)

# --- Re - CF Process ---
# Fit the curve
# Generate fitted curve for plotting
Re_array = np.array(Re_array)
CF_array = np.array(CF_array)

sort_idx = np.argsort(Re_array)
Re_array_sorted = Re_array[sort_idx]
CF_array = CF_array[sort_idx]
print(Re_array_sorted)
print(CF_array)
indices_to_exclude = [0, 1, 4]  # e.g., 1st, 6th, and 12th data points

# Create a mask for values NOT in the exclusion list
mask = np.ones(len(Re_array_sorted), dtype=bool)
mask[indices_to_exclude] = False

# Fit using only the filtered data
params, covariance = curve_fit(power_law, Re_array_sorted[mask], CF_array[mask])
a, b = params

Re_fit = np.linspace(0, max(Re_array_sorted), 100)
CF_fit = power_law(Re_fit, a, b)
CF_theory = power_law(Re_fit, 0.074, -0.2)

# Calculate C_f

# --- Plots ---
# Plotting Viscous Drag Coeff
fig1, ax1 = plt.subplots(1, 1, figsize=(10, 6), constrained_layout=True)
ax1.plot(
    Re_array_sorted[mask],
    CF_array[mask],
    "d",
    label="Experimental data (fitted)",
    color="black",
)
ax1.plot(
    Re_array_sorted[~mask],
    CF_array[~mask],
    "x",
    color="red",
    markersize=10,
    label="Excluded data",
)
ax1.plot(
    Re_fit,
    CF_fit,
    "--",
    color="black",
    label=f"Curve Fit: $C_F = {a:.2f} \\cdot Re^{{{b:.3f}}}$",
)
ax1.plot(
    Re_fit,
    CF_theory,
    "-",
    color="purple",
    label="Theory: $C_F = 0.074 \\cdot Re^{-0.2}$",
)
ax1.set_xlabel("Reynolds Number $Re$", fontsize=20)
ax1.set_ylabel("Viscous Drag Coefficient $C_F$", fontsize=20)
ax1.set_title("$C_F$ over $Re$", fontsize=22)
ax1.legend(fontsize=16)
# print(f"Power-law fit: CF = {a:.6e} * Re^({b:.4f})")
ax1.grid(True)
plt.savefig("ReCF.png", dpi=600)

# Sort C_f_array using the same indices as Re_array
C_f_array_sorted = np.array(C_f_array)[sort_idx]
# Fit using the same mask (excluding the same Reynolds numbers)
params_cf, covariance_cf = curve_fit(
    power_law, Re_array_sorted[mask], C_f_array_sorted[mask]
)
a_cf, b_cf = params_cf

Re_fit = np.linspace(200000, max(Re_array_sorted), 100)

# Generate fitted curve
Cf_fit = power_law(Re_fit, a_cf, b_cf)
Cf_theory = power_law(Re_fit, 0.059, -0.2)  # Blasius theory for turbulent BL

# Replace the existing fig2 plot with:
fig2, ax2 = plt.subplots(1, 1, figsize=(10, 6), constrained_layout=True)
ax2.plot(
    Re_array_sorted[mask],
    C_f_array_sorted[mask],
    "p",
    color="blue",
    label="Experimental data (fitted)",
)
ax2.plot(
    Re_fit,
    Cf_fit,
    "--",
    color="blue",
    label=f"Fit: $C_f = {a_cf:.2e} \\cdot Re^{{{b_cf:.3f}}}$",
)
ax2.plot(
    Re_fit,
    Cf_theory,
    "-",
    color="darkgreen",
    label="Theory: $C_f = 0.059 \\cdot Re^{-0.2}$",
)
ax2.set_xlabel("$Re$", fontsize=20)
ax2.set_ylabel("Local friction coefficient $C_f$", fontsize=20)
ax2.set_title("$C_f$ over $Re$", fontsize=22)
ax2.legend(fontsize=16)
ax2.grid(True)
plt.savefig("CfRe.png", dpi=600)

# Plotting velocity profiles
fig3, ax3 = plt.subplots(1, 1, figsize=(12, 10), constrained_layout=True)
# In dimensional coordinates
sorted_indices = np.argsort(Re_array)
for index in sorted_indices:
    bench = index // 2 + 1
    ax3.plot(
        u_ms_array[index],
        y_mm_array[index],
        "o",
        label=f"$Re$ {np.round(Re_array[index] / 10**3, 0)}k",
        linestyle="--",
    )

ax3.set_xlabel("Velocity (m/s)", fontsize=20)
ax3.set_ylabel("Distance from wall (mm)", fontsize=20)
ax3.set_title(
    "Boundary Layer Velocity Profiles of one side of a flat plate", fontsize=22
)
plt.legend(loc="best", fontsize=12)
ax3.grid(True)
plt.savefig("profiles.png", dpi=600)

# In non-dimensional coordinates
"""for index in range(len(u_ms_array)):
    bench = index // 2 + 1
    ax3[1].plot(u_ms_array[index]/Uinf_array[index],
             y_mm_array[index]/delta_star_array[index], 's',
             label=f'Ue {np.round(Ue_array[index],2 )} m/s')

ax3[1].set_xlabel("$u/U_{\infty}$")
ax3[1].set_ylabel('$y/\delta_{99}$')
ax3[1].set_title('Velocity Profiles (Non-Dimensional)')
ax3[1].grid(True)
ax3[1].grid(True)
ax3[1].legend()
# plt.legend()
ax3[1].grid(True)"""

plt.show()
