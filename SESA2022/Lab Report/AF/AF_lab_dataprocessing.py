import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from numpy import cos, sin
from scipy.integrate import simpson

plt.rcParams["font.family"] = "Helvetica Neue"


def Aerofoil(x):
    """This function takes the x coordinate and returns the positive y
    coordinate of the aerofoil"""
    t = 0.20
    y = (5 * t) * (
        (0.2969 * x**0.5)
        - (0.1260 * x)
        - (0.3516 * x**2)
        + (0.2843 * x**3)
        - (0.1015 * x**4)
    )
    return y


def find_vel(total_manometer, static_manometer):
    U = np.sqrt(
        2
        * 0.001
        * (total_manometer - static_manometer)
        * np.cos(theta * np.pi / 180)
        * rho_water
        * g
        / rho
    )
    return U


def find_cn_ct(x_c_u, cp_u, x_c_l, cp_l):
    ct = simpson((cp_u * upper_slope), x_c_u) - simpson((cp_l * lower_slope), x_c_l)

    cn = simpson((cp_l), x_c_l) - simpson((cp_u), x_c_u)

    return cn, ct


def gradient(x):
    dydx = (
        5
        * thickness
        * (
            0.2969 * 0.5 * x ** (-0.5)
            - 0.1260
            - 0.3516 * 2 * x
            + 0.2843 * 3 * x**2
            - 0.1015 * 4 * x**3
        )
    )
    return dydx


def find_cp(manometer_reading, total_manometer, static_manometer):
    cp = (manometer_reading - static_manometer) / (total_manometer - static_manometer)
    return cp


def process_data(data):
    aa, ab = data.columns
    # var = np.array(aa)
    reading = np.array(data[ab])
    p_total = reading[0]
    p_static = reading[1]
    p_lower = reading[2:8]
    p_upper = reading[8:]

    vel_inf = find_vel(p_total, p_static)

    Re = vel_inf * chord / nu

    vel_lower = find_vel(p_total - p_static, p_lower - p_static)
    vel_upper = find_vel(p_total - p_static, p_upper - p_static)

    cp_lower = find_cp(p_lower, p_total, p_static)
    cp_upper = find_cp(p_upper, p_total, p_static)

    return vel_lower, vel_upper, cp_lower, cp_upper, Re


def xflr_cp(data):
    cols = data.columns.tolist()
    cp_xflr = []
    for col in cols:
        if col == "X":
            xc = np.array(data[col])
        else:
            cp = np.array(data[col])
            cp_xflr.append(cp)
    return xc, cp_xflr


def xflr_cl(data):
    aa, ab = data.columns.tolist()
    alpha = np.array(data[aa])
    cl = np.array(data[ab])
    return alpha, cl


# Print graph?
gr_airfoil = False
gr_cl = True
gr_cp = True

# --- Env Data ---
Temperature_deg = 22.5  # degrees C
Pressure = 100100  # Pa

Temperature_kelvin = Temperature_deg + 273  # K

rho = Pressure / (287 * Temperature_kelvin)  # Density of air in kg/m^3

mu = 1.81 * 10**-5  # Pa*s dynamic viscosity
nu = 1.81 * 10**-5 / rho  # m^2/s - kinematic viscosity of air

rho_water = 1000  # kg/m^3
g = 9.81  # m/s^2

theta = 0  # deg - manometer angle

thickness = 0.2  # thickness in terms of chord length
chord = 0.063

# This is set up so if the aerofoil is turned clockwise the upper surface will
# have the tapping positons below otherwise
# just switch upper and lower surface around
x = [1, 2, 4.5, 7.5, 11, 14.5, 20, 26, 32.1, 38, 44, 50]
xc_upper = []
xc_lower = []

for i in range(12):
    if i % 2 == 0:
        xc_lower.append(x[i] / 63)
    if i % 2 == 1:
        xc_upper.append(x[i] / 63)

# Airfoil profile graph
if gr_airfoil:
    fig0, ax0 = plt.subplots(1, 1, figsize=[8, 6])
    ax0.plot(
        [Aerofoil(i) for i in xc_upper], xc_upper, "ro", label="Upper Surface Tappings"
    )
    ax0.plot(
        [Aerofoil(i) for i in np.linspace(0, 1, 100)], np.linspace(0, 1, 100), "g-"
    )
    ax0.plot(
        [-Aerofoil(i) for i in xc_lower], xc_lower, "bo", label="Lower Surface Tappings"
    )
    ax0.plot(
        [-Aerofoil(i) for i in np.linspace(0, 1, 100)], np.linspace(0, 1, 100), "g-"
    )
    ax0.invert_yaxis()
    ax0.legend(loc="center left", bbox_to_anchor=(1.02, 0.5))
    # ax0.tight_layout()
    ax0.grid()
    ax0.axis("scaled")
    plt.savefig("airfoil.png", dpi=300)
    plt.show()


xloc1 = 0.265  # distance in m of the loc of measurement from the leading edge

# --- Data setup ---
angles = []
cl_values = []

# --- Data processing ---

files = [
    {"file": "-8_deg1.csv", "angle": -8},
    {"file": "-8_deg2.csv", "angle": -8},
    {"file": "-4_deg1.csv", "angle": -4},
    {"file": "-4_deg2.csv", "angle": -4},
    {"file": "-2_deg1.csv", "angle": -2},
    {"file": "-2_deg2.csv", "angle": -2},
    {"file": "0_deg1.csv", "angle": 0},
    {"file": "0_deg2.csv", "angle": 0},
    {"file": "2_deg1.csv", "angle": 2},
    {"file": "2_deg2.csv", "angle": 2},
    {"file": "4_deg1.csv", "angle": 4},
    {"file": "4_deg2.csv", "angle": 4},
    {"file": "8_deg1.csv", "angle": 8},
    {"file": "8_deg2.csv", "angle": 8},
    {"file": "12_deg1.csv", "angle": 12},
    {"file": "12_deg2.csv", "angle": 12},
    {"file": "14_deg1.csv", "angle": 14},
    {"file": "14_deg2.csv", "angle": 14},
]


# -- Reading and processing xflr output
# For Cp
cp_data = pd.read_csv("/Users/dfps16/Desktop/LabReportSESA2022/Cp_Graph.csv")
angle_to_col = {-8: 0, -4: 1, -2: 2, 0: 3, 2: 4, 4: 5, 8: 6, 12: 7, 14: 8}
xc_xflr, cp_xflr = xflr_cp(cp_data)

# For Cl
cl_data = pd.read_csv("/Users/dfps16/Desktop/LabReportSESA2022/xflr_cl.csv")
alpha_xflr, cl_xflr = xflr_cl(cl_data)
alpha_xflr = alpha_xflr * np.pi / 180  # Converting AoA to rad

data_id = 0
unique_angles = sorted(set([d["angle"] for d in files]))
n_angles = len(unique_angles)

selected_angles = [-4, 2, 12]  # ordered whitelist of AoA to plot

if gr_cp:
    fig1, axes = plt.subplots(
        1, len(selected_angles), figsize=(18, 6), constrained_layout=True
    )
    axes = np.atleast_1d(axes).flatten()
    angle_to_subplot = {angle: i for i, angle in enumerate(selected_angles)}
    plotted_angles = set()

for index, data in enumerate(files):
    data_id += 1
    profile = pd.read_csv(f"/Users/dfps16/Desktop/LabReportSESA2022/AF/{data['file']}")
    u_l, u_u, cp_l, cp_u, Re = process_data(profile)

    # Plotting cp data - only plot first occurrence of each angle
    # only plot these angles

    # inside the loop, change the plotting condition to also check membership
    if (
        gr_cp
        and data["angle"] not in plotted_angles
        and data["angle"] in angle_to_subplot
    ):
        subplot_idx = angle_to_subplot[data["angle"]]
        ax = axes[subplot_idx]

        ax.plot(
            xc_lower,
            -cp_l,
            marker="d",
            color="grey",
            linestyle="--",
            label="Lower Surface",
        )
        ax.plot(
            xc_upper,
            -cp_u,
            marker="o",
            color="brown",
            linestyle="--",
            label="Upper Surface",
        )
        ax.plot(
            xc_xflr,
            -cp_xflr[angle_to_col[data["angle"]]],
            color="black",
            linestyle="-",
            label="XFOIL",
            alpha=0.75,
        )
        ax.set_xlabel("x/c", fontsize=20)
        ax.set_ylabel("-$C_p$", fontsize=20)
        ax.set_title(f"{data['angle']}° AoA", fontsize=22)
        ax.legend(fontsize=15)
        ax.grid(True)

        plotted_angles.add(data["angle"])

    upper_slope = gradient(np.array(xc_upper))
    lower_slope = -gradient(np.array(xc_lower))

    c_n, c_t = find_cn_ct(xc_upper, cp_u, xc_lower, cp_l)

    cl = (c_n * cos(data["angle"] * np.pi / 180)) - c_t * sin(
        data["angle"] * np.pi / 180
    )
    # print(data["angle"], cl, Re)
    angles.append(data["angle"])
    cl_values.append(cl)
plt.savefig("CpDistrib.png", dpi=600)
angles = np.array(angles) * np.pi / 180
cl_values = np.array(cl_values)

coeffs = np.polyfit(angles, cl_values, 1)
lift_slope = coeffs[0]
zero_lift_angle = coeffs[1]

cxflr = np.polyfit(alpha_xflr, cl_xflr, 1)
ls_xflr = cxflr[0]
zla_xflr = cxflr[1]

if gr_cl:
    fig2, ax2 = plt.subplots(figsize=(10, 6), constrained_layout=True)
    ax2.scatter(
        angles,
        cl_values,
        label="Experimental Data",
        color="orange",
    )
    ax2.plot(
        angles,
        np.polyval(coeffs, angles),
        color="orange",
        linestyle="--",
        label=f"Linear Fit (slope={lift_slope:.3f})",
    )

    ax2.plot(
        angles,
        np.polyval([2 * np.pi, 0], angles),
        color="black",
        linestyle="-",
        label="Thin airfoil theory (slope=2$\pi$)",
    )
    ax2.plot(
        alpha_xflr,
        cl_xflr,
        color="purple",
        linestyle="-",
        label="XFOIL",
    )
    ax2.plot(
        alpha_xflr,
        np.polyval(cxflr, alpha_xflr),
        color="purple",
        linestyle="--",
        label=f"XFOIL Best Fit (slope={ls_xflr:.3f})",
    )

    ax2.set_xlabel("Angle of Attack (rad)", fontsize=20)
    ax2.set_ylabel("$C_l$", fontsize=20)
    ax2.set_title("Lift Slope of a NACA0020 2D airfoil at 100k $Re$", fontsize=22)
    ax2.grid(True)
    # plt.tight_layout()
    plt.legend(fontsize=16)
    plt.savefig("LiftSlope2D.png", dpi=600)

plt.show()
