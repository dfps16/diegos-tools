import numpy as np
import matplotlib.pyplot as plt
from numpy import cos
from numpy import sin
import pandas as pd
from scipy.optimize import curve_fit
from pathlib import Path

plt.rcParams['font.family'] = 'Helvetica Neue'


def drag_polar(CL, CD0, k):
    """Drag polar equation: CD = CD0 + k*CL^2"""
    return CD0 + k * CL**2


def read_env(filepath):
    env = pd.read_csv(filepath, nrows=5, header=None, usecols=[0, 1])
    aa, ab = env.columns
    values = np.array(env[ab])
    u_inf = values[0]
    # p_dyn = values[1]
    T_K = values[2] + 273  # Amb temp in K
    p_a = values[3] * 100  # Atm pressure in Pa
    rho = p_a / (287 * T_K)
    return u_inf, rho


def read_data(filepath):
    data = pd.read_csv(filepath, skiprows=6)
    cols = data.columns.tolist()
    readings = []
    for col in cols:
        array = np.array(data[col])
        readings.append(array)

    # Assigning readings to arrays
    lift_reading = readings[0]
    side_reading = readings[1]
    drag_reading = readings[2]
    roll_reading = readings[3]
    pitch_reading = readings[4]
    # yaw_nm = readings[5]
    yaw_deg = readings[6]

    # Work out mean readings
    lift_mean = np.mean(lift_reading)
    side_mean = np.mean(side_reading)
    drag_mean = np.mean(drag_reading)
    yaw_mean = np.mean(yaw_deg)
    return side_mean, drag_mean, yaw_mean


def find_csvs(folder_path):
    folder = Path(folder_path)
    csv_files = list(folder.glob('*.csv'))
    return csv_files


def read_xflr_fw(filepath):
    intro = pd.read_csv(filepath, nrows=2, skiprows=3, usecols=[1], header=None)
    analysis_name = str(intro.iloc[0, 0])
    data = pd.read_csv(filepath, skiprows=5)

    if 'LLT' in analysis_name:
        target_dict = xflr_300k_LLT if '11.0' in analysis_name else xflr_600k_LLT
    elif 'VLM1' in analysis_name:
        target_dict = xflr_300k_VLM1 if '11.0' in analysis_name else xflr_600k_VLM1
    elif 'VLM2' in analysis_name:
        target_dict = xflr_300k_VLM2 if '11.0' in analysis_name else xflr_600k_VLM1
    elif 'Panel' in analysis_name:
        target_dict = xflr_300k_Panel if '11.0' in analysis_name else xflr_600k_Panel

    # Extracting the data and populating the targeted dictionary
    target_dict['Cl'] = np.array(data[' CL'].values)
    target_dict['Cd'] = np.array(data[' CD'].values)
    target_dict['aoa'] = np.array(data['alpha'].values)

    sort_idx = np.argsort(target_dict['aoa'])
    target_dict['aoa'] = np.array(data['alpha'].values)[sort_idx]
    target_dict['Cl'] = np.array(data[' CL'].values)[sort_idx]
    target_dict['Cd'] = np.array(data[' CD'].values)[sort_idx]
    target_dict['Max L/D'] = np.max(target_dict['Cl'] / target_dict['Cd'])


    # Work out lift slope & zero lift AoA
    coefficients = np.polyfit(target_dict['aoa'] * np.pi / 180,
                              target_dict['Cl'], 1)
    target_dict['Cl slope'] = coefficients[0]
    target_dict['Zero lift angle'] = coefficients[1] * 180 / np.pi

    return target_dict


def find_oswald(k):
    delta = k * np.pi * ar - 1
    e = 1 / (1 + delta)
    return e


def create_summary(target_dict, Re):
    summary = (f'Summary for: {Re} Re\n'
               f'Lift slope (1/rad): {np.round(target_dict['Cl slope'], 2)}\n'
               f'Zero lift angle (rad): '
               f'{np.round(target_dict['Zero lift angle'], 5)}\n'
               f'Cd0: {np.round(target_dict['Cd0'], 4)}\n'
               f'k: {np.round(target_dict['k'], 4)}\n'
               f'e: {np.round(target_dict['e'], 4)}\n'
               f'Max L/D: {np.round(target_dict['Max L/D'], 4)}\n'
               )
    return summary


def compare(target_dict, xflr_dict, Re, model):
    comparison = (f'Comparing exp. data to XFLR5 w/ {model} model.\n'
                  f'For: {Re}\n'
                  f'Exp: Lift slope (1/rad):'
                  f' {np.round(target_dict['Cl slope'], 2)}\n'
                  f'XFLR: Lift slope (1/rad):'
                  f' {np.round(xflr_dict['Cl slope'], 2)}\n'
                  f'Exp: Zero lift angle (rad): '
                  f'{np.round(target_dict['Zero lift angle'], 5)}\n'
                  f'XFLR: Zero lift angle (rad): '
                  f'{np.round(xflr_dict['Zero lift angle'] * np.pi / 180, 5)}\n'
                  f'Exp: Cd0: {np.round(target_dict['Cd0'], 4)}\n'
                  f'XFLR: Cd0: {np.round(xflr_dict['Cd0'], 4)}\n'
                  f'Exp: k: {np.round(target_dict['k'], 4)}\n'
                  f'XFLR: k: {np.round(xflr_dict['k'], 4)}\n'
                  f'Exp: e: {np.round(target_dict['e'], 4)}\n'
                  f'XFLR: e: {np.round(xflr_dict['e'], 4)}\n'
                  f'Exp: Max L/D: {np.round(target_dict['Max L/D'], 4)}\n'
                  f'XFLR: Max L/D: {np.round(xflr_dict['Max L/D'], 4)}\n'
                  )
    print(comparison)
    for val in target_dict:
        if val in ['Cl', 'Cd', 'aoa']:  # Skip these keys
            continue
        abs_difference = np.abs(target_dict[val] - xflr_dict[val])
        prc_difference = 100 * (np.abs(target_dict[val] - xflr_dict[val]) /
                           xflr_dict[val])
        print(f'Absolute difference for {val}: {abs_difference}')
        print(f'% difference for {val}: {np.round(prc_difference, 2)}%')

    return comparison
# Air properties
mu = 1.81 * 10**-5  # Pa*s dynamic viscosity


# --- Wing Profile ---
b = 1.11
c = 0.4

bf = 2 * b  # Full wingspan
s = bf * c  # Total wing area
ar = bf ** 2 / s  # Aspect ratio

# --- Processing data ---
data_300k = {'Cl': [], 'Cd': [], 'aoa': [],
             'Cl slope': 0, 'Zero lift angle': 0,
             'Cd0': 0, 'k': 0, 'e': 0, 'Max L/D': 0}
data_600k = {'Cl': [], 'Cd': [], 'aoa': [],
             'Cl slope': 0, 'Zero lift angle': 0,
             'Cd0': 0, 'k': 0, 'e': 0, 'Max L/D': 0}
xflr_300k_LLT = {'Cl': [], 'Cd': [], 'aoa': [],
                 'Cl slope': 0, 'Zero lift angle': 0,
                 'Cd0': 0, 'k': 0, 'e': 0, 'Max L/D': 0}
xflr_300k_VLM1 = {'Cl': [], 'Cd': [], 'aoa': [],
                 'Cl slope': 0, 'Zero lift angle': 0,
                 'Cd0': 0, 'k': 0, 'e': 0, 'Max L/D': 0}
xflr_300k_VLM2 = {'Cl': [], 'Cd': [], 'aoa': [],
                 'Cl slope': 0, 'Zero lift angle': 0,
                 'Cd0': 0, 'k': 0, 'e': 0, 'Max L/D': 0}
xflr_300k_Panel = {'Cl': [], 'Cd': [], 'aoa': [],
                 'Cl slope': 0, 'Zero lift angle': 0,
                 'Cd0': 0, 'k': 0, 'e': 0, 'Max L/D': 0}
xflr_600k_LLT = {'Cl': [], 'Cd': [], 'aoa': [],
                 'Cl slope': 0, 'Zero lift angle': 0,
                 'Cd0': 0, 'k': 0, 'e': 0, 'Max L/D': 0}
xflr_600k_VLM1 = {'Cl': [], 'Cd': [], 'aoa': [],
                 'Cl slope': 0, 'Zero lift angle': 0,
                 'Cd0': 0, 'k': 0, 'e': 0, 'Max L/D': 0}
xflr_600k_VLM2 = {'Cl': [], 'Cd': [], 'aoa': [],
                 'Cl slope': 0, 'Zero lift angle': 0,
                 'Cd0': 0, 'k': 0, 'e': 0, 'Max L/D': 0}
xflr_600k_Panel = {'Cl': [], 'Cd': [], 'aoa': [],
                 'Cl slope': 0, 'Zero lift angle': 0,
                 'Cd0': 0, 'k': 0, 'e': 0, 'Max L/D': 0}

files = find_csvs('/Users/dfps16/Desktop/LabReportSESA2022/GP16')

for file in files:
    u_inf, rho = read_env(file)
    # print(rho)
    # print(u_inf)
    nu = mu / rho  # m^2/s - kinematic viscosity of air
    Re = u_inf * c / nu
    F_n, F_t, alpha = read_data(file)
    L = F_n * cos(alpha * np.pi / 180) - F_t * sin(alpha * np.pi / 180)
    D = F_n * sin(alpha * np.pi / 180) + F_t * cos(alpha * np.pi / 180)
    Cl = 2 * L / (0.5 * rho * u_inf ** 2 * s)
    Cd = 2 * D / (0.5 * rho * u_inf ** 2 * s)

    target = data_300k if Re < 350000 else data_600k
    target['Cl'].append(Cl)
    target['Cd'].append(Cd)
    target['aoa'].append(alpha)

for data in [data_300k, data_600k]:
    for key in data:
        data[key] = np.array(data[key])
    sort_idx = np.argsort(data['aoa'])
    data['aoa'] = data['aoa'][sort_idx]
    data['Cl'] = data['Cl'][sort_idx]
    data['Cd'] = data['Cd'][sort_idx]
    data['Max L/D'] = np.max(data['Cl'] / data['Cd'])

# --- Extracting xflr results ---
xflr_results = find_csvs(
    '/Users/dfps16/Desktop/LabReportSESA2022/XFLR_FW_Results')

for file in xflr_results:
    res = read_xflr_fw(file)

# --- Creating Line fits by working out lift slope and 0 lift angle ---
for data in [data_300k, data_600k]:
    # Find indices of the two largest AoA (stalled)
    stalled = np.argsort(data['aoa'][:-3])

    # Mask the indices
    mask = np.ones(len(data['aoa']), dtype=bool)
    mask[stalled] = False

    # Apply mask to AoA and Cl
    angles = - data['aoa'][mask] * np.pi / 180  # Converting to radians
    # print(angles)
    cl = - data['Cl'][mask]
    # print(cl)

    # Fit linear model to find the zero lift angle and lift slope
    coeffs = np.polyfit(angles, cl, 1)
    data['Cl slope'] = coeffs[0]
    data['Zero lift angle'] = coeffs[1]
    # print(data['Cl slope'])
    # print(data['Zero lift angle'])

# Create best fit lines for plot
aoa_range_300k = np.linspace(min(-data_300k['aoa']), max(-data_300k['aoa']), 100)
aoa_range_600k = np.linspace(min(-data_600k['aoa']), max(-data_600k['aoa']), 100)

# Convert slopes to degrees and create fit line
cl_fit_300k = (data_300k['Cl slope'] * np.pi / 180 * aoa_range_300k +
               data_300k['Zero lift angle'] * np.pi / 180)
cl_fit_600k = (data_600k['Cl slope'] * np.pi / 180 * (aoa_range_600k) +
               data_600k['Zero lift angle'] * np.pi / 180)

# Create best fit lines for drag polar (quadratic fit) and work out Cd0 and k
for data in [data_300k, data_600k]:
    # Find indices of the two largest AoA (stalled)
    stalled = np.argsort(data['aoa'][-2:])

    # Mask the indices
    mask = np.ones(len(data['aoa']), dtype=bool)
    mask[stalled] = False

    # Apply mask to Cd and Cl
    cd = data['Cd'][mask]
    # print(angles)
    cl = - data['Cl'][mask]
    # print(cl)
    consts, dummy = curve_fit(drag_polar, cl, cd)
    data['Cd0'] = consts[0]
    data['k'] = consts[1]
    data['e'] = find_oswald(data['k'])


for data in [xflr_300k_LLT, xflr_600k_LLT]:
    cd = data['Cd']
    # print(angles)
    cl = data['Cl']
    # print(cl)
    consts, dummy = curve_fit(drag_polar, cl, cd)
    data['Cd0'] = consts[0]
    data['k'] = consts[1]
    data['e'] = find_oswald(data['k'])

# ---- Plots ----

# Lift vs AoA plot
fig1, ax1 = plt.subplots(1, 1, figsize=(10, 6),
                         constrained_layout=True)

# Plotting experimental data
'''ax1.plot(-data_300k['aoa'], -data_300k['Cl'],
            marker='x',
            color='g',
            linestyle='None',
            label='Experimental data - 300k Re'
            )'''
ax1.plot(-data_600k['aoa'], -data_600k['Cl'],
            marker='x',
            color='darkblue',
            linestyle='None',
            label='Experimental data - 600k Re'
            )

# Plotting best fit lines
'''ax1.plot(aoa_range_300k, cl_fit_300k,
         color='g', linestyle='--', alpha=0.7, linewidth=2,
         label=f'300k Re fit (slope: {data_300k['Cl slope']:.2f})')'''

ax1.plot(aoa_range_600k, cl_fit_600k,
         color='darkblue', linestyle='--', alpha=0.7, linewidth=2,
         label=f'Experimental data fit (slope={data_600k['Cl slope']:.2f})')

# Plotting XFLR results
'''ax1.plot(xflr_300k_LLT['aoa'], xflr_300k_LLT['Cl'],
         linestyle='-',
         color='r',
         label=f'LLT - 300k Re, slope: '
               f'{np.round(xflr_300k_LLT['Cl slope'], 2)}'
         )'''
ax1.plot(xflr_600k_LLT['aoa'], xflr_600k_LLT['Cl'],
         linestyle='-',
         color='red',
        alpha=0.8,
         label=f'XFLR LLT Analysis, (slope='
               f'{np.round(xflr_600k_LLT['Cl slope'], 2)})'
         )

ax1.set_xlabel('Angle of Attack (degrees)', fontsize=20)
ax1.set_ylabel('$C_L$', fontsize=20)
ax1.set_title(f'Lift Slope of a NACA0020 3D finite wing at 600k $Re$',
              fontsize=22)
ax1.legend(fontsize=16)
ax1.grid(True)
plt.savefig('LiftSlope3D.png', dpi=600)

# Cl vs Cd plot
fig2, ax2 = plt.subplots(1, 1, figsize=(10, 6),
                         constrained_layout=True)

'''ax2.plot(data_300k['Cd'], -data_300k['Cl'],
            marker='o',
            color='purple',
            label='300k Re',
            linestyle='None',
            )'''
ax2.plot(data_600k['Cd'], -data_600k['Cl'],
            marker='o',
            color='brown',
            label='Experimental data',
            linestyle='None',
            )
ax2.plot(xflr_600k_LLT['Cd'], xflr_600k_LLT['Cl'],
            # marker='o',
            color='darkgreen',
            label='XFLR5 - LLT ($C_{D0}$ = '
                  f'{np.round(xflr_600k_LLT['Cd0'], 4)}, '
                  f'k = {np.round(xflr_600k_LLT['k'], 4)})',
            linestyle='-',
            )
'''ax2.plot(xflr_300k_LLT['Cd'], xflr_300k_LLT['Cl'],
            # marker='o',
            color='g',
            label='LLT - 300k Re',
            linestyle='-',
            )'''

'''ax2.plot(xflr_300k_LLT['Cd'], xflr_300k_LLT['Cl'],
            # marker='o',
            color='g',
            label='LLT - 300k Re',
            linestyle='-',
            )'''
cl_range = np.linspace(min(-data_600k['Cl']), max(-data_600k['Cl']), num=100)

ax2.plot(drag_polar(cl_range, data_600k['Cd0'], data_600k['k']),
         cl_range,
            # marker='o',
            color='brown',
            label='Curve fit ($C_{D0}$ = '
                  f'{np.round(data_600k['Cd0'], 4)}, '
                  f'k = {np.round(data_600k['k'], 4)})',
            linestyle='--',
            )

# Recreate the mask used for the drag-polar fit on data_600k
stalled = np.argsort(data_600k['aoa'][-2:])
mask = np.ones(len(data_600k['aoa']), dtype=bool)
mask[stalled] = False

# Excluded (stalled) points
excluded_cd = data_600k['Cd'][~mask]
excluded_cl = -data_600k['Cl'][~mask]
# Add these to the existing Cl vs Cd plot
ax2.scatter(excluded_cd, excluded_cl,
            marker='x',
            color='black',
            s=80,
            label='Excluded (stalled)',
            zorder=5)

ax2.set_xlabel('$C_D$', fontsize=16)
ax2.set_ylabel('$C_L$', fontsize=16)
ax2.set_title(f'$C_L$ over $C_D$ for a NACA0020 3D finite wing at 600k $Re$',
              fontsize=22)
ax2.legend(fontsize=14)
ax2.grid(True)
plt.savefig('clcd.png', dpi=600)


'''fig3, ax3 = plt.subplots(1, 1, figsize=(10, 6),
                         constrained_layout=True)
ax3.plot(data_600k['aoa'], - data_600k['Cl'] / data_600k['Cd'])
'''
plt.show()

print(create_summary(data_600k, '600k'))

pff = compare(data_600k, xflr_600k_LLT, '600 000', 'LLT')

