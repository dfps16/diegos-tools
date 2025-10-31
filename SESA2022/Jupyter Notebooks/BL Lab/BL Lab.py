import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd

# General Data
Temperature_deg = 22 #Degrees C
Pressure= 100100 #Pa

Temperature_kelvin = Temperature_deg + 273 #K
rho = Pressure / (287 * Temperature_kelvin) #Kg/m^3 (You need to enter the formula to determine density from pressure and temperature with gas constant)

xloc1 = 0.265  # distance in m of the location of measurement from the leading edge.

# --- Data collection ---
# Enter the micrometer reading when the probe is just touching the wall
# This is value in mm
wall_offset = 10.19  # mm

nu = 1.5e-5  # m^2/s - kinematic viscosity of air (not that an accurate value can be computed using the measured temperature)
rho_water = 1000  # kg/m^3
g = 9.81  # m/s^2
theta = 60  # Angle of inclination of the manometer with the vertical in degrees.

# Reynolds number 1
total_manometer_1 = 188  # Manometer reading in mm for total pressure
static_manometer_1 = 146  # Manometer reading in mm for static pressure
print(f"Total p reading: {total_manometer_1}, Static p reading: {static_manometer_1}"
      )

# Reynolds number 2
total_manometer_2 = 224  # Manometer reading in mm for total pressure
static_manometer_2 = 154  # Manometer reading in mm for static pressure

# Need to convert height to m and then convert m of water in to velocity.
# The below equation is just converting total and static pressure to velocity
# The cos(theta) is there to take the inclination of the manometer in to account - Part 1 Thermofluids!!

Ue1 = np.sqrt(2 * 0.001 * (total_manometer_1 - static_manometer_1) * np.cos(
    theta * np.pi / 180) * rho_water * g / rho)
Ue2 = np.sqrt(2 * 0.001 * (total_manometer_2 - static_manometer_2) * np.cos(
    theta * np.pi / 180) * rho_water * g / rho)

ReL1 = Ue1 * xloc1 / nu
ReL2 = Ue2 * xloc1 / nu

# --- Re Case 1 ---
#Probe should be ~10mm of wall (These are approximate distances)
Micrometerreading_1 = 20.25
Manometerreading_1 = 188

#Probe should be ~8mm of wall (These are approximate distances)
Micrometerreading_2 =18
Manometerreading_2 = 187

#Probe should be ~7mm of wall (These are approximate distances)
Micrometerreading_3 = 17
Manometerreading_3 = 186

#Probe should be ~6mm of wall (These are approximate distances)
Micrometerreading_4 = 16
Manometerreading_4 = 185

#Probe should be ~5mm of wall (These are approximate distances)
Micrometerreading_5 = 15
Manometerreading_5 = 185

#Probe should be ~4mm of wall (These are approximate distances)
Micrometerreading_6 = 14
Manometerreading_6 = 182

#Probe should be ~3mm of wall (These are approximate distances)
Micrometerreading_7 = 13
Manometerreading_7 = 178

#Probe should be ~2mm of wall (These are approximate distances)
Micrometerreading_8 = 12
Manometerreading_8 = 172

#Probe should be ~1.8mm of wall (These are approximate distances)
Micrometerreading_9 = 11.8
Manometerreading_9 = 170

#Probe should be ~1.4mm of wall (These are approximate distances)
Micrometerreading_10 = 11.4
Manometerreading_10 = 169

#Probe should be ~1mm of wall (These are approximate distances)
Micrometerreading_11 = 11
Manometerreading_11 = 165

#Probe should be ~0.6mm of wall (These are approximate distances)
Micrometerreading_12 = 10.6
Manometerreading_12 = 161

#Probe should be ~0.2mm of wall(These are approximate distances)
Micrometerreading_13 = 10.2
Manometerreading_13 = 156

# Enter the micrometer reading when the probe is just touching the wall
#This is value in mm
wall_offset = 9.69

# We are going to take the above manometer values, compile that into an array and convert the pressure readings to velocity.
# Note that the numbers are entered in to the array in reverse order so that the data appears from the wall to the freestream
#This is easier to plot and also write to file and carry out curve fits during analysis

micrometer = np.array([Micrometerreading_13,Micrometerreading_12,Micrometerreading_11,Micrometerreading_10,Micrometerreading_9,
                        Micrometerreading_8,Micrometerreading_7,Micrometerreading_6,Micrometerreading_5,Micrometerreading_4,
                        Micrometerreading_3,Micrometerreading_2,Micrometerreading_1])

manometer = np.array([Manometerreading_13,Manometerreading_12,Manometerreading_11,Manometerreading_10,Manometerreading_9,
                        Manometerreading_8,Manometerreading_7,Manometerreading_6,Manometerreading_5,Manometerreading_4,
                        Manometerreading_3,Manometerreading_2,Manometerreading_1])


#This subtracts the wall offset and the static manometer readings and converts the manometer to pressure
micrometer = micrometer - wall_offset

#Converting the pressure to velocity
U = np.sqrt(2*0.001*(manometer - static_manometer_1)*np.cos(theta*np.pi/180)*rho_water*g/rho)

#The below uses Pandas to write the data to a csv file that you can later on read for analysis
filename = 'Bench_6_0.5_throttle.csv'

plt.plot(U, micrometer, label='Re Case 1')
plt.title('Re Case 1')
plt.xlabel('U (m/s)')
plt.ylabel('x (mm)')
plt.show()

data = pd.DataFrame({"Wall position (mm)" : micrometer, "Manometer reading (mm)" : manometer})
data.to_csv(filename, index=False)

# --- Re Case 2 ---

#Probe should be ~10mm of wall (These are approximate distances)
Micrometerreading_1 = 20
Manometerreading_1 = 224

#Probe should be ~8mm of wall (These are approximate distances)
Micrometerreading_2 = 18
Manometerreading_2 = 224

#Probe should be ~7mm of wall (These are approximate distances)
Micrometerreading_3 = 17
Manometerreading_3 = 223

#Probe should be ~6mm of wall (These are approximate distances)
Micrometerreading_4 = 16
Manometerreading_4 = 223

#Probe should be ~5mm of wall (These are approximate distances)
Micrometerreading_5 = 15
Manometerreading_5 = 223

#Probe should be ~4mm of wall (These are approximate distances)
Micrometerreading_6 = 14
Manometerreading_6 = 220

#Probe should be ~3mm of wall (These are approximate distances)
Micrometerreading_7 = 13
Manometerreading_7 = 216

#Probe should be ~2mm of wall (These are approximate distances)
Micrometerreading_8 = 12
Manometerreading_8 = 208

#Probe should be ~1.8mm of wall (These are approximate distances)
Micrometerreading_9 = 11.8
Manometerreading_9 = 206

#Probe should be ~1.4mm of wall (These are approximate distances)
Micrometerreading_10 = 11.4
Manometerreading_10 = 200

#Probe should be ~1mm of wall (These are approximate distances)
Micrometerreading_11 = 11
Manometerreading_11 = 194

#Probe should be ~0.6mm of wall (These are approximate distances)
Micrometerreading_12 = 10.6
Manometerreading_12 = 192

#Probe should be ~0.2mm of wall(These are approximate distances)
Micrometerreading_13 = 10.2
Manometerreading_13 = 188

# Enter the micrometer reading when the probe is just touching the wall
#This is value in mm
wall_offset = 9.7

# We are going to take the above manometer values, compile that into an array and convert the pressure readings to velocity.
# Note that the numbers are entered in to the array in reverse order so that the data appears from the wall to the freestream
#This is easier to plot and also write to file and carry out curve fits during analysis

micrometer = np.array([Micrometerreading_13,Micrometerreading_12,Micrometerreading_11,Micrometerreading_10,Micrometerreading_9,
                        Micrometerreading_8,Micrometerreading_7,Micrometerreading_6,Micrometerreading_5,Micrometerreading_4,
                        Micrometerreading_3,Micrometerreading_2,Micrometerreading_1])

manometer = np.array([Manometerreading_13,Manometerreading_12,Manometerreading_11,Manometerreading_10,Manometerreading_9,
                        Manometerreading_8,Manometerreading_7,Manometerreading_6,Manometerreading_5,Manometerreading_4,
                        Manometerreading_3,Manometerreading_2,Manometerreading_1])


#This subtracts the wall offset and the static manometer readings and converts the manometer to pressure
micrometer = micrometer - wall_offset

#Converting the pressure to velocity
U = np.sqrt(2*0.001*(manometer - static_manometer_2)*np.cos(theta*np.pi/180)*rho_water*g/rho)

#The below uses Pandas to write the data to a csv file that you can later on read for analysis
filename = 'Bench_6_0.75_throttle.csv'
plt.plot(U, micrometer, label='Re Case 2', color='red')
plt.title('Re Case 2')
plt.title('Re Case 2')
plt.xlabel('U (m/s)')
plt.ylabel('x (mm)')
plt.show()

data = pd.DataFrame({"Wall position (mm)" : micrometer, "Manometer reading (mm)" : manometer})
data.to_csv(filename, index=False)