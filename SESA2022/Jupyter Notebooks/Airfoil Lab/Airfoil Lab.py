import numpy as np
import matplotlib.pyplot as plt
from numpy import cos
from numpy import sin
import os
import pandas as pd

def Aerofoil(x):
    '''This function takes the x coordinate and returns the positive y coordinate of the aerofoil'''
    t=0.20
    y = (5*t)*((0.2969*x**0.5)-(0.1260*x)-(0.3516*x**2)+(0.2843*x**3)-(0.1015*x**4))
    return(y)

# --- General Info ---
Temperature_deg = 22.5 #degrees C
Pressure= 100100 #Pa

Temperature_kelvin = Temperature_deg + 273 #K

rho = Pressure / (287 * Temperature_kelvin) #Density of air in kg/m^3

# This is set up so if the aerofoil is turned clockwise the upper surface with have the tapping positons below otherwise
# just switch upper and lower surface around
x = [1, 2, 4.5, 7.5, 11, 14.5, 20, 26, 32.1, 38, 44, 50]
xc_upper = []
xc_lower = []
for i in range(12):
    if i % 2 == 0:
        xc_lower.append(x[i] / 63)
    if i % 2 == 1:
        xc_upper.append(x[i] / 63)

plt.plot([Aerofoil(i) for i in xc_upper], xc_upper, 'ro',
         label='Upper Surface Tappings')
plt.plot([Aerofoil(i) for i in np.linspace(0, 1, 100)], np.linspace(0, 1, 100),
         'g-')
plt.plot([-Aerofoil(i) for i in xc_lower], xc_lower, 'bo',
         label='Lower Surface Tappings')
plt.plot([-Aerofoil(i) for i in np.linspace(0, 1, 100)],
         np.linspace(0, 1, 100), 'g-')
plt.gca().invert_yaxis()
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.grid()
plt.axis('scaled')
plt.show()

# density of air (fill with your own value depending on today's weather)
rho_w=1000 # density of water
U=25
g=9.81
theta=0 # inclination of manometer to the vertical in radians
delta_h = 0.5 * (rho/rho_w) * (1/(g*np.cos(theta))) * U**2
print("Required delta h: " + str(delta_h*1000) + "mm")

Angle_1 = 2      #Enter your angle of attack in degrees

Total_deg = 230    #This is the total pressure at infinity and is also used in computing Cp
Static_deg = 193   #This is the static pressure at infinity and is also used in computing Cp

#The below readings for one left hand side of the foil (if you turn the foil clockwise - this will be lower side)
deg_tapping1 = 210
deg_tapping2 = 176
deg_tapping3 = 168
deg_tapping4 = 170
deg_tapping5 = 176
deg_tapping6 = 180
#The below readings are right hand side of the foil (if you turn the foil clockwise - this will be the upper side)
deg_tapping7 = 172
deg_tapping8 = 156
deg_tapping9 = 158
deg_tapping10 = 166
deg_tapping11 = 172
deg_tapping12 = 184

deg1 = [Total_deg,Static_deg,deg_tapping1,deg_tapping2,deg_tapping3,deg_tapping4,
           deg_tapping5,deg_tapping6,deg_tapping7,deg_tapping8,deg_tapping9,deg_tapping10,
           deg_tapping11,deg_tapping12]


Index_deg1 = ['P_tot', 'P_static', 'LWR x/c = 0.0159', 'LWR x/c = 0.0714', 'LWR x/c = 0.175', 'LWR x/c = 0.317',
              'LWR x/c = 0.510','LWR x/c = 0.698','UPPER x/c = 0.0317', 'UPPER x/c = 0.119', 'UPPER x/c = 0.230',
              'UPPER x/c = 0.413', 'UPPER x/c = 0.603', 'UPPER x/c = 0.794']

d = {'Variable': Index_deg1,'Manometer Reading (mm)':deg1}
df = pd.DataFrame(d)

filename = str(Angle_1) + "_deg.csv"
#YOU NEED TO CHANGE THE FILE NAME TO MATCH THE ANGLE OF ATTACK THAT YOU USE
df.to_csv(filename, index=False)

Angle_2 = -4      #Enter your angle of attack in degrees

Total_deg = 231
Static_deg = 193   #This is the static pressure at infinity and is also used in computing Cp

#The below readings for one left hand side of the foil (if you turn the foil clockwise - this will be lower side)
deg_tapping1 = 155
deg_tapping2 = 135
deg_tapping3 = 141
deg_tapping4 = 152
deg_tapping5 = 168
deg_tapping6 = 178
#The below readings are right hand side of the foil (if you turn the foil clockwise - this will be the upper side)
deg_tapping7 = 213
deg_tapping8 = 185
deg_tapping9 = 178
deg_tapping10 = 178
deg_tapping11 = 181
deg_tapping12 = 185

deg2 = [Total_deg,Static_deg,deg_tapping1,deg_tapping2,deg_tapping3,deg_tapping4,
           deg_tapping5,deg_tapping6,deg_tapping7,deg_tapping8,deg_tapping9,deg_tapping10,
           deg_tapping11,deg_tapping12]
Index_deg2 = ['P_tot', 'P_static', 'LWR x/c = 0.0159', 'LWR x/c = 0.0714', 'LWR x/c = 0.175', 'LWR x/c = 0.317',
              'LWR x/c = 0.510','LWR x/c = 0.698','UPPER x/c = 0.0317', 'UPPER x/c = 0.119', 'UPPER x/c = 0.230',
              'UPPER x/c = 0.413', 'UPPER x/c = 0.603', 'UPPER x/c = 0.794']

d = {'Variable': Index_deg2,'Manometer Reading (mm)':deg2}
df = pd.DataFrame(d)

filename = str(Angle_2) + "_deg.csv"
#YOU NEED TO CHANGE THE FILE NAME TO MATCH THE ANGLE OF ATTACK THAT YOU USE
df.to_csv(filename, index=False)

Angle_3 = 12      #Enter your angle of attack in degrees

Total_deg = 232
Static_deg = 196   #This is the static pressure at infinity and is also used in computing Cp

#The below readings for one left hand side of the foil (if you turn the foil clockwise - this will be lower side)
deg_tapping1 = 234
deg_tapping2 = 218
deg_tapping3 = 200
deg_tapping4 = 192
deg_tapping5 = 190
deg_tapping6 = 190
#The below readings are right hand side of the foil (if you turn the foil clockwise - this will be the upper side)
deg_tapping7 = 82
deg_tapping8 = 108
deg_tapping9 = 124
deg_tapping10 = 156
deg_tapping11 = 170
deg_tapping12 = 182

deg3 = [Total_deg,Static_deg,deg_tapping1,deg_tapping2,deg_tapping3,deg_tapping4,
           deg_tapping5,deg_tapping6,deg_tapping7,deg_tapping8,deg_tapping9,deg_tapping10,
           deg_tapping11,deg_tapping12]
Index_deg3 = ['P_tot', 'P_static', 'LWR x/c = 0.0159', 'LWR x/c = 0.0714', 'LWR x/c = 0.175', 'LWR x/c = 0.317',
              'LWR x/c = 0.510','LWR x/c = 0.698','UPPER x/c = 0.0317', 'UPPER x/c = 0.119', 'UPPER x/c = 0.230',
              'UPPER x/c = 0.413', 'UPPER x/c = 0.603', 'UPPER x/c = 0.794']

d = {'Variable': Index_deg3,'Manometer Reading (mm)':deg3}
df = pd.DataFrame(d)

filename = str(Angle_3) + "_deg.csv"
#YOU NEED TO CHANGE THE FILE NAME TO MATCH THE ANGLE OF ATTACK THAT YOU USE
df.to_csv(filename, index=False)