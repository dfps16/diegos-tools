# -*- coding: utf-8 -*-
"""
Created on Mon Apr 21 10:56:55 2025
Generate a wing AR, TR and NACA 4-digit airfoil section from a student ID
with predictions from thin aerofoil and lifting line theory
@author: nds9
"""

import numpy as np


# ------------------------------------------------------------------------------
# Function to solve thin aerofoil theory to find the zero lift angle
# ------------------------------------------------------------------------------
def zero_lift(m, p, n):
    dtheta = np.pi / n
    intlast = 0.0
    sum0 = 0.0
    for i in range(n):
        theta = (i + 1) * dtheta
        x = 0.5 * (1.0 - np.cos(theta))
        # use NACA 4-digit camber line, defined in two parts x<p and x>p
        if x < p:
            intnext = (
                2.0 * m / p**2 * (p - x) * (1.0 - np.cos(theta))
            )  # i.e. (dz/dx)*(1-cos(theta))
        else:
            intnext = 2.0 * m / (1.0 - p) ** 2 * (p - x) * (1.0 - np.cos(theta))
        sum0 = (
            sum0 + 0.5 * dtheta * (intlast + intnext) / np.pi
        )  # trapezoid integration
        intlast = intnext
    alpha_zero_lift = sum0 * 180.0 / np.pi  # final result in degrees
    return alpha_zero_lift


# ------------------------------------------------------------------------------
# Function to solve lifting line theory (basic formulation, no twist)
# ------------------------------------------------------------------------------
def lifting_line(alpha, AR, tr, a0, alpha_zero, n):
    A = np.zeros((n, n))  # initialise the A matrix
    rhs = np.zeros(n)
    coeffs = np.zeros(n)
    for i in range(n):
        theta = (
            np.pi / 2.0 * np.real(2 * (i + 1) - 1) / np.real(2 * n)
        )  # collocation points for theta
        c = (2.0 * (1.0 + (tr - 1.0) * np.cos(theta))) / (
            1.0 + tr
        )  # chord for tapered wing
        rhs[i] = alpha - alpha_zero  # set rhs vector
        for j in range(n):
            k = 2 * (j + 1) - 1
            A[i, j] = np.sin(k * theta) * (
                4.0 * AR / (a0 * c) + k / np.sin(theta)
            )  # set matrix elements
    # solve matrix system for Fourier coefficients
    coeffs = np.linalg.solve(A, rhs)
    # find aerodynamic coefficients
    CL = coeffs[0] * np.pi * AR  # wing lift coefficient
    sum0 = 0.0
    for j in range(n):
        k = 2 * (j + 1) - 1
        sum0 = sum0 + k * (coeffs[j] / coeffs[0]) ** 2
    delta = sum0 - 1.0  # induced drag factor
    CDi = np.pi * AR * coeffs[0] ** 2 * sum0  # induced drag
    return CL, CDi, delta


# ------------------------------------------------------------------------------
# Start of main program
# ------------------------------------------------------------------------------
# First get the airfoil properties from ID
ID = input("\nEnter student ID (must be 8 chararacters): ")
if len(ID) == 8:
    NACA1 = int(1 + 2 * int(ID[2]) / 5)
    NACA2 = int(3 + 2 * int(ID[3]) / 5)
    NACA4 = int(int(ID[4]) / 2)
    NACA = str(NACA1) + str(NACA2) + "1" + str(NACA4)
    print("\nNACA airfoil = ", NACA)
    ARhalf = str(int(3 + 1 * int(ID[6]) / 2))
    print("Half-wing aspect ratio = ", ARhalf)
    TR = "0." + str(int(3 + int(ID[7]) / 2))
    print("Taper ratio = ", TR)
    Velocity = 35 + int(ID[1]) + int(ID[5])
    print("Velocity = ", Velocity, " m/s")

    # Find the zero lift angle using thin aerofoil theory (see SESA2022 for theory)
    nthinaerofoil = (
        64  # number of points to use for thin aerofoil evaluation of zero lift angle
    )
    alpha_zero_lift = zero_lift(
        np.real(NACA1) / 100.0, np.real(NACA2) / 10.0, nthinaerofoil
    )
    print("\nThin aerofoil theory results:")
    print("   Zero lift angle (degrees) = ", round(alpha_zero_lift, 3))

    # Solve the lifting line equation for a range of angles of incidence alpha (see SESA2022 for theory)
    a0 = 2.0 * np.pi  # lift slope from 2D thin aerofoil theory
    n_lift_line = 64  # number of spanwise locations to use for lifting line theory

    print("\nLifting line theory results (using full span wing AR):")
    AR = 2.0 * float(ARhalf)
    dalpha = 0.01
    (CL0, CDi, delta) = lifting_line(
        0.0, AR, float(TR), a0, alpha_zero_lift * np.pi / 180.0, n_lift_line
    )
    (CL1, CDi, delta) = lifting_line(
        dalpha, AR, float(TR), a0, alpha_zero_lift * np.pi / 180.0, n_lift_line
    )
    a = (CL1 - CL0) / dalpha
    tau = np.pi * AR / a0 * ((a0 / a) - 1.0) - 1.0
    print("   Lift slope a = ", round(a, 3))
    print("   Lift correction tau = ", round(tau, 3))
    print("   Induced drag correction delta = ", round(delta, 3))

    print("\nalpha(deg)  CL        CDi")
    for ialpha in range(9):
        alpha = (-6.0 + ialpha * 2.0) * np.pi / 180.0
        (CL, CDi, delta) = lifting_line(
            alpha, AR, float(TR), a0, alpha_zero_lift * np.pi / 180.0, n_lift_line
        )
        # print(round(alpha*180/np.pi,2),'      ',round(CL,4),'  ',round(CDi,5))
        print(
            "{:5.2f}".format(alpha * 180 / np.pi),
            "{:11.4f}".format(CL),
            "{:10.5f}".format(CDi),
        )

else:
    print("Not 8 characters, try again")
# ------------------------------------------------------------------------------
# End of program
# ------------------------------------------------------------------------------
