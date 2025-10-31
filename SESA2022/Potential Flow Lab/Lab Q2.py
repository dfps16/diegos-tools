import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
from scipy.interpolate import RectBivariateSpline as interp_2D
from scipy.optimize import fsolve


def create_grid(Domain_x_min,Domain_x_max,Domain_z_min,Domain_z_max,dx,dz):
    """This function creates a 2D grid between xmin and xmax as well zmin and zmax
    with spacing dx and dz. The function also returns a vector (1D array) x and z coordinates"""

    #This command makes 1-D array (or vectors) for x and z axis
    #If you ever want to know just the x-coordinates or the z-coordinate then just
    # use xa and za or that
    xa = np.arange(Domain_x_min,Domain_x_max+dx,dx)
    za = np.arange(Domain_z_min,Domain_z_max+dx,dz)

    #This command makes a mesh (2-D array) using xa and za, This is now our domain
    x,z = np.meshgrid(xa,za)
    return[xa,za,x,z]


 #The input for this function is the streamfunction
def findvelocities(psi):
    """This command computes the w and u components of velocity as gradient(psi). The gradient along the columns is taken in first and then the rows. The function gradient computes the first derivative of a variable along a given axis."""
    u,w = np.gradient(psi)

    u = u/dz
    w = -1*w/dx

    return [u,w]


def mask_box(xms,xmf,zms,zmf,psi,u,w):
    """This function takes the four corners of a rectangular box and masks the data within the box"""
    psi_m = ma.masked_where((x>xms)&(x<=xmf)&(z>zms)&(z<=zmf),psi)
    u_m = ma.masked_where((x>xms)&(x<=xmf)&(z>zms)&(z<=zmf),u)
    w_m = ma.masked_where((x>xms)&(x<=xmf)&(z>zms)&(z<=zmf),w)
    return[psi_m,u_m,w_m]


def find_index(axis_vec,loc):
    """This function finds the nearest index in a vector that meets your condition"""
    #This finds x-index for the location
    difference_array = np.absolute(axis_vec-loc)
    # find the index of minimum element from the array
    index = difference_array.argmin()
    return index


def interpolate_values(xa,za,func,xl,zl):
    """This function takes the 1D array (or vector) of x and z coordinates as well as the 2D array of data that you want to interpolate. It porvide the value of the function at xl and zl points specified by you."""
    fun =  interp_2D(za, xa, func) #You dont need to use 2D array of x and y for this.
    #Just the 1D array with coordinates of y and x are needed.
    #You also need to pass y coordinates first and then the x coordinates
    interp_value = fun(zl,xl)[0] #This interpolates the values of the function f in the new set of points.
    #Note that you need to pass y coodinate first and then the x-coordinate
    return interp_value

# ------- Typical Flows -------

def uniform_flow(Uinf,alpha):
    """This command is the stream function for Uniform flow. The input are Uinf and angle alpha (alpha is in degrees)"""
    psiu = Uinf*np.cos(np.pi*alpha/180)*z - Uinf*np.sin(np.pi*alpha/180)*x

    return psiu


def source(Q,x0,z0):
     """The input to this function is the strength and the location (x0,z0) of source/sink. Positive Q is source and negative Q is sink."""
     #This command finds the radial coordinate r from the cartesian coordinate
     #relative to this specific flow element

     r = np.sqrt((x-x0)**2+(z-z0)**2)

     #This command finds the tangential angle theta from the cartesian coordinate
     #relative to this specific flow element

     theta = np.arctan2(z-z0,x-x0)

     #This command is the stream function for source/sink
     psis = Q*theta/(2*np.pi)

     return psis


def doublet(Kappa,x0,z0):
      """The input to this function is the strength and the location (x0,z0) of the doublet. Positive doublet is source on the left and sink on the right. Negative doublet is sink on the left and source on ther right."""
      #This command finds the radial coordinate r from the cartesian coordinate
      #relative to this specific flow element
      r = np.sqrt((x-x0)**2+(z-z0)**2)

      #This command finds the tangential angle theta from the cartesian coordinate
      #relative to this specific flow element
      theta = np.arctan2(z-z0,x-x0)

      #This command is the stream function for source/sink
      psid = -Kappa*(z-z0)/(2*np.pi*r*r)

      return psid


# -------- Plotting Flows --------
# --------------------------------
Domain_x_min = -2
Domain_x_max = 2
Domain_z_min = -2
Domain_z_max = 2
dx = dz = 0.22

# We first create grid to compute the flow
[xa,za,x,z] = create_grid(Domain_x_min,Domain_x_max,Domain_z_min,Domain_z_max,dx,dz)

Uinf = 10
alpha = 0

# Stream function input
psi = doublet(2 * 3.14, - 0.5, 0) - doublet(2 * 3.14, 0.5, 0)
# ---------------------

w = np.gradient(psi) / dx

c_p = 1 - w**2 / (2 * 3.14 / 2 * np.pi() * 0.5 ** 2)
#

plt.figure(figsize=(8,8), dpi=100)
plt.plot(z, w)

[u,w] = findvelocities(psi)

plt.figure(figsize=(5,5),dpi=300)# make the plot
plt.contour(x,z,psi,50)       # streamlines
#Note that you can change the number of contours

plt.colorbar()
plt.quiver(x,z,u,w, scale=300)
plt.axis('equal')
plt.xlabel('x')
plt.ylabel('z')
plt.show()