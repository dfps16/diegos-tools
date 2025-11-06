import numpy as np
import matplotlib.pyplot as plt

# For masking the data
import numpy.ma as ma
from numpy.random import uniform

# Interpolation method
from scipy.interpolate import RectBivariateSpline as interp_2D

# For root-finding
from scipy.optimize import fsolve


def create_grid(Domain_x_min,Domain_x_max,Domain_z_min,Domain_z_max,dx,dz):
    """This function creates a 2D grid between xmin and xmax as well zmin and
    zmax with spacing dx and dz. The function also returns a vector (1D array)
    x and z coordinates"""

    # This command makes 1-D array (or vectors) for x and z axis
    # If you ever want to know just the x-coordinates or the z-coordinate
    # then jus t use xa and za or that
    xa= np.arange(Domain_x_min,Domain_x_max+dx,dx)
    za= np.arange(Domain_z_min,Domain_z_max+dx,dz)
    # This command makes a mesh (2-D array) using xa and za, This is now our
    # domain
    x,z = np.meshgrid(xa,za)
    return[xa,za,x,z]


def findvelocities(psi):
    """This command computes the w and u components of velocity as gradient
    (psi). The gradient along the columns is taken in first and then the rows.
    The function gradient computes the first derivative of a variable along a
    given axis."""
    u,w = np.gradient(psi)
    u = u/dz
    w = -1*w/dx

    return [u,w]


def findpressure(u,w,Uref):
    """This command computes the pressure coefficient from the velocities and a
    reference speed"""
    Cp = 1-((u**2+w**2)/Uref**2)
    return Cp


def mask_box(xms,xmf,zms,zmf,psi,u,w):
    """This function takes the four corners of a rectangular box and masks the
    data within the box"""
    psi_m = ma.masked_where((x > xms) & (x <= xmf) & (z > zms) & (z <= zmf), psi)
    u_m = ma.masked_where((x > xms) & (x <= xmf) & (z > zms) & (z <= zmf), u)
    w_m = ma.masked_where((x > xms) & (x <= xmf) & (z > zms) & (z <= zmf), w)
    return[psi_m,u_m,w_m]


def mask_data(val,val_min,val_max,psi,u,w):
    """This function masks data between a min and max value for a given parameter.
    We use 1 = psi, 2 = u and 3 = w"""
    if (val == 0):
        valf = psi
    if (val == 1):
        valf = u
    if (val == 2):
        valf = w

    psi_m = ma.masked_where((valf > val_min) & (valf < val_max), psi)
    u_m = ma.masked_where((valf > val_min) & (valf < val_max), u)
    w_m = ma.masked_where((valf > val_min) & (valf < val_max), w)
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

# ------ Flows ------

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


def vortex(Gamma,x0,z0):
  """The input to this function is the circulation and the location (x0,z0) of the vortex. Positive Gamma is rotation in clockwise direction. Negative Gamma is rotation in counter-clockwise direction."""

  #This command finds the radial coordinate r from the cartesian coordinate
  #relative to this specific flow element
  r = np.sqrt((x-x0)**2+(z-z0)**2)

  #This command finds the tangential angle theta from the cartesian coordinate
  #relative to this specific flow element
  theta = np.arctan2(z-z0,x-x0)

  #This command is the stream function for a vortex
  psiv = Gamma*np.log(r)/(2*np.pi)

  return psiv


def find_corner_max(u, x, z, corner='top-right', frac=0.2):
    """
    Find largest u inside one corner region.
    corner: 'top-right', 'top-left', 'bottom-right', 'bottom-left'
    frac: fraction of domain to treat as the corner (0 < frac <= 1)
    Returns: (max_val, x_loc, z_loc, (i,j)) or None if corner empty
    """
    xa_min, xa_max = x.min(), x.max()
    za_min, za_max = z.min(), z.max()
    dx = xa_max - xa_min
    dz = za_max - za_min

    if corner == 'top-right':
        mask = (x >= xa_max - frac*dx) & (z >= za_max - frac*dz)
    elif corner == 'top-left':
        mask = (x <= xa_min + frac*dx) & (z >= za_max - frac*dz)
    elif corner == 'bottom-right':
        mask = (x >= xa_max - frac*dx) & (z <= za_min + frac*dz)
    elif corner == 'bottom-left':
        mask = (x <= xa_min + frac*dx) & (z <= za_min + frac*dz)
    else:
        raise ValueError("corner must be one of: top-right, top-left, bottom-right, bottom-left")

    inds = np.where(mask)
    if inds[0].size == 0:
        return None

    vals = u[inds]
    arg = np.argmax(vals)
    i = inds[0][arg]
    j = inds[1][arg]
    return u[i, j], x[i, j], z[i, j], (i, j)


# ------- CODE INPUT -------

# YOUR CODE HERE TO PLAY AROUND WITH THE VARIOUS FLOW ELEMENTS!!!
Domain_x_min = -4
Domain_x_max = 4
Domain_z_min = -4
Domain_z_max = 4
dx = dz = 0.01

# We first create grid to compute the flow
[xa,za,x,z] = create_grid(Domain_x_min,Domain_x_max,Domain_z_min,Domain_z_max,dx,dz)

# HERE I AM DOING A UNIFORM FLOW + VORTEX + SOURCE
Uinf = 2.4
alpha = 0
Q = - 2*np.pi
Gamma = 3.3

psi_flow = uniform_flow(Uinf, 0) + vortex(3.3, 0, 0)
wall = 0
psi = psi_flow + wall

[u, w] = findvelocities(psi)

xms = 1
xmf = 1
zms = 1
zmf = 1
# psi,u,w = mask_box(-0.5*dx,0.5*dx,-0.5*dz,0.5*dz,psi,u,w)

max_u = np.max(u)
print(max_u)

# Find streamline value at point
x1 = 2
x2 = 1
z1 = 4
z2 = 3.2
midpoint_x = (x1+x2)/2
midpoint_z = (z1+z1)/2
st_x_index = find_index(xa,midpoint_x)
st_z_index = find_index(za,midpoint_z)
vol_flow = psi_flow[st_x_index,st_z_index]
m_flow = vol_flow * 1
print(np.round(m_flow, 2))

Cp = findpressure(u, w, Uinf)

x_axis = find_index(za, 0)
z_axis = find_index(xa, 0)
Cp_slice_x = Cp[x_axis, :]
Cp_slice_z = Cp[:, z_axis]

min_idx = np.argmin(Cp_slice_x)
max_idx = np.argmax(Cp_slice_x)

min_idz = np.argmin(Cp_slice_z)
max_idz = np.argmax(Cp_slice_z)

plt.figure(figsize=(6, 5), dpi=120)
plt.plot(xa, Cp_slice_x)

# Plot min point
plt.plot(xa[min_idx], Cp_slice_x[min_idx], 'go', markersize=8,
         label=f'Min = {Cp_slice_x[min_idx]:.2f}, loc = {xa[min_idx]:.2f}')

# Plot max point
plt.plot(xa[max_idx], Cp_slice_x[max_idx], 'ro', markersize=8,
         label=f'Max = {Cp_slice_x[max_idx]:.2f}, Loc = {xa[max_idx]:.2f}')
plt.legend(loc='best')

plt.title('Cp variation along x')
plt.xlabel('x')
plt.ylabel('Cp')

plt.figure(figsize=(6, 5), dpi=120)
plt.plot(za, Cp_slice_z)
# Plot min point
plt.plot(za[min_idz], Cp_slice_z[min_idz], 'go', markersize=8,
         label=f'Min = {Cp_slice_z[min_idz]:.2f}, loc = {za[min_idz]:.2f}')

# Plot max point
plt.plot(za[max_idz], Cp_slice_z[max_idz], 'ro', markersize=8,
         label=f'Max = {Cp_slice_z[max_idz]:.2f}, Loc = {za[max_idz]:.2f}')
plt.legend(loc='best')

plt.title('Cp variation along z')
plt.xlabel('z')
plt.ylabel('Cp')



plt.figure(figsize=(5, 5), dpi=120)
plt.contour(x, z, psi, 50)
plt.quiver(x, z, u, w)
plt.colorbar()

plt.show()

#plt.figure(figsize=(5,5),dpi=120)
#plt.contourf(x,z,Cp,cmap='bwr')
#plt.colorbar()


