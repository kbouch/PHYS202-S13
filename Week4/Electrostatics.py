import numpy as np

k = 8.987552e9 # units: Nm^2/C^2
e = 1.60217649e-19 # units: Coulumbs

def pointPotential(x,y,q,xo=0.0,yo=0.0,k=8.987552e9):
    """Takes two numbers,lists, tuples, or two arrays with the same shape, such as two 2D arrays x and y 
    created from np.meshgrid(). Also take charge of the point particle as a float including sign.
    Ouputs an array of the same shape with the scalar electric potential values evaluated at each array element.
    xo and yo are the coordinates of the charge in the xy plane, both default to zero. k is the coulumb force constant 
    1/4pi(epsilon_not)defualts to the SI value 8.987552e9 with units (Newton)(meter^2)/(Coulumb^2).
    Note: this module contains k and the fundamental charge e as constants.
    """
    if type(x) in (list,tuple):
        V = x
        for i in range(len(V)):
            for j in range(len(V[0])):
                V[i,j] = (k*q)/((x[i,j]-xo)**2.0 + (y[i,j]-yo)**2)**0.5
    else:
        V = V = (k*q)/((x-xo)**2.0 + (y-yo)**2)**0.5
    return V

def dipolePotential(x,y,q,d,a,k=8.987552e9,ang=0.0,xc=0.0,yc=0.0):
    """x, y are equally shaped ndarrays or lists, tuples, or coordinates. Returns a new array
    with the same shape, returning the result of evaluating the dipole potential for the respective x and y elements.
    Also needs k,q,d specified.
    k is the coulumb force constant 1/4pi(epsilon_not)defualts to the SI value 8.987552e9 with units (Newton)(meter^2)/(Coulumb^2).
    q is the charge magnitude of both separated point charges. If q > 0, then the positive charge goes on the positive axis direction.
    d is the distance each charge is from the origin along the y axis. The positive charge is on the positive y axis.
    a is the axis designation: takes 'x', 'y', or 'a'. If a second character is added after the letter, both charges will have the same sign.
    The optional arguments ang, xc, yc orient the dipole axis at an angle to the x axis. xc, yc change the central location 
    """
    if a[0] == 'x':
        xo = d
        yo = 0.0
    elif a[0] == 'y':
        xo = 0.0
        yo = d
    elif a[0] == 'a':
        xo = d*np.cos(ang)
        yo = d*np.sin(ang)
    else:
        raise 'Error: dipole axis invalid arg. 3rd arg must be x or y'
    if len(a) > 1:
        s = -1
    else:
        s = 1
    V = pointPotential(x,y,q,xo+xc,yo+yc,k) + pointPotential(x,y,-q * s,-xo+xc,-yo+yc)
    return V
