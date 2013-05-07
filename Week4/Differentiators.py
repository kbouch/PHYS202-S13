import numpy as np

def finiteDifference(x,y):
    """Finds the numerical derivative of a correlation given by x and y, 1D arrays of the same shape,
    where y = f(x). Uses centered finite difference differentiation for interior data points,
    forward or backward difference for the endpoints.
    """
    import numpy as np
    dxdy = np.zeros(x.shape,dtype='float')
    dxdy[1:-1] = 1.0*(y[2:] - y[:-2])/(x[2:] - x[:-2])
    dxdy[0] = 1.0*(y[1] - y[0])/(x[1]-x[0])
    dydx[-1] = 1.0*(y[-1] - y[-2])/(x[-1]-x[-2])
    return dydx


