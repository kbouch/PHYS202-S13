def CDFinv(z,F,x):
    """ Takes an array z of randomly sampled floats on the interval 
    [0.,1.) such as with numpy.random.random(10000). Also takes an array of 
    values of the CDF transform function F evaluated at each value in
    an array x. x is also passed. F and x have the same shape. F is the CDF
    transform of a probability distribution funcion f = probx
    Returns the x values by taking the inverse transform of the z
    (random floats from 'continous uniform' distribution which are spread 
    over F (cumulative probability density of f)). Returns x_samples
    that have probability given by f.
    """
    x_indicies = np.searchsorted(F,z)
    x_samples = []
    for ind in x_indicies:
        x_samples.append(x[ind])
    return np.array(x_samples,dtype=float)


def distribute1D(x,probx,N,zseed=None):
    """Takes an array x of values of a quantity and a second array of the
    same length giving the relative probability of corresponding
    values in x to occur. Generates N random vaules with same distribution
    given by probx, with outputs landing only on values given in x
    (so space x values closely enough depending on the accuray you want)."""
    
    from scipy import integrate
    import numpy.random as rand
    
    A = 1.0/integrate.trapz(probx,x)
    
    F = np.ones(x.shape[-1],dtype=float)
    
    for i in range(x.shape[-1]):
        F[i] = A*integrate.trapz(probx[:i+1],x[:i+1])
    
    rand.seed(zseed)
    z = rand.random_sample(N)
    
    return CDFinv(z,F,x)

def woods_saxon(r,rho0=1.0,R0=3.,a=0.5,w=0.0):
    """ Radial density of nucleon matter in a nucleus.
    rho0 is the density at the center (r=0). R0 is the average
    radius. If w=0, R0 is the radius at which the density reaches
    half rho0. w and a are other parameters
    """
    import numpy as np
    return rho0*(1 + w*(r*1.0/R0)**2)/(1.0 + np.exp(r-R0*1.0/a))

def MCdistribute1D(probx,probxmax,xmin,xmax,N=1,seed=None):
    """Returns an array with N random values following
    the probability distribution function pro. Uses Monte
    Carlo hit-or-miss method. probx function should be callable
    requiring one argument.
    """
    import numpy as np
    x_samples_L = []
    while len(x_samples_L) < N:
        shotx = (xmax - xmin)*np.random.random_sample() + xmin
        height = probxmax**np.random.random_sample()
        if heigth <= probx(shotx):
            x_samples_L.append(shotx)
    return np.array(x_samples_L,dtype=float)

def rough_find_max(xmin,xmax,f,step):
    import numpy as np
    xtest = np.arange(xmin,xmax+step,step)
    return np.max(f(xtest))
