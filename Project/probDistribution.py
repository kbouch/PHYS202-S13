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