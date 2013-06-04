def linear_least_squares(x,y):
    """ Takes two arrays which pair into points (x,y) for linearly varying data.
    Uses linear regression to find the slope and y-intercept of the best fit line.
    Returns a 4-tuple (slope,intercept,slope uncertainty,intercept uncertainty)
    """
    
    import numpy as np
    
    x_avg = np.sum(x)*1.0/x.shape[0]
    x2_avg = np.sum(x**2)*1.0/x.shape[0]
    xy_avg = np.sum(x*y)*1.0/x.shape[0]
    y_avg = np.sum(y)*1.0/y.shape[0]
    
    denom = 1.0*(x2_avg - x_avg**2)
    
    m = (xy_avg - (x_avg*y_avg))/denom
    b = (y_avg*x2_avg - x_avg*xy_avg)/denom
    
    dev = y - (m*x + b)
    dev2_avg = np.sum(dev**2)*1.0/dev.shape[0]
    
    m_unc = np.sqrt(dev2_avg/denom/(x.shape[0] - 2))
    b_unc = np.sqrt(dev2_avg*x2_avg/denom/(x.shape[0] - 2))
    
    return (m,b,m_unc,b_unc)

def weighted_linear_least_squares(x,y,yerr):
    """ Takes two arrays which pair into points (x,y) for linearly varying data.
    Also takes an array of the uncertainty in y at each point.
    Uses weighted linear regression to find the slope and y-intercept of the best fit line.
    Returns a 4-tuple (slope,intercept,slope uncertainty,intercept uncertainty)
    """
    
    import numpy as np
    
    if all(yerr == yerr[0]):
        w = np.ones(yerr.shape,dtype=float)
    else:
        w = np.divide(np.ones(yerr.shape,dtype=float),yerr**2)
    
    w_sum = np.sum(w) 
    wx_sum = np.sum(w*x)
    wx2_sum = np.sum(w*x**2)
    wxy_sum = np.sum(w*x*y)
    wy_sum = np.sum(w*y)
    
    wdenom = 1.0*((w_sum*wx2_sum) - wx_sum**2)
    
    m = ((w_sum*wxy_sum) - (wx_sum*wy_sum))/wdenom
    b = (wy_sum*wx2_sum - wx_sum*wxy_sum)/wdenom
    
    m_unc = np.sqrt(w_sum/wdenom)
    b_unc = np.sqrt(wx2_sum/wdenom)
    
    return (m,b,m_unc,b_unc)
