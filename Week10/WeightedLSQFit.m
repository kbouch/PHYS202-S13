function [m b m_unc b_unc] = WeightedLSQFit(x,y,yerr)
%WeightedLSQFit(x,y,yerr)
%    Takes two arrays which pair into points (x,y) for linearly varying 
%    data. Also takes an array of the uncertainty in y at each point.
%
%    Uses weighted linear regression to find the slope and y-intercept of 
%    the best fit line. Returns a 1x4 array
%    [slope  intercept  slope_uncertainty  intercept_uncertainty]

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