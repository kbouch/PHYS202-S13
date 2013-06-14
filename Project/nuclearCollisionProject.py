def nuclei_log(symb):
    """Takes the symbol (string) for one of the nuclei of atoms that are in nuclei_dic.
    returns a tuple with (A,Z,N,rho0,R0,a,w).
    A is atomic mass number, Z is number of protons, N is number of neutrons.
    denisty of nuclear matter rho0 = 0.169 [u/fm**3] 
    u is atomic mass units, fm = femptometers.
    """
    nuclei_dict = {'C' : ( 12,  6,   6, 0.169, 2.47 , 0.   , 0.    ),# a should not be zero for Carbon?
                   'O' : ( 16,  8,   8, 0.169, 2.608, 0.513,-0.051 ),
                   'Al': ( 27, 13,  14, 0.169, 3.07 , 0.519, 0.    ),
                   'S' : ( 32, 16,  16, 0.169, 3.458, 0.61 , 0.    ),
                   'Ca': ( 40, 20,  20, 0.169, 3.76 , 0.586,-0.161 ),
                   'Ni': ( 58, 28,  30, 0.169, 4.309, 5.16 ,-0.1308),
                   'Cu': ( 63, 29,  34, 0.169, 4.2  , 0.869, 0.    ),
                   'W' : (186, 74, 112, 0.169, 6.51 , 0.535, 0.    ),
                   'Au': (197, 79, 118, 0.169, 6.38 , 0.535, 0.    ),
                   'Pb': (208, 82, 126, 0.169, 6.68 , 0.546, 0.    ),
                   'U' : (238, 92, 146, 0.169, 6.68 , 0.6  , 0.    )
                   }
    return nuclei_dict(symb)

def gen_polar_direc(N=1):
    """ generates N random polar coordinate direction
    (azimuthal,polar) angles. 0<=azimuth<=2*pi, 0<=polar<=pi.
    """
    import numpy.random.random_sample
    import numpy.pi
    out = np.zeros((N,2),dtype=float)
    for i in range(N):
        out[i,0] = (2*pi)*random_sample()
        out[i,1] = (pi)*random_sample()
    return out

def gen_nucleus2D(symb,rpob='CDF'):
    """ Generates random (x,y) coordinates for nucleons, up to the number
    of nucleons in the nucleus of the atom given by its chemical symbol.
    The origin is at the center of the nucleus. Allows nucleons to fall
    within 1.1*R0 femptometers of the center. Indicate CDF of MC to choose
    cumulative distribution function transformation (numerical integration)
    for randomly choosing radial position of nucleons. MC uses monte carlo
    hit-or-miss method instead.
    """
    import probDistributions
    import numpy as np
    
    elem_nucl = nuclei_log(symb)
    A = elem_nucl[0]
    R0 = elem_nucl[4]
    
    xypos = np.zeros((A,2),dtype=float)
    
    azimupos = gen_polar_direc(A)[:,0]
    if rprob == 'CDF':
        rpos = ditsribute1D(np.arange(0,1.101*R0,R0/1000)\
                            ,lambda r: woods_saxon(r,*elem_nucl[3:])\
                            ,A)
    elif rprob == 'MC':
        probx = lambda r: woods_saxon(r,*elem_nucl[3:])
        probxmax = rough_find_max(0,1.1*R0,probx,R0/1000)
        rpos = MCdistribute1D(probx,probxmax,0.,1.1*R0,A)
    
    xypos[:0] = rpos*np.cos(azimupos)
    xypos[:1] = rpos*np.sin(azimupos)
    
    return xypos

def interaction_detect(xypos_n1,xypos_n2,crosssec):
    """Finds all inelastic collisions between nucleons in nucleus 1
    with position xypos_n1 and nuceons in nucleus 2 with position 
    xypos_n2. Returns two arrays with the integer number of collisions
    each nucleon experienced. The total number of collisions is the 
    sum of elements in one of these arrays (both have same total).
    The number of participants in a given nucleus is the number of 
    nonzero elements in one coll_count array. Total number of participants
    is the number on nonzero elements in both output arrays.
    """
    import numpy as np
    d = np.sqrt(crosssec/np.pi)
    coll_count1 = np.zeros(xypos_1.shape[0],dtype=int)
    coll_count2 = np.zeros(xypos_2.shape[0],dtype=int)
    for j in range(xypos_n1.size):
        for k in range(xypos_n2.size):
            if np.sqrt((xypos_n1[j,0]-xypos_n2[k,0])**2\
                      +(xypos_n1[j,1]-xypos_n2[k,1])**2) < d:
                coll_count1[j] += 1
                coll_count2[k] += 1
    return coll_count1, coll_count2
