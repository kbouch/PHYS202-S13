def nuclei_log(symb):
    """Takes the symbol (string) for one of the nuclei of atoms that are in nuclei_dic.
    returns a tuple with (A,Z,N,mass [MeV/c**2],rho0 ,R0 [fm],a [fm**(-1)],w).
    A is atomic mass number, Z is number of protons, N is number of neutrons.
    denisty of nuclear matter rho0 = 0.169 [u/fm**3] 
    u is atomic mass units, fm = femptometers.
    Masses are 0.9315 GeV/(u*c**2) * (atomic mass [u]), The atomic masses for specific
    isotopes are from my physics textbook.
    """
    nuclei_dict = {'C' : ( 12,  6,   6, .9315*12        , 0.169, 2.47 , 0.   , 0.    ),# a should not be zero for Carbon?
                   'O' : ( 16,  8,   8, .9315*15.994915 , 0.169, 2.608, 0.513,-0.051 ),
                   'Al': ( 27, 13,  14, .9315*26.981539 , 0.169, 3.07 , 0.519, 0.    ),
                   'S' : ( 32, 16,  16, .9315*31.97207  , 0.169, 3.458, 0.61 , 0.    ),
                   'Ca': ( 40, 20,  20, .9315*39.962591 , 0.169, 3.76 , 0.586,-0.161 ),
                   'Ni': ( 58, 28,  30, .9315*57.935346 , 0.169, 4.309, 5.16 ,-0.1308),
                   'Cu': ( 63, 29,  34, .9315*62.939598 , 0.169, 4.2  , 0.869, 0.    ),
                   'W' : (186, 74, 112, .9315*185.954357, 0.169, 6.51 , 0.535, 0.    ),
                   'Au': (197, 79, 118, .9315*196.966543, 0.169, 6.38 , 0.535, 0.    ),
                   'Pb': (208, 82, 126, .9315*207.976627, 0.169, 6.68 , 0.546, 0.    ),
                   'U' : (238, 92, 146, .9315*238.050941, 0.169, 6.68 , 0.6  , 0.    )
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

def collision_events_sim(symb1,beamEcm,N,symb2=None):
    """ Runs the interaction with nuclei given by symb1, symb2 N times.
    Choose the beam center of mass energy beamEcm. Impact parameter is 
    random.
    """
    import numpy as np
    if symb2 == None:
        symb2 = symb1
    elem1 = nuclei_log(symb1)
    A1 = elem1[0]
    R01 = elem1[5]
    elem2 = nuclei_log(symb2)
    A2 = elem2[0]
    R02 = elem2[5]
        
    b_max = 1.1*(R01 + R02)
    crosssec,popt_tcs,popt_ecs = inelastic_cs(beamEcm)
    
    out_coords1 = np.empty((N,A1,2),dtype=float)
    out_coords2 = np.empty((N,A2,2),dtype=float)
    out_counts1 = np.empty((N,A1),dtype=int)
    out_counts2 = np.empty((N,A2),dtype=int)
        
    for rr in range(N):        
        azimu = gen_polar_direc()[0]
        b = b_max*np.random.random_sample()
        
        xypos_n1 = gen_nucleus2D(symb1)
        
        temp_xypos =gen_nucleus2D(symb2)
        
        xypos_n2[0] = b*np.cos(azimu)*temp_xypos[0]
        xypos_n2[1] = b*np.cos(azimu)*temp_xypos[1]
        
        out_coords1[rr,:,:] = xypos_n1
        out_coords2[rr,:,:] = xypos_n2
        
        out_counts1[rr,:], out_counts2[rr,:] = interaction_detect(xypos_n1,xypos_n2,crosssec)
        
    N_part = np.empty(N,dtype=int)
    N_coll = np.empty(N,dtype=int)
    
    for rr in range(N):
        total_part = 0
        total_coll = 0
        for NN in outcounts1[rr]:
            total_coll += NN
            if NN != 0:
                total_part += 1
        for NN in outcounts2[rr]:
            if NN != 0:
                total_part += 1
                
        N_part[rr] = total_part
        N_coll[rr] = total_coll
        
    return N_part,N_coll,out_coords1,out_coords2,out_counts1,out_counts2

def hadron_cs_model(sqrt_s,Z,Y1,Y2,m1,m2):
    """ Hadron scattering cross section versus system center of mass
    total energy (system is two particles with mass m1 and m2 colliding).
    Z,Y1,Y2 are parameter to be fit to data and should
    m1, m2 should be in GeV
    """
    import numpy as np
    eta1 = 0.462 #unitless
    eta2 = 0.550 # unitless
    M = 2.15 # GeV/c**2 (but let c = speed of light = 1)
    sM = (m1 + m2 + M)**2 
    B = np.pi*(9.390800635203e-22/M)**2  # hbar = 9.3908e-16 eV = 9.3908e-22 MeV
    
    return Z + B*(np.log((sqrt_s**2)*1.0/sM))**2 + Y1*(sM*1.0/(sqrt_s**2))**eta1 - Y2*(sM*1.0/(sqrt_s**2))**eta2

p_mass = 1.007276*0.9315 #proton mass u * (0.9315 GeV/(u*c**2))

def pp_cs_model(sqrt_s,Z,Y1,Y2):
    m1 = 1.007276*0.9315 #proton mass u * (0.9315 GeV/(u*c**2))
    m2 = m1
    return hadron_cs_model(sqrt_s,Z,Y1,Y2,m1,m2)

def convert_plab_Ecm(plab_arr):
    import numpy as np
    p_mass = 1.007276*0.9315 #proton mass u * (0.9315 GeV/(u*c**2))
    return np.sqrt(2*(p_mass)**2 + 2*p_mass*np.sqrt(plab_arr**2 + (p_mass)**2))


def inelastic_cs(Ecm):
    tcs_Ecm = convert_plab_Ecm(tcs_points[:,1])
    ecs_Ecm = convert_plab_Ecm(ecs_points[:,1])
    
    from scipy.optimize import curve_fit
    
    mask5GeV = ecs_Ecm > 5
    #         Z   Y1  Y2
    ecs_p0 = [8., 12.,7.]
    popt_ecs,pcov_ecs = curve_fit(pp_cs_model,ecs_Ecm[mask5GeV],ecs_points[:,4][mask5GeV] , ecs_p0\
                                 ,ecs_points[:,6][mask5GeV]+ecs_points[:,5][mask5GeV])
    
    mask5GeV_tcs = tcs_Ecm > 5
    tcs_p0 = [34., 12.,7.]
    popt_tcs,pcov_tcs = curve_fit(pp_cs_model,tcs_Ecm[mask5GeV_tcs],tcs_points[:,4][mask5GeV_tcs] , tcs_p0\
                                 ,tcs_points[:,6][mask5GeV_tcs]+tcs_points[:,5][mask5GeV_tcs])
    
    return 0.1*(pp_cs_model(Ecm,*popt_tcs)-pp_cs_model(Ecm,*popt_ecs)),popt_tcs,popt_ecs
