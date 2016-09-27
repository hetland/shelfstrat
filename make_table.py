import numpy as np
import netCDF4

def polyfit_omega(n=6):
    'fit an order-n polynomial to the maximum growth rate as a function of delta, the slope parameter.'
    delta, mu = np.mgrid[-1.2:2.2:1001j, 0:4.2:1001j]
    
    tmu = np.tanh(mu)
    omega = np.sqrt( (1.0+delta)*(mu - tmu)/tmu 
                     -0.25*(delta/tmu + mu)**2 ).real
    omega2 = (1.0+delta)*(mu - tmu)/tmu - 0.25*(delta/tmu + mu)**2
    
    omega = np.ma.masked_where(np.isnan(omega), omega)
    
    omega_max = omega.max(axis=1)
    idx = np.where(~omega_max.mask)
    
    omega_max = omega_max[idx]
    delta = delta[:, 0][idx]
    
    p = np.polyfit(delta, omega_max, n)
    
    return p

omega_poly = polyfit_omega()

deltas = []
Ris = []
Ss = []
M2s = []
fs = []

N2 = 1e-4
alpha = 1e-3

for S in [0.1, 0.2, 0.3, 0.5]:
    for Ri in [1, 2, 3, 5, 10]:
        Ss.append(S)
        Ris.append(Ri)
        
        delta = S*np.sqrt(Ri)
        deltas.append(delta)
        
        f = np.sqrt(N2)*alpha/S
        fs.append(f)
        
        M2 = np.sqrt(N2*f**2/Ri)
        M2s.append(M2)

for delta in [0.1, 0.2, 0.3, 0.5]:
    for Ri in [1, 2, 3, 5, 10]:
        deltas.append(delta)
        Ris.append(Ri)
        
        S = delta/np.sqrt(Ri)
        Ss.append(S)
        
        f = np.sqrt(N2)*alpha/S
        fs.append(f)
        
        M2 = np.sqrt(N2*f**2/Ri)
        M2s.append(M2)
        
tab = np.vstack((M2s, fs, Ris, deltas, Ss)).T
_, idx = np.unique(tab.sum(axis=-1), return_index=True)
idx.sort()
tab = tab[idx]

# shelfstrat_M2_3.33e-07_N2_1.00e-04_f_5.77e-05

times = []
omegas = []
Tis = []
Tfs = []

for case in tab:
    # print('%5.2e, %5.2e, %5.2f, %5.2f, %5.2f' % tuple(case))
    filename = 'simulations/shelfstrat_M2_%5.2e_N2_1.00e-04_f_%5.2e/shelfstrat_his.nc' % (case[0], case[1])
    print filename
    nc = netCDF4.Dataset(filename)
    times.append(nc.variables['ocean_time'][-1]/86400)
    
    S = case[4]
    delta = case[3]
    Ri = case[2]
    f = case[1]
    
    omega = np.polyval(omega_poly, delta) # non-dim
    omega_dim = 86400.0 * omega * f / np.sqrt(Ri) # rad/days
    
    omegas.append(omega_dim)  
    
    Tis.append(np.sqrt(S) / omega_dim)
    
    ###################################
    diafilename = 'simulations/shelfstrat_M2_%5.2e_N2_1.00e-04_f_%5.2e/shelfstrat_dia.nc' % (case[0], case[1])
    nd = netCDF4.Dataset(diafilename)
    
    ubar_bstr = nd.variables['ubar_bstr'][:, 1:51, :]
    stress_mag = np.sqrt(ubar_bstr**2)
    stress = np.sqrt((stress_mag**2).mean(axis=-1).mean(axis=-1))
    
    ubar = Ri**(-1./2.) * np.sqrt(1e-4) * 50.0 / 2.0
    
    Tfs.append(ubar/stress[1]/86400.0)




tab = np.vstack((tab.T, times, omegas, Tis, Tfs)).T

np.savetxt('table.dat', tab, fmt='%5.2e & %5.2e & %5.2f & %5.2f & %5.2f & %6.1f & %5.3f & %5.3f & %5.3f \\\\')



