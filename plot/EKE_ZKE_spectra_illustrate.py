import matplotlib
font = {'size': 18,
        'family':'sans-serif','sans-serif':['Helvetica']}

matplotlib.rc('font', **font)
matplotlib.rc('lines', linewidth=2)
matplotlib.rc('text', usetex=True)

import numpy as np
import matplotlib.pyplot as plt
k_beta = 15
k_Rh = 4
k_z  = 40
k_f  = 60
k_nu = 180
epsilon = 1.0
beta = (2*k_beta)**(5./3)
ks_all  = np.arange(0,240).astype(float)
EKE = 6*epsilon**(2./3)*ks_all**(-5./3) * np.tanh(ks_all/k_Rh*1.5)**6
ZKE = 0.5*beta**2 * ks_all**(-5.) * np.tanh(ks_all/k_Rh)**8
ks_z_lg = np.arange(1,k_Rh+2)
ks_z_certain = np.arange(k_Rh+1,k_z)
ks_z_small = np.arange(k_z,k_z+10)
ks_eddy_lg = np.arange(1,k_Rh+5)
ks_eddy_certain = np.arange(k_Rh+5,k_f)
ks_smallscales = ks_all[k_f:].astype(int)
EKE[k_f:]  = EKE[k_f]*((ks_smallscales/1.0/k_f)**-3)
EKE[k_nu:] = EKE[k_nu]*(np.arange(k_nu, 240)/1.0/k_nu)**-10
plt.loglog(ks_all[ks_z_lg], ZKE[ks_z_lg],'--r')
plt.loglog(ks_all[ks_z_certain], ZKE[ks_z_certain],'r',
           label=r'Zonal kinetic energy')
plt.loglog(ks_all[ks_z_small], ZKE[ks_z_small],'--r')
plt.loglog(ks_all[ks_eddy_lg], EKE[ks_eddy_lg],'--b')
plt.loglog(ks_all[ks_eddy_certain], EKE[ks_eddy_certain], 'b',label=r'Eddy kinetic energy')
plt.loglog(ks_all[k_f:k_nu], EKE[k_f:k_nu],'b')
plt.loglog(ks_all[k_nu:], EKE[k_nu:], '--b')
plt.plot([k_Rh, k_Rh], [ZKE[k_Rh], 1e-3], '--k')
plt.plot([k_beta-0.82, k_beta-0.82], [EKE[k_beta-1], 1e-3], '--k')
plt.plot([k_f, k_f], [EKE[k_f], 1e-3], '--k')

fontsize_text = 16
fontsize_legend = 18
plt.text(k_Rh-0.5, 5e-4, r'$k_{r}\approx k_{Rh}$',fontsize=fontsize_text)
plt.text(k_beta, 5e-4, r'$k_\varepsilon$',fontsize=fontsize_text)
plt.text(k_f, 5e-4, r'$k_{f}$',fontsize=fontsize_text)
plt.text(8, 8e-1, r'$\beta^{2}k^{-5}$', color='red',fontsize=fontsize_text, rotation=-60)
plt.text(22, 4e-2, r'$\varepsilon^{2/3}k^{-5/3}$', color='blue',fontsize=fontsize_text, rotation=-30)
plt.text(82, 3e-3, r'$\eta^{2/3}k^{-3}$', color='blue',fontsize=fontsize_text, rotation=-50)
plt.legend(loc='upper right',frameon=False,fontsize=fontsize_text)
plt.xlabel('$k$')
plt.ylabel(r'$\mathcal{E}(k)$')
plt.setp( plt.gca().get_xticklabels(), visible=False)
plt.setp( plt.gca().get_yticklabels(), visible=False)
ax = plt.gca()
ax.annotate(r'$\varepsilon$', xy=(k_f-32, EKE[k_f-32]/1.8),
            xytext=(k_f-15, EKE[k_f-15]/3), arrowprops=dict(facecolor='gray', shrink=0.05, color='gray'),
            color='gray')
ax.annotate(r'$\eta$', xy=(k_f+60, EKE[k_f+60]/2),
            xytext=(k_f+10, EKE[k_f+10]/2.5), arrowprops=dict(facecolor='gray', shrink=0.05, color='gray'),
            color='gray')
plt.tight_layout()
plt.savefig('EKE_ZKE_spectra_illustrate.eps')
#plt.show()
