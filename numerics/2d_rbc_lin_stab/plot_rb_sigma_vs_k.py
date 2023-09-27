"""
Usage:
  plot_sigma_vs_k.py [--Pr=<prandtl> --nRaT=<nRaT> --nkx=<nkx> --nz=<nz>]

Options:
  --Pr=<prandtl>            Prandtl number [default: 1]
  --nz=<nz>                 Chebyshev discretization in z [default: 32]
"""


"""
read data and plot neutral stability curves 

"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from pylab import *
import h5py
from docopt import docopt

##########################################
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 30}
#matplotlib.rc('font', **font)
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
mpl.rcParams.update({'font.size': 20})

font = {'family' : 'monospace',
        'weight' : 'bold',
        'size'   : 20}

from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Lucida Grande']})
# for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
mpl.rcParams.update({'font.size': 20})
##########################################
logger = logging.getLogger(__name__)
args=docopt(__doc__)

Pr = float(args['--Pr'])
nz = int(args['--nz'])

RaTskip = 1 # for plotting, plot sigma vs k from Repoints with RaTskip gap in between
spurious = -999

root_name = "rb_sigma_vs_kx_Pr_{:.3e}".format(Pr) 
sigma_vs_k_file = 'results/'+root_name+'.h5'
hf = h5py.File(sigma_vs_k_file, 'r')
RaTpoints = hf.get("RaTpoints")
sigma_vs_k = hf.get("sigma_vs_k")


RaTpoints = np.array(RaTpoints)
sigma_vs_k = np.array(sigma_vs_k)

nRaT = len(RaTpoints)
nk = len(sigma_vs_k[0, :, 0])

# sigma_vs_k_filtered = np.zeros((nRe, nk, 3))
# remove spurious eigenvals
cutoff = 1e2

for i in range(0, RaTskip, nRaT):
        for j in range(0, nk):
                if (abs(sigma_vs_k[i, j, 1]) > cutoff):
                        sigma_vs_k[i, j, 1] = spurious


fig = plt.figure(0)
ax = fig.gca()

for i in range(0, len(RaTpoints)):
        RaT = RaTpoints[i]
        kx = sigma_vs_k[i, :, 0]
        sigma_r = sigma_vs_k[i, :, 1]
        idx = []
        for j in range(0, nk):
                if(sigma_r[j] == spurious):
                        idx.append(j)
        kx = np.delete(kx, [idx])
        sigma_r = np.delete(sigma_r, [idx])

        ax.plot(kx, sigma_r, linewidth=2, label=r'$Ra_T = \,$' + "{:.3e}".format(RaT))
        ax.axhline(y=0, linewidth = 2, color='k')
# ax.set_ylim([-0.001, 0.001])
# ax.set_xlim([0, 7])
ax.legend(loc=4)
ax.set_xlabel(r'$k_{x}$',fontsize=24)
ax.set_ylabel(r'$\sigma_r$',fontsize=24)
plt.tight_layout()
fig.savefig("figs/" + root_name + ".png")
##########################################