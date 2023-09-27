#!/bin/bash

# make sure dedalus env is active
#------------------------------------
# params
nproc=4
Pr=1

nz=32
nRaT=4
nkx=20
#------------------------------------
printf '############################################################################# \n' >> logfiles/rb_sigma_vs_k_Pr_${Pr}.log
python3 rb_sigma_vs_k.py --Pr=${Pr} --nRaT=${nRaT} --nkx=${nkx} --nz=${nz} --RaT_range='1700, 1712' --kx_range='2.5, 3.5' >> logfiles/rb_sigma_vs_k_Pr_${Pr}.log
printf '############################################################################# \n' >> logfiles/rb_sigma_vs_k_Pr_${Pr}.log
#------------------------------------


#------------------------------------

