#!/bin/bash

# make sure dedalus env is active
#------------------------------------
# params
nproc=4
Pr=1000

nz=32
nRaT=20
nkx=20
#------------------------------------
printf '############################################################################# \n' >> logfiles/rb_neutral_Pr_${Pr}.log
mpiexec -np ${nproc} python3 rb_stability.py --Pr=${Pr} --nRaT=${nRaT} --nkx=${nkx} --nz=${nz} --RaT_range='1000, 3000' --kx_range='2, 4' >> logfiles/rb_neutral_Pr_${Pr}.log
printf '############################################################################# \n' >> logfiles/rb_neutral_Pr_${Pr}.log
#------------------------------------


#------------------------------------

