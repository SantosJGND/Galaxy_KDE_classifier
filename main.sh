#!/bin/bash

python -u Kernel_mPLib3D_FM36_Galaxy.py \
--CHR 1 \
--geno Simulation_test/test_.05_.2/test.geno \
--fam Simulation_test/test_.05_.2/test.fam \
--bim Simulation_test/test_.05_.2/test.bim \
--ref Simulation_test/refs_sim.txt \
--admx Simulation_test/admx_sim.txt \
--out test_.05_.2 \
-w 20 \
--overlap 15 \
--proc 5 \
--dr PCA \
--MSprint