./plot_fukaid_2D.py \
--bhns \
-f Example_id/converged_BHNS_ECC_RED.togashi.35.0.6.0.52.3.6.q0.487603.0.1.11.dat \
--vars cPsi cLapsePsi --extent -50 50 -50 50 \
--log --vmin -10 --vmax -2 --npts 256 \
--qvar shift --qpts 64 --qscale 0.05 \
-qf Example_id/converged_BHNS_ECC_RED.togashi.35.0.6.0.52.3.6.q0.487603.0.1.11.dat --cmap viridis --qcmap coolwarm \
--cbar

