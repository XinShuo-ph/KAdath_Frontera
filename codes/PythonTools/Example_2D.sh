f=${HOME_KADATH}/codes/PythonTools/Example_id/converged_BHNS_ECC_RED.togashi.35.0.6.0.52.3.6.q0.487603.0.1.11.dat
qf=$f
./plot_fukaid_2D.py \
--bhns \
-f $f \
--vars cPsi --extent -45 45 0 45 \
--log --vmin -10 --vmax -2 --npts 256 \
--qvar shift --qpts 64 --qscale 0.05 \
-qf $qf --cmap viridis --qcmap coolwarm \
--cbar --landscape

