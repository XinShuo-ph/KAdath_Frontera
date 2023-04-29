f=/home/tootle/lib/kadath/codes/FUKAv2_Solvers/BHNS/mass_ratio/BHNS_ECC_RED.togashi.134.4.0.6.0.8.22.4.q0.0666667.0.0.13.dat
qf=$f
./plot_fukaid_2D.py \
--bhns \
-f $f \
--vars cPsi --extent -82.5 -42.5 -20 20 \
--log --vmin -10 --vmax -2 --npts 256 \
--qvar shift --qpts 64 --qscale 0.05 \
-qf $qf --cmap viridis --qcmap coolwarm \
--cbar

