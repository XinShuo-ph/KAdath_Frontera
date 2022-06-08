import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from bin_props import *
from glob import glob
import os
import collections

def makehash():
    return collections.defaultdict(makehash)

def rel_diff(a, b):
  return 2.0 * np.abs(a - b) / (np.abs(a) + np.abs(b))

#data_dir = "../bin/Release/data/"
data_dir = "./"

print "Reading BNS ID..."

path = "seq"
dsets = makehash()

#for res_dir in glob(os.path.join(path,'*')):
for res_dir in glob(os.path.join(path,'11')):
  res = os.path.basename(res_dir)

  for config_dir in glob(os.path.join(res_dir,'*')):
    config = os.path.basename(config_dir)

    for dset in glob(os.path.join(config_dir,'spinning*.dat')):
      print "reading " + dset

      bns_id = bns_reader(dset)
      dsets[config][res]["qe"][bns_id.vars["dist"]] = bns_id

    for dset in glob(os.path.join(config_dir,'ecc*.dat')):
      print "reading " + dset

      bns_id = bns_reader(dset)
      dsets[config][res]["pn"][bns_id.vars["dist"]] = bns_id

def irrot_PN(nu,x):
  EboM = 1 \
       + (-3. / 4. - 1. / 12. * nu) * x \
       + (- 27. / 8. + 19. / 8. * nu - 1. / 24. * nu**2) * x**2 \
       + (- 675. / 64. + (34445. / 576. - 205. / 96. * np.pi**2) * nu - 155. / 96. * nu**2 - 35. / 5184. * nu**3) * x**3 \
       + (- 3969. / 128. + 448. / 15. * nu * np.log(x) \
        + (- 123671. / 5760. + 9037. / 1536. * np.pi**2 + 1792. / 15. * np.log(2.) + 896. / 15. * np.e) * nu \
        + (- 498449. / 3456. + 3157. / 576. * np.pi**2) * nu**2 \
        + 301. / 1728. * nu**3 \
        + 77. / 31104. * nu**4) * x**4

  return - nu * x / 2. * EboM


#fig, axs = plt.subplots(3, 2, sharex='col', sharey='row')
fig, axs = plt.subplots(2, 2, sharex='col', sharey='row')

for i,cfg in enumerate(dsets):
  axs[i, 0].set_title('quasi-equilibrium ' + cfg)
  axs[i, 1].set_title('post-newtonian ' + cfg)

  for res in dsets[cfg]:
    for j,t in enumerate(dsets[cfg][res]):
      mO = []
      diff = []

      for dist in sorted(dsets[cfg][res][t]):
        bns_id = dsets[cfg][res][t][dist]

        M1 = bns_id.vars["Madm1"]
        M2 = bns_id.vars["Madm2"]
        Minf = bns_id.vars["Minf"]

        mOmega = bns_id.vars["mOmega"]

        EboM = bns_id.vars["EboM"]

        nu = M1 * M2 / Minf**2
        x = np.power(mOmega, 2./3.)
        EboM_pn = irrot_PN(nu,x)

        mO.append(mOmega)
        diff.append(rel_diff(EboM_pn, EboM))
        #diff.append(EboM)

      axs[i, j].scatter(mO, diff, label=res)
      axs[i, j].set_yscale('log')
      axs[i, j].set_ylim([1e-5, 1])
      axs[i, j].legend()

#plt.legend()
#plt.show()
plt.savefig("eb.pdf", bbox_inches="tight", dpi=300)