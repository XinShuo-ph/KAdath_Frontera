import numpy as np
import matplotlib.pyplot as plt
import matplotlib, numpy as np
from kadath_readers import *

def rel_diff(a, b):
  return 2.0 * np.abs(a - b) / (np.abs(a) + np.abs(b))

#data_dir = "../bin/Release/data/"
data_dir = "./"

print("Reading BBH ID...")

path = data_dir

files = [#"bns.9.45.0.6.0.6.3.q1.dat",
#         "sly_q2_py.dat",
#         "sly_q2_py_tot_1.dat"]
#         "sly_q2_nopy.dat"]
#         "spinning.H.33.8524.0.0.520048.4.15861.q2.03878.dat"
#'converged_BNS_TOTAL.33.8524.-0.45.0.45.2.74.q1.09.dat',
#'converged_BNS_TOTAL.33.8524.-0.4.0.4.2.74.q1.13.dat'
#'initbin.dat'
'converged_BBH_TOTAL_BC.20.0.0.2.q1.09.dat'
#'converged_BNS_ECC_RED.33.8524.-0.4.0.4.2.74.q1.13.dat'
#'bns.dat',
]

bns_ids = [bbh_reader(path + f) for f in files]
#bns_id = bns_reader()

print("Constructing cartesian grid...")

num_points = 10000
dist = 200
x_mid = dist /2.
delta = 50

x_range = np.linspace(x_mid, x_mid + delta, num=num_points)
y_range = np.linspace(- delta, delta, num=num_points)
z_range = np.linspace(- delta, delta, num=num_points)

coords_lst = [[x, 0, 0] for x in x_range]
#coords_lst = [[x_mid, y, 0] for y in y_range]
#coords_lst = [[x_mid, 0, z] for z in z_range]

print("Getting data at {} points...".format(num_points))
#quants = ['h','eps','rho','press']
#quants = ['rho']
#quants = ['delta']
#quants = ['h']
#quants = ['drhodx']
#quants = ['dHdx']
#quants = ['delta']
#quants = ['cP','cNP'] #,'c2','c3','c4']
#quants = ['c0','c1','c3']
#quants = ['c5','c6']
quants = ['conf']

data = [{q:np.array(bns_id.getFieldValues(q, coords_lst, -1)) for q in quants} for i,bns_id in enumerate(bns_ids)]
central_data = [{q:bns_id.getFieldValues(q, [[x_mid,0,0]], -1) for q in quants} for i,bns_id in enumerate(bns_ids)]

print("Plotting...")

for i,d in enumerate(data):
#  h = 1 + d['eps'] + d['press'] / d['rho']
#  err = rel_diff(h,d['h'])
#  plt.plot(x_range, err, label=files[i])

  for q in d:
    plt.plot(x_range, d[q], label=files[i] + " " + q)
#    plt.plot(x_range, np.log10(np.abs(d[q])), label=files[i] + " " + q)
#    plt.plot(x_range, np.log10(np.abs(d[q] - central_data[i][q])), label=files[i] + " " + q)
for k in bns_ids[0].vars.keys():
  if 'BH_R' in k:
    vplus = x_mid+float(bns_ids[0].vars[k])
    vminus = x_mid-float(bns_ids[0].vars[k])
    print(vplus, vminus)
    print(x_mid + delta, x_mid - delta)
    if vplus < x_mid + delta:
      plt.vlines(vplus,-5, 5)
    if vminus > x_mid - delta:
      plt.vlines(vminus,-5, 5)
plt.legend()

plt.show()
