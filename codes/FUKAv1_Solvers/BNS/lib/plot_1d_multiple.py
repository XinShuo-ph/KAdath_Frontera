import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from kadath_readers import *

def rel_diff(a, b):
  return 2.0 * np.abs(a - b) / (np.abs(a) + np.abs(b))

#data_dir = "../bin/Release/data/"
data_dir = "./"

print("Reading BNS ID...")

path = data_dir

files = [#"bns.9.45.0.6.0.6.3.q1.dat",
#         "sly_q2_py.dat",
#         "sly_q2_py_tot_1.dat"]
#         "sly_q2_nopy.dat"]
#         "spinning.H.33.8524.0.0.520048.4.15861.q2.03878.dat"
'bns.dat']

bns_ids = [bns_reader(path + f) for f in files]

print("Constructing cartesian grid...")

num_points = 10000

x_mid = 16.
delta = 15

x_range = np.linspace(x_mid - delta, x_mid + delta, num=num_points)
y_range = np.linspace(- delta, delta, num=num_points)
z_range = np.linspace(- delta, delta, num=num_points)

coords_lst = [[x, 0, 0] for x in x_range]
#coords_lst = [[x_mid, y, 0] for y in y_range]
#coords_lst = [[x_mid, 0, z] for z in z_range]

print("Getting data at {} points...".format(num_points))
#quants = ['h','eps','rho','press']
quants = ['rho']
quants = ['h']
#quants = ['drhodx']
#quants = ['dHdx']
#quants = ['delta']
#quants = ['c0','c1','c2','c3','c4']
#quants = ['c0','c1','c3']
#quants = ['c5','c6']
#quants = ['firstint']

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

plt.legend()

plt.show()
