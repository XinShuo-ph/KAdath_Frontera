import numpy as np, sys
import matplotlib.pyplot as plt
import matplotlib
from itertools import cycle
from kadath_readers import *
sys.path.append('/home/user/lib/codes/Python_Codes/mylibs')
import analysis_tools as at

def rel_diff(a, b):
  return 2.0 * np.abs(a - b) / (np.abs(a) + np.abs(b))

data_dir = "./"

path = data_dir

files = [
#'bbh.dat',
'converged_BBH_TOTAL_BC.20.0.5.0.1.q0.0526316.09.dat',
#'converged_BBH_TOTAL_BC.20.0.5.0.1.q0.0526316.09-fixed_omega.dat',
'converged_BBH_TOTAL_BC.20.0.5.0.1.q0.0526316.09-oldbounds.dat',
#'converged_BBH_TOTAL_BC.20.0.5.0.1.q0.0526316.11.dat',
#'converged_BBH_TOTAL_BC.20.0.5.0.1.q0.0526316.13.dat',
]

bbh_ids = [bbh_reader(path + f) for f in files]
bconfigs = [at.get_bconfig(path+f[0:len(f)-4]+".info") for f in files]
gomega = float(bconfigs[len(bconfigs)-1]['binary']['global_omega'])
for conf in bconfigs:
#    print(conf['binary']['global_omega'])
    print(rel_diff(float(conf['binary']['global_omega']), gomega))

#determine most number of shells to draw in plot
draw_lines = False
r_id = 0
cnt_max = 0
for i in range(len(bbh_ids)):
    cid = bbh_ids[i]
    rcnt = 0
    for key in cid.vars.keys():
        if key[0:5] == 'BH2_R':
            rcnt += 1
    if rcnt > cnt_max:
        cnt_max = rcnt
        r_id = i

num_points = 10000

x_mid = float(bbh_ids[0].vars['xc2'])
delta = 5
if draw_lines == True:
    x_mid = 0
    delta = 20
#x_mid = 0
print(x_mid)

x_range = np.linspace(x_mid - delta, x_mid + delta, num=num_points)
y_range = np.linspace(-delta, delta, num=num_points)
z_range = np.linspace(-delta, delta, num=num_points)

coords_lst = [[x, 0, 0] for x in x_range]
#coords_lst = [[x_mid, y, 0] for y in y_range]
#coords_lst = [[x_mid, 0, z] for z in z_range]

#print "Getting data at {} points...".format(num_points)
quants = [
  'cP',
#  'cNP'
]
data = [{q:np.array(bbh_id.getFieldValues(q, coords_lst, -1)) for q in quants} for i,bbh_id in enumerate(bbh_ids)]

print("Plotting...")
maxv=0

plt.clf()
linestyles = cycle(['-','--','-.'])
for i,d in enumerate(data):
  ls = linestyles.__next__() 
  for j,q in enumerate(d):
    plt.plot(x_range, d[q], label="{}_{}".format(q,i), linestyle=ls)

def add_lines(x_mid, idx, shells):
    global bbh_ids
    plt.axvline(x_mid, color='red', linestyle='--')
    for i in range(1,shells+1):
        plt.axvline(x_mid+bbh_ids[idx].vars['BH2_R'+str(i)], color='black')
        plt.axvline(x_mid-bbh_ids[idx].vars['BH2_R'+str(i)], color='black')

draw_lines=True
if draw_lines:
    add_lines(x_mid, r_id, cnt_max)
plt.legend()
plt.show()
