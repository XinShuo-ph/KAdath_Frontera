import numpy as np, sys
import matplotlib.pyplot as plt
import matplotlib
from kadath_readers import *
sys.path.append('/home/user/lib/codes/Python_Codes/mylibs')
import analysis_tools as at

id_file='bbh.dat'

bbh_id = bbh_reader(id_file)
bconfig = at.get_bconfig(id_file[0:len(id_file)-4]+".info")
gomega = float(bconfig['binary']['global_omega'])

print("Constructing cartesian grid...")

num_points = 20000

x_mid = float(bbh_id.vars['xc2'])
delta = 0.5
print(x_mid)

minusdelta = delta
plusdelta = 0
x_range = np.linspace(x_mid - minusdelta, x_mid+plusdelta, num=num_points)
y_range = np.linspace(-delta, delta, num=num_points)
z_range = np.linspace(-delta, delta, num=num_points)

coords_lst = [[x, 0, 0] for x in x_range]
#coords_lst = [[x_mid, y, 0] for y in y_range]
#coords_lst = [[x_mid, 0, z] for z in z_range]

print("Getting data at {} points...".format(num_points))
quants = ['cP','cNP']

data = {q:np.array(bbh_id.getFieldValues(q, coords_lst, -1)) for q in quants}

plt.clf()
for j,q in enumerate(data):
  plt.plot(x_range, data[q], label=q)

def add_lines(x_mid, shells):
    global bbh_id
    
    plt.axvline(x_mid, color='red', linestyle='--')
    for i in range(1,shells+1):
        plt.axvline(x_mid+bbh_id.vars['BH2_R'+str(i)], color='black')
        plt.axvline(x_mid-bbh_id.vars['BH2_R'+str(i)], color='black')

draw_lines=True
if draw_lines == True:
  add_lines(x_mid, int(bconfig['bh2']['nshells'])+3)

plt.xlim((x_mid-minusdelta)*0.98, (x_mid+plusdelta)*1.02)
plt.legend()
plt.show()
