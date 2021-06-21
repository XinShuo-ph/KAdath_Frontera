import numpy as np
from bin_props import *
from glob import glob
import os
import collections
import pickle
from multiprocessing import Pool
import itertools
import time

def interp_data(params):
  file, coords = params
  coords = coords.tolist()

  bns_id = bns_reader(file)
  rho = np.array(bns_id.getEOSValues(coords))[:,0]

  return rho

data_dir = "./"
data_dir = "./tog_ecc_reuc/total_17/"

data = ["tog_highq_irrot",
        "tog_highq_spin"]

data = ["spinning.H.20.0.0.3.q1",
        "spinning.H.20.0.45.0.45.3.q1"]

data = ["spinning.H.30.0.0.3.73461.q2.2",
        "spinning.H.30.0.0.6.3.62791.q2.2"]

files = [data_dir + d + ".dat" for d in data]

num_points = 5000
num_procs = 40

x_mid = 0
x_delta = 30
y_delta = 11
z_delta = 11

x_range = np.linspace(x_mid - x_delta, x_mid + x_delta, num=num_points)
y_range = np.linspace(-y_delta, y_delta, num=num_points)
z_range = np.linspace(-z_delta, z_delta, num=num_points)

ranges = [x_range,y_range,z_range]
np.savez_compressed('data_ranges.npz', *ranges)

coords_lst_xy = [[x, y, 0] for x in x_range for y in y_range]
coords_lst_xz = [[x, 0, z] for x in x_range for z in z_range]

coords_lst_xy_split = np.array_split(np.array(coords_lst_xy), num_procs);
coords_lst_xz_split = np.array_split(np.array(coords_lst_xz), num_procs);

data_xy = []
data_xz = []

p = Pool(num_procs)

print "Starting interpolation..."

for file in files:
  input_xy = [[file,coords_xy] for coords_xy in coords_lst_xy_split]
  input_xz = [[file,coords_xz] for coords_xz in coords_lst_xz_split]

  print file,"xy"
  start = time.time()
  data_xy.append(np.concatenate(p.map(interp_data, input_xy)))
  end = time.time()
  print "This took ", end - start, " s"

  print file, "xz"
  start = time.time()
  data_xz.append(np.concatenate(p.map(interp_data, input_xz)))
  end = time.time()
  print "This took ", end - start, " s"

data_grids_xy = [dset.reshape(len(x_range),len(y_range)).swapaxes(0,1) for dset in data_xy]
data_grids_xz = [dset.reshape(len(x_range),len(z_range)).swapaxes(0,1) for dset in data_xz]

np.savez_compressed('data_xy.npz', *data_grids_xy)
np.savez_compressed('data_xz.npz', *data_grids_xz)