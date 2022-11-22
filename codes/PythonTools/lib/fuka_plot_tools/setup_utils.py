import sys, os
import matplotlib.pyplot as plt
pyFUKA_libspath = os.getenv('HOME_KADATH')+'/codes/PythonTools/lib/'
sys.path.append(pyFUKA_libspath)

# Set default matplotlib settings
plt.rcParams.update({
    'figure.figsize'    : [8.0, 8.0],
    'text.usetex'       : True,
    'font.family'       : "sans-serif",
    'font.serif'        : "cm",
    'xtick.major.size'  : 6,
    'xtick.minor.size'  : 3,
    'ytick.major.size'  : 6,
    'ytick.minor.size'  : 3,
})
LabelSize=15
TickSize=15
CbarLabelSize=15

def extract_path_filename(f):
  # extract filename and path
  idx = f.rfind('/')
  if idx != -1:
    path = f[0:idx+1]
    f = f[idx+1:]
  else:
    path = './'
  return path, f

def gen_cbarlabel(var_name,sq=False,inv=False,log=False):
  cbarlabel = var_name
  if sq:
    cbarlabel += r'^2'
  if inv:
    cbarlabel = r'1 / ' + cbarlabel
  if log:
    cbarlabel = r'\log \left( '+cbarlabel+r' \right)'
  cbarlabel = r'$'+cbarlabel+r'$'
  return cbarlabel

def check_ID_filename(f):
  path, f = extract_path_filename(f)
  
  # make sure we have the correct ext
  ext = f.split('.')[-1].lower()
  ispickle = (ext == 'pickle')
  
  if ext != 'dat' and not ispickle:
    ext = '.' + f.split('.')[-1]
    f = f.strip(ext) + '.dat'
  
  return path+f, ispickle

def get_reader(args, filepathabs):
  if args.bhns:
    from fukaID_readers.bhns import bhns_reader
    return bhns_reader(filepathabs)
  elif args.bbh:
    from fukaID_readers.bbh import bbh_reader
    return bbh_reader(filepathabs)
  elif args.bh:
    from fukaID_readers.bh import bh_reader
    return bh_reader(filepathabs)
  elif args.bns:
    from fukaID_readers.bns import bns_reader
    return bns_reader(filepathabs)
  elif args.ns:
    from fukaID_readers.ns import ns_reader
    return ns_reader(filepathabs)

def parse_var(var):
  inv = False
  sq =  False
  plotz=False
  if var[-2:] == "sq":
    sq = True
    var = var[:-2]
  if var[-3:] == "inv":
    inv = True
    var = var[:-3]
  if var[-2:] == "-z":
    plotz=True
    var = var[:-2]
  return inv, sq, plotz, var

name_dict = {
  'conf' : r'\Psi',
  'lapse' : r'\alpha',
  'shift' : r'\beta',
  'cPsi' : r'\mathcal{C}_{\Psi}',
  'cLapsePsi' : r'\mathcal{C}_{\alpha\Psi}',
  'cShift_1' : r'\mathcal{C}_{\beta}^1',
  'cShift_2' : r'\mathcal{C}_{\beta}^2',
  'cShift_3' : r'\mathcal{C}_{\beta}^3',
}

def var_name(var):
  if var in name_dict.keys():
    return name_dict[var]
  else: 
    return var

def extract_data(
  reader, 
  var, 
  x_coords,
  y_coords,
  plotz=False,
  logscale=False,
  square=False,
  inverse=False):
  
  import numpy as np

  xlen = len(x_coords)
  ylen = len(y_coords)
  
  if not plotz:
    coords_lst = [[x, y, 0] for y in y_coords for x in x_coords]
  else:
    coords_lst = [[x, 0, z] for z in y_coords for x in x_coords]
  
  # FIXME add excision filling extrapolation
  # if excision_extraction is not None:
    # if excision_extraction == "spherical":
    #   data = bh.getExFieldVal_sphere_extrap(
    #       reader, 
    #       coords_lst, 
    #       interpolation_order,
    #       interpolation_offset_r,
    #       relative_spacing_dr,
    #       -1)
    # elif excision_extraction == "oldaxial":
    #   data = bh.getExFieldVal_axial_extrap(
    #       reader, 
    #       coords_lst, 
    #       interpolation_order,
    #       interpolation_offset_r,
    #       relative_spacing_dr,
    #       -1)
    # elif excision_extraction == "axial":
    #   data = bh.getExFieldVal_spectral_axial_extrap(
    #       reader, 
    #       coords_lst, 
    #       interpolation_order,
    #       interpolation_offset_r,
    #       relative_spacing_dr,
    #       -1)

    # else:
  data = reader.getFieldValues(var, coords_lst, -1)
  data = np.array(data)
  
  if square:
    data *= data

  if inverse:
    data = 1. / data

  if logscale:
    data = np.array(np.log10(np.abs(data)+1e-15))
  
  if xlen == ylen:
    data=data.reshape(xlen,ylen)
  return data

def get_quiver_vars(var, plane):
  qvars = list()
  for c in plane:
    if c.lower() == "x" or c.lower() == '1':
      qvars.append(var[0]+"_1")
    elif c.lower() == "y" or c.lower() == '2':
      qvars.append(var[0]+"_2")
    elif c.lower() == "z" or c.lower() == '3':
      qvars.append(var[0]+"_3")
    else:
       ValueError("Quiver Components: Quiver only works with xyz components")
  if len(qvars) > 2:
    ValueError("Quiver Components: Quiver only works with 2D components")
  return qvars