#!/usr/bin/python3 -u
import re
import numpy as np
import pickle

# trajs = sorted(set(range(1920, 2260, 10)) - set([1020, 1060, 1100, 1340, 1620, 1830, 1880, 1920, 1980]))   # [1620, 1830, 1920, 1980] have one time slice missing
# trajs = sorted(set(range(1500, 1550, 10)) - set([1020, 1060, 1100, 1340, 1620, 1830, 1880, 1920, 1980]))   # [1620, 1830, 1920, 1980] have one time slice missing
# trajs = sorted(set(range(1010, 1200, 10)) - set([1020, 1060, 1100, 1340, 1620, 1830, 1880, 1920, 1980]))   # [1620, 1830, 1920, 1980] have one time slice missing
# trajs = sorted(set(range(1010, 1550, 10)) - set([1090, 1130, 1300, 1430]) - set([1020, 1060, 1100, 1340, 1620, 1830, 1880, 1920, 1980]))   # [1620, 1830, 1920, 1980] have one time slice missing
# trajs = sorted(set(range(1010, 2260, 10)) - set([1020, 1060, 1100, 1340, 1620, 1830, 1880, 1920, 1980]))   # [1620, 1830, 1920, 1980] have one time slice missing
# trajs = sorted(set(range(1870, 2260, 10)) - set([1020, 1060, 1100, 1340, 1620, 1830, 1880, 1920, 1980]))   # [1620, 1830, 1920, 1980] have one time slice missing
trajs = sorted((set(range(1410, 1710, 10)) | set(range(1870, 2260, 10))) - set([1520, 1020, 1060, 1100, 1340, 1620, 1830, 1880, 1920, 1980]))   # [1620, 1830, 1920, 1980] have one time slice missing
# trajs = [1450, 1460]
print(trajs)
# assert 0

##################################################################

# for typeI and typeII, tseps = [6, 8, 10, 12, 14]
# for type III, type IV, type V, I calculated all tseps
# for cs diagrams, tseps = [12, 14, 16, 18]

T = 64
tsep_max = 20  # for type III, IV and V, only keep  - tsep_max <= tK <= 0

def get_array(text, keyword, var_name):
  rst = re.findall(keyword + r'.*: \[(.*)\]', text)
  assert len(rst) == 1 and len(rst[0]) > 0, f"keyword: {keyword}. Not found."
  rst = rst[0]
  rst = re.findall(r'[\d\-.e+]+', rst)
  rst = np.array(rst, dtype=float)
  rst = rst[::2] + 1j * rst[1::2]
  if var_name.startswith('cs_'):
    rst = rst.reshape((4, 19, 21))  # (N_tseps, NUM_xt, num_R)
  else:
    rst = rst.reshape((-1, 64))

  # if var_name.startswith(('type', 'sBar_d')) and rst.shape[0] == T:  # for the diagrams to calculate <JJ Hw K> (type III, type IV, type V) and <JJ sBar d K> # convert table2d from (T, T) to (tsep_max, T) to save disk space 
  #     indices = [-x for x in range(tsep_max)]
  #     rst = rst[indices]

  return rst

all_data = {
    'control_R_new_strategy/typeIII/output_new_strategy/typeIII_D1a.out': [('table2d_Q1', 'typeIII_D1aQ1'), ('table2d_Q2', 'typeIII_D1aQ2')],  
    'control_R_new_strategy/typeIII/output_new_strategy/typeIII_D1b.out': [('table2d_Q1', 'typeIII_D1bQ1'), ('table2d_Q2', 'typeIII_D1bQ2')],  
    'control_R_new_strategy/typeIV_a/output_new_strategy/typeIV_D1a.out': [('table2d_Q1', 'typeIV_D1aQ1'), ('table2d_Q2', 'typeIV_D1aQ2'), ('table2d sBar_d T2 diagram a', 'sBar_d_T2D1a')],  

    # 'control_R_new_strategy/typeIV_b/output_new_strategy/typeIV_D1b.out': [('table2d_Q1', 'typeIV_D1bQ1'), ('table2d_Q2', 'typeIV_D1bQ2')],  # , ('table2d sBar_d T2 diagram b', 'sBar_d_T2D1b')

    'control_R_new_strategy/typeV/output_new_strategy/typeV.out': [('table2d_D1Q1', 'typeV_D1Q1'), ('table2d_D1Q2', 'typeV_D1Q2'), ('table2d_D2Q1', 'typeV_D2Q1'), ('table2d_D2Q2', 'typeV_D2Q2'), ('table2d sBar_d_D1', 'sBar_d_T3D1'), ('table2d sBar_d_D2', 'sBar_d_T3D2') ], 

    'control_R_new_strategy/c_s_eta/output_new_strategy/c_s_eta.out': [('table3d sBar_d_T1D1', 'cs_sBar_d_T1D1'), ('table3d sBar_d_T1D2', 'cs_sBar_d_T1D2'), ('table3d sBar_d_T2D1', 'cs_sBar_d_T2D1'), ('table3d sBar_d_T2D2', 'cs_sBar_d_T2D2'), ('table3d Hw_T1D1Q1', 'cs_Hw_T1D1Q1'), ('table3d Hw_T1D1Q2', 'cs_Hw_T1D1Q2'), ('table3d Hw_T2D1Q1', 'cs_Hw_T2D1Q1'), ('table3d Hw_T2D1Q2', 'cs_Hw_T2D1Q2'), ('table3d Hw_T2D2Q1', 'cs_Hw_T2D2Q1'), ('table3d Hw_T2D2Q2', 'cs_Hw_T2D2Q2'), ('table3d Hw_T3D1Q1', 'cs_Hw_T3D1Q1'), ('table3d Hw_T3D1Q2', 'cs_Hw_T3D1Q2'), ('table3d Hw_T3D2Q1', 'cs_Hw_T3D2Q1'), ('table3d Hw_T3D2Q2', 'cs_Hw_T3D2Q2')],
    }


for traj in trajs:
  for file_prefix in all_data:
    fname = file_prefix + '.' + str(traj) # 'typeI/output/typeI_D1a.out.2000'
    # print(fname)
    with open(fname) as f:
      fcontent = f.read()
      for keyword, var_name in all_data[file_prefix]:
        print(traj, var_name)
        exec(f"{var_name} = get_array(fcontent, keyword, var_name)")

  # print(cs_Hw_T2D1Q2.shape)  

  arrays = [pair[1] for l in all_data.values() for pair in l]
  print(arrays)
  rst = {array: eval(array) for array in arrays}  # map from variable name to variable

  pickle.dump(rst, open(f'all_results_control_R/{traj}.pkl', 'wb'))
  # pickle.dump(rst, open(f'all_results_typeV_control_R/{traj}.pkl', 'wb'))

