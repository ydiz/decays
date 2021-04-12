#!/usr/bin/python3 -u
import re
import numpy as np
import pickle

# traj = 2000
# trajs = range(2100, 2260, 10)
trajs = range(2000, 2100, 10)

##################################################################

T = 64
main_diagram_tseps = [6, 8, 10, 12, 14]

def get_array(text, keyword, var_name):
  rst = re.findall(keyword + r': \[(.*)\]', text)
  assert len(rst) == 1 and len(rst[0]) > 0, f"keyword: {keyword}. Not found."
  rst = rst[0]
  rst = re.findall(r'[\d\-.e+]+', rst)
  rst = np.array(rst, dtype=float)
  rst = rst[::2] + 1j * rst[1::2]

  if keyword != "JJ_pion":
    rst = rst.reshape((-1, 64))

  if var_name.startswith(('type', 'sBar_d')):  # for the diagrams to calculate <JJ Hw K> or <JJ sBar d K>
    if rst.shape[0] == T:     # convert table2d from (T, T) to (tsep, T) so that all tables have the same shape
      indices = [-x for x in main_diagram_tseps]
      rst = rst[indices]

  return rst

all_data = {
    'typeI/output/typeI_D1a.out': [('table2d_Q1', 'typeI_D1aQ1'), ('table2d_Q2', 'typeI_D1aQ2')], # filename -> (keyword, variable name)
    'typeI/output/typeI_D1b.out': [('table2d_Q1', 'typeI_D1bQ1'), ('table2d_Q2', 'typeI_D1bQ2')], 
    'typeI/output/typeI_D2a.out': [('table2d_Q1', 'typeI_D2aQ1'), ('table2d_Q2', 'typeI_D2aQ2')],
    'typeI/output/typeI_D2b.out': [('table2d_Q1', 'typeI_D2bQ1'), ('table2d_Q2', 'typeI_D2bQ2')],

    'typeII/output/typeII_D1a.out': [('table2d_Q1', 'typeII_D1aQ1'), ('table2d_Q2', 'typeII_D1aQ2'), ('table2d_sBar_d T1D1', 'sBar_d_T1D1a')],
    'typeII/output/typeII_D1b.out': [('table2d_Q1', 'typeII_D1bQ1'), ('table2d_Q2', 'typeII_D1bQ2'), ('table2d_sBar_d T1D1', 'sBar_d_T1D1b')],

    'typeIII/output/typeIII_D1a.out': [('table2d_Q1', 'typeIII_D1aQ1'), ('table2d_Q2', 'typeIII_D1aQ2')],
    'typeIII/output/typeIII_D1b.out': [('table2d_Q1', 'typeIII_D1bQ1'), ('table2d_Q2', 'typeIII_D1bQ2')],

    'typeIV_a/output/typeIV_D1a.out': [('table2d_Q1', 'typeIV_D1aQ1'), ('table2d_Q2', 'typeIV_D1aQ2'), ('table2d sBar_d T2 diagram a', 'sBar_d_T2D1a')], 
    'typeIV_a/output/typeIV_D2a.out': [('table2d_Q1', 'typeIV_D2aQ1'), ('table2d_Q2', 'typeIV_D2aQ2'), ('table2d sBar_d T2 diagram a', 'sBar_d_T2D2a')],  
    'typeIV_b/output/typeIV_D1b.out': [('table2d_Q1', 'typeIV_D1bQ1'), ('table2d_Q2', 'typeIV_D1bQ2'), ('table2d sBar_d T2 diagram b', 'sBar_d_T2D1b')],  
    'typeIV_b/output/typeIV_D2b.out': [('table2d_Q1', 'typeIV_D2bQ1'), ('table2d_Q2', 'typeIV_D2bQ2'), ('table2d sBar_d T2 diagram b', 'sBar_d_T2D2b')],  

    'typeV/output/typeV.out': [('table2d_D1Q1', 'typeV_D1Q1'), ('table2d_D1Q2', 'typeV_D1Q2'), ('table2d_D2Q1', 'typeV_D2Q1'), ('table2d_D2Q2', 'typeV_D2Q2'), ('table2d sBar_d_D1', 'sBar_d_T3D1'), ('table2d sBar_d_D2', 'sBar_d_T3D2')],  

    'JJ_pion/output/JJ_pi.out': [('JJ_pion', 'JJ_pion')],

    'c_s_eta/output/c_s_eta.out': [('table2d sBar_d_T1D1', 'cs_sBar_d_T1D1'), ('table2d sBar_d_T1D2', 'cs_sBar_d_T1D2'), ('table2d sBar_d_T2D1', 'cs_sBar_d_T2D1'), ('table2d sBar_d_T2D2', 'cs_sBar_d_T2D2'), ('table2d Hw_T1D1Q1', 'cs_Hw_T1D1Q1'), ('table2d Hw_T1D1Q2', 'cs_Hw_T1D1Q2'), ('table2d Hw_T2D1Q1', 'cs_Hw_T2D1Q1'), ('table2d Hw_T2D1Q2', 'cs_Hw_T2D1Q2'), ('table2d Hw_T2D2Q1', 'cs_Hw_T2D2Q1'), ('table2d Hw_T2D2Q2', 'cs_Hw_T2D2Q2'), ('table2d Hw_T3D1Q1', 'cs_Hw_T3D1Q1'), ('table2d Hw_T3D1Q2', 'cs_Hw_T3D1Q2'), ('table2d Hw_T3D2Q1', 'cs_Hw_T3D2Q1'), ('table2d Hw_T3D2Q2', 'cs_Hw_T3D2Q2')],

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

  # print(typeI_D1aQ2)
  # print(typeII_D1aQ2)
  # print(sBar_d_T1D1a)
  # print(sBar_d_T2D1a)  # ?? Some sites  very large; 1e15 
  # print(sBar_d_T2D2b.shape)  
  # print(sBar_d_T3D2.shape)  
  # print(JJ_pion)  
  # print(cs_Hw_T2D1Q2.shape)  


  arrays = [pair[1] for l in all_data.values() for pair in l]
  print(arrays)
  rst = {array: eval(array) for array in arrays}  # map from variable name to variable

  pickle.dump(rst, open(f'all_results/{traj}.pkl', 'wb'))

