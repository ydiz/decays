import os

root = os.getcwd()
main_dirs = ['typeI', 'typeII', 'typeIII', 'typeIV_a', 'typeIV_b', 'typeV', 'c_s_eta', 'JJ_pion']

for d in main_dirs:
  os.chdir(os.path.join(root, d))
  print(os.getcwd())
  cmd = "python submit_all_diagrams.py"
  print(cmd)
  # os.sytem(cmd)

