#!/usr/bin/python3 -u

import os
import time

def run_cmd(cmd):
  print(cmd)
  os.system(cmd)
  time.sleep(3)

root = os.getcwd()
# main_dirs = ['typeI', 'typeII', 'typeIII', 'typeIV_a', 'typeIV_b', 'typeV', 'c_s_eta', 'JJ_pion']
main_dirs = ['typeI', 'typeII', 'typeIV_a', 'typeIV_b', 'typeV', 'c_s_eta', 'JJ_pion']

trajs = [2000]

for d in main_dirs:
  os.chdir(os.path.join(root, d))
  print(os.getcwd())

  for traj in trajs:
    cmd = "python submit.py " + str(traj)
    run_cmd(cmd)

