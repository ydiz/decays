#!/usr/bin/python3 -u

import os
import time

def run_cmd(cmd):
  print(cmd)
  os.system(cmd)
  # time.sleep(3)

root = os.getcwd()
dirs = ['typeI', 'typeII', 'typeIII', 'typeIV_a', 'typeIV_b', 'typeV', 'c_s_eta', 'JJ_pion']

for d in dirs:
  os.chdir(os.path.join(root, d))
  print(os.getcwd())

  cmd = "nohup make &"
  run_cmd(cmd)

