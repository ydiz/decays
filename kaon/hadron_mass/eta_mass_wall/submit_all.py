#!/usr/bin/python3 -u

import os
import time
import re

def run_cmd(cmd):
  print(cmd)
  os.system(cmd)
  time.sleep(3)

program = "eta_mass_wall"

for traj_start in range(1000, 2260, 100):
  print(traj_start)
  with open("submit.sh", 'r+') as f:
    content = f.read()
    content = re.sub(r'#SBATCH -J (.*)', '#SBATCH -J ' + f'{program}.{traj_start}', content)
    content = re.sub(r'program=(.+)', f'program={program}', content)
    content = re.sub(r'traj_start=(.+)', f'traj_start={traj_start}', content)
    f.truncate(0)   # erase all content of the file
    f.seek(0)       # truncate does not move file pointer
    f.write(content)

  time.sleep(1)

  run_cmd("sbatch submit.sh")
