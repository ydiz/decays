#!/usr/bin/python3 -u

import re
import os
import time
import glob
import math

def run_cmd(cmd):
  print(cmd)
  os.system(cmd)
  time.sleep(3)


trajs = [2020]
# trajs = range(2010, 2070, 10)

#################################

for f in glob.glob('slurm-*.out'):
  os.remove(f)

prefix = "/global/cfs/cdirs/mp13/ydzhao/24ID/sequential"
TOTAL_NUM_POINTS = 512
POINTS_PER_JOB = 100  # number of propagators to calculate in one job
point_l_path = "/global/cscratch1/sd/ydzhao/point_l"

for traj in trajs:
  if not os.path.isdir(os.path.join(point_l_path, str(traj))):
    print(f'Trajectory {traj} does not have point_l; skipping.')

  with open("submit.sh", 'r+') as f:
    content = f.read()
    content = re.sub(r'traj=(.+)', 'traj='+str(traj), content)
    f.truncate(0)   # erase all content of the file
    f.seek(0)       # truncate does not move file pointer
    f.write(content)
  time.sleep(1)

  path = os.path.join(prefix, str(traj))
  if os.path.exists(path):
    N_files = len(os.listdir(path))   # number of points that have already been calculated
  else:
    N_files = 0

  n_jobs = math.ceil( (TOTAL_NUM_POINTS - N_files) / POINTS_PER_JOB)  # submit n_jobs

  print(f'traj: {traj}. Submitting {n_jobs} jobs')
  for _ in range(n_jobs):
    run_cmd("sbatch submit.sh")

