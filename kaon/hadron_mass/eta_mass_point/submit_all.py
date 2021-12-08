#!/usr/bin/python3 -u

import os
import time
import re

def run_cmd(cmd):
  print(cmd)
  os.system(cmd)
  time.sleep(3)

program = "eta_mass_point"

point_l_dir = "/global/cfs/cdirs/mp13/ydzhao/24ID/point_l"
point_s_dir = "/global/cfs/cdirs/mp13/ydzhao/24ID/point_s"

# for traj in range(1000, 2260, 10):
# for traj in range(1300, 1320, 10):
for traj in [1380]:
  print(traj)

  if not os.path.isdir(point_l_dir + '/' + str(traj)):
    print(f"traj {traj} does not have point_l propagators; skipping this trajectory")
    continue
  if not os.path.isdir(point_s_dir + '/' + str(traj)):
    print(f"traj {traj} does not have point_s propagators; skipping this trajectory")
    continue

  output_file = f'./eta_mass_point_output/{program}.out.{traj}'
  if os.path.exists(output_file):
# if os.stat(output_file).st_size < 10000:  # Less than 10K
# 	print(f"{d}/{output_file} is too small; will submit a job to redo this")
# else:
    print(f"eta_mass_point_output/{output_file} already exists; skipping this job")
    continue




  with open("submit.sh", 'r+') as f:
    content = f.read()
    content = re.sub(r'#SBATCH -J (.*)', '#SBATCH -J ' + f'{program}.{traj}', content)
    content = re.sub(r'program=(.+)', f'program={program}', content)
    content = re.sub(r'traj_start=(.+)', f'traj_start={traj}', content)
    f.truncate(0)   # erase all content of the file
    f.seek(0)       # truncate does not move file pointer
    f.write(content)

  time.sleep(1)

  run_cmd("sbatch submit.sh")
