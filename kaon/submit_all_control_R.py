#!/usr/bin/python3 -u

import os
import time
import glob
import re

def run_cmd(cmd):
  print(cmd)
  os.system(cmd)
  time.sleep(3)

# # trajs = range(1900, 2010, 10)
# # trajs = range(1010, 1170, 10)
# trajs = range(2000, 2260, 10)
# trajs = range(1800, 1870, 10)
# trajs = range(1500, 1550, 10)
# trajs = range(1200, 1550, 10)
# trajs = range(1100, 1200, 10)
# trajs = range(1000, 1550, 10)
# trajs = range(1550, 1800, 10)
trajs = [1090, 1130]
# trajs = [2250]
# trajs = [1910]


########################  directories and diagrams

# root = os.getcwd()
root = os.getcwd() + '/control_R_new_strategy'
# dirs = ['typeV']
# dirs = ['typeIII', 'typeIV_a', 'typeIV_b']
# dirs = ['typeIII', 'typeIV_a']
dirs = ['typeIV_b']

diagrams_dir = {
    # 'typeI': ["typeI_D1a", "typeI_D1b", "typeI_D2a", "typeI_D2b"],
    # 'typeII': ["typeII_D1a", "typeII_D1b"],
    'typeIII': ["typeIII_D1a", "typeIII_D1b"],
    'typeIV_a': ["typeIV_D1a"],
    'typeIV_b': ["typeIV_D1b"],
    'typeV': ["typeV"],
    # 'c_s_eta': ['c_s_eta'],
    # 'JJ_pion': ['JJ_pi']
    }

#################### Remove trajectories that cannot be calculated or have already been calculated.


sequential_dir = "/global/cfs/cdirs/mp13/ydzhao/24ID/sequential"
point_l_dir = "/global/cfs/cdirs/mp13/ydzhao/24ID/point_l"
point_s_dir = "/global/cscratch1/sd/ydzhao/point_s"
Lxx_dir = "/global/cfs/cdirs/mp13/ydzhao/24ID/Lxx"

doable_trajs = []
for traj in trajs:
  if 'typeIII' in dirs:
    if not os.path.isdir(sequential_dir + '/' + str(traj)):   # !!! FIXME: uncomment this if running typeIII
      print(f"traj {traj} does not have sequential propagators; skipping this trajectory")
      continue
  if not os.path.isdir(point_l_dir + '/' + str(traj)):
    print(f"traj {traj} does not have point_l propagators; skipping this trajectory")
    continue
  # if not os.path.isdir(point_s_dir + '/' + str(traj)):
  #   print(f"traj {traj} does not have point_s propagators; skipping this trajectory")
  #   continue
  if not os.path.isfile(Lxx_dir + '/Lxx.' + str(traj)):
    print(f"traj {traj} does not have Lxx; skipping this trajectory")
    continue

  doable_trajs.append(traj)

trajs = sorted(doable_trajs)

print(f"Will submit jobs for these trajectories: {trajs}")


################### submit jobs

for d in dirs:
  os.chdir(os.path.join(root, d))
  print('=' * 70)
  print("Current directory:", os.getcwd())

  for f in glob.glob('slurm-*.out'):
    os.remove(f)

  diagrams = diagrams_dir[d]

  for traj in trajs:
    for diagram in diagrams:

      output_file = f'./output_new_strategy/{diagram}.out.{traj}'
      if os.path.exists(output_file):
        if os.stat(output_file).st_size < 10000:  # Less than 10K
          print(f"{d}/{output_file} is too small; will submit a job to redo this")
        else:
          print(f"{d}/{output_file} already exists; skipping this job")
          continue

      with open("submit.sh", 'r+') as f:
        content = f.read()
        content = re.sub(r'#SBATCH -J (.*)', '#SBATCH -J ' + f'{diagram}.{traj}', content)
        content = re.sub(r'traj=(.+)', 'traj='+str(traj), content)
        content = re.sub(r'diagram=(.+)', 'diagram='+diagram, content)
        f.truncate(0)   # erase all content of the file
        f.seek(0)       # truncate does not move file pointer
        f.write(content)

      time.sleep(1)

      run_cmd("sbatch submit.sh")

