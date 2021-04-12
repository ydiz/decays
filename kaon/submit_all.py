#!/usr/bin/python3 -u

import os
import time
import glob
import re

def run_cmd(cmd):
  print(cmd)
  os.system(cmd)
  time.sleep(3)

# trajs = [1990]
# trajs = [2000]
# trajs = range(2100, 2150, 10)
# trajs = [2000] + list(range(2100, 2150, 10))
# trajs = range(2150, 2260, 10)
# trajs = [2010, 2020, 2070, 2080, 2090]
trajs = [2030, 2040, 2050, 2060]


########################  directories and diagrams

root = os.getcwd()
dirs = ['typeI', 'typeII', 'typeIII', 'typeIV_a', 'typeIV_b', 'typeV', 'c_s_eta', 'JJ_pion']
# dirs = ['typeI', 'typeII', 'JJ_pion']

diagrams_dir = {
    'typeI': ["typeI_D1a", "typeI_D1b", "typeI_D2a", "typeI_D2b"],
    'typeII': ["typeII_D1a", "typeII_D1b"],
    'typeIII': ["typeIII_D1a", "typeIII_D1b"],
    'typeIV_a': ["typeIV_D1a", "typeIV_D2a"],
    'typeIV_b': ["typeIV_D1b", "typeIV_D2b"],
    'typeV': ["typeV"],
    'c_s_eta': ['c_s_eta'],
    'JJ_pion': ['JJ_pi']
    }

#################### Remove trajectories that cannot be calculated or have already been calculated.


sequential_dir = "/global/cfs/cdirs/mp13/ydzhao/24ID/sequential"
point_l_dir = "/global/cfs/cdirs/mp13/ydzhao/24ID/point_l"
point_s_dir = "/global/cscratch1/sd/ydzhao/point_s"
Lxx_dir = "/global/cfs/cdirs/mp13/ydzhao/24ID/Lxx"

doable_trajs = []
for traj in trajs:
  if not os.path.isdir(sequential_dir + '/' + str(traj)):
    print(f"traj {traj} does not have sequential propagators; skipping this trajectory")
    continue
  if not os.path.isdir(point_l_dir + '/' + str(traj)):
    print(f"traj {traj} does not have point_l propagators; skipping this trajectory")
    continue
  if not os.path.isdir(point_s_dir + '/' + str(traj)):
    print(f"traj {traj} does not have point_s propagators; skipping this trajectory")
    continue
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

      output_file = f'./output/{diagram}.out.{traj}'
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

