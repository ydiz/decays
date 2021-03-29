# python3 submit.py [traj]
import re
import os
import time
import sys
import glob

def run_cmd(cmd):
  print(cmd)
  os.system(cmd)
  time.sleep(3)

for f in glob.glob('slurm-*.out'):
  os.remove(f)

if len(sys.argv) != 2:
  print("There must be exactly one argument -- the trajectory number")
  sys.exit(0)

traj = sys.argv[1]

diagrams = ["typeIII_D1a", "typeIII_D1b"]

for diagram in diagrams:
  with open("submit.sh", 'r+') as f:
    content = f.read()
    content = re.sub(r'traj=(.+)', 'traj='+str(traj), content)
    content = re.sub(r'diagram=(.+)', 'diagram='+diagram, content)
    f.truncate(0)   # erase all content of the file
    f.seek(0)       # truncate does not move file pointer
    f.write(content)

  time.sleep(1)

  run_cmd("sbatch submit.sh")
