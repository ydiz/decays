# python3 -u submit_all_diagrams.py

import re
import os
import time

diagrams = ["typeV"]

for diagram in diagrams:

  with open("submit.sh", 'r') as f:
    content = f.read()
    new_content = re.sub(r'diagram=(.+)', 'diagram='+diagram, content)
  with open('submit.sh', 'w') as f:
    f.write(new_content)
    time.sleep(1)

  cmd = "sbatch submit.sh"
  print(cmd)
  os.system(cmd)
  time.sleep(3)
