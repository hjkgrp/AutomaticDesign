import os
import subprocess
## check with the queue manager
base_name = gen
ll = subprocess.check_output('qstat -j ' + str(base_name),shell=True)
print(ll)  
 
