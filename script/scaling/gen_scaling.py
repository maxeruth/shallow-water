import math
import subprocess

def main():
  proc_dimx = range( 1,  4 )
  tstep     = range( 1,  3 )
  str_res   = 1000

  # generate submission scripts
  for dimx in proc_dimx:
    for step in tstep:
      nproc = dimx * dimx
      nodes = int(math.ceil( nproc / 12.0 ))
      tasks = int(math.ceil( float(nproc) / nodes ))

      filename = "strong_scaling_{dimx}x{dimx}_{step}.sub".format( dimx = dimx, step = step )
      with open( filename, 'w' ) as f:
        script = """#!/bin/bash

#SBATCH -J shallow_strong_{dimx}x{dimx}_{step}
#SBATCH -o /home/%u/shallow-water/script/scaling/strong_{dimx}x{dimx}_{step}.out
#SBATCH -e /home/%u/shallow-water/script/scaling/strong_{dimx}x{dimx}_{step}.err
#SBATCH --nodes={nodes}
#SBATCH --ntasks={nproc}
#SBATCH --tasks-per-node={tasks}
#SBATCH --cpus-per-task=1
#SBATCH --get-user-env
#SBATCH -t 00:20:00
#SBATCH --mem-per-cpu=1000
#SBATCH --partition=cs5220

source /etc/profile.d/modules.sh
module load openmpi-4.0.0
cd $HOME/shallow-water
mpirun -np {nproc} src/lshallow tests.lua dam {str_res} {dimx} {dimx} {step}
"""
        script = script.format( dimx = dimx, step = step, nodes = nodes,
                                nproc = nproc, tasks = tasks, str_res = str_res )
        f.write( script )
        f.close

      weak_res = dimx * 200
      filename = "weak_scaling_{dimx}x{dimx}_{step}.sub".format( dimx = dimx, step = step )
      with open( filename, 'w' ) as f:
        script = """#!/bin/bash

#SBATCH -J shallow_weak_{dimx}x{dimx}_{step}
#SBATCH -o /home/%u/shallow-water/script/scaling/weak_{dimx}x{dimx}_{step}.out
#SBATCH -e /home/%u/shallow-water/script/scaling/weak_{dimx}x{dimx}_{step}.err
#SBATCH --nodes={nodes}
#SBATCH --ntasks={nproc}
#SBATCH --tasks-per-node={tasks}
#SBATCH --cpus-per-task=1
#SBATCH --get-user-env
#SBATCH -t 00:20:00
#SBATCH --mem-per-cpu=1000
#SBATCH --partition=cs5220

source /etc/profile.d/modules.sh
module load openmpi-4.0.0
cd $HOME/shallow-water
mpirun -np {nproc} src/lshallow tests.lua dam {weak_res} {dimx} {dimx} {step}
"""
        script = script.format( dimx = dimx, step = step, nodes = nodes,
                                nproc = nproc, tasks = tasks, weak_res = weak_res )
        f.write( script )
        f.close


  # submit and delete submission scripts
  for dimx in proc_dimx:
    for step in tstep:
      filename = "strong_scaling_{dimx}x{dimx}_{step}.sub".format( dimx = dimx, step = step )
      subprocess.call( [ "sbatch", filename ] )
      subprocess.call( [ "rm",     filename ] )

      filename = "weak_scaling_{dimx}x{dimx}_{step}.sub".format( dimx = dimx, step = step )
      subprocess.call( [ "sbatch", filename ] )
      subprocess.call( [ "rm",     filename ] )


main()
