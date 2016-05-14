cat > $tjobname <<EOF 
#!/bin/bash
#
# -p regular
#SBATCH -p $jtype
#SBATCH --qos=$qos
#SBATCH -N $nnode
#SBATCH -t ${whr}:${wmin}:00
#SBATCH -e ${logdir}/${exename}_%j.err
#SBATCH -o ${logdir}/${exename}_%j.out
#SBATCH --ccm
#SBATCH --mail-type=all

###=============== start running bash script ==================

ulimit -s unlimited
export OMP_STACKSIZE=2g # 1GB for thread stack size
export OMP_NUM_THREADS=$nomp

srun -n ${nmpi_tasks} -N ${nnode}  -c $nomp ./bin/${exename} ${parsname} > ${logdir}/${exename}_${nside}${j0}${nj}.log

#-S ${ncpu_numa}
#-S number of MPI per numa node (numa node contains 6 cores, and a hopper node contains 4 numa nodes)
#-N number of MPI jobs per hopper node
#-n total number of MPI jobs
echo ""
echo ""
echo "srun -n ${nmpi_tasks} -c $nomp ./bin/$exename ${parsname}"

EOF

echo "srun -n ${nmpi_tasks} -N ${nnode} -c $nomp ./bin/$exename ${parsname}"

echo '***** To submit job do ***'
echo "sbatch $tjobname"