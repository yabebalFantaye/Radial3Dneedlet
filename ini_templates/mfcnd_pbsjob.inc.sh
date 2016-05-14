cat > $tjobname <<EOF 
#!/bin/bash
#
#Set a job name. This will allow you to identify your job in PBS.
#PBS -e ${logdir}/job.\${PBS_JOBID}.isim${nsim_start}_fsim${nnsim}.$disktype.err
#PBS -o ${logdir}/job.\${PBS_JOBID}.isim${nsim_start}_fsim${nnsim}.$disktype.log 
#
#PBS -q ${jtype}
#PBS -l mppwidth=${ncoretotal}

#PBS -lwalltime=${whr}:00:00
#PBS -V
# send me mail on job abort, begin, end
#PBS -m abe



###=============== start running bash script ==================

cd $PBS_O_WORKDIR
cd /global/homes/f/fantaye/planck

aprun -n ${ncpu_total} -N ${ncpu_hopper} -S ${ncpu_numa} ./bin/mfpol_cnd_planck ${parsname} > ${logdir}/${disktype}_isim${nsim_start}_fsim${nnsim}.log

#-S number of MPI per numa node (numa node contains 6 cores, and a hopper node contains 4 numa nodes)
#-N number of MPI jobs per hopper node
#-n total number of MPI jobs
echo ""
echo ""
echo "aprun -n ${nnode}  ./bin/mfpol_cnd_planck ${parsname}"

EOF
