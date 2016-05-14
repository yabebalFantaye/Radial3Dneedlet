#!/bin/bash


NOW=$(date +"%m-%d-%Y_%T")

workdir=`pwd`

exename="ball2almn"
echo "narg=$#" 

if (( $#>0 )); then
	exename=$1
fi
echo "excutable file is: $exename"

nomp_in=-1
if (( $#>1 )); then
	nomp_in=$2
fi

keyname="SurfDensMap"
dirname="fits4096/${keyname}"

whr=1

ifirstshell=0
ilastshell=255
nshell=`echo $ilastshell+1 | bc`

#needlet params
j0=1
nj=10
bb=2.0
glnpow=2
wav=0 #standard needlet

nside=1024
lmax=2000
nnmax=`echo "$nshell/2" | bc`
nside_jmap=$nside

#include machine specific variables
if [ $NERSC_HOST == "edison" ]; then
	nsockets=2
	ncore_per_socket=12
	ncore_per_node=24
	nomp=12      #-c
fi

if [ $NERSC_HOST == "cori" ]; then
	nsockets=2
	ncore_per_socket=16 #two 16 haswell
	ncore_per_node=32
	nomp=16      #-c
fi


min_number() {
    printf "%s\n" "$@" | sort -g | head -n1
}

max_number() {
    printf "%s\n" "$@" | sort -g | tail -n1
}


if [ $exename == "almn2beta" ]; then
	nomp=`echo "$nsockets*$nomp" | bc`
fi
if (( $nomp_in>0 )); then
	nomp=$nomp_in
fi


#based on variable to MPI parallize, set key params
nmpi_tasks=$nj #`echo "($nnsim+$nomp-1)/$nomp" | bc` #-n
if [ $glnpow == 2 ]; then
	nmpi_tasks=1 
fi
ntask_per_node=`echo "($ncore_per_node+$nomp-1)/$nomp" | bc` #ceiling 
#take the min of ntask_per_node and nmpi_tasks
ntask_per_node="$(min_number ${ntask_per_node} ${nmpi_tasks})" 
nnode=`echo "(${nmpi_tasks}+${ntask_per_node}-1)/$ntask_per_node" | bc` #-N


echo "-----------------"
echo "Using $nnode $NERSC_HOST nodes: mpi_tasks=$nmpi_tasks with $ntask_per_node MPI and $nomp OMP tasks per node"
echo "-----------------"


#QOS 
#qos="normal"
qos="premium"
jtype="regular"


odir="${SCRATCH}/seedlets/${dirname}"
pdir="${odir}/parsers"
pbsdir="${odir}/pbs"
logdir="jobs/$dirname"

balltype=1
iwidth=1
if [ $exename == "ball2almn" ]; then
	tempdir="${odir}/b2a/"
	osubdir="${tempdir}/almnout/"

	inputfile="data/3Dneedlets/HealPixMap_nside1024_NoRnd_Radius0_std/SurfDensMap_snap_049."
	outputfile="${osubdir}/almn_all.unf"
	finext="fits"
fi

if [ $exename == "almn2ball" ]; then
	tempdir="${odir}/a2b/"
	osubdir="${tempdir}/maps/"

	inputfile="${odir}/b2a/almnout/almn_all.unf"
	outputfile="${osubdir}/map.unf"
	finext="unf"
fi

if [ $exename == "almn2beta" ]; then
	tempdir="${odir}/a2beta/"
	osubdir="${tempdir}/maps/"

	inputfile="${odir}/b2a/almnout/almn_all.unf"
	outputfile="${osubdir}/map.unf"
	finext="unf"
fi

if [ $exename == "beta2ball" ]; then
	tempdir="${odir}/beta2b/"
	osubdir="${tempdir}/maps/"
	inputfile="${odir}/a2beta/maps/map.unf"
	outputfile="${osubdir}/map.unf"
	finext="unf"
fi

mkdir -p $odir $osubdir $pdir $logdir $pbsdir $tempdir
#================  param file

parsname="$pdir/parser_${exename}.ini"

echo "parfile: $parsname"

cat > $parsname <<EOF
feedback=1

j0=$j0
nj=$nj
bb=$bb
glnpow=$glnpow
wav=$wav
nside_jmap=$nside_jmap

nshell=$nshell
nside=$nside
lmax=$lmax
nnmax=$nnmax
balltype=$balltype
ifirstshell=$ifirstshell
ilastshell=$ilastshell
iwidth=$iwidth
inputfile=$inputfile
outputfile=$outputfile
tempdir=$tempdir
finext=$finext

EOF

#=========================================


#===========================

tjobname="${pbsdir}/job_${exename}.pbs"

echo "job script : $tjobname"

cat > $tjobname <<EOF 
#!/bin/bash
#
#SBATCH -p regular
#SBATCH -N $nnode
#SBATCH -t ${whr}:00:00
#SBATCH -e ${logdir}/{exename}_%j.err
#SBATCH -o ${logdir}/{exename}_%j.out
#SBATCH --ccm
#SBATCH --mail-type=all

###=============== start running bash script ==================

ulimit -s unlimited
export OMP_STACKSIZE=2g # 1GB for thread stack size
export OMP_NUM_THREADS=$nomp

srun -n ${nmpi_tasks} -N ${nnode}  -c $nomp ./bin/${exename} ${parsname} > ${logdir}/${exename}_%j.log

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
#sbatch $tjobname


