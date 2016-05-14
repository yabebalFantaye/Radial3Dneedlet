#!/bin/bash

##--------------------------------------------------
#Author: 
#    Yabebal Fantaye
#     May 12, 2016
#      
#Purpose: 
#   This script creates parameter file and cluster job 
#   files based on scripts found in ./ini_template/ folder.
#
#Optional Input args:
#   1 - excutable name in ../bin   
#   2 - OMP_NUM_THREADS  
#   3 - Nside
#   4 - j0 (Starting Needlet j)
#   4 - nj (Number of j from j0)
#   6 - bb (Needlet B)
#
#Example:
#  ./qsub_3dneedlet.sh ball2almn 16 1024 1 11 2
#--------------------------------------------------

scheduler='slurm' #or pbs

NOW=$(date +"%m-%d-%Y_%T")

workdir=`pwd`

echo "------------------------"
echo "number of passed arguments =$#" 
echo "------------------------"

#which excutable to use
exename="ball2almn"
if (( $#>0 )); then
	exename=$1
fi
echo "excutable file is: $exename"

#key names
keyname="SurfDensMap"
dirname="fits1024/${keyname}"

#walltime hr and min
whr=1
wmin=0
echo "walltime: $whr:$wmin"

#----------------------------------------------------------
#------------------ set parameter file variables  ---------
#----------------------------------------------------------
#number of shell 
ifirstshell=0
ilastshell=255
nshell=`echo $ilastshell+1 | bc`

#needlet params
j0=1
nj=11
bb=2.0
glnpow=1
wav=0 #standard needlet

nside=1024
lmax=3000
nnmax=`echo "$nshell/2" | bc`
nside_jmap=$nside


#openMP variabl: value for OMP_NUM_THREADS
nomp_in=-1
if (( $#>1 )); then
	nomp_in=$2
fi

if (( $#>2 )); then
	nside=$3
fi
if (( $#>3 )); then
	j0=$4
fi
if (( $#>4 )); then
	nj=$5
fi
if (( $#>5 )); then
	bb=$6
fi


#----------------------------------------------------------
#--------------- set cluster job file parameters ---------
#----------------------------------------------------------
#include machine specific variables
source ini_templates/cluster_spec.inc.sh

#QOS 
#qos="normal"
qos="premium"
jtype="regular"


min_number() {
    printf "%s\n" "$@" | sort -g | head -n1
}

#
#-----based on variable to MPI parallize, set key params
#
nmpi_tasks=$nj #`echo "($nnsim+$nomp-1)/$nomp" | bc` #-n
if [ $exename == "ball2almn" ] || [ $exename == "almn2ball" ] ; then
	nmpi_tasks=`echo "($nshell+$nomp-1)/$nomp" | bc` #ceiling
fi
#if [ $glnpow == 2 ]; then
#	nmpi_tasks=1 
#fi
ntask_per_node=`echo "($ncore_per_node+$nomp-1)/$nomp" | bc` #ceiling 
#take the min of ntask_per_node and nmpi_tasks
ntask_per_node="$(min_number ${ntask_per_node} ${nmpi_tasks})" 
nnode=`echo "(${nmpi_tasks}+${ntask_per_node}-1)/$ntask_per_node" | bc` #-N

echo "-----------------"
echo "Using $nnode $NERSC_HOST nodes: mpi_tasks=$nmpi_tasks with $ntask_per_node MPI and $nomp OMP tasks per node"
echo "Que type: partition=$jtype and QOS=$qos"
echo "-----------------"

#----------------------------------------------------------
#------------- define directory and IO file names ---------
#----------------------------------------------------------
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

#----------------------------------------------------------
#------------- By now all needed variables are defined ----
#------------------ produce *.ini and *.job files ---------
#----------------------------------------------------------

#*.ini
parsname="$pdir/parser_${exename}.ini"
echo "parfile: $parsname"
#include the bash parameter file template 
source ini_templates/r3dneedlet_params.inc.sh

#*.job
#===========================
tjobname="${pbsdir}/job_${exename}.pbs"
echo "job script : $tjobname"
#include the appropriate job script
source ini_templates/${scheduler}_job.inc.sh

#sbatch $tjobname


