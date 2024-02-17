#!/bin/bash

########################################################################################################################################

## This script is used to run MD simulation of membrane proteins in AMBER using AMBERFF19SB/LIPID21/GAFF2 force field. Protein is expected to be inserted into the custom lipid membrane prior to this step. The system is first minimized and then heated up to 303K. The protein heavy atoms are restrained before production run to relax the system. Furthermore, after 100 ns of production run, the timestep is modified from 0.002 ps to 0.004 ps with hydrogen mass repartitioning to accelerate the sampling.

## Written by Nidhin Thomas (dr.thomas@nidhin-thomas.com)

########################################################################################################################################

## sbatch run_MD.sh --prmtop system.prmtop --name system --hmass system.hmass.prmtop --dest ./

########################################################################################################################################

#SBATCH --time=2:00:00
#SBATCH -N 1 -n 1
#SBATCH -p gpu
#SBATCH --gres=gpu:1

module load amber/22

usage() { echo -e "\nThis script is used to run MD simulation of membrane proteins in AMBER using AMBERFF19SB/LIPID21/GAFF2 force field. Protein is expected to be inserted into the custom lipid membrane prior to this step. The system is first minimized and then heated up to 303K. The protein heavy atoms are restrained before production run to relax the system. Furthermore, after 100 ns of production run, the timestep is modified from 0.002 ps to 0.004 ps with hydrogen mass repartitioning to accelerate the sampling.\n";

          echo -e "sbatch run_MD.sh [--prmtop <system.prmtop>] [--name <system>] [--hmass <system.hmass.prmtop>] [--dest <destination directory>] [--help]\n" 1>&2;

          echo -e "Example: sbatch path/to/run_MD.sh --prmtop system.prmtop --name system --hmass system.hmass.prmtop --dest ./ \n";

          exit 1; 
        }

for arg in "$@"; do
  shift
  case "$arg" in
    '--help')      set -- "$@" '-h'   ;;
    '--prmtop')   set -- "$@" '-p'    ;;
    '--name')    set -- "$@" '-n'     ;;
    '--hmass')      set -- "$@" '-m'  ;;
    '--dest')      set -- "$@" '-d'   ;;
    '--script')    set -- "$@" '-s'   ;;
    *)             set -- "$@" "$arg" ;;
    
  esac
done

OPTIND=1
while getopts "hp:n:m:d:s:" opt
do
  case "$opt" in
    'h') usage; exit 0 ;;
    'p') prmtop=$OPTARG ;;
    'n') name=$OPTARG;;
    'm') hmass=$OPTARG  ;;
    'd') dest=$OPTARG  ;;
    's') script=$OPTARG  ;;
    '?') usage >&2; exit 1 ;;
  esac
done
shift $(expr $OPTIND - 1)

script_dir=${script}
echo "${script_dir}"
input_dir=${script_dir}/md_inputs

echo "${input_dir}"

AMBERHOME='path/to/amber/22'

## Minimise ##
$AMBERHOME/bin/pmemd \
-O \
-i ${input_dir}/01_Min.in \
-o ${dest}/01_Min.out \
-p ${dest}/$prmtop \
-c ${dest}/${name}.inpcrd \
-r ${dest}/01_Min.rst 

$AMBERHOME/bin/pmemd.cuda \
-O \
-i ${input_dir}/02_Min2.in \
-o ${dest}/02_Min2.out \
-p ${dest}/$prmtop \
-c ${dest}/01_Min.rst \
-r ${dest}/02_Min2.rst 

### Equilibration ##
$AMBERHOME/bin/pmemd.cuda \
-O \
-i ${input_dir}/03_Heat.in \
-o ${dest}/03_Heat.out \
-p ${dest}/$prmtop \
-c ${dest}/02_Min2.rst \
-r ${dest}/03_Heat.rst \
-x ${dest}/03_Heat.nc \
-ref ${dest}/02_Min2.rst

$AMBERHOME/bin/pmemd.cuda \
-O \
-i ${input_dir}/04_Heat2.in \
-o ${dest}/04_Heat2.out \
-p ${dest}/$prmtop \
-c ${dest}/03_Heat.rst \
-r ${dest}/04_Heat2.rst \
-x ${dest}/04_Heat2.nc \
-ref ${dest}/03_Heat.rst

## Peptide backbone restrained
$AMBERHOME/bin/pmemd.cuda \
-O \
-i ${input_dir}/05_Back.in \
-o ${dest}/05_Back.out \
-p ${dest}/$prmtop \
-c ${dest}/04_Heat2.rst \
-r ${dest}/05_Back.rst \
-x ${dest}/05_Back.nc \
-ref ${dest}/04_Heat2.rst \
-inf ${dest}/05_Back.mdinfo

## Peptide C-alpha atoms only restrained
$AMBERHOME/bin/pmemd.cuda \
-O \
-i ${input_dir}/06_Calpha.in \
-o ${dest}/06_Calpha.out \
-p ${dest}/$prmtop \
-c ${dest}/05_Back.rst \
-r ${dest}/06_Calpha.rst \
-x ${dest}/06_Calpha.nc \
-ref ${dest}/05_Back.rst \
-inf ${dest}/06_Calpha.mdinfo

## 100ns NPT run ##
$AMBERHOME/bin/pmemd.cuda \
-O \
-i ${input_dir}/07_Prod.in \
-o ${dest}/07_Prod_$name.out \
-p ${dest}/$prmtop \
-c ${dest}/06_Calpha.rst \
-r ${dest}/07_Prod_$name.rst \
-x ${dest}/07_Prod_$name.nc \
-inf ${dest}/07_Prod_$name.mdinfo

## 500ns HMASS ##
$AMBERHOME/bin/pmemd.cuda \
-O \
-i ${input_dir}/08_Long.in \
-o ${dest}/08_Long_$name.out \
-p ${dest}/$hmass \
-c ${dest}/07_Prod_$name.rst \
-r ${dest}/08_Long_$name.rst \
-x ${dest}/08_Long_$name.nc \
-inf ${dest}/08_Long_$name.mdinfo \
