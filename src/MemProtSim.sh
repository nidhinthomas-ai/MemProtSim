#!/bin/bash

###############################################################################################################################################
## This script is used to prepare AMBER topology for the membrane-protein system. In this workflow, we use amber19FFSB FF for protein, Lipid21 
## for POPC, TIP3P for water, ionsjc_tip3p for ions and GAFF2 for small molecules. A system_build.leap input file is created to run in "tleap"
## program. 

## Usage:

##    bash MemProtSim.sh --input input.pdb --dest ./ --chain AB --disulfide ACYX96 ACYX176 ACYX413 ACYX416

## Written by Nidhin Thomas (dr.thomas@nidhin-thomas.com)

###############################################################################################################################################

usage() { echo -e "\nThis script is used to prepare AMBER topology for the membrane-protein system.\n";

          echo -e "bash MemProtSim.sh [--input <input.pdb>] [--dest <destination directory>] [--chain <chains to be considered for insertion>] [--disulfide <list of residues in disulfide bonds>] [--help]\n" 1>&2;

          echo -e "Example: bash MemProtSim.sh --input input.pdb --dest ./ --chain AB --disulfide ACYX96 ACYX176 ACYX413 ACYX416\n";

          exit 1; 
        }

for arg in "$@"; do
  shift
  case "$arg" in
    '--help')      set -- "$@" '-h'   ;;
    '--input')     set -- "$@" '-f'   ;;
    '--dest')      set -- "$@" '-d'   ;;
    '--chain')     set -- "$@" '-c'   ;;
    '--disulfide') set -- "$@" '-s'   ;;
    *)             set -- "$@" "$arg" ;;
    
  esac
done

# Default behavior
dest="./";

OPTIND=1
while getopts "hf:d:c:s:" opt
do
  case "$opt" in
    'h') usage; exit 0 ;;
    'f') input=$OPTARG ;;
    'd') dest=$OPTARG ;;
    'c') chain=$OPTARG ;;
    's') disulfide=$OPTARG ;;    
    '?') usage >&2; exit 1 ;;
  esac
done
shift $(expr $OPTIND - 1)

script_dir=$(dirname "$(realpath "$0")")


source ~/mambaforge/etc/profile.d/conda.sh
conda activate MemProtSim
    
module load amber/22
module load gcc

export AMBERHOME=path/to/amber/22

function pqr2pdb () {
inputPDB=$1
script_dir=$2
destination=$3
pdbfilename=$(basename "$inputPDB")  
pdbname=${pdbfilename%.*}
pdb2pqr30 --with-ph 7.4 --ff=AMBER --titration-state-method propka --ffout=AMBER --pdb-output ${destination}/${pdbname}_pka7.pdb ${inputPDB} ${destination}/${pdbname}_pka7.pqr
propka3 ${inputPDB}
mv ${pdbname}.pka ${destination}/${pdbname}.pka
python ${script_dir}/pdbGMXprep.py --pkafile ${destination}/${pdbname}.pka --pH 7.4 --mutatedPDB ${destination}/protein.pdb --pdb ${destination}/${pdbname}_pka7.pdb
}

function ppmOrient () {
inputPDB=$1
destination=$2
script_dir=$3
${chain}=$4
python ${script_dir}/prepPPMinput.py --pdb ${inputPDB} --dest ${destination} --ppm ${script_dir}/ppm --chain ${chain}
${script_dir}/ppm/immers<${destination}/protein.in>${destination}/protein.out
sed '/DUM/d' ${destination}/proteinout1.pdb > ${destination}/protein_ppm3.pdb
packmol-memgen --pdb ${destination}/protein_ppm3.pdb --lipids POPC --ratio 1 --preoriented --salt --salt_c Na+ --saltcon 0.15 --dist 10 --dist_wat 15 --notprotonate --nottrim
}

##Protonate the protein properly using PROPKA.
pqr2pdb ${input} ${script_dir} ${dest}

##Runs the PPM local version to insert the protein into the lipid bilayer.
ppmOrient ${dest}/protein.pdb ${dest} ${script_dir} ${chain}

## Prepare input files for system parameter generation and hydrogen mass repartitioning
python ${script_dir}/prepLeapHmass.py --prot ${dest}/bilayer_protein_ppm3.pdb --leap ${dest}/system.leap --hmass ${dest}/hmass.in --box ${boxsize} --disu ${disulfide}

## TLEAP for system parametrization
tleap -f ${dest}/system.leap

## HMASS for hydrogen mass repartitioning
parmed -i ${dest}/hmass.in

## Final MD simulation
# python ${script_dir}/membProtMD.py --input ${dest}/system.pdb --json ${json} --dest ${dest}

sbatch ${script_dir}/run_MD.sh --prmtop ${dest}/system.prmtop --name ${dest}/system --hmass ${dest}/system.hmass.prmtop --dest ${dest}