# MembProtSim

A workflow to automate MD simulations of membrane-protein complexes using AMBER

## Conducting MD Simulations of Membrane-Proteins Using AMBER

This guide outlines the process for conducting equilibrium molecular dynamics (MD) simulations of proteins within lipid bilayers, assuming the user possesses a basic understanding of molecular dynamics and structural chemistry. The process is segmented into five key phases:

1. Preparation of initial protein input via PROPKA+PDB2PQR
2. Aligning proteins within lipid membranes correctly with the PPM local edition
3. Embedding proteins into selected lipid membranes using PACKMOL-MEMGEN
4. Creating system topology with tleap
5. Performing equilibrium MD simulation of the complex system using AMBER PMEMD

### 1. Preparation of Initial Protein Input via PROPKA+PDB2PQR
Issues such as incorrect protonation states of residues, termini, and definitions of disulfide bonds can arise from the original method of structure generation, necessitating thorough pre-modeling cleaning of the input protein structure. For initial cleaning, it is advised to use simple PDB editing tools like "pdb-tools" or Parmed. Structures originating from Rosetta might include residue-level energy values at the end of PDB files, which should be removed at this stage for optimal structure integrity. Similar to the GROMACS procedure, the input PDB structure needs to be examined for correct protonation. PROPKA and PDB2PQR can suggest modifications if the protonation states differ from the original files, with pdbGMXprep.py adjusting the PDB files based on PROPKA's recommendations. Unlike GROMACS, which identifies salt bridges during "pdb2gmx", AMBER's MD package identifies disulfide bonds early on. Accurate identification of chain IDs, residues, and termini is crucial for the successful execution of these preparatory steps.

<pre>
function pqr2pdb () {
inputPDB=$1
script_dir=$2
destination=$3
pdbfilename=$(basename "$inputPDB")  
pdbname=${pdbfilename%.*}
pdb2pqr30 --with-ph 7.4 --ff=AMBER --titration-state-method propka --ffout=AMBER \ 
    --pdb-output ${destination}/${pdbname}_pka7.pdb ${inputPDB} ${destination}/${pdbname}_pka7.pqr
propka3 ${inputPDB}
mv ${pdbname}.pka ${destination}/${pdbname}.pka
python ${script_dir}/pdbGMXprep.py --pkafile ${destination}/${pdbname}.pka --pH 7.4 \
    --mutatedPDB ${destination}/${pdbname}_protonated.pdb --pdb ${destination}/${pdbname}_pka7.pdb
}  
</pre>

### 2. Orienting Proteins within Lipid Membranes Using PPM3 Local Version
To accurately place transmembrane proteins in lipid membranes, orientations are often derived from the OPM database or the PPM server. However, a local variant of this tool is included in the "PPM" directory of our pipeline. By using the local PPM version, we ensure the protein aligns correctly with a DOPC bilayer membrane. Further information on this protocol and its outcomes is available at the provided link. Detailed instructions by the original developers are located in the PPM directory, named "ppm3_instructions.docx". PPM3 necessitates the use of gfortran, and the format for the input file (*.inp) can be found in the "example" directory. Below is an example bash script for this process, with all required files located in the working directory.

<pre>
module load gcc
./immers<2membrane.inp>output_file
</pre>

The generated output PDB file includes the adjusted protein structure and dummy beads to depict the membrane-water interface, typically phosphorus atoms in DOPC or the Carbonyl group at the interface. Prior to the next step, it's necessary to remove these "HETATM" lines using the "sed" command.
<pre>
sed '/HETATM\|DUM/d' input.pdb > output.pdb
</pre>

### 3. Embedding Proteins into Lipid Membranes with PACKMOL-MEMGEN
Within the AMBER module, PACKMOL-MEMGEN comes pre-installed. With AMBER activated, users can employ PACKMOL-MEMGEN to create systems of membranes with or without embedded proteins. The necessary arguments for this phase are detailed in the commands below. Since PPM3 local version has already determined protein orientation, the next steps focus on embedding into the lipid membrane.

**--pdb pdbname_from_ppm3.pdb**: Specifies the input transmembrane protein PDB
**--lipids POPC**: Constructs a POPC lipid membrane (lipid composition is adjustable)
**--ratio 1**: Sets the POPC ratio to 1 (experiment with different compositions)
**--preoriented**: Indicates the protein is oriented along the z-axis as per PPM coordinates
**--salt --salt_c Na+ --saltcon 0.15**: Adds a 0.15 M NaCl salt concentration to the aqueous layer
**--dist 10**: Sets the minimum distance from the protein to the box boundary in x, y, z directions
**--dist_wat 15**: Specifies a 15 Ångström thickness for the water layer
**--notprotonate --nottrim**: Omits further processing of the input protein PDB file

<pre>
packmol-memgen --help packmol-memgen --pdb destination/pdbname_from_ppm3.pdb --lipids POPC --ratio 1 --preoriented --salt --salt_c Na+ --saltcon 0.15 --dist 10 --dist_wat 15 --notprotonate --nottrim
</pre>
    
### 4. Creating System Topology with tleap
System topology files, containing necessary force field parameters for MD simulations, can be generated in AMBER using various methods. Here, tleap is used for topology construction, which can be done either interactively or in a single step. For novel systems, it's recommended to utilize tleap's interactive mode via the command line. For parallelized runs, a single-step topology construction is preferred. To initiate tleap on Lovelace and Marjorie, load the AMBER package.

For an interactive session:

Enter each line in topol_build.leap one by one. Modify the lines if user has additional molecules in the system or  user wants to change the system force field parameters.  

For execution of topology building in single step:

<pre>
module load amber/22
tleap
</pre>

Follow the instructions in **topol_build.leap**, modifying as needed for additional molecules or to alter force field parameters.

For single-step execution:

<pre>
tleap -f system_build.leap
</pre>

Each line in the topol_build.leap is explained below.  

**source leaprc.protein.ff19SB**: AMBERFF19SB force field is used for protein-in-solution  
**source leaprc.water.tip3p**: TIP3P water molecules are added to solvate the system  
**source leaprc.gaff2**: If small molecules or ligands are added to the system, they need be to parametrized. In AMBER, we can use GAFF2 to parametrize these molecules using AM1-BCC theory.  
**source leaprc.lipid21**: Latest lipid force field parametes used by AMBER. The choice of the lipid force field is dependent on the combination of FF being used in the system.  
**loadamberparams frcmod.ionsjc_tip3p**: Monovalent ion parameters for Ewald and TIP3P  
**water.receptor = loadpdb protein.pdb**: Adding the transmembrane protein to the FF  
**membrane = loadpdb POPC_amber.pdb**: Adding the membrane lipid to the FF  
**rism_wat = loadpdb rism_sele.pdb**: Adding the initial water molecules to the FF  
**system = combine{ receptor membrane rism_wat }**: Combining these FF parametes together  
**set system box {XX.XXX XX.XXX XX.XXX }**: System box size is added  
**saveamberparm system protein_membrane.prmtop protein_membrane.inpcrd**: Saving the parameters  
**savepdb system protein_membrane.pdb**: Save the PDB files after building the topology  
**quit**: Quit tleap after saving the files  

**Hydrogen Mass Repartitioning** : Optionally user can also increase the simulation timestep to 4fs from 2fs with the help of hydrogen mass repartitioning. This can be done with the help of the Parmed package.  

<pre>
parmed -i hmass_parmed.in
</pre>

### 5. Equilibrium MD Simulation of the Complex System in AMBER PMEMD
Equilibrium MD simulation encompasses initial energy minimization, NVT equilibration, NPT equilibration, and production runs, executable as a unified package. These steps are optimized for GPU execution (with CPU-based minimization). This phase can integrate into the overall workflow or function independently. For initial explorations or system setups, the standalone bash script "run_MD.sh" is recommended. This script allows for flexible management of input files and simulation parameters, like timestep, thermostat, and barostat settings. The simulation steps and corresponding input files are outlined below, with execution possible via the command line: 

**01_Min.in**: very short minimization on CPU with pmemd. This is advised for membrane systems, which often have bad clashes of lipid chains to resolve. The CPU code is more robust in dealing with these than the GPU  
**02_Min2.in**: longer minimization on GPU  
**03_Heat.in, 04_Heat2.in**: heating to 303 K, with restraints on protein and membrane lipids  
**05_Back.in**: run 1 ns NPT with restraints on receptor backbone atoms  
**06_Calpha.in**: run 1 ns NPT with restraints on receptor carbon-alpha atoms  
**07_Prod.in**: run 100 ns NPT equilibration, all restraints removed  
**08_Long.in**: run 500 ns NPT production, with Monte Carlo barostat and hydrogen mass repartitioning  

The script can be run from command line:  

<pre>
sbatch run_MD.sh --prmtop system.prmtop --name system --hmass system.hmass.prmtop --dest ./
</pre>

Overall, the entire script can be run as given below:

<pre>
bash MemProtSim.sh [--input <input.pdb>] [--dest <destination directory>] [--chain <chains to be considered for insertion>] [--disulfide <list of residues in disulfide bonds>] [--help]
</pre>
