import argparse
import re

def create_disulfide_bond(residue1, residue2, leapfile):
    with open(leapfile, 'a') as f:
        f.write(f"bond system.{residue1}.SG system.{residue2}.SG\n")
    
def create_disulfide_bond(residue1, residue2, leapfile):
    if residue1[1:4] != 'CYX' or residue2[1:4] != 'CYX':
        raise ValueError("Disulfide bond can only be created between CYS residues.")
    elif residue1[4:] == residue2[4:]:
        raise ValueError("CYX residues must be in pairs for disulfide bond.")
    with open(leapfile, 'a') as f:
        f.write(f"bond system.{residue1[0]}.{residue1[1:]}.SG system.{residue2[0]}.{residue2[1:]}.SG\n")

def write_leap_file(pdb, box, leapfile, disulfide_bonds=[]):
    with open(leapfile, 'w') as f:
        f.write('source leaprc.protein.ff19SB\n')
        f.write('source leaprc.water.tip3p\n')
        f.write('source leaprc.gaff2\n')
        f.write('source leaprc.lipid21\n')
        f.write('loadamberparams frcmod.ionsjc_tip3p\n')
        f.write(f'system = loadpdb {pdb}\n')
        f.write(f'set system box {box}\n')
        
        # Create disulfide bonds
        for bond in disulfide_bonds:
            create_disulfide_bond(bond[0], bond[1], leapfile)
        
        f.write('saveamberparm system system.prmtop system.inpcrd\n')
        f.write('savepdb system membrane_protein.pdb\n')
        f.write('quit\n')

def write_hmass_file(input_pdb, hmassfile):
    with open(hmassfile, 'w') as f:
        f.write('parm system.prmtop\n')
        f.write('HMassRepartition\n')
        f.write('outparm system.hmass.prmtop \n')
        f.write('quit\n')

def main():
    parser = argparse.ArgumentParser(description='Generate a .leap file.')
    parser.add_argument('--prot', type=str, help='path to protein_membrane PDB file')
    parser.add_argument('--leap', type=str, help='output path to leap file')
    parser.add_argument('--hmass', type=str, help='output path to hmass input file')
    parser.add_argument('--disu', type=str, nargs='+', default=[], help='Residue pairs to create disulfide bond; e.g.: CYX1 CYX2 CYX3 CYX4')
    parser.add_argument('--log', type=str, help='path to packmol-memgen.log file')
    args = parser.parse_args()
    
    print(args.disu)
    
    # Read box dimensions from packmol-memgen.log file
    with open(args.log, 'r') as f:
        log_data = f.read()
        x_len = float(re.search(r'x_len\s+=\s+([\d\.]+)', log_data).group(1))
        y_len = float(re.search(r'y_len\s+=\s+([\d\.]+)', log_data).group(1))
        z_len = float(re.search(r'z_len\s+=\s+([\d\.]+)', log_data).group(1))
    
    box = f'{{x_len:.3f} {y_len:.3f} {z_len:.3f}}'

    disulfide_bonds = []
    for i in range(0, len(args.disu), 2):
        disulfide_bonds.append((args.disu[i], args.disu[i+1]))
    
    write_leap_file(args.prot, box, args.leap, disulfide_bonds)
    write_hmass_file("membrane_protein.pdb", args.hmass)

if __name__ == '__main__':
    main()