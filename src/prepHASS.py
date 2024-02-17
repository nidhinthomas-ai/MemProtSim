import argparse

def write_hmass_file(input_pdb, output_parmed):
    with open('hmass.in', 'w') as f:
        f.write('parm system.prmtop\n')
        f.write('HMassRepartition\n')
        f.write('outparm system.hmass.prmtop \n')
        f.write('quit\n')

def main():
    parser = argparse.ArgumentParser(description='Generate a .in file.')
    parser.add_argument('input_prmtop', type=str, help='path to input prmtop file')
    parser.add_argument('output_prmtop', type=str, help='path to output prmtop file')
    args = parser.parse_args()

    write_hmass_file(args.input_prmtop, args.output_prmtop)

if __name__ == '__main__':
    main()
