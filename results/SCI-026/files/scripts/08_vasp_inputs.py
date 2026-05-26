from pathlib import Path
from textwrap import dedent

ROOT = Path(__file__).resolve().parents[1]
INPUT_DIR = ROOT / 'inputs'
INPUT_DIR.mkdir(exist_ok=True)

files = {
    'INCAR_relax': dedent('''
        SYSTEM = Li6PS5Cl_LiCoO2_interface_relax
        ENCUT = 520
        EDIFF = 1E-5
        EDIFFG = -0.02
        IBRION = 2
        ISIF = 2
        NSW = 150
        ISPIN = 2
        PREC = Accurate
        ALGO = Normal
        ISMEAR = 0
        SIGMA = 0.05
        LASPH = .TRUE.
        LREAL = Auto
        LCHARG = .FALSE.
        LWAVE = .FALSE.
    ''').strip() + '\n',
    'INCAR_neb': dedent('''
        SYSTEM = Li_migration_interface_CI_NEB
        ENCUT = 520
        EDIFF = 1E-5
        IBRION = 3
        POTIM = 0
        NSW = 200
        IMAGES = 5
        SPRING = -5
        LCLIMB = .TRUE.
        IOPT = 1
        ISMEAR = 0
        SIGMA = 0.05
        LCHARG = .FALSE.
        LWAVE = .FALSE.
    ''').strip() + '\n',
    'INCAR_aimd': dedent('''
        SYSTEM = Li6PS5Cl_LiCoO2_interface_AIMD
        ENCUT = 450
        EDIFF = 1E-4
        IBRION = 0
        NSW = 5000
        POTIM = 2.0
        TEBEG = 600
        TEEND = 600
        SMASS = 0
        MDALGO = 2
        ISMEAR = 1
        SIGMA = 0.1
        PREC = Normal
        LWAVE = .FALSE.
        LCHARG = .FALSE.
    ''').strip() + '\n',
    'KPOINTS': dedent('''
        Interface supercell
        0
        Gamma
        3 3 1
        0 0 0
    ''').strip() + '\n',
    'lammps_interface.in': dedent('''
        units           metal
        dimension       3
        boundary        p p p
        atom_style      charge

        read_data       interface.data
        pair_style      buck/coul/long 10.0
        kspace_style    pppm 1.0e-5
        pair_coeff      * * 0.0 1.0 0.0

        neighbor        2.0 bin
        neigh_modify    every 1 delay 0 check yes
        timestep        0.001

        thermo          200
        thermo_style    custom step temp pe ke etotal press
        min_style       cg
        minimize        1.0e-8 1.0e-10 2000 10000

        velocity        all create 600.0 4928459 rot yes dist gaussian
        fix             1 all nvt temp 600.0 600.0 0.1
        run             20000
        unfix           1

        write_data      interface_equilibrated.data
    ''').strip() + '\n',
}

for name, content in files.items():
    path = INPUT_DIR / name
    path.write_text(content, encoding='utf-8')

print('Generated input files:')
for name in files:
    print(f' - {INPUT_DIR / name}')
