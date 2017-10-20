import os
from . import helpers


def input_block(config, section):
    """Return the input block as a string."""
    block = ''

    if section not in config:
        return ''

    for key in config[section]:
        value = config[section][key].strip()

        if value:
            block += '{} {}\n'.format(key, value)
        else:
            block += key + '\n'

    return block


def create_link_files(config, dirname='.'):
    """Create a 'link-files' executable script."""
    file = os.path.join(dirname, 'link-files')

    with open(file, 'a') as f:
        f.write('#!/bin/bash\n'
                'QE="../../qe"\n'
                '\n'
                'ln -sf $QE/2-wfn/WFN 1-epsilon/WFN\n'
                'ln -sf $QE/3-wfnq/WFN 1-epsilon/WFNq\n'
                'ln -sf $QE/4-wfn_co/RHO 2.1-sigma/RHO\n'
                'ln -sf $QE/4-wfn_co/WFN 2.1-sigma/WFN_inner\n'
                'ln -sf $QE/4-wfn_co/WFN 2.2-inteqp/WFN_co\n'
                'ln -sf $QE/4-wfn_co/WFN 3-kernel/WFN_co\n'
                'ln -sf $QE/4-wfn_co/WFN 4-absorption/WFN_co\n'
                'ln -sf $QE/4-wfn_co/vxc.dat 2.1-sigma/vxc.dat\n'
                'ln -sf $QE/5-wfn_fi/WFN 4-absorption/WFN_fi\n'
                'ln -sf $QE/5-wfn_fi/WFN 2.2-inteqp/WFN_fi\n'
                'ln -sf $QE/6-wfnq_fi/WFN 4-absorption/WFNq_fi\n'
                'ln -sf $QE/6-wfnq_fi/WFN 2.2-inteqp/WFNq_fi\n'
                '\n'
                'ln -sf ../1-epsilon/eps0mat 2.1-sigma\n'
                'ln -sf ../1-epsilon/eps0mat 3-kernel\n'
                'ln -sf ../1-epsilon/eps0mat 4-absorption\n'
                'ln -sf ../1-epsilon/epsmat 2.1-sigma\n'
                'ln -sf ../1-epsilon/epsmat 3-kernel\n'
                'ln -sf ../1-epsilon/epsmat 4-absorption\n'
                'ln -sf ../2.1-sigma/eqp1.dat ./2.2-inteqp/eqp_co.dat\n'
                'ln -sf ../2.1-sigma/eqp1.dat 4-absorption\n'
                'ln -sf ../3-kernel/bsedmat 4-absorption\n'
                'ln -sf ../3-kernel/bsexmat 4-absorption\n'
                )

    helpers.make_executable(file)


def create_kgrid_in(config, dirname='.'):
    """Create a kgrid.in input file following the config."""
    file_in = os.path.join(dirname, 'kgrid.in')
    nat = str(helpers.num_lines(config['ATOMIC_POSITIONS']['value']))
    k_points = config['K_POINTS']['value'].strip().split()
    nk = ' '.join(k_points[:3])
    dk = ' '.join([str(float(k) * 0.5) for k in k_points[3:]])
    atomic_species = config['ATOMIC_SPECIES']['value'].strip().splitlines()
    elements = [s.split()[0] for s in atomic_species]
    positions = ''

    for atom in config['ATOMIC_POSITIONS']['value'].strip().splitlines():
        atom_split = atom.split()
        element = atom_split[0]
        positions += '{:d} {}\n'.format(elements.index(element) + 1,
                                        ' '.join(atom_split[1:]))

    with open(file_in, 'a') as f:
        f.write(nk + '\n')
        f.write(dk + '\n')
        f.write(config['kgrid']['q-shift'].strip() + '\n\n')
        f.write(config['kgrid']['cell'].strip() + '\n')
        f.write(nat + '\n')
        f.write(positions)
        f.write('20 20 20\n')
        f.write('.false.\n')


def create_epsilon(config, dirname='.'):
    """Create 1-epsilon directory and its input files."""
    dirpath = os.path.join(dirname, '1-epsilon')
    get_kgrid = os.path.join(dirpath, 'get-kgrid')
    inp = os.path.join(dirpath, 'epsilon.inp')
    clean = os.path.join(dirpath, 'clean')
    k_points = config['K_POINTS']['value'].strip().split()
    override = {
        'epsilon': {
            'number_bands': str(int(config['&system']['nbnd']) - 1),
        },
    }
    kgrid_override = {
        'K_POINTS': {
            'value': '{} 0 0 0'.format(' '.join(k_points[:3]))
        },
        'kgrid': {
            'q-shift': '0.0 0.0 0.0',
        },
    }

    os.makedirs(dirpath)

    create_kgrid_in(helpers.deep_merge(config, kgrid_override), dirpath)

    with open(get_kgrid, 'w') as f:
        f.write('#!/bin/bash\n'
                'kgrid.x kgrid.in kgrid.out kgrid.log\n'
                )

    with open(inp, 'a') as f:
        f.write(input_block(helpers.deep_merge(config, override), 'epsilon'))
        f.write('\nbegin qpoints\n')
        f.write(config['kgrid']['q-shift'] + ' 1.0 1\n')
        f.write('end\n')

    with open(clean, 'w') as f:
        f.write('#!/bin/bash\n'
                'rm chi_converge.dat eps0mat epsmat {kgrid,epsilon}.log *.out '
                '2> /dev/null\n'
                )

    helpers.make_executable(get_kgrid)
    helpers.make_executable(clean)


def create_sigma(config, dirname='.'):
    """Create 2.1-sigma directory and its input files."""
    dirpath = os.path.join(dirname, '2.1-sigma')
    inp = os.path.join(dirpath, 'sigma.inp')
    clean = os.path.join(dirpath, 'clean')
    override = {
        'sigma': {
            'number_bands': str(int(config['&system']['nbnd']) - 1),
            'band_index_min': config['pp_in']['vxc_diag_nmin'],
            'band_index_max': config['pp_in']['vxc_diag_nmax'],
        },
    }

    os.makedirs(dirpath)

    with open(inp, 'a') as f:
        f.write(input_block(helpers.deep_merge(config, override), 'sigma'))
        f.write('\nbegin kpoints\n')
        f.write('end\n')

    with open(clean, 'w') as f:
        f.write('#!/bin/bash\n'
                'rm eqp*.dat chi*.dat x.dat *.{log,out} ch_converge.dat 2> '
                '/dev/null\n'
                )

    helpers.make_executable(clean)


def create_inteqp(config, dirname='.'):
    """Create 2.2-inteqp directory and its input files."""
    dirpath = os.path.join(dirname, '2.2-inteqp')
    inp = os.path.join(dirpath, 'inteqp.inp')
    clean = os.path.join(dirpath, 'clean')

    os.makedirs(dirpath)

    with open(inp, 'a') as f:
        f.write(input_block(config, 'inteqp'))

    with open(clean, 'w') as f:
        f.write('#!/bin/bash\n'
                'rm chi*.dat x.dat *.{log,out} ch_converge.dat 2> /dev/null\n'
                )

    helpers.make_executable(clean)


def create_kernel(config, dirname='.'):
    """Create 3-kernel directory and its input files."""
    dirpath = os.path.join(dirname, '3-kernel')
    inp = os.path.join(dirpath, 'kernel.inp')
    clean = os.path.join(dirpath, 'clean')

    os.makedirs(dirpath)

    with open(inp, 'a') as f:
        f.write(input_block(config, 'kernel'))

    with open(clean, 'w') as f:
        f.write('#!/bin/bash\n'
                'rm bse* *.log *.out 2> /dev/null\n'
                )

    helpers.make_executable(clean)


def create_absorption(config, dirname='.'):
    """Create 4-absorption directory and its input files."""
    dirpath = os.path.join(dirname, '4-absorption')
    inp = os.path.join(dirpath, 'absorption.inp')
    clean = os.path.join(dirpath, 'clean')
    override = {
        'absorption': {
            'number_val_bands_coarse':
                config['inteqp']['number_val_bands_coarse'],
            'number_cond_bands_coarse':
                config['inteqp']['number_cond_bands_coarse'],
            'number_val_bands_fine':
                config['inteqp']['number_val_bands_fine'],
            'number_cond_bands_fine':
            config['inteqp']['number_cond_bands_fine'],
        },
    }

    os.makedirs(dirpath)

    with open(inp, 'a') as f:
        f.write(input_block(helpers.deep_merge(config, override),
                            'absorption'))

    with open(clean, 'w') as f:
        f.write('#!/bin/bash\n'
                'rm *.log *.out absorption_noeh.dat bandstructure.dat '
                'dtmat* \\\n'
                '  eigenvalues_noeh.dat dcmat_norm.dat dvmat_norm.dat '
                'epsdiag.dat vmtxel \\\n'
                '  eqp.dat eqp_q.dat absorption_eh.dat eigenvalues.dat '
                'eigenvectors 2> /dev/null\n'
                )

    helpers.make_executable(clean)


def create_bgw(config, dirname='.'):
    dirpath = os.path.join(dirname, 'bgw')

    os.makedirs(dirpath)
    create_link_files(config, dirpath)
    create_epsilon(config, dirpath)
    create_sigma(config, dirpath)
    create_inteqp(config, dirpath)
    create_kernel(config, dirpath)
    create_absorption(config, dirpath)
