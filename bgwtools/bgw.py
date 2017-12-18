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
                'ln -sf $QE/4-wfn_co/RHO 2-sigma/RHO\n'
                'ln -sf $QE/4-wfn_co/WFN 2-sigma/WFN_inner\n'
                'ln -sf $QE/4-wfn_co/WFN 3-kernel/WFN_co\n'
                'ln -sf $QE/4-wfn_co/WFN 4-absorption/WFN_co\n'
                'ln -sf $QE/4-wfn_co/vxc.dat 2-sigma/vxc.dat\n'
                'ln -sf $QE/5-wfn_fi/WFN 4-absorption/WFN_fi\n'
                'ln -sf $QE/6-wfnq_fi/WFN 4-absorption/WFNq_fi\n'
                '\n'
                'ln -sf ../1-epsilon/eps0mat 2-sigma\n'
                'ln -sf ../1-epsilon/eps0mat 3-kernel\n'
                'ln -sf ../1-epsilon/eps0mat 4-absorption\n'
                'ln -sf ../1-epsilon/epsmat 2-sigma\n'
                'ln -sf ../1-epsilon/epsmat 3-kernel\n'
                'ln -sf ../1-epsilon/epsmat 4-absorption\n'
                'ln -sf ../2-sigma/eqp1.dat 4-absorption/eqp_co.dat\n'
                'ln -sf ../3-kernel/bsedmat 4-absorption\n'
                'ln -sf ../3-kernel/bsexmat 4-absorption\n'
                )

    helpers.make_executable(file)


def create_epsilon(config, dirname='.'):
    """Create 1-epsilon directory and its input files."""
    dirpath = os.path.join(dirname, '1-epsilon')
    inp = os.path.join(dirpath, 'epsilon.inp')
    clean = os.path.join(dirpath, 'clean')

    os.makedirs(dirpath)

    with open(inp, 'a') as f:
        f.write(input_block(config, 'epsilon'))
        f.write('\nbegin qpoints\n')
        f.write(config['kgrid']['q-shift'] + ' 1.0 1\n')
        f.write('end\n')

    with open(clean, 'w') as f:
        f.write('#!/bin/bash\n'
                'rm chi_converge.dat eps0mat epsmat {kgrid,epsilon}.log *.out '
                '2> /dev/null\n'
                )

    helpers.make_executable(clean)


def create_sigma(config, dirname='.'):
    """Create 2-sigma directory and its input files."""
    dirpath = os.path.join(dirname, '2-sigma')
    inp = os.path.join(dirpath, 'sigma.inp')
    clean = os.path.join(dirpath, 'clean')

    os.makedirs(dirpath)

    with open(inp, 'a') as f:
        f.write(input_block(config, 'sigma'))
        f.write('\nbegin kpoints\n')
        f.write('end\n')

    with open(clean, 'w') as f:
        f.write('#!/bin/bash\n'
                'rm eqp*.dat chi*.dat x.dat *.{log,out} ch_converge.dat 2> '
                '/dev/null\n'
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
    create_kernel(config, dirpath)
    create_absorption(config, dirpath)
