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
                'QE="../../1-qe"\n'
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
    setup = os.path.join(dirpath, '0-setup.sh')
    qpoints = os.path.join(dirpath, 'qpoints')

    os.makedirs(dirpath)

    with open(setup, 'a') as f:
        f.writelines([
            '#!/bin/bash\n'
            'num_kp=$(cat qpoints | wc -l)\n',
            '\n',
            '# Create epsilon.inp for every kpoints inside the qpoints file\n',
            'for i in $(seq 1 $num_kp); do\n',
            '	dir="eps$(seq -f "%02g" $i $i)"\n',
            '\n',
            '	if [[ -z $1 ]]; then\n',
            '		mkdir ${dir}\n',
            '		cd $dir\n',
            '\n',
            '		cat > epsilon.inp <<- EOM\n',
        ])
        f.writelines([
            '\t\t\t{}\n'.format(x) for x in input_block(
                config, 'epsilon').split('\n') if x.strip()
        ])
        f.writelines([
            '\n',
            '			begin qpoints\n',
            '			$(sed -n ${i}p ../qpoints)\n',
            '			end\n',
            '		EOM\n',
            '\n',
            '		ln -s ../WFN .\n',
            '		ln -s ../WFNq .\n',
            '\n',
            '		cd ..\n',
            '	elif [[ $1 == "clean" ]]; then\n',
            '		rm -rf ${dir}\n',
            '	fi\n',
            'done\n',
            '\n',
            '# Create an epsmat merge folder\n',
            'if [[ -z $1 ]]; then\n',
            '	nkp=$((num_kp-1))\n',
            '\n',
            '	mkdir merge\n',
            '	cd merge\n',
            '\n',
            '	echo "{} $nkp" > epsmat_merge.inp\n'.format(
                config['epsilon']['epsilon_cutoff']),
            '\n',
            '	for i in $(seq 2 $num_kp); do\n',
            '		kpoint=$(sed -n ${i}p ../qpoints)\n',
            '		echo "${kpoint%?}" >> epsmat_merge.inp\n',
            '\n',
            '		dir="eps$(seq -f "%02g" $i $i)"\n',
            '		epsmat="epsmat$(seq -f "%02g" $i $i)"\n',
            '		ln -s ../$dir/epsmat $epsmat\n',
            '	done\n',
            '\n',
            '	echo "$nkp" >> epsmat_merge.inp\n',
            '\n',
            '	for i in $(seq 2 $num_kp); do\n',
            '		epsmat="epsmat$(seq -f "%02g" $i $i)"\n',
            '		echo "$epsmat 1" >> epsmat_merge.inp\n',
            '	done\n',
            '	cd ..\n',
            'elif [[ $1 == "clean" ]]; then\n',
            '	rm -rf merge\n',
            'fi\n',
        ])

    with open(qpoints, 'a') as f:
        f.write('# Replace this file with all the qpoints for epsilon.inp\n')

    helpers.make_executable(setup)


def create_sigma(config, dirname='.'):
    """Create 2-sigma directory and its input files."""
    dirpath = os.path.join(dirname, '2-sigma')
    setup = os.path.join(dirpath, '0-setup.sh')
    kpoints = os.path.join(dirpath, 'kpoints')
    merge = os.path.join(dirpath, '2-merge.sh')
    override = {
        'sigma': {
            'band_index_min': config['pp_in']['vxc_diag_nmin'],
            'band_index_max': config['pp_in']['vxc_diag_nmax'],
        },
    }
    config = helpers.deep_merge(config, override)

    os.makedirs(dirpath)

    with open(setup, 'a') as f:
        f.writelines([
            '#!/bin/bash\n'
            'num_kp=$(cat kpoints | wc -l)\n',
            '\n',
            '# Create sigma.inp for every kpoints inside the kpoints file\n',
            'for i in $(seq 1 $num_kp); do\n',
            '	dir="sig$(seq -f "%02g" $i $i)"\n',
            '\n',
            '	if [[ -z $1 ]]; then\n',
            '		mkdir ${dir}\n',
            '		cd $dir\n',
            '\n',
            '		cat > sigma.inp <<- EOM\n',
        ])
        f.writelines([
            '\t\t\t{}\n'.format(x) for x in input_block(
                config, 'sigma').split('\n') if x.strip()
        ])
        f.writelines([
            '\n',
            '			begin kpoints\n',
            '			$(sed -n ${i}p ../kpoints)\n',
            '			end\n',
            '		EOM\n',
            '\n',
            '		ln -s ../RHO .\n',
            '		ln -s ../WFN_inner .\n',
            '		ln -s ../eps0mat .\n',
            '		ln -s ../epsmat .\n',
            '		ln -s ../vxc.dat .\n',
            '\n',
            '		cd ..\n',
            '	elif [[ $1 == "clean" ]]; then\n',
            '		rm -rf ${dir}\n',
            '	fi\n',
            'done\n',
        ])

    with open(kpoints, 'a') as f:
        f.write('# Replace this file with all the kpoints for sigma.inp\n')

    with open(merge, 'a') as f:
        f.writelines([
            '#!/bin/bash\n',
            'num_kp=$(cat kgrid | wc -l)\n',
            '\n',
            'for i in $(seq 1 $num_kp); do\n',
            '	dir="sig$(seq -f "%02g" $i $i)"\n',
            '	cat $dir/eqp0.dat >> eqp0.dat\n',
            '	cat $dir/eqp1.dat >> eqp1.dat\n',
            'done\n',
        ])

    helpers.make_executable(setup)
    helpers.make_executable(merge)


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
            'number_val_bands_coarse': config['kernel']['number_val_bands'],
            'number_cond_bands_coarse': config['kernel']['number_cond_bands'],
        },
    }
    config = helpers.deep_merge(config, override)

    os.makedirs(dirpath)

    with open(inp, 'a') as f:
        f.write(input_block(config, 'absorption'))

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


def create_inteqp(config, dirname='.'):
    """Create 5-inteqp directory and its input files."""
    dirpath = os.path.join(dirname, '5-inteqp')
    inp = os.path.join(dirpath, 'inteqp.inp')
    clean = os.path.join(dirpath, 'clean')

    os.makedirs(dirpath)

    with open(inp, 'a') as f:
        f.write(input_block(config, 'inteqp'))

    with open(clean, 'w') as f:
        f.write('#!/bin/bash\n'
                'rm *.log *.out bandstructure.dat dtmat* dcmat_norm.dat '
                'dvmat_norm.dat \\\n'
                '  eqp.dat eqp_q.dat 2> /dev/null\n'
                )

    helpers.make_executable(clean)


def create_bgw(config, dirname='.'):
    """Create a new directory '2-bgw' and all its directories."""
    dirpath = os.path.join(dirname, '2-bgw')

    os.makedirs(dirpath)
    create_link_files(config, dirpath)
    create_epsilon(config, dirpath)
    create_sigma(config, dirpath)
    create_kernel(config, dirpath)
    create_absorption(config, dirpath)
    create_inteqp(config, dirpath)
