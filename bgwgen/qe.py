import os
from . import helpers


def namelist_block(config, namelist):
    """Return the namelist block as a string."""
    if namelist not in config:
        return ''

    block = namelist + '\n'

    for key in config[namelist]:
        if config[namelist][key]:
            value = config[namelist][key].strip()
            block += '{} = {}\n'.format(key, value)

    block += '/\n'

    return block


def card_block(config, card):
    """Return the card block as a string."""
    if card not in config:
        return ''

    block = card

    if 'option' in config[card]:
        block += ' ' + config[card]['option']

    block += '\n'
    value = config[card]['value'].strip()

    if value:
        block += value + '\n'

    return block


def create_link_files(config, dirname='.'):
    """Create a 'link-files' executable script."""
    file = os.path.join(dirname, 'link-files')

    with open(file, 'a') as f:
        f.write('#!/bin/bash\n')
        f.write('SEED={}\n\n'.format(config['&control']['prefix']))

        f.write('for d in {2..7}-*/; do\n'
                '  cd $d\n'
                '  rm -rf ${SEED}.save\n'
                '  mkdir -p ${SEED}.save\n'
                '  ln -sf ../../1-scf/${SEED}.save/data-file.xml '
                '${SEED}.save\n'
                '  ln -sf ../../1-scf/${SEED}.save/charge-density.dat '
                '${SEED}.save\n'
                '  cd ..\n'
                'done\n\n'
                'for d in {3..7}-*/; do\n'
                '  cd $d\n'
                '  ln -sf ../2-wfn/clean clean\n'
                '  ln -sf ../2-wfn/get-kgrid get-kgrid\n'
                '  cd ..\n'
                'done\n\n'
                'rm 7-bands/get-kgrid\n')

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


def create_in(config, dirname='.'):
    """Create an 'in' input file following the config."""
    file_in = os.path.join(dirname, 'in')
    nat = str(helpers.num_lines(config['ATOMIC_POSITIONS']['value']))
    ntyp = str(helpers.num_lines(config['ATOMIC_SPECIES']['value']))

    if '&system' in config:
        config['&system']['nat'] = nat
        config['&system']['ntyp'] = ntyp

    with open(file_in, 'a') as f:
        f.write(namelist_block(config, '&control'))
        f.write(namelist_block(config, '&system'))
        f.write(namelist_block(config, '&electrons'))
        f.write(card_block(config, 'CELL_PARAMETERS'))
        f.write(card_block(config, 'ATOMIC_SPECIES'))
        f.write(card_block(config, 'ATOMIC_POSITIONS'))
        f.write(card_block(config, 'K_POINTS'))


def create_pp_in(config, dirname='.', wfng_kgrid=True):
    """Create an 'pp_in' input file following the config."""
    file_in = os.path.join(dirname, 'pp_in')
    k_points = config['K_POINTS']['value'].strip().split()
    q_shift = config['kgrid']['q-shift'].strip().split()
    nk = k_points[:3]
    dk = [float(k) * 0.5 for k in k_points[3:]]
    dk = [str(dk[i] + float(q_shift[i]) * float(nk[i])) for i in range(3)]

    with open(file_in, 'a') as f:
        f.write('&input_pw2bgw\n')
        f.write('prefix = {}\n'.format(config['&control']['prefix']))
        f.write('wfng_flag = .true.\n')

        if wfng_kgrid:
            f.write('wfng_kgrid = .true.\n')
            f.write('wfng_nk1 = {}\n'.format(nk[0]))
            f.write('wfng_nk2 = {}\n'.format(nk[1]))
            f.write('wfng_nk3 = {}\n'.format(nk[2]))
            f.write('wfng_dk1 = {}\n'.format(dk[0]))
            f.write('wfng_dk2 = {}\n'.format(dk[1]))
            f.write('wfng_dk3 = {}\n'.format(dk[2]))

        if 'pp_in' in config:
            for key in config['pp_in']:
                if config['pp_in'][key]:
                    f.write('{} = {}\n'.format(key, config['pp_in'][key]))

        f.write('/\n')


def create_scf(config, dirname='.'):
    """Create 1-scf directory and its input files."""
    dirpath = os.path.join(dirname, '1-scf')
    clean = os.path.join(dirpath, 'clean')
    override = {
        '&control': {
            'calculation': '\'scf\'',
        },
        '&system': {
            'nbnd': '',
        },
        'K_POINTS': {
            'option': 'automatic',
        },
    }

    os.makedirs(dirpath)
    create_in(helpers.deep_merge(config, override), dirpath)

    with open(clean, 'w') as f:
        f.write('#!/bin/bash\n'
                'rm -rf CRASH *out* charge* *.{igk,mix,wfc,save}* 2> '
                '/dev/null\n'
                )

    helpers.make_executable(clean)


def create_wfn(config, dirname='.'):
    """Create 2-wfn directory and its input files."""
    dirpath = os.path.join(dirname, '2-wfn')
    get_kgrid = os.path.join(dirpath, 'get-kgrid')
    clean = os.path.join(dirpath, 'clean')
    override = {
        '&control': {
            'calculation': '\'bands\'',
            'wf_collect': '.true.',
        },
        'K_POINTS': {
            'option': 'crystal',
            'value': '',
        },
    }
    kgrid_override = {
        'kgrid': {
            'q-shift': '0.0 0.0 0.0',
        },
    }
    pp_in_config = helpers.deep_merge(config, kgrid_override)
    pp_in_config.pop('pp_in', None)

    os.makedirs(dirpath)
    create_kgrid_in(helpers.deep_merge(config, kgrid_override), dirpath)
    create_pp_in(pp_in_config, dirpath)
    create_in(helpers.deep_merge(config, override), dirpath)

    with open(get_kgrid, 'w') as f:
        f.write('#!/bin/bash\n'
                'kgrid.x kgrid.in kgrid.out kgrid.log\n'
                'sed -n 2,99999p kgrid.out >> in\n')

    with open(clean, 'w') as f:
        f.write('#!/bin/bash\n'
                'rm -rf CRASH *{paw,igk,wfc,log,out,.xmgr}* MoS2.save/K00* '
                'MoS2.save/*.UPF RHO \\\n'
                '  WFN vxc.dat *.ps bands.dat* 2> /dev/null\n'
                )

    helpers.make_executable(get_kgrid)
    helpers.make_executable(clean)


def create_wfnq(config, dirname='.'):
    """Create 3-wfnq directory and its input files."""
    dirpath = os.path.join(dirname, '3-wfnq')
    override = {
        '&control': {
            'calculation': '\'bands\'',
            'wf_collect': '.true.',
        },
        '&system': {
            'nbnd': '',
        },
        'K_POINTS': {
            'option': 'crystal',
            'value': '',
        },
    }
    pp_in_config = helpers.deep_merge(config, {})
    pp_in_config.pop('pp_in', None)

    os.makedirs(dirpath)
    create_kgrid_in(config, dirpath)
    create_pp_in(pp_in_config, dirpath)
    create_in(helpers.deep_merge(config, override), dirpath)


def create_wfn_co(config, dirname='.'):
    """Create 4-wfn_co directory and its input files."""
    dirpath = os.path.join(dirname, '4-wfn_co')
    k_points = config['K_POINTS']['value'].strip().split()[:3]
    override = {
        '&control': {
            'calculation': '\'bands\'',
            'wf_collect': '.true.',
        },
        'K_POINTS': {
            'option': 'crystal',
            'value': '',
        },
    }
    kgrid_override = {
        'K_POINTS': {
            'value': '{} 0 0 0'.format(' '.join([str(k) for k in k_points]))
        },
        'kgrid': {
            'q-shift': '0.0 0.0 0.0',
        },
        'pp_in': {
            'rhog_flag': '.true.',
            'vxc_flag': '.true.',
        },
    }

    os.makedirs(dirpath)
    create_kgrid_in(helpers.deep_merge(config, kgrid_override), dirpath)
    create_pp_in(helpers.deep_merge(config, kgrid_override), dirpath)
    create_in(helpers.deep_merge(config, override), dirpath)


def create_wfn_fi(config, dirname='.'):
    """Create 5-wfn_fi directory and its input files."""
    dirpath = os.path.join(dirname, '5-wfn_fi')
    k_points = config['K_POINTS']['value'].strip().split()
    nk = [int(int(k_points[i]) * 2) if int(k_points[i]) != 1 else
          int(k_points[i]) for i in range(3)]
    override = {
        '&control': {
            'calculation': '\'bands\'',
            'wf_collect': '.true.',
        },
        'K_POINTS': {
            'option': 'crystal',
            'value': '',
        },
    }
    kgrid_override = {
        'K_POINTS': {
            'value': '{} 0 0 0'.format(' '.join([str(k) for k in nk]))
        },
        'kgrid': {
            'q-shift': '0.0 0.0 0.0',
        },
    }
    pp_in_config = helpers.deep_merge(config, kgrid_override)
    pp_in_config.pop('pp_in', None)

    os.makedirs(dirpath)
    create_kgrid_in(helpers.deep_merge(config, kgrid_override), dirpath)
    create_pp_in(pp_in_config, dirpath)
    create_in(helpers.deep_merge(config, override), dirpath)


def create_wfnq_fi(config, dirname='.'):
    """Create 6-wfnq_fi directory and its input files."""
    dirpath = os.path.join(dirname, '6-wfnq_fi')
    k_points = config['K_POINTS']['value'].strip().split()
    q_shift = config['kgrid']['q-shift'].strip().split()

    for q in q_shift:
        if float(q) > 0:
            q_shift = float(q)
            break

    nk = [int(int(k_points[i]) * 2) if int(k_points[i]) != 1 else
          int(k_points[i]) for i in range(3)]
    override = {
        '&control': {
            'calculation': '\'bands\'',
            'wf_collect': '.true.',
        },
        '&system': {
            'nbnd': '',
        },
        'K_POINTS': {
            'option': 'crystal',
            'value': '',
        },
    }
    kgrid_override = {
        'K_POINTS': {
            'value': '{} {}'.format(
                ' '.join([str(k) for k in nk]),
                ' '.join([str(2 * k * q_shift) if k != 1 else
                          '0' for k in nk]))
        },
        'kgrid': {
            'q-shift': '0.0 0.0 0.0',
        },
    }
    pp_in_config = helpers.deep_merge(config, kgrid_override)
    pp_in_config.pop('pp_in', None)

    os.makedirs(dirpath)
    create_kgrid_in(helpers.deep_merge(config, kgrid_override), dirpath)
    create_pp_in(pp_in_config, dirpath)
    create_in(helpers.deep_merge(config, override), dirpath)


def create_bands(config, dirname='.'):
    """Create 7-bands directory and its input files."""
    dirpath = os.path.join(dirname, '7-bands')
    override = {
        '&control': {
            'calculation': '\'bands\'',
            'wf_collect': '.true.',
        },
        'K_POINTS': config['K_POINTS_bands'],
    }
    in_config = helpers.deep_merge(config, override)
    in_config.pop('pp_in', None)

    os.makedirs(dirpath)
    create_pp_in(in_config, dirpath, False)
    create_in(in_config, dirpath)


def create_qe(config, dirname='.'):
    """Create a new directory '1-qe' and all its directories."""
    dirpath = os.path.join(dirname, '1-qe')

    os.makedirs(dirpath)
    create_link_files(config, dirpath)
    create_scf(config, dirpath)
    create_wfn(config, dirpath)
    create_wfnq(config, dirpath)
    create_wfn_co(config, dirpath)
    create_wfn_fi(config, dirpath)
    create_wfnq_fi(config, dirpath)
    create_bands(config, dirpath)
