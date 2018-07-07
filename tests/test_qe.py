import copy
import os
from bgwgen import qe

fixtures_dir = os.path.join('tests', 'fixtures', '1-qe')

config = {
    '&control': {
        'prefix': '\'MoS2\'',
        'pseudo_dir': '\'../pseudo\'',
    },
    '&system': {
        'ibrav': '4',
        'celldm(1)': '3.169',
        'celldm(3)': '20.0',
        'ecutwfc': '45.0',
        'nbnd': '100',
    },
    'ATOMIC_POSITIONS': {
        'option': 'angstrom',
        'value': ('\n'
                  'Mo  1.5845  0.9148  3.0810\n'
                  'S   0.0000  1.8296  1.5158\n'
                  'S   0.0000  1.8296  4.6461\n'
                  )
    },
    'ATOMIC_SPECIES': {
        'value': ('\n'
                  'Mo   95.95  Mo.UPF\n'
                  'S    32.06  S.UPF'
                  )
    },
    'K_POINTS': {
        'value': '12 12 1 1 1 0',
    },
    'K_POINTS_bands': {
        'option': 'crystal_b',
        'value': ('\n'
                  '4\n'
                  '0.000 0.000 0.000  50\n'
                  '0.500 0.000 0.000  50\n'
                  '0.333 0.333 0.000  50\n'
                  '0.000 0.000 0.000  0'
                  )
    },
    'kgrid': {
        'q-shift': '0.001 0.0 0.0',
        'cell': ('1.0 0.0 0.0\n'
                 '0.0 1.0 0.0\n'
                 '0.0 0.0 1.0'
                 )
    },
    'pp_in': {
        'vxc_diag_nmin': '1',
        'vxc_diag_nmax': '44',
    },
}


def test_namelist_block():
    """Returns the correct string for the control namelist."""
    expect = ('&control\n'
              'prefix = \'MoS2\'\n'
              'pseudo_dir = \'../pseudo\'\n'
              '/\n'
              )
    assert qe.namelist_block(config, '&control') == expect


def test_card_block():
    """Returns the correct string for the specified card or an empty string if
    card does not exists in the config."""
    card = 'ATOMIC_POSITIONS'
    expected = '{} angstrom\n{}\n'.format(card, config[card]['value'].strip())

    assert qe.card_block(config, card) == expected
    assert qe.card_block({}, 'foo') == ''


def test_create_link_files(tmpdir):
    """Creates a 'link-files' executable bash script."""
    qe.create_link_files(config, tmpdir.realpath())

    with open(os.path.join(fixtures_dir, 'link-files.expected'), 'r') as f:
        assert tmpdir.join('link-files').read() == f.read()


def test_create_kgrid_in(tmpdir):
    """Creates an 'kgrid.in' input file with the correct config."""
    qe.create_kgrid_in(config, tmpdir.realpath())

    with open(os.path.join(fixtures_dir, 'kgrid.in.expected'), 'r') as f:
        assert tmpdir.join('kgrid.in').read() == f.read()


def test_create_in(tmpdir):
    """Creates an 'in' input file with the correct config."""
    c = copy.deepcopy(config)
    c['K_POINTS']['option'] = 'automatic'

    qe.create_in(c, tmpdir.realpath())

    with open(os.path.join(fixtures_dir, 'create_in.expected'), 'r') as f:
        assert tmpdir.join('in').read() == f.read()


def test_create_pp_in(tmpdir):
    """Creates an 'pp_in' input file with the correct config."""
    qe.create_pp_in(config, tmpdir.realpath())

    with open(os.path.join(fixtures_dir, 'create_pp_in.expected'), 'r') as f:
        assert tmpdir.join('pp_in').read() == f.read()


def test_create_scf(tmpdir):
    """Creates a new directory '1-scf' and all its input files."""
    dirname = '1-scf'
    d = tmpdir.join(dirname)
    expected_dir = os.path.join(fixtures_dir, dirname)

    qe.create_scf(config, tmpdir.realpath())

    with open(os.path.join(expected_dir, 'in.expected'), 'r') as f:
        assert d.join('in').read() == f.read()

    with open(os.path.join(expected_dir, 'clean.expected'), 'r') as f:
        assert d.join('clean').read() == f.read()


def test_create_wfn(tmpdir):
    """Creates a new directory '2-wfn' and all its input files."""
    dirname = '2-wfn'
    d = tmpdir.join(dirname)
    expected_dir = os.path.join(fixtures_dir, dirname)

    qe.create_wfn(config, tmpdir.realpath())

    with open(os.path.join(expected_dir, 'kgrid.in.expected'), 'r') as f:
        assert d.join('kgrid.in').read() == f.read()

    with open(os.path.join(expected_dir, 'in.expected'), 'r') as f:
        assert d.join('in').read() == f.read()

    with open(os.path.join(expected_dir, 'pp_in.expected'), 'r') as f:
        assert d.join('pp_in').read() == f.read()

    with open(os.path.join(expected_dir, 'get-kgrid.expected'), 'r') as f:
        assert d.join('get-kgrid').read() == f.read()

    with open(os.path.join(expected_dir, 'clean.expected'), 'r') as f:
        assert d.join('clean').read() == f.read()


def test_create_wfnq(tmpdir):
    """Creates a new directory '3-wfnq' and all its input files."""
    dirname = '3-wfnq'
    d = tmpdir.join(dirname)
    expected_dir = os.path.join(fixtures_dir, dirname)

    qe.create_wfnq(config, tmpdir.realpath())

    with open(os.path.join(expected_dir, 'kgrid.in.expected'), 'r') as f:
        assert d.join('kgrid.in').read() == f.read()

    with open(os.path.join(expected_dir, 'in.expected'), 'r') as f:
        assert d.join('in').read() == f.read()

    with open(os.path.join(expected_dir, 'pp_in.expected'), 'r') as f:
        assert d.join('pp_in').read() == f.read()


def test_create_wfn_co(tmpdir):
    """Creates a new directory '4-wfn_co' and all its input files."""
    dirname = '4-wfn_co'
    d = tmpdir.join(dirname)
    expected_dir = os.path.join(fixtures_dir, dirname)

    qe.create_wfn_co(config, tmpdir.realpath())

    with open(os.path.join(expected_dir, 'kgrid.in.expected'), 'r') as f:
        assert d.join('kgrid.in').read() == f.read()

    with open(os.path.join(expected_dir, 'in.expected'), 'r') as f:
        assert d.join('in').read() == f.read()

    with open(os.path.join(expected_dir, 'pp_in.expected'), 'r') as f:
        assert d.join('pp_in').read() == f.read()


def test_create_wfn_fi(tmpdir):
    """Creates a new directory '5-wfn_fi' and all its input files."""
    dirname = '5-wfn_fi'
    d = tmpdir.join(dirname)
    expected_dir = os.path.join(fixtures_dir, dirname)

    qe.create_wfn_fi(config, tmpdir.realpath())

    with open(os.path.join(expected_dir, 'kgrid.in.expected'), 'r') as f:
        assert d.join('kgrid.in').read() == f.read()

    with open(os.path.join(expected_dir, 'in.expected'), 'r') as f:
        assert d.join('in').read() == f.read()

    with open(os.path.join(expected_dir, 'pp_in.expected'), 'r') as f:
        assert d.join('pp_in').read() == f.read()


def test_create_wfnq_fi(tmpdir):
    """Creates a new directory '6-wfnq_fi' and all its input files."""
    dirname = '6-wfnq_fi'
    d = tmpdir.join(dirname)
    expected_dir = os.path.join(fixtures_dir, dirname)

    qe.create_wfnq_fi(config, tmpdir.realpath())

    with open(os.path.join(expected_dir, 'kgrid.in.expected'), 'r') as f:
        assert d.join('kgrid.in').read() == f.read()

    with open(os.path.join(expected_dir, 'in.expected'), 'r') as f:
        assert d.join('in').read() == f.read()

    with open(os.path.join(expected_dir, 'pp_in.expected'), 'r') as f:
        assert d.join('pp_in').read() == f.read()


def test_create_bands(tmpdir):
    """Creates a new directory '7-bands' and all its input files."""
    dirname = '7-bands'
    d = tmpdir.join(dirname)
    expected_dir = os.path.join(fixtures_dir, dirname)

    qe.create_bands(config, tmpdir.realpath())

    with open(os.path.join(expected_dir, 'in.expected'), 'r') as f:
        assert d.join('in').read() == f.read()

    with open(os.path.join(expected_dir, 'pp_in.expected'), 'r') as f:
        assert d.join('pp_in').read() == f.read()


def test_create_qe(tmpdir):
    """Create a new directory '1-qe' and all its directories."""
    qedir = tmpdir.join('1-qe')

    qe.create_qe(config, tmpdir.realpath())

    assert os.path.isdir(qedir)
    assert os.path.isfile(qedir.join('link-files'))
    assert os.path.isdir(qedir)
    assert os.path.isdir(qedir.join('1-scf'))
    assert os.path.isdir(qedir.join('2-wfn'))
    assert os.path.isdir(qedir.join('3-wfnq'))
    assert os.path.isdir(qedir.join('4-wfn_co'))
    assert os.path.isdir(qedir.join('5-wfn_fi'))
    assert os.path.isdir(qedir.join('6-wfnq_fi'))
    assert os.path.isdir(qedir.join('7-bands'))
