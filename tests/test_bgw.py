import os
from bgwgen import bgw

fixtures_dir = os.path.join('tests', 'fixtures', 'bgw')

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
    'epsilon': {
        'epsilon_cutoff': '5.0',
        'cell_slab_truncation': '',
    },
    'sigma': {
        'screened_coulomb_cutoff': '5.0',
        'bare_coulomb_cutoff': '45.0',
        'cell_slab_truncation': '',
        'screening_semiconductor': '',
    },
    'kernel': {
        'screened_coulomb_cutoff': '35.0',
        'bare_coulomb_cutoff': '100.0',
        'number_val_bands': '14',
        'number_cond_bands': '14',
        'cell_slab_truncation': '',
        'use_symmetries_coarse_grid': '',
        'screening_semiconductor': '',
    },
    'absorption': {
        'number_val_bands_fine': '2',
        'number_cond_bands_fine': '2',
        'diagonalization': '',
        'use_symmetries_coarse_grid': '',
        'use_symmetries_fine_grid': '',
        'no_symmetries_shifted_grid': '',
        'screening_semiconductor': '',
        'use_velocity': '',
        'cell_slab_truncation': '',
        'gaussian_broadening': '',
        'energy_resolution': '0.05',
        'eqp_co_corrections': '',
        'write_eigenvectors': '10',
    },
    'inteqp': {
        'number_val_bands_coarse': '14',
        'number_cond_bands_coarse': '14',
        'number_val_bands_fine': '2',
        'number_cond_bands_fine': '2',
        'use_symmetries_coarse_grid': '',
        'no_symmetries_fine_grid': '',
        'no_symmetries_shifted_grid': '',
    },
}


def test_input_block():
    """Returns the correct string for the epsilon input."""
    expected = ('epsilon_cutoff 5.0\n'
                'cell_slab_truncation\n'
                )
    assert bgw.input_block(config, 'epsilon') == expected
    assert bgw.input_block(config, 'foo') == ''


def test_create_link_files(tmpdir):
    """Creates a 'link-files' executable bash script."""
    bgw.create_link_files(config, tmpdir.realpath())

    with open(os.path.join(fixtures_dir, 'link-files.expected'), 'r') as f:
        assert tmpdir.join('link-files').read() == f.read()


def test_create_epsilon(tmpdir):
    """Creates a new directory '1-epsilon' and all its input files."""
    dirname = '1-epsilon'
    d = tmpdir.join(dirname)
    expected_dir = os.path.join(fixtures_dir, dirname)

    bgw.create_epsilon(config, tmpdir.realpath())

    with open(os.path.join(expected_dir, 'epsilon.inp.expected'), 'r') as f:
        assert d.join('epsilon.inp').read() == f.read()

    with open(os.path.join(expected_dir, 'clean.expected'), 'r') as f:
        assert d.join('clean').read() == f.read()


def test_create_sigma(tmpdir):
    """Creates a new directory '2-sigma' and all its input files."""
    dirname = '2-sigma'
    d = tmpdir.join(dirname)
    expected_dir = os.path.join(fixtures_dir, dirname)

    bgw.create_sigma(config, tmpdir.realpath())

    with open(os.path.join(expected_dir, 'sigma.inp.expected'), 'r') as f:
        assert d.join('sigma.inp').read() == f.read()

    with open(os.path.join(expected_dir, 'clean.expected'), 'r') as f:
        assert d.join('clean').read() == f.read()


def test_create_kernel(tmpdir):
    """Creates a new directory '3-kernel' and all its input files."""
    dirname = '3-kernel'
    d = tmpdir.join(dirname)
    expected_dir = os.path.join(fixtures_dir, dirname)

    bgw.create_kernel(config, tmpdir.realpath())

    with open(os.path.join(expected_dir, 'kernel.inp.expected'), 'r') as f:
        assert d.join('kernel.inp').read() == f.read()

    with open(os.path.join(expected_dir, 'clean.expected'), 'r') as f:
        assert d.join('clean').read() == f.read()


def test_create_absorption(tmpdir):
    """Creates a new directory '4-absorption' and all its input files."""
    dirname = '4-absorption'
    d = tmpdir.join(dirname)
    expected_dir = os.path.join(fixtures_dir, dirname)

    bgw.create_absorption(config, tmpdir.realpath())

    with open(os.path.join(expected_dir, 'absorption.inp.expected'), 'r') as f:
        assert d.join('absorption.inp').read() == f.read()

    with open(os.path.join(expected_dir, 'clean.expected'), 'r') as f:
        assert d.join('clean').read() == f.read()


def test_create_inteqp(tmpdir):
    """Creates a new directory '5-inteqp' and all its input files."""
    dirname = '5-inteqp'
    d = tmpdir.join(dirname)
    expected_dir = os.path.join(fixtures_dir, dirname)

    bgw.create_inteqp(config, tmpdir.realpath())

    with open(os.path.join(expected_dir, 'inteqp.inp.expected'), 'r') as f:
        assert d.join('inteqp.inp').read() == f.read()

    with open(os.path.join(expected_dir, 'clean.expected'), 'r') as f:
        assert d.join('clean').read() == f.read()


def test_create_bgw(tmpdir):
    """Create a new directory '2-bgw' and all its directories."""
    bgwdir = tmpdir.join('2-bgw')

    bgw.create_bgw(config, tmpdir.realpath())

    assert os.path.isdir(bgwdir)
    assert os.path.isfile(bgwdir.join('link-files'))
    assert os.path.isdir(bgwdir)
    assert os.path.isdir(bgwdir.join('1-epsilon'))
    assert os.path.isdir(bgwdir.join('2-sigma'))
    assert os.path.isdir(bgwdir.join('3-kernel'))
    assert os.path.isdir(bgwdir.join('4-absorption'))
    assert os.path.isdir(bgwdir.join('5-inteqp'))
