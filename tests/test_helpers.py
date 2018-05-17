import os
import stat
from bgwgen import helpers


def test_deep_merge():
    """deep_merge function returns a single deep merged dictionary."""
    original = {
        '&control': {
            'calculation': 'bands',
            'prefix': 'MoS2',
            'pseudo_dir': '../pseudo',
            'wf_collect': '.true.',
        },
        'ATOMIC_SPECIES': {
            'type': 'angstrom',
            'value': '\nMo   95.95  Mo.UPF\nS    32.06  S.UPF',
        },
    }
    override = {
        '&control': {
            'calculation': 'scf',
            'tprnfor': '.true.',
            'tstress': '.true.',
        }
    }
    expected = {
        '&control': {
            'calculation': 'scf',
            'prefix': 'MoS2',
            'pseudo_dir': '../pseudo',
            'wf_collect': '.true.',
            'tprnfor': '.true.',
            'tstress': '.true.',
        },
        'ATOMIC_SPECIES': {
            'type': 'angstrom',
            'value': '\nMo   95.95  Mo.UPF\nS    32.06  S.UPF',
        },
    }

    assert helpers.deep_merge(original, override) == expected


def test_num_lines():
    """num_lines function returns the number of lines in a multiline string."""
    text = ('\nMo   0.914811501530962   1.5845   3.081\n' +
            'Mo   3.659246006123848   0.0000   3.081000\n' +
            'S    1.829623003061924   0.0000   1.515852\n' +
            'S    1.829623003061924   0.0000   4.646148\n' +
            'S    4.574057507654810   1.5845   1.515852\n' +
            'S    4.574057507654810   1.5845   4.646148\n')
    assert helpers.num_lines(text) == 6


def test_make_executable(tmpdir):
    file = tmpdir.join('executable')
    with open(file, 'a') as f:
        f.write('#!/bin/bash')
    helpers.make_executable(file)
    assert stat.S_IXUSR & os.stat(file)[stat.ST_MODE]
