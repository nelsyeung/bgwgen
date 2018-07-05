import os
from bgwgen import parser

fixtures_dir = os.path.join('tests', 'fixtures')


def test_parse():
    """parse function parses all settings."""
    input_file = os.path.join(fixtures_dir, 'bgwgen.ini')
    config = parser.parse(input_file)
    expected = {
        '&control': {
            'prefix': 'MoS2',
            'pseudo_dir': '../pseudo',
        },
        '&system': {
            'ibrav': '4',
            'celldm(1)': '3.169',
            'celldm(3)': '20.0',
            'ecutwfc': '45.0',
        },
        'ATOMIC_SPECIES': {
            'value': ('\nMo   95.95  Mo.UPF\n'
                      'S    32.06  S.UPF'
                      ),
        },
        'epsilon': {
            'cell_slab_truncation': '',
        },
    }

    assert config == expected
    for key in config:
        for c_key in config[key]:
            assert config[key][c_key] == expected[key][c_key]
