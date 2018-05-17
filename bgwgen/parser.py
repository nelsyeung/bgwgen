import configparser


def parse(file):
    config = configparser.ConfigParser()
    config.read(file)
    c = {}

    for section in config.sections():
        c[section] = dict(config.items(section))

    return c
