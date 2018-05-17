import copy
import os
import stat


def deep_merge(original, update):
    """Recursively update a dictionary. Subdictionaries will not be overwritten
    but also updated.

    Keyword arguments:
    original -- the original dictionary to be overwritten
    update   -- the dictionary to overwrite the original dictionary
    """
    for key in original:
        if key not in update:
            update[key] = original[key]
        elif isinstance(original[key], dict):
            deep_merge(original[key], update[key])

    return copy.deepcopy(update)


def num_lines(text):
    """Return the number of lines for a multiline text string."""
    return len(text.strip().splitlines())


def make_executable(file):
    os.chmod(file, os.stat(file).st_mode | stat.S_IEXEC)
