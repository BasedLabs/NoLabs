import os
import time
import json
from pathlib import Path


def wait_until_file_generation(path, interval_sec=1, maximum_sec=None) -> bool:
    """Function waits until a given file is created or a maximum duration is reached. Returns "True" if file
       got created and False if maximum duration is exceeded. Potentially runs forever if "maximum_sec" is "None"."""
    counter = 0
    while not os.path.exists(path):
        time.sleep(interval_sec)
        counter = counter + 1
        if maximum_sec is not None and counter * interval_sec >= maximum_sec:
            break
    if os.path.exists(path):
        return True
    else:
        return False


def move_up_directory(path, n=1):
    """Function, to move up 'n' directories for a given "path"."""
    # add +1 to take file into account
    if os.path.isfile(path):
        n += 1
    for _ in range(n):
        path = os.path.dirname(os.path.abspath(path))
    return path


def attach_root_path(path):
    """Function to attach the root path of the module for a given "path"."""
    ROOT_DIR = move_up_directory(os.path.abspath(__file__), n=2)
    return os.path.join(ROOT_DIR, path)


def lines_in_file(path):
    with open(path) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def dict_from_json_file(path):
    with open(path, 'r') as f:
        return json.load(f)


def any_in_file(path, strings):
    if isinstance(strings, str):
        strings = [strings]
    if os.path.isfile(path):
        with open(path, 'r') as f:
            file_raw = f.readlines()
            for string in strings:
                if any(string in line for line in file_raw):
                    return True
            return False
    else:
        return False


def generate_folder_structure(filepath: str):
    folder_path = os.path.dirname(filepath)
    Path(folder_path).mkdir(parents=True, exist_ok=True)
