import argparse


def str2bool(inp):
    if isinstance(inp, bool):
        return inp
    if inp.lower() in ("yes", "true", 't', 'y', '1'):
        return True
    elif inp.lower() in ("no", "false", 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError("Expected castable string or boolean value as input.")
