import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--nch', type=int, default=16)

args = parser.parse_args()
print(args.nch)
