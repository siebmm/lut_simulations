"""Validate geometry, acquisition metadata, and every velocity-history CSV."""
import argparse
import json
from pathlib import Path
from dataset_io import load_dataset


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('dataset', type=Path, help='Directory containing dataset.json')
    args = parser.parse_args()
    try:
        *_, report = load_dataset(args.dataset)
    except (ValueError, KeyError, OSError, StopIteration) as error:
        parser.exit(1, f'Invalid dataset: {error}\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
