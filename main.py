#!/usr/bin/env python
import utils
from parameterize import *
from argparse import ArgumentParser


def main():
    parser = ArgumentParser()
    parser.add_argument('-c', '--config', help='The YAML config file to use.', default='config.yaml')
    args = parser.parse_args()

    config = utils.load_config(args.config)

    check_cutoff_biotypes_internal(config)


if __name__ == '__main__':
    main()
