#!/usr/bin/env python3
"""Script to demonstrate vaspy.incar functionality."""
import argparse
from logging import DEBUG, Formatter, StreamHandler, getLogger
from pathlib import Path

import vaspy
import vaspy.incar

LOGLEVEL = DEBUG
logger = getLogger(__name__)
fmt = "%(asctime)s %(levelname)s %(name)s :%(message)s"
formatter = Formatter(fmt)
handler = StreamHandler()
handler.setLevel(LOGLEVEL)
logger.setLevel(LOGLEVEL)
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.propagate = False


parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument(
    "-r",
    help="""Show reformatted INCAR (Use -i if edit in place)""",
    action="store_true",
)
parser.add_argument("-i", help="""Edit the INCAR file in place""", action="store_true")
parser.add_argument(
    "--lint",
    help="""Tyny and private version of code checker for vasp""",
    action="store_true",
)
parser.add_argument("incar_file", metavar="INCAR_file", nargs=1)
args = parser.parse_args()
assert not (
    args.lint and (args.i or args.r)
), "Lint option and re-format option (-i, -r) is exclusive."
logger.debug(f"args: {args}")
incar: vaspy.incar.Incar = vaspy.load(args.incar_file[0])

if args.i:
    with Path(args.incar_file[0]).open(mode="w") as incar_file:
        incar_file.write(incar.__str__())
if args.r:
    print(incar)

if args.lint:
    lint_msg = incar.lint_all()
    if lint_msg:  # if python 3.8 lint_msg:= incar.lint_all() can be used...
        print(lint_msg)
    else:
        print("ALL OK. Submit the job!!")
