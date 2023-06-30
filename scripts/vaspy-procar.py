#!/usr/bin/env python
"""Script to use (demonstrate) vaspy.procar module."""

import argparse
import functools as ft
import re
from itertools import chain
from logging import INFO, Formatter, StreamHandler, getLogger

from vaspy import procar, tools
from vaspy.outcar import OUTCAR

# logger
LOGLEVEL = INFO
logger = getLogger(__name__)
fmt = "%(asctime)s %(levelname)s %(name)s :%(message)s"
formatter = Formatter(fmt)
handler = StreamHandler()
handler.setLevel(LOGLEVEL)
logger.setLevel(LOGLEVEL)
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.propagate = True

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument(
    "--outcar",
    metavar="outcar_file",
    help="""Use "OUTCAR" for the Fermi level correction.
outcar_file must be specified.
NOTE:  E-fermi of OUTCAR generated in
Band-calculation may NOT be reliable.""",
)
parser.add_argument(
    "--fermi",
    metavar="value",
    type=float,
    help="""Fermi level correction
Energy shifts by this value
if --outcar is set, this option is ignored""",
)
parser.add_argument(
    "--site",
    metavar="atom_indices",
    dest="atomindex",
    action="append",
    type=tools.atom_selection_to_list,
    help="""atom index specifed with range.
Use "-" or ","
    (ex.) --site 1,2,7-9""",
)
parser.add_argument(
    "--as",
    metavar="name",
    nargs="+",
    dest="atomsetname",
    action="append",
    help="""the name of the sites identified
by --site option
the name is used in the title of the column""",
)
parser.add_argument(
    "--orbital",
    metavar="orbitals",
    action="append",
    type=ft.partial(re.split, r"[,:]"),
    help="""orbital name deliminated by ":" or ",".
orbital names are:
s, p, pxpy, pz, d, dxy, dyz, dz2, dxz, dx2
    (ex.) --orbital s:pxpy:d""",
)
parser.add_argument("procar", metavar="PROCAR_file", help="""PROCAR file""")

args = parser.parse_args()
logger.debug("Debugging...")
recvec = None
if args.outcar is not None:
    outcar = OUTCAR(args.outcar)
    fermi = outcar.fermi
    logger.debug(f"Fermi: {fermi}")
    recvec = [[v * 6.283185307179586 for v in recv] for recv in outcar.recvec]
    logger.debug(f"recvec: {recvec}")
elif args.fermi is not None:
    fermi = args.fermi
else:
    fermi = 0.0
procar = procar.PROCAR(args.procar)
procar.fermi_correction(fermi)
if recvec:
    logger.debug(f"procar.k_vectors: {procar.k_vectors}")
    procar.to_physical_kvector(recvec)
    logger.debug(f"procar.k_vectors in (AA): {procar.k_vectors}")

#
assert (
    len(args.atomindex) == len(args.orbital) == len(args.atomsetname)
), "--atom, --as and --orbital are mismatched."
#
sitenames = tuple(chain.from_iterable(args.atomsetname))
flat_orbitals = tuple(chain.from_iterable(args.orbital))
#
# As atomindex used here begins with "1", but siteindex used in procar.py
# internally begins with "0".  (VASP is written with fortran!)
siteindex = [[i - 1 for i in internal] for internal in args.atomindex]
#
logger.debug(f"siteindex: {siteindex}")
logger.debug(f"sitenames: {sitenames}")
logger.debug(f"args.orbitals: {args.orbital}")
logger.debug(f"orbitals flat: {flat_orbitals}")
#
for sites, name in zip(siteindex, sitenames):
    procar.append_sumsite(tuple(sites), name)
for orb in tuple(set(flat_orbitals)):
    procar.append_sumorbital(procar.orbital_index(orb), orb)
#
logger.debug("label['site'] is {}".format(procar.label["site"]))
logger.debug("label['orbital'] is {}".format(procar.label["orbital"]))
#
site_indexes: list[int] = []
orbtal_indexes_sets = []
tmp: list[str] | tuple[str, ...]
for site in sitenames:
    site_indexes.append(procar.label["site"].index(site))
for orbitals in args.orbital:
    tmp = []
    for orbital_in_site in orbitals:
        tmp.append(procar.label["orbital"].index(orbital_in_site))
    tmp = tuple(tmp)
    orbtal_indexes_sets.append(tmp)
orbital_indexes_sets = tuple(orbtal_indexes_sets)
#
logger.debug(f"site_indexes: {site_indexes}")
logger.debug(f"orbital_indexes_sets: {orbital_indexes_sets}")
#
print(
    procar.text_sheet(
        site_indexes=site_indexes, orbital_indexes_sets=orbital_indexes_sets,
    ),
)
