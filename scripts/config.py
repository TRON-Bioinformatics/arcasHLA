"""
Centralized configuration.
"""

import json

from os.path import dirname, realpath

root_dir = dirname(realpath(__file__)) + "/../"

# Static, repository-checked-in assets. These are independent of any
# specific built arcasHLA reference: `reference build` copies them
# unchanged into every reference it produces, and `genes`/`populations`
# are used to validate CLI arguments before a reference directory is
# even resolved. The reference data itself (indexes, conversion tables,
# etc.) lives in a built reference directory resolved dynamically via
# `ref_paths` (see `arcas_utilities.get_ref_dir`/`ref_path`).
info_dir = root_dir + "dat/info/"
parameters_json = info_dir + "parameters.json"

with open(parameters_json, "r") as param_file:
    genes, populations, versions = json.load(param_file)
    genes = set(genes)
    populations = set(populations)
