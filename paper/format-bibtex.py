#!/usr/bin/env python

"""
Format a bibtex file.

1. Organize alphabetically by ID.
2. Only include the following fields: author, year, title, journal,
   volumne, number, pages, doi.

Usage:
  python format-bibtex.py file.bib

Note that it overwrites the original bibtex file.
"""

import os
import sys
import bibtexparser as bib

bibtex_input = sys.argv[1]
assert os.path.isfile(bibtex_input), "bibtex file exists"

with open(bibtex_input) as bibtex_file:
    bib_database = bib.load(bibtex_file)

bib_dict = bib_database.entries_dict
keys = bib_dict.keys()
keys = list(keys)
keys.sort()

# Overwrite original file
bibtex_output = open(bibtex_input, "w")

for k in keys:
    entry = bib_dict[k]
    out = "@article{%s,\n"%(entry["ID"])
    for field in ["author", "year", "title", "journal", "volume", "number",
                  "pages", "doi"]:
        out = out + "%s = {%s},\n"%(field, entry[field])
    out = out.rstrip(",\n") + "\n}\n"
    bibtex_output.write(out)

bibtex_output.close()

