#!/usr/bin/env python

"""
Format a bibtex file.

1. Organize alphabetically by ID.
2. Only include the following fields: author, year, title, journal,
   volume, number, pages, doi.

Usage:
  python format-bibtex.py file.bib

Note that it overwrites the original bibtex file, but only if the
script completes successfully. If the script fails due to an error,
the original file is unmodified and a temporary file of the same name
is available in /tmp. The last entry in this temporary file is the
last bibtex entry that was successfully formated.
"""

import os
import sys
import bibtexparser as bib
import shutil

bibtex_input = sys.argv[1]
assert os.path.isfile(bibtex_input), "bibtex file does not exist"

with open(bibtex_input) as bibtex_file:
    bib_database = bib.load(bibtex_file)

bib_dict = bib_database.entries_dict
keys = bib_dict.keys()
keys = list(keys)
keys.sort()


bibtex_input_base = os.path.basename(bibtex_input)
bibtex_output = open("/tmp/" + bibtex_input_base, "w")

for k in keys:
    entry = bib_dict[k]
    entrytype = entry["ENTRYTYPE"]
    out = "@%s{%s,\n"%(entrytype, entry["ID"])
    if entrytype == "article":
        for field in ["author", "year", "title", "journal", "volume", "number",
                      "pages", "doi"]:
            # For consortia, encase author name in double brackets so that the
            # full name is used instead of parsing it into a first and last name
            if field == "author" and "," not in entry["author"]:
                out = out + "%s = {{%s}},\n"%(field, entry[field])
            else:
                out = out + "%s = {%s},\n"%(field, entry[field])
    else:
        fields_all = list(entry.keys())
        fields_all.sort()
        for field in fields_all:
            if not field.islower():
                continue
            out = out + "%s = {%s},\n"%(field, entry[field])

    out = out.rstrip(",\n") + "\n}\n"
    bibtex_output.write(out)

bibtex_output.close()

# Overwrite original file
shutil.move("/tmp/" + bibtex_input_base, bibtex_input)
