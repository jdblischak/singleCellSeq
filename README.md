# SingleCellSeq

## Important links

See the [site][] to view the results.

Read the [paper][].

Download the raw FASTQ files at GEO record [GSE77288][geo].

See the [contributing guidelines][contrib] to add a new analysis.

[site]: http://jdblischak.github.io/singleCellSeq/analysis
[paper]: http://www.nature.com/articles/srep39921
[geo]: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77288
[contrib]: https://github.com/jdblischak/singleCellSeq/blob/master/CONTRIBUTING.md

## Project organization

*  `analysis/` contains the R Markdown files documenting the results
*  `paper/` contains the R Markdown files used to write the paper
*  `code/` contains command-line scripts
*  `data/` contains various summarized data files
*  `tests/` contains some sanity checks

The majority of the generated content, e.g. figures and HTML files, are in the `analysis/` subdirectory of the [gh-pages][] branch.

[gh-pages]: https://github.com/jdblischak/singleCellSeq/tree/gh-pages/analysis

## Useful data files

*  [reads.txt][reads] - tab-delimited file contains read counts for the 19,027 Ensembl protein-coding genes with at least one observed read in at least one of the 864 single cell samples
*  [molecules.txt][molecules] - same as `reads.txt`, but for the molecule counts
*  [reads-filter.txt][reads-filter] - tab-delimited file contains read counts for the 13,106 Ensembl protein-coding genes that passed our expression level filter for the 564 single cells that passed our quality control filters
*  [molecules-filter.txt][molecules-filter] - same as `reads-filter.txt`, but for the molecule counts
*  [molecules-final.txt][molecules-final] - tab-delimited file contains the final gene expression values after all of our processing steps


[reads]: https://github.com/jdblischak/singleCellSeq/blob/gh-pages/data/reads.txt
[molecules]: https://github.com/jdblischak/singleCellSeq/blob/gh-pages/data/molecules.txt
[reads-filter]: https://github.com/jdblischak/singleCellSeq/blob/gh-pages/data/reads-filter.txt
[molecules-filter]: https://github.com/jdblischak/singleCellSeq/blob/gh-pages/data/molecules-filter.txt
[molecules-final]: https://github.com/jdblischak/singleCellSeq/blob/gh-pages/data/molecules-final.txt
