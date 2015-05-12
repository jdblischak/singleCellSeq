# How to contribute

## Outline

*  [Running RStudio Server](#running-rstudio-server)
*  [Creating a new analysis](#creating-a-new-analysis)
*  [Style guide](#style-guide)
*  [Building the site](#building-the-site)

## Running RStudio Server

After ssh'ing into spudhead, request an interactive session with sufficient memory.

```bash
ql 8g
```

Next start an instance of RStudio Server running in the background:

```bash
/mnt/lustre/data/tools/rstudio_new/bin/rserver &
```

Then on your local computer, run the following:

```bash
ssh -N -f -L localhost:8787:spudling##:8787 user-name@pps-gateway.uchicago.edu
```

replacing spudling## with the name of the spudling where you started the RStudio Server instance, e.g. `spudling87` or `bigmem01`, and user-name with your login ID.

Open a browser to the address http://127.0.0.1:8787/ to access your RStudio instance.
From here, you can choose "Open Project" and select `singleCellSeq.Rproj`.

When you're finished, you can close the browser tab, and then kill the RStudio instance.
Find the 5-digit process id (PID) for rserver and rsession using the command `ps`.
End these processes by running `kill #####`.

## Creating a new analysis

Here are the steps for creating a new analysis:

*  Open RStudio project `singleCellSeq.Rproj`.
*  Set working directory to `analysis` with `setwd` or "Session" -> "Set Working Directory" -> "To Files Pane Location".
*  Create a copy of [template.Rmd][].
*  Change the author, title, and date at the top of the file.
*  Add the analysis code.
*  Use the RStudio button "Preview HTML" to view the results.
*  Add the analysis to the list in [index.Rmd][].
*  Add, commit, and push the new analysis.

```bash
cd analysis
git add new-analysis.Rmd index.Rmd
git commit -m "Add new analysis on..."
git push origin master
```

[template.Rmd]: https://raw.githubusercontent.com/jdblischak/singleCellSeq/master/analysis/template.Rmd
[index.Rmd]: https://raw.githubusercontent.com/jdblischak/singleCellSeq/master/analysis/index.Rmd

## Style guide

For consistency, we'll use the following conventions:

*  Name variables using snakecase, e.g. `gene_exp_mat`.
*  Name files with dashes, e.g. `this-is-a-long-filename.txt`.
*  Name directories with camelCase, e.g `data`, `rawData`.
*  Use `<-` for assignment.
*  Surround binary operators with spaces, e.g. `1 + 1`, not `1+1`.
*  Use two spaces for indentation.

When in doubt, use the style indicated either in [Google's R Style Guide][google-style] or [Hadley's R Style Guide][hadley-style].

[google-style]: https://google-styleguide.googlecode.com/svn/trunk/Rguide.xml
[hadley-style]: http://r-pkgs.had.co.nz/style.html

## Building the site

```bash
git checkout gh-pages
git merge master
cd analysis
make
git add -f *html figure/*
git commit -m "Build site."
git push origin gh-pages
```
