# README

## Build the document

This directory contains the source files for creating the manuscript.
It is compiled in two steps.
First, the R Markdown files are converted to Markdown using knitr.
Second, pandoc converts the markdown files to HTML, PDF, and Word.
This is all automated with a Makefile.
Use one of the following commands to build the paper.

```bash
# Build paper.docx
$ make word
# Build paper.html
$ make html
# Build paper.pdf
$ make pdf
# Build all 3
$ make
```

## Manual formatting steps

The final submission is the Word document.
Unfortunately not all the formatting steps can be automated using pandoc.
The following steps must be performed manually to complete the formatting of the Word document:

*  Convert the list the authors and their affilations to the Author Style
*  Add the necessary blank lines in the title page
*  Insert page break so that the Abstract starts on the second page
*  Insert page numbers on the bottom right of the page
*  Add line numbers starting after the title page (Page Layout -> Line Numbers)
*  Resize the figures to have a width of 6 inches (after you resize the first figure, then you only have to select the remaining figures individually and press F4 - [directions][f4]
*  Crop the supplementary figures
*  Extract the supplementary figures and tables into separate document
*  Convert supplement to PDF

[f4]: http://answers.microsoft.com/en-us/office/forum/office_2010-word/select-multiple-images-in-word-to-resize-all/2061eed7-0522-4127-9b84-f490c02e2d81

## Additional manual formatting steps for Windows

If you're running Microsoft Word on Windows, a few of the automated pandoc formats are not instituted.
Perform the following additional manual steps:

*  Select the references and assign the style Bibiography
*  Set the side margins to be 0.7 inches
