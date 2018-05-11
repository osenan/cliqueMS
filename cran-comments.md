## Test environments

1. Ubuntu 16.04, R 3.4.3
1. win-builder R-devel

## R CMD check results

### round 1

No ERRORs and WARNINGs, two NOTEs. One for the first release
and another for installed size, but it is needed external raw
data for the vignette.

### round 2

Corrected the "Depends" on R version with patch level 0.
Included "cran-comments.md" in ".Rbuildignore".

### round 3

Changed "\dontrun" for "\donttest". Added CITATION. Modified
DESCRIPTION, including adding a reference
(current reference is a conference, an article has been
submitted and is under peer reviewing process). 