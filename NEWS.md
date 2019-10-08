# cliqueMS 0.99.3

* Changes in NEWS.md format.
* Change in 'anClique' class from S3 to S4,
new documentation together with this change.
* Internal changes in functions to work with S4 class.


# cliqueMS 0.99.0

* First Submission to Biocoductor. Changes in format the of
the vignette, and functions.


# cliqueMS 0.3.2

* Added functionality in 'getAnnotation' function,
now reported scores can be normalized and then
the comparison between annotations in different
groups is easier.

## BUG FIXES

* A bug in the score report was corrected.

# cliqueMS 0.3.1

* Vignette modified due to new 'xcsm'
package, inside new R version 3.6.

# cliqueMS 0.3.0

* New methods for functions 'anClique' and 'createNetwork'.
Methods 'anClique.XCMSnExp' and 'createNetwork.XCMSnExp'
developed for mz processed objects of class 'XCMSnExp'
from xcms package.
* In addition, a new internal function
has been created to compute cosine similarity. 
Therefore, now the installation of package 'CAMERA'
is in 'Suggests', because is only needed for
analysing processed mzdata with the older 'xcmsSet'
class of 'xcms'.

# cliqueMS 0.2.4

* Function "getAnnotation" now it is more
accurate, because for isotopic features,
it uses charge information from "getIsotopes"
function.

# cliqueMS 0.2.3

* Function "computeCliques" corrected, there was a
small bug in the "step" variable.
Function "getIsotopes" corrected, bug when
"filterIso" leave no isotopes in the group.

## BUG FIXES

* In function "filterIso", small bug corrected.
maxGrade parameter was not working correctly, fixed
for function "getIsotopes".

# cliqueMS 0.2.2

## BUG FIXES

* Function "getIsotopes" updated: there was a
bug when the sample had no isotopes, now corrected.
* Function "getAnnotation" updated: there was a
bug when the sample had no isotopes, now corrected.

# cliqueMS 0.2.1

*First version of cliqueMS