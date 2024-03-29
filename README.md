# IV Checklist
Code for implementing the IV checklist in Felton and Stewart (2023). The abstract for the paper is as follows:

> Instrumental variables (IV) analysis is a powerful, but fragile, tool for drawing causal inferences from observational data. Sociologists increasingly turn to this strategy in settings where unmeasured confounding between the treatment and outcome is likely. This paper reviews the assumptions required for IV and the consequences of violating them, focusing on sociological applications. We highlight three methodological problems IV faces: (i) identification bias, an asymptotic bias from assumption violations; (ii) estimation bias, a finite-sample bias that persists even when assumptions hold; and (iii) type-M error, the exaggeration of effects given statistical significance. In each case, we emphasize how weak instruments exacerbate these problems and make results sensitive to minor violations of assumptions. We survey IV papers from top sociology journals, showing that assumptions often go unstated and robust uncertainty measures are rarely used. We provide a practical checklist to show how IV, despite its fragility, can still be useful when handled with care.

Comments are very welcome! You can e-mail us at christopher_felton@gse.harvard.edu and bms4@princeton.edu.

# Obtaining the Data

Obtaining the data requires an IPUMS account. Go to https://usa.ipums.org/usa-action/variables/group to create an extract using the 1990 Census 5% microdata. We list the exact variables used in cg2006_IPUMS_summary.pdf. We exported the data as a .csv file and named it ``census.csv``. If you use a different name, be sure to change it in ``cg2006tidying.R``. 

# Running the Code

The easiest way to run everything will be to place ``IV_checklist.Rproj`` into a directory with subfolders "data" and "code." You can download the "code" folder here. Create a new folder called "data" and follow the directions above to obtain ``census.csv``. Then put it into "data" and launch ``IV_checklist.Rproj``. The scripts use the ``here`` package, so you won't have to specify a working directory.

``cg2006tidying.R`` prepares the data for the main IV analysis, and ``iv_checklist.R`` carries out the procedures we describe in the checklist (not including things like stating assumptions or clarifying who the compliers are). We were unable to replicate Conley and Glauber's (2006) results exactly, but our treatment effect and standard error estimates are very close to theirs. 
