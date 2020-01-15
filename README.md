# nicerpy_xrayanalysis
Personal suite of programs for X-ray analysis of NICER Data 

v1 - January 14 2019 

- To create light curves, pulse profiles, color diagrams, and power spectra from NICER data

v2 - June 3 2019 

- To create light curves, pulse profiles, color diagrams, and power spectra from data processed through NICERsoft. Also edited Lv2_lc, Lv2_ps, Lv2_phase, Lv2_color to allow NICERsoft data.

v2.1 - July 18 2019

- Creating averaged power spectra from either one observation ID, or a list of observation IDs! One has to take care of the naming scheme for the merged event files though. 

v2.2 - August 23 2019

- Generalized Lv2_phase.py, where I've added in higher order correction terms to the phase, i.e., including second order and third order corrections from the Taylor series. This changes the input in Lv3_main.py a bit, but it was straightforward to adapt.
- Also added in routines particularly for NGC300 ULX. Could be generalized in the future, but for now, looking at just NGC300.

v3.0 - January 9 2020

- Significant revision to many scripts! Look at Evernote from first week of January 2020. The main revision is to use event files as inputs in most of the functions/modules, rather than ObsIDs, to generalize the algorithms a little. This is in anticipation for other missions (RXTE, for example), and for NICER quicklook stuff. It's mostly tested. Documentation hasn't been updated yet, but will do so on 1/10 or 1/11. 

- Most of the scripts were tested; most of the Lv1 and Lv2 scripts have been tested function by function, and as a whole. Some Lv3 scripts have been tested at this stage as well. The more substantial ones like Lv3_average_ps_main.py, Lv3_main.py, and Lv3_presto_main.py, are well-tested. The only exception is the "prepfold" and "ps2pdf" functions in Lv3_presto_main.py, because my personal PRESTO installation is screwing up. Will check in with this once my issue with the PRESTO installation is resolved.

v3.1 - January 11 2020

- Added Lv2_efsearch (for epoch folding) and Lv3_quicklook.py, a script that is intended to easily cycle through different basic search techniques (power spectrum, PRESTO, sideband searches, and epoch folding) as new data come out. 

- There is also testslider.py which is a work in progress - was thinking of adding a slider for the power spectra, but this might not be so necessary.

Will make a more comprehensive README in the future, but this is the current state... Can refer to the documentation! 

v3.2 - January 15 2020

- Adding Lv2_TBOs_method, which takes care of searching for burst candidates and adding contour/colormap plots. Edited Lv3_TBOs.py to make it more an executive script than a methods script.
