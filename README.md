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

Will make a more comprehensive README in the future, but this is the current state... Can refer to the documentation! 
