1. High rate data procesing
Run_Read_Plot_HighRateCASESdata_routines.m 
Calls routines: 
Fn_ReadHighRate_CASESdata.m to read CASES high rate data and
Fn_Plot_HighRate_CASESdata.m to plot CASES high rate data in several ways (all scintillating events, S4 and sigma_phi etc.). 

Needs to input the correct path name and folder name for the CASES data (where the daily .log files are). 

See Fn_ReadHighRate_CASESdata.m for details on signal type. Generally, it is 0 or 2 for CASES. Thus, i=1 or 3. Also, different plot types that can be created with the plotting routines are mentioned in the same code. Select what is necessary, or just plot everything.

The rest of the .m files are just functions required by these codes.


**********************
If you have any trouble running these codes or if you have any comments/ suggestions, please send me an email at kshitija@vt.edu
Thanks, Kshitija Deshpande, Virginia Tech.
**********************

