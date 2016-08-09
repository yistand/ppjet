//numbers from http://people.physics.tamu.edu/sakuma/star/jets/c101121_event_selection/s0150_mclist_001/web.php

#ifndef __PPCROSSSECPERPT
#define __PPCROSSSECPERPT

const int NUMBEROFPT = 11;
const char *PTBINS[NUMBEROFPT]={"3_4","4_5","5_7","7_9","9_11","11_15","15_25","25_35","35_45","45_55","55_65"};
float XSEC[NUMBEROFPT] = {1.3e9, 3.15e8, 1.37e8, 2.3e7, 5.53e6, 2.2e6, 3.9e5, 1.02e4, 5.01e2, 2.86e1, 1.45};
float NUMBEROFEVENT[NUMBEROFPT] = {686000, 500000, 398000, 420000, 414307, 420000, 397200, 400000, 110000, 118000, 120000};

#endif // __PPCROSSSECPERPT
