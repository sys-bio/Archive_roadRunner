// == PREAMBLE ================================================

// * Licensed under the Apache License, Version 2.0; see README

// == FILEDOC =================================================

/** @file ss-solver-test.c
* @author JKM
* @date 08/04/2015
* @copyright Apache License, Version 2.0
* @brief Tests the RoadRunner steady state solver from C
**/

#include "rrc_api.h"
#include "rrc_logging_api.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int main(int argc, char* argv[])
{
    int n;
    // enable logging
//     setLogLevel("debug");
//     setLogLevel("warning");
//     enableLoggingToFile();
//     {
//         char* t = getLogFileName();
//         fprintf(stderr,"Enabling logging to %s\n", t);
//         freeText(t);
//     }

    // print version
    {
        char* t = getVersionStr();
        fprintf(stderr,"RoadRunner version %s\n", t);
        freeText(t);
    }

    fprintf(stderr,"Initializing RoadRunner...\n");
    RRHandle rr = createRRInstance();

    // Dante's 9th circle
    {
        RRListPtr syms = getAvailableTimeCourseSymbols(rr);
        if(syms) {
            // Dante's 8th circle
            int n, k;
            for(n=0; n<getListLength(syms); ++n) {
                // Dante's 7th circle
                if(isListItemList(getListItem(syms, n))) {
                    // Almost there, my son...
                    RRListPtr cat = getList(getListItem(syms, n));
                    for(k=0; k<getListLength(cat); ++k) {
                        fprintf(stderr, "  %s\n", getStringListItem(getListItem(cat, k)));
                    }
                } else
                    fprintf(stderr, "  ???\n");
            }
        }
        freeRRList(syms);
    }

    freeRRInstance(rr);
	
    return 0;
}



