#pragma once
#include "global_var.h"
#include "general_tools.h"

void generate_grid();


void generate_grid(){

    register int i;

    dx=(xf-x0)/Nelem;

    for(i=0; i<Nfaces; i++) {

        X[i]   = dx * (i)  + x0 ;  // node 0, element i
    }

    for (i=0; i<Nelem; i++) {

        h_j[i]= X[i+1]-X[i];
    }

    return;
}
