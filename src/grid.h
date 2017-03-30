

void generate_grid();


generate_grid(){

    for(i=0; i<Nfaces; i++) {

        X[i]   = (xf-x0) * (i)   / Nfaces + x0 ;  // node 0, element i
    }

    for (i=0; i<Nelem; i++) {

        h_j[i]= X[i+1]-X[i];
    }

    return;
}
