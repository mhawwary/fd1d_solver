#include"solver_tools.h"


void test_cyclic_tridiag(const int& n){

    register int i;
    double *a=nullptr,*b=nullptr,*c=nullptr,*d=nullptr,*x=nullptr;
    x = new double[n];
    construct_tridiag_example(n,a,b,c,d);
    a[0]=2.0;
    c[n-1]=2.0;

    cyclic_tridiag_solve_mh(n,a,b,c,d,x);

    printf("\nChecking x:\n------------------------------------\n");
    for(i=0; i<n; i++)
        printf("%1.4f\n",x[i]);

    emptyarray(a);
    emptyarray(b);
    emptyarray(c);
    emptyarray(d);
    emptyarray(x);

    return;
}

void test_tridaig(const int& n){

    register int i;
    double *a=nullptr,*b=nullptr,*c=nullptr,*d=nullptr,*x=nullptr;
    x = new double[n];
    construct_tridiag_example(n,a,b,c,d);
    tridiag_solve_mh(n,a,b,c,d,x);

    printf("\nChecking x:\n------------------------------------\n");
    for(i=0; i<n; i++)
        printf("%1.4f\n",x[i]);

    emptyarray(a);
    emptyarray(b);
    emptyarray(c);
    emptyarray(d);
    emptyarray(x);

    return;
}

void construct_tridiag_example(const int n, double *&a, double *&b, double *&c
                               , double *&d){
    register int i;
    a = new double[n];
    b = new double[n];
    c = new double[n];
    d = new double[n];

    for(i=0; i<n; i++){
        a[i]=-1.0;
        b[i]=3.0;
        c[i]=1.0;
        d[i]=1.0;
    }

    d[0]=2.0;
    d[n-1] =2.0;

    return;
}

void tridiag_solve_mh(const int n, const double *a, const double *b
                     , const double *c, const double *d, double *&x){
    register int i;
    double *gam=nullptr;
    gam = new double[n];

    double bet = b[0];
    if(bet == 0)
        FatalErrorST("warning, b[0] == 0, division by zero is expected");

    x[0] = d[0]/bet;

    for(i=1; i<n; i++){       // forward elimination
        gam[i] = c[i-1]/bet;
        bet = b[i] - a[i] * gam[i];

        if(bet == 0){
            char temp_str[50];
            sprintf(temp_str,"warning, b[%d] == 0, division by zero is expected",i);
            FatalErrorST(temp_str);
        }
        x[i] = (d[i] - a[i] * x[i-1])/ bet;
    }

    for (i=n-2; i>-1; i--){  // backward substitution
        x[i] = x[i] - gam[i+1] * x[i+1];
    }

    emptyarray(gam);

    return;
}

void cyclic_tridiag_solve_mh(const int n, const double *a, const double *b
                     , const double *c, const double *d, double *&x){
    register int i;

    double *u=nullptr, *z=nullptr, *bb=nullptr;
    u = new double[n];
    z = new double[n];
    bb = new double[n];

    double beta = a[0];
    double alpha = c[n-1];
    double gam = - b[0];
    bb[0] = b[0] - gam;
    bb[n-1] = b[n-1] - alpha * beta/gam;
    u[0] = gam;
    u[n-1] = alpha;

    for(i=1; i<n-1; i++){
        u[i]=0.0;
        bb[i]=b[i];
    }

    tridiag_solve_mh(n,a,bb,c,d,x);
    tridiag_solve_mh(n,a,bb,c,u,z);

    double fact = (x[0] + (beta * x[n-1] / gam))
            / (1.0 + z[0] + (beta * z[n-1]/gam));

    for(i=0; i<n; i++)
        x[i] = x[i] - fact*z[i];

    emptyarray(u);
    emptyarray(z);
    emptyarray(bb);

    return;
}

