* Create a SAS data set named distance;
* Convert miles to kilometres;
DATA distance;
   Miles = 26.22;
   Kilometers = 1.61 * Miles;
RUN;
* Print the results;
PROC PRINT DATA = distance;
RUN;

* Start IML Session;
proc IML;
    start approx(x);
    y=1;
    do until(w<1e-3);
        z=y;
        y=.5#(z+x/z);
        w=abs(y-z);
    end;
    return(y);
    finish approx;
    t=approx({3,5,7,9});
    print t;

    z={1 2, 3 4, 5 6};
    y=z'*z;
    print y;
    quit;
