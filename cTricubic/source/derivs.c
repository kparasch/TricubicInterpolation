#include "derivs.h"

double* tricubic_finite_diff(int shape1, int shape2, double* A, int ix, int iy, int iz)
{
    double* b = (double*)malloc(64*sizeof(double));

    double  p0 = A[(((ix  ) * shape1) + (iy  )) * shape2 + (iz  )];
    double  p1 = A[(((ix+1) * shape1) + (iy  )) * shape2 + (iz  )];
    double  p2 = A[(((ix  ) * shape1) + (iy+1)) * shape2 + (iz  )];
    double  p3 = A[(((ix+1) * shape1) + (iy+1)) * shape2 + (iz  )];
    double  p4 = A[(((ix  ) * shape1) + (iy  )) * shape2 + (iz+1)];
    double  p5 = A[(((ix+1) * shape1) + (iy  )) * shape2 + (iz+1)];
    double  p6 = A[(((ix  ) * shape1) + (iy+1)) * shape2 + (iz+1)];
    double  p7 = A[(((ix+1) * shape1) + (iy+1)) * shape2 + (iz+1)];
    double  p8 = A[(((ix-1) * shape1) + (iy  )) * shape2 + (iz  )];
    double  p9 = A[(((ix+2) * shape1) + (iy  )) * shape2 + (iz  )];
    double p10 = A[(((ix-1) * shape1) + (iy+1)) * shape2 + (iz  )];
    double p11 = A[(((ix+2) * shape1) + (iy+1)) * shape2 + (iz  )];
    double p12 = A[(((ix-1) * shape1) + (iy  )) * shape2 + (iz+1)];
    double p13 = A[(((ix+2) * shape1) + (iy  )) * shape2 + (iz+1)];
    double p14 = A[(((ix-1) * shape1) + (iy+1)) * shape2 + (iz+1)];
    double p15 = A[(((ix+2) * shape1) + (iy+1)) * shape2 + (iz+1)];
    double p16 = A[(((ix  ) * shape1) + (iy-1)) * shape2 + (iz  )];
    double p17 = A[(((ix+1) * shape1) + (iy-1)) * shape2 + (iz  )];
    double p18 = A[(((ix  ) * shape1) + (iy+2)) * shape2 + (iz  )];
    double p19 = A[(((ix+1) * shape1) + (iy+2)) * shape2 + (iz  )];
    double p20 = A[(((ix  ) * shape1) + (iy-1)) * shape2 + (iz+1)];
    double p21 = A[(((ix+1) * shape1) + (iy-1)) * shape2 + (iz+1)];
    double p22 = A[(((ix  ) * shape1) + (iy+2)) * shape2 + (iz+1)];
    double p23 = A[(((ix+1) * shape1) + (iy+2)) * shape2 + (iz+1)];
    double p24 = A[(((ix  ) * shape1) + (iy  )) * shape2 + (iz-1)];
    double p25 = A[(((ix+1) * shape1) + (iy  )) * shape2 + (iz-1)];
    double p26 = A[(((ix  ) * shape1) + (iy+1)) * shape2 + (iz-1)];
    double p27 = A[(((ix+1) * shape1) + (iy+1)) * shape2 + (iz-1)];
    double p28 = A[(((ix  ) * shape1) + (iy  )) * shape2 + (iz+2)];
    double p29 = A[(((ix+1) * shape1) + (iy  )) * shape2 + (iz+2)];
    double p30 = A[(((ix  ) * shape1) + (iy+1)) * shape2 + (iz+2)];
    double p31 = A[(((ix+1) * shape1) + (iy+1)) * shape2 + (iz+2)];
    double p32 = A[(((ix-1) * shape1) + (iy-1)) * shape2 + (iz  )];
    double p33 = A[(((ix+2) * shape1) + (iy-1)) * shape2 + (iz  )];
    double p34 = A[(((ix-1) * shape1) + (iy+2)) * shape2 + (iz  )];
    double p35 = A[(((ix+2) * shape1) + (iy+2)) * shape2 + (iz  )];
    double p36 = A[(((ix-1) * shape1) + (iy-1)) * shape2 + (iz+1)];
    double p37 = A[(((ix+2) * shape1) + (iy-1)) * shape2 + (iz+1)];
    double p38 = A[(((ix-1) * shape1) + (iy+2)) * shape2 + (iz+1)];
    double p39 = A[(((ix+2) * shape1) + (iy+2)) * shape2 + (iz+1)];
    double p40 = A[(((ix-1) * shape1) + (iy  )) * shape2 + (iz-1)];
    double p41 = A[(((ix+2) * shape1) + (iy  )) * shape2 + (iz-1)];
    double p42 = A[(((ix-1) * shape1) + (iy+1)) * shape2 + (iz-1)];
    double p43 = A[(((ix+2) * shape1) + (iy+1)) * shape2 + (iz-1)];
    double p44 = A[(((ix-1) * shape1) + (iy  )) * shape2 + (iz+2)];
    double p45 = A[(((ix+2) * shape1) + (iy  )) * shape2 + (iz+2)];
    double p46 = A[(((ix-1) * shape1) + (iy+1)) * shape2 + (iz+2)];
    double p47 = A[(((ix+2) * shape1) + (iy+1)) * shape2 + (iz+2)];
    double p48 = A[(((ix  ) * shape1) + (iy-1)) * shape2 + (iz-1)];
    double p49 = A[(((ix+1) * shape1) + (iy-1)) * shape2 + (iz-1)];
    double p50 = A[(((ix  ) * shape1) + (iy+2)) * shape2 + (iz-1)];
    double p51 = A[(((ix+1) * shape1) + (iy+2)) * shape2 + (iz-1)];
    double p52 = A[(((ix  ) * shape1) + (iy-1)) * shape2 + (iz+2)];
    double p53 = A[(((ix+1) * shape1) + (iy-1)) * shape2 + (iz+2)];
    double p54 = A[(((ix  ) * shape1) + (iy+2)) * shape2 + (iz+2)];
    double p55 = A[(((ix+1) * shape1) + (iy+2)) * shape2 + (iz+2)];
    double p56 = A[(((ix-1) * shape1) + (iy-1)) * shape2 + (iz-1)];
    double p57 = A[(((ix+2) * shape1) + (iy-1)) * shape2 + (iz-1)];
    double p58 = A[(((ix-1) * shape1) + (iy+2)) * shape2 + (iz-1)];
    double p59 = A[(((ix+2) * shape1) + (iy+2)) * shape2 + (iz-1)];
    double p60 = A[(((ix-1) * shape1) + (iy-1)) * shape2 + (iz+2)];
    double p61 = A[(((ix+2) * shape1) + (iy-1)) * shape2 + (iz+2)];
    double p62 = A[(((ix-1) * shape1) + (iy+2)) * shape2 + (iz+2)];
    double p63 = A[(((ix+2) * shape1) + (iy+2)) * shape2 + (iz+2)];



    b[ 0] = p0;
    b[ 1] = p1;
    b[ 2] = p2;
    b[ 3] = p3;
    b[ 4] = p4;
    b[ 5] = p5;
    b[ 6] = p6;
    b[ 7] = p7;

    b[ 8] = 0.5*(p1 - p8);
    b[ 9] = 0.5*(p9 - p0);
    b[10] = 0.5*(p3 - p10);
    b[11] = 0.5*(p11 - p2);
    b[12] = 0.5*(p5 - p12);
    b[13] = 0.5*(p13 - p4);
    b[14] = 0.5*(p7 - p14);
    b[15] = 0.5*(p15 - p6);

    b[16] = 0.5*(p2 - p16);
    b[17] = 0.5*(p3 - p17);
    b[18] = 0.5*(p18 - p0);
    b[19] = 0.5*(p19 - p1);
    b[20] = 0.5*(p6 - p20);
    b[21] = 0.5*(p7 - p21);
    b[22] = 0.5*(p22 - p4);
    b[23] = 0.5*(p23 - p5);

    b[24] = 0.5*(p4 - p24);
    b[25] = 0.5*(p5 - p25);
    b[26] = 0.5*(p6 - p26);
    b[27] = 0.5*(p7 - p27);
    b[28] = 0.5*(p28 - p0);
    b[29] = 0.5*(p29 - p1);
    b[30] = 0.5*(p30 - p2);
    b[31] = 0.5*(p31 - p3);

    b[32] = 0.25*((p3 - p10) + (- p17 + p32));
    b[33] = 0.25*((p11 - p2) + (- p33 + p16));
    b[34] = 0.25*((p19 - p34) + (- p1 + p8));
    b[35] = 0.25*((p35 - p18) + (- p9 + p0));
    b[36] = 0.25*((p7 - p14) + (- p21 + p36));
    b[37] = 0.25*((p15 - p6) + (- p37 + p20));
    b[38] = 0.25*((p23 - p38) + (- p5 + p12));
    b[39] = 0.25*((p39 - p22) + (- p13 + p4));

    b[40] = 0.25*((p5 - p12) + (- p25 + p40));
    b[41] = 0.25*((p13 - p4) + (- p41 + p24));
    b[42] = 0.25*((p7 - p14) + (- p27 + p42));
    b[43] = 0.25*((p15 - p6) + (- p43 + p26));
    b[44] = 0.25*((p29 - p44) + (- p1 + p8));
    b[45] = 0.25*((p45 - p28) + (- p9 + p0));
    b[46] = 0.25*((p31 - p46) + (- p3 + p10));
    b[47] = 0.25*((p47 - p30) + (- p11 + p2));

    b[48] = 0.25*((p6 - p20) + (- p26 + p48));
    b[49] = 0.25*((p7 - p21) + (- p27 + p49));
    b[50] = 0.25*((p22 - p4) + (- p50 + p24));
    b[51] = 0.25*((p23 - p5) + (- p51 + p25));
    b[52] = 0.25*((p30 - p52) + (- p2 + p16));
    b[53] = 0.25*((p31 - p53) + (- p3 + p17));
    b[54] = 0.25*((p54 - p28) + (- p18 + p0));
    b[55] = 0.25*((p55 - p29) + (- p19 + p1));

    b[56] = 0.125*(((p7 - p14) + (- p21 + p36)) + ((- p27 + p42) + (+ p49 - p56)));
    b[57] = 0.125*(((p15 - p6) + (- p37 + p20)) + ((- p43 + p26) + (+ p57 - p48)));
    b[58] = 0.125*(((p23 - p38) + (- p5 + p12)) + ((- p51 + p58) + (+ p25 - p40)));
    b[59] = 0.125*(((p39 - p22) + (- p13 + p4)) + ((- p59 + p50) + (+ p41 - p24)));
    b[60] = 0.125*(((p31 - p46) + (- p53 + p60)) + ((- p3 + p10) + (+ p17 - p32)));
    b[61] = 0.125*(((p47 - p30) + (- p61 + p52)) + ((- p11 + p2) + (+ p33 - p16)));
    b[62] = 0.125*(((p55 - p62) + (- p29 + p44)) + ((- p19 + p34) + (+ p1 - p8)));
    b[63] = 0.125*(((p63 - p54) + (- p45 + p28)) + ((- p35 + p18) + (+ p9 - p0)));

    return b;
}

double* tricubic_exact_diff(int shape1, int shape2, int shape3, double *A, int ix, int iy, int iz, double dx, double dy, double dz)
{
    double scale[8] = {1., dx, dy, dz, dx * dy, dx * dz, dy * dz, (dx * dy) * dz};

    double* b = (double*)malloc(64*sizeof(double));



    int p0 = (((ix  ) * shape1 + (iy  )) * shape2 + iz  ) * shape3; 
    int p1 = (((ix+1) * shape1 + (iy  )) * shape2 + iz  ) * shape3;  
    int p2 = (((ix  ) * shape1 + (iy+1)) * shape2 + iz  ) * shape3;  
    int p3 = (((ix+1) * shape1 + (iy+1)) * shape2 + iz  ) * shape3;  
    int p4 = (((ix  ) * shape1 + (iy  )) * shape2 + iz+1) * shape3;
    int p5 = (((ix+1) * shape1 + (iy  )) * shape2 + iz+1) * shape3;
    int p6 = (((ix  ) * shape1 + (iy+1)) * shape2 + iz+1) * shape3;
    int p7 = (((ix+1) * shape1 + (iy+1)) * shape2 + iz+1) * shape3;

    for(int l = 0 ; l < 8; l++)
    {
        b[8*l+0] = A[p0 + l] * scale[l];
        b[8*l+1] = A[p1 + l] * scale[l];
        b[8*l+2] = A[p2 + l] * scale[l];
        b[8*l+3] = A[p3 + l] * scale[l];
        b[8*l+4] = A[p4 + l] * scale[l];
        b[8*l+5] = A[p5 + l] * scale[l];
        b[8*l+6] = A[p6 + l] * scale[l];
        b[8*l+7] = A[p7 + l] * scale[l];
    }

    return b;
}
