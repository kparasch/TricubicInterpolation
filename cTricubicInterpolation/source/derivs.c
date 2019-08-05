#include "derivs.h"

double* finite_diff(int shape1, int shape2, double A[][shape1][shape2], int ix, int iy, int iz)
{
    double* b = (double*)malloc(64*sizeof(double));

    b[ 0] = A[ix  ][iy  ][iz  ];
    b[ 1] = A[ix+1][iy  ][iz  ];
    b[ 2] = A[ix  ][iy+1][iz  ];
    b[ 3] = A[ix+1][iy+1][iz  ];
    b[ 4] = A[ix  ][iy  ][iz+1];
    b[ 5] = A[ix+1][iy  ][iz+1];
    b[ 6] = A[ix  ][iy+1][iz+1];
    b[ 7] = A[ix+1][iy+1][iz+1];

    b[ 8] = 0.5*(A[ix+1][iy  ][iz  ] - A[ix-1][iy  ][iz  ]);
    b[ 9] = 0.5*(A[ix+2][iy  ][iz  ] - A[ix  ][iy  ][iz  ]);
    b[10] = 0.5*(A[ix+1][iy+1][iz  ] - A[ix-1][iy+1][iz  ]);
    b[11] = 0.5*(A[ix+2][iy+1][iz  ] - A[ix  ][iy+1][iz  ]);
    b[12] = 0.5*(A[ix+1][iy  ][iz+1] - A[ix-1][iy  ][iz+1]);
    b[13] = 0.5*(A[ix+2][iy  ][iz+1] - A[ix  ][iy  ][iz+1]);
    b[14] = 0.5*(A[ix+1][iy+1][iz+1] - A[ix-1][iy+1][iz+1]);
    b[15] = 0.5*(A[ix+2][iy+1][iz+1] - A[ix  ][iy+1][iz+1]);

    b[16] = 0.5*(A[ix  ][iy+1][iz  ] - A[ix  ][iy-1][iz  ]);
    b[17] = 0.5*(A[ix+1][iy+1][iz  ] - A[ix+1][iy-1][iz  ]);
    b[18] = 0.5*(A[ix  ][iy+2][iz  ] - A[ix  ][iy  ][iz  ]);
    b[19] = 0.5*(A[ix+1][iy+2][iz  ] - A[ix+1][iy  ][iz  ]);
    b[20] = 0.5*(A[ix  ][iy+1][iz+1] - A[ix  ][iy-1][iz+1]);
    b[21] = 0.5*(A[ix+1][iy+1][iz+1] - A[ix+1][iy-1][iz+1]);
    b[22] = 0.5*(A[ix  ][iy+2][iz+1] - A[ix  ][iy  ][iz+1]);
    b[23] = 0.5*(A[ix+1][iy+2][iz+1] - A[ix+1][iy  ][iz+1]);

    b[24] = 0.5*(A[ix  ][iy  ][iz+1] - A[ix  ][iy  ][iz-1]);
    b[25] = 0.5*(A[ix+1][iy  ][iz+1] - A[ix+1][iy  ][iz-1]);
    b[26] = 0.5*(A[ix  ][iy+1][iz+1] - A[ix  ][iy+1][iz-1]);
    b[27] = 0.5*(A[ix+1][iy+1][iz+1] - A[ix+1][iy+1][iz-1]);
    b[28] = 0.5*(A[ix  ][iy  ][iz+2] - A[ix  ][iy  ][iz  ]);
    b[29] = 0.5*(A[ix+1][iy  ][iz+2] - A[ix+1][iy  ][iz  ]);
    b[30] = 0.5*(A[ix  ][iy+1][iz+2] - A[ix  ][iy+1][iz  ]);
    b[31] = 0.5*(A[ix+1][iy+1][iz+2] - A[ix+1][iy+1][iz  ]);

    b[32] = 0.25*((A[ix+1][iy+1][iz  ] - A[ix-1][iy+1][iz  ]) + (- A[ix+1][iy-1][iz  ] + A[ix-1][iy-1][iz  ]));
    b[33] = 0.25*((A[ix+2][iy+1][iz  ] - A[ix  ][iy+1][iz  ]) + (- A[ix+2][iy-1][iz  ] + A[ix  ][iy-1][iz  ]));
    b[34] = 0.25*((A[ix+1][iy+2][iz  ] - A[ix-1][iy+2][iz  ]) + (- A[ix+1][iy  ][iz  ] + A[ix-1][iy  ][iz  ]));
    b[35] = 0.25*((A[ix+2][iy+2][iz  ] - A[ix  ][iy+2][iz  ]) + (- A[ix+2][iy  ][iz  ] + A[ix  ][iy  ][iz  ]));
    b[36] = 0.25*((A[ix+1][iy+1][iz+1] - A[ix-1][iy+1][iz+1]) + (- A[ix+1][iy-1][iz+1] + A[ix-1][iy-1][iz+1]));
    b[37] = 0.25*((A[ix+2][iy+1][iz+1] - A[ix  ][iy+1][iz+1]) + (- A[ix+2][iy-1][iz+1] + A[ix  ][iy-1][iz+1]));
    b[38] = 0.25*((A[ix+1][iy+2][iz+1] - A[ix-1][iy+2][iz+1]) + (- A[ix+1][iy  ][iz+1] + A[ix-1][iy  ][iz+1]));
    b[39] = 0.25*((A[ix+2][iy+2][iz+1] - A[ix  ][iy+2][iz+1]) + (- A[ix+2][iy  ][iz+1] + A[ix  ][iy  ][iz+1]));

    b[40] = 0.25*((A[ix+1][iy  ][iz+1] - A[ix-1][iy  ][iz+1]) + (- A[ix+1][iy  ][iz-1] + A[ix-1][iy  ][iz-1]));
    b[41] = 0.25*((A[ix+2][iy  ][iz+1] - A[ix  ][iy  ][iz+1]) + (- A[ix+2][iy  ][iz-1] + A[ix  ][iy  ][iz-1]));
    b[42] = 0.25*((A[ix+1][iy+1][iz+1] - A[ix-1][iy+1][iz+1]) + (- A[ix+1][iy+1][iz-1] + A[ix-1][iy+1][iz-1]));
    b[43] = 0.25*((A[ix+2][iy+1][iz+1] - A[ix  ][iy+1][iz+1]) + (- A[ix+2][iy+1][iz-1] + A[ix  ][iy+1][iz-1]));
    b[44] = 0.25*((A[ix+1][iy  ][iz+2] - A[ix-1][iy  ][iz+2]) + (- A[ix+1][iy  ][iz  ] + A[ix-1][iy  ][iz  ]));
    b[45] = 0.25*((A[ix+2][iy  ][iz+2] - A[ix  ][iy  ][iz+2]) + (- A[ix+2][iy  ][iz  ] + A[ix  ][iy  ][iz  ]));
    b[46] = 0.25*((A[ix+1][iy+1][iz+2] - A[ix-1][iy+1][iz+2]) + (- A[ix+1][iy+1][iz  ] + A[ix-1][iy+1][iz  ]));
    b[47] = 0.25*((A[ix+2][iy+1][iz+2] - A[ix  ][iy+1][iz+2]) + (- A[ix+2][iy+1][iz  ] + A[ix  ][iy+1][iz  ]));

    b[48] = 0.25*((A[ix  ][iy+1][iz+1] - A[ix  ][iy-1][iz+1]) + (- A[ix  ][iy+1][iz-1] + A[ix  ][iy-1][iz-1]));
    b[49] = 0.25*((A[ix+1][iy+1][iz+1] - A[ix+1][iy-1][iz+1]) + (- A[ix+1][iy+1][iz-1] + A[ix+1][iy-1][iz-1]));
    b[50] = 0.25*((A[ix  ][iy+2][iz+1] - A[ix  ][iy  ][iz+1]) + (- A[ix  ][iy+2][iz-1] + A[ix  ][iy  ][iz-1]));
    b[51] = 0.25*((A[ix+1][iy+2][iz+1] - A[ix+1][iy  ][iz+1]) + (- A[ix+1][iy+2][iz-1] + A[ix+1][iy  ][iz-1]));
    b[52] = 0.25*((A[ix  ][iy+1][iz+2] - A[ix  ][iy-1][iz+2]) + (- A[ix  ][iy+1][iz  ] + A[ix  ][iy-1][iz  ]));
    b[53] = 0.25*((A[ix+1][iy+1][iz+2] - A[ix+1][iy-1][iz+2]) + (- A[ix+1][iy+1][iz  ] + A[ix+1][iy-1][iz  ]));
    b[54] = 0.25*((A[ix  ][iy+2][iz+2] - A[ix  ][iy  ][iz+2]) + (- A[ix  ][iy+2][iz  ] + A[ix  ][iy  ][iz  ]));
    b[55] = 0.25*((A[ix+1][iy+2][iz+2] - A[ix+1][iy  ][iz+2]) + (- A[ix+1][iy+2][iz  ] + A[ix+1][iy  ][iz  ]));

    b[56] = 0.125*(((A[ix+1][iy+1][iz+1] - A[ix-1][iy+1][iz+1]) + (- A[ix+1][iy-1][iz+1] + A[ix-1][iy-1][iz+1])) + ((- A[ix+1][iy+1][iz-1] + A[ix-1][iy+1][iz-1]) + (+ A[ix+1][iy-1][iz-1] - A[ix-1][iy-1][iz-1])));
    b[57] = 0.125*(((A[ix+2][iy+1][iz+1] - A[ix  ][iy+1][iz+1]) + (- A[ix+2][iy-1][iz+1] + A[ix  ][iy-1][iz+1])) + ((- A[ix+2][iy+1][iz-1] + A[ix  ][iy+1][iz-1]) + (+ A[ix+2][iy-1][iz-1] - A[ix  ][iy-1][iz-1])));
    b[58] = 0.125*(((A[ix+1][iy+2][iz+1] - A[ix-1][iy+2][iz+1]) + (- A[ix+1][iy  ][iz+1] + A[ix-1][iy  ][iz+1])) + ((- A[ix+1][iy+2][iz-1] + A[ix-1][iy+2][iz-1]) + (+ A[ix+1][iy  ][iz-1] - A[ix-1][iy  ][iz-1])));
    b[59] = 0.125*(((A[ix+2][iy+2][iz+1] - A[ix  ][iy+2][iz+1]) + (- A[ix+2][iy  ][iz+1] + A[ix  ][iy  ][iz+1])) + ((- A[ix+2][iy+2][iz-1] + A[ix  ][iy+2][iz-1]) + (+ A[ix+2][iy  ][iz-1] - A[ix  ][iy  ][iz-1])));
    b[60] = 0.125*(((A[ix+1][iy+1][iz+2] - A[ix-1][iy+1][iz+2]) + (- A[ix+1][iy-1][iz+2] + A[ix-1][iy-1][iz+2])) + ((- A[ix+1][iy+1][iz  ] + A[ix-1][iy+1][iz  ]) + (+ A[ix+1][iy-1][iz  ] - A[ix-1][iy-1][iz  ])));
    b[61] = 0.125*(((A[ix+2][iy+1][iz+2] - A[ix  ][iy+1][iz+2]) + (- A[ix+2][iy-1][iz+2] + A[ix  ][iy-1][iz+2])) + ((- A[ix+2][iy+1][iz  ] + A[ix  ][iy+1][iz  ]) + (+ A[ix+2][iy-1][iz  ] - A[ix  ][iy-1][iz  ])));
    b[62] = 0.125*(((A[ix+1][iy+2][iz+2] - A[ix-1][iy+2][iz+2]) + (- A[ix+1][iy  ][iz+2] + A[ix-1][iy  ][iz+2])) + ((- A[ix+1][iy+2][iz  ] + A[ix-1][iy+2][iz  ]) + (+ A[ix+1][iy  ][iz  ] - A[ix-1][iy  ][iz  ])));
    b[63] = 0.125*(((A[ix+2][iy+2][iz+2] - A[ix  ][iy+2][iz+2]) + (- A[ix+2][iy  ][iz+2] + A[ix  ][iy  ][iz+2])) + ((- A[ix+2][iy+2][iz  ] + A[ix  ][iy+2][iz  ]) + (+ A[ix+2][iy  ][iz  ] - A[ix  ][iy  ][iz  ])));

    return b;
}

double* exact_diff(int shape1, int shape2, int shape3, double A[][shape1][shape2][shape3], int ix, int iy, int iz, double dx, double dy, double dz)
{
    double scale[8] = {1., dx, dy, dz, dx * dy, dx * dz, dy * dz, (dx * dy) * dz};

    double* b = (double*)malloc(64*sizeof(double));

    for(int l = 0 ; l < 8; l++)
    {
        b[8*l+0] = A[ix  ][iy  ][iz  ][l] * scale[l];
        b[8*l+1] = A[ix+1][iy  ][iz  ][l] * scale[l];
        b[8*l+2] = A[ix  ][iy+1][iz  ][l] * scale[l];
        b[8*l+3] = A[ix+1][iy+1][iz  ][l] * scale[l];
        b[8*l+4] = A[ix  ][iy  ][iz+1][l] * scale[l];
        b[8*l+5] = A[ix+1][iy  ][iz+1][l] * scale[l];
        b[8*l+6] = A[ix  ][iy+1][iz+1][l] * scale[l];
        b[8*l+7] = A[ix+1][iy+1][iz+1][l] * scale[l];
    }

    return b;
}
