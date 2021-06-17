//
//  designspacetest.c
//
//
//  Created by Miguel Valderrama on 10/23/18.
//
//
#include <stdio.h>
#include <glpk.h>
#include <string.h>

#include "glpk_test.h"

int main(void) {

    glp_prob *lp;
    int ia[1+1000], ja[1+1000], i;
    double ar[1+1000];
    lp = glp_create_prob();
    glp_set_obj_dir(lp, GLP_MIN);
    glp_add_rows(lp, 5);

    char * strings[12] = {"\0"};
    strings[0] = strdup("X1");
    strings[1] = strdup("X2");
    strings[2] = strdup("k2");
    strings[3] = strdup("v1");
    strings[4] = strdup("v2");
    strings[5] = strdup("vs");
    strings[6] = strdup("ks");
    strings[7] = strdup("k12");
    strings[8] = strdup("s");
    strings[9] = strdup("k1");
    strings[10] = strdup("k2");
    strings[11] = strdup("z1");

    //set row boundaries.
    // The first three rows (corresponding to the boundaries) are set to have a boundary of GLP_UP of zero.
    for (i=1; i<4; i++){
        glp_set_row_bnds(lp, i, GLP_UP, 0.0, 0.0);
    }
    // correspond to the fourth condition, where the value of X2 is defined.
    glp_set_row_bnds(lp, 4, GLP_FX, 0.0, 0.0);

    // correspond to the fifth condition, the knife edge condition.
    glp_set_row_bnds(lp, 5, GLP_FR, 0, 0.0); // (lp, 5, GLP_FX, 0, 0.0),

    // correspond to lower than knife edge
//    glp_set_row_bnds(lp, 5, GLP_UP, 0, 0.0);
//
//    glp_set_row_bnds(lp, 5, GLP_LO, 0, 0.0);


    // Add 12 columns corresponding to the variables.
    glp_add_cols(lp, 12);

    // now set boundaries for columns
    for(i=1; i<12; i++){
        glp_set_col_bnds(lp, i, GLP_DB, -3, 3);
    }

    // bounds for X1. Set blow up or blow down.
    glp_set_col_bnds(lp, 1, GLP_FX, -12, -12);

    // constraint z
    glp_set_col_bnds(lp, 12, GLP_LO, -1, 0);
    glp_set_obj_coef(lp, 12, 1.0);

    // fill values for matrices

    ia[1] = 1, ja[1] = 2, ar[1] = 2;
    ia[2] = 1, ja[2] = 7, ar[2] = -2;
    ia[3] = 1, ja[3] = 8, ar[3] = -2;
    ia[4] = 1, ja[4] = 9, ar[4] = 2;
    ia[5] = 1, ja[5] = 12, ar[5] = -1;

    ia[6] = 2, ja[6] = 1, ar[6] = -1;
    ia[7] = 2, ja[7] = 10, ar[7] = 1;
    ia[8] = 2, ja[8] = 12, ar[8] = -1;


    ia[9] = 3, ja[9] = 2, ar[9] = 1;
    ia[10] = 3, ja[10] = 11, ar[10] = -1;
    ia[11] = 3, ja[11] = 12, ar[11] = -1;


    ia[12] = 4, ja[12] = 2, ar[12] = -1;
    ia[13] = 4, ja[13] = 3, ar[13] = 1;
    ia[14] = 4, ja[14] = 4, ar[14] = 1;
    ia[15] = 4, ja[15] = 5, ar[15] = -1;


    ia[16] = 5, ja[16] = 4, ar[16] = -1;
    ia[17] = 5, ja[17] = 5, ar[17] = 2;
    ia[18] = 5, ja[18] = 6, ar[18] = -1;
    ia[19] = 5, ja[19] = 7, ar[19] = 2;
    ia[20] = 5, ja[20] = 8, ar[20] = 2;
    ia[21] = 5, ja[21] = 9, ar[21] = -2;
    ia[22] = 5, ja[22] = 11, ar[22] = -2;

//    ia[23] = 5, ja[23] = 12, ar[23] = 1;

    glp_load_matrix(lp, 22, ia, ja, ar);

    printf("matrix loaded \n");
    glp_simplex(lp, NULL);

    printf("the parameter values are: \n");
    for(i=1; i<13; i++){

        printf("%s = %g \n",strings[i-1], glp_get_col_prim(lp, i) );


    }

    if(glp_get_obj_val(lp) <= -1E-14 && glp_get_prim_stat(lp) == GLP_FEAS){
        printf("conditions are feasible \n");
     }
     else{
        printf("conditions are not feasible \n");

     }

    // now let's print the results.


    printf("the value of left-hand side of the first condition is: %f \n", 2*glp_get_col_prim(lp, 7)
                                                                        -2*glp_get_col_prim(lp, 9)
                                                                        +2*glp_get_col_prim(lp, 8)
                                                                        -2*glp_get_col_prim(lp, 2));

    printf("the value of left-hand side of the second condition is: %f \n", -glp_get_col_prim(lp, 10)
                                                                        +glp_get_col_prim(lp, 1));

    printf("the value of left-hand side of the third condition is: %f \n", glp_get_col_prim(lp, 11)
                                                                        -glp_get_col_prim(lp, 2));

    printf("the value of left-hand side of the fourth condition (eq. for X2) is: %f \n", glp_get_col_prim(lp, 2)
                                                                        -glp_get_col_prim(lp, 3)
                                                                        -glp_get_col_prim(lp, 4)
                                                                        +glp_get_col_prim(lp, 5));

    printf("the value of left-hand side of the fifth condition (knife-edge) is: %f \n", glp_get_col_prim(lp, 4)
                                                                        -glp_get_col_prim(lp, 6)
                                                                        +2*glp_get_col_prim(lp, 7)
                                                                        +2*glp_get_col_prim(lp, 8)
                                                                        -2*glp_get_col_prim(lp, 9)
                                                                        +2*glp_get_col_prim(lp, 5)
                                                                        -2*glp_get_col_prim(lp, 11)
                                                                        -2*glp_get_col_prim(lp, 4));


        return 0;
}


