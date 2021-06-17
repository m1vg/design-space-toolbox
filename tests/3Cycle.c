//
//  designspacetest.c
//
//
//  Created by Miguel Valderrama on 08/27/18.
//
//
#include <stdio.h>
#include <string.h>
#include <designspace/DSStd.h>
#include <gsl/gsl_matrix_double.h>

int main(int argc, const char ** argv) {
        int i;
        char * strings[3] = {"\0"};
        strings[0] = strdup("x1. = a11 + 2*b31*x3 - b11*x2 - 2*b12*(x1^2)");
        strings[1] = strdup("x2. = b12*(x1^2) - b23*x2 - b22*x2");
        strings[2] = strdup("x3. = a31 + b23*x2 - b31*x3 - b33*x3");
        printf("the variable now is: %s\n", *strings);
    
        DSDesignSpace * ds;
        DSExpression ** expr = NULL;
        DSGMASystem *gma = NULL;

        
        ds = DSDesignSpaceByParsingStrings(strings, NULL, 3);
        DSDesignSpaceCalculateCyclicalCases(ds);
        printf("Number of valid cases is: %i\n", DSDesignSpaceNumberOfValidCases(ds));
        printf("DSDesignSpaceNumberOfValidCases passed!\n");
        
//        DSDesignSpaceFree(ds);
        printf("DSDesignSpaceFree passed!\n");


        // construct gma structure
        gma = DSGMASystemByParsingStrings(strings, NULL, 3);

        // print a Alpha
        printf("Showing Alpha: \n ");
        int rows, columns;
        for (int rows=0; rows < gma->alpha->rows; rows++)
        {
            for(int columns=0; columns < gma->alpha->columns; columns++)
            {
            printf("%f     ", gsl_matrix_get(gma->alpha->mat, rows, columns));
            }
        printf("\n");
        }

         // beta
         printf("Showing Beta: \n ");
        for (int rows=0; rows < gma->beta->rows; rows++)
        {
            for(int columns=0; columns < gma->beta->columns; columns++)
            {
            printf("%f     ", gsl_matrix_get(gma->beta->mat, rows, columns));
            }
        printf("\n");
        }

        // printing array of matrixes.

        printf("Printing Array of Matrizes for Gd \n ");

        for (int m=0; m < DSMatrixArrayNumberOfMatrices(gma->Gd); m++)
        {
        DSMatrix Matrix = *DSMatrixArrayMatrix(gma->Gd, m);
        printf("Printing Gd. Matrix Nr. %i \n", m);

        for (int rows=0; rows < Matrix.rows; rows++)
        {
            for(int columns=0; columns < Matrix.columns; columns++)
            {
            printf("%f     ", gsl_matrix_get(Matrix.mat, rows, columns));
            }
        printf("\n");
        }

        }

        // printing array of matrixes.

        printf("Printing Array of Matrizes for Gi \n ");

        for (int m=0; m < DSMatrixArrayNumberOfMatrices(gma->Gi); m++)
        {
        DSMatrix Matrix = *DSMatrixArrayMatrix(gma->Gi, m);
        printf("Printing Gi. Matrix Nr. %i \n", m);

        for (int rows=0; rows < Matrix.rows; rows++)
        {
            for(int columns=0; columns < Matrix.columns; columns++)
            {
            printf("%f     ", gsl_matrix_get(Matrix.mat, rows, columns));
            }
        printf("\n");
        }

        }

       // printing array of matrixes.

        printf("Printing Array of Matrizes for Hd \n ");

        for (int m=0; m < DSMatrixArrayNumberOfMatrices(gma->Hd); m++)
        {
        DSMatrix Matrix = *DSMatrixArrayMatrix(gma->Hd, m);
        printf("Printing Hd. Matrix Nr. %i \n", m);

        for (int rows=0; rows < Matrix.rows; rows++)
        {
            for(int columns=0; columns < Matrix.columns; columns++)
            {
            printf("%f     ", gsl_matrix_get(Matrix.mat, rows, columns));
            }
        printf("\n");
        }

        }

               // printing array of matrixes.

        printf("Printing Array of Matrices for Hi \n ");

        for (int m=0; m < DSMatrixArrayNumberOfMatrices(gma->Hi); m++)
        {
        DSMatrix Matrix = *DSMatrixArrayMatrix(gma->Hi, m);
        printf("Printing Hi. Matrix Nr. %i \n", m);

        for (int rows=0; rows < Matrix.rows; rows++)
        {
            for(int columns=0; columns < Matrix.columns; columns++)
            {
            printf("%f     ", gsl_matrix_get(Matrix.mat, rows, columns));
            }
        printf("\n");
        }

        }

        // Let's print list of variables. First Xd.

        printf(" Printing Xd. Dependent Variables of the model: \n");

        for(int j = 0; j < gma->Xd->numberOfVariables; j++)
        {
        printf("Name %s ; Value %f  \n", gma->Xd->variables[j]->name, gma->Xd->variables[j]->value );
        }

         // Let's print list of variables. First Xd_a

        printf(" Printing Xd_a. Algebraic Dependent Variables of the model: \n");

        for(int j = 0; j < gma->Xd_a->numberOfVariables; j++)
        {
        printf("Name %s ; Value %f  \n", gma->Xd_a->variables[j]->name, gma->Xd_a->variables[j]->value );
        }

        // Let's print list of variables. First Xd_t

        printf(" Printing Xd_a. Algebraic Dependent Variables of the model: \n");

        for(int j = 0; j < gma->Xd_t->numberOfVariables; j++)
        {
        printf("Name %s ; Value %f  \n", gma->Xd_t->variables[j]->name, gma->Xd_t->variables[j]->value );
        }

        // Let's print list of variables. First Xi

        printf(" Printing Xd_a. Algebraic Dependent Variables of the model: \n");

        for(DSUInteger j = 0; j < gma->Xi->numberOfVariables; j++)
        {
        printf("Name %s ; Value %f  \n", gma->Xi->variables[j]->name, gma->Xi->variables[j]->value );
        }

        // Print a DSDictionary containing valid Cases

        printf("Printing field validCases of the DSDesignSpace : \n");
        printf("number of cases %u: \n", ds->validCases->count);

        for(int w = 0; w<ds->validCases->count; w++)
        {
        printf( "names are: %s \n", ds->validCases->names[w]);
        }


        // No let's print the cyclical cases.

        printf("Printing field cyclicalCases of the DSDesignSpace : \n");
        printf("number of cases %u: \n", ds->cyclicalCases->count);

        for(int w = 0; w<ds->cyclicalCases->count; w++)
        {
        printf( "names are: %s \n", ds->cyclicalCases->names[w]);
        }


        printf("Printing CasePrefix: %c \n", ds->modifierFlags);
        printf("caca %u \n",3);

        DSDesignSpaceFree(ds);




        return 0;
}
