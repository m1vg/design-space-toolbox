//
//  designspacetest.c
//
//
//  Created by Miguel Valderrama on 10/23/18.
//
//
#include <stdio.h>
#include <string.h>
#include <designspace/DSStd.h>
#include "designspacegtest.h"

int main(int argc, const char ** argv) {
        int i;

        // simplest case
//        char * strings[2] = {"\0"};
//        strings[0] = strdup("X1. = Vo - Vm*(D^-1)*K^(-1)*X1");
//        strings[1] = strdup("D = 1 + (K^-1)*X1");


        // Motif1
//        char * strings[4] = {"\0"};
//        strings[0] = strdup("R. = (100^-1)*(DR^-1) + (DR^-1)*(KA^-3)*(A^3) - 1.0000002*R");
//        strings[1] = strdup("A. = (DA^-1) + (KAA^-3)*(A^3)*(DA^-1) + (100^-1)*(KR^-3)*(R^3)*(DA^-1)  - 1.0000001*A");
//        strings[2] = strdup("DR = 1 + (KA^-3)*(A^3)");
//        strings[3] = strdup("DA = 1 + (KAA^-3)*(A^3) + (KR^-3)*(R^3)");

//        // positive feedback, two pools.
        char * strings[5] = {"\0"};
        strings[0] = strdup("X1. = Vs*(KS^-2)*(S^2)*(K12^-2)*(X2^2)*(Ds^-1) - V1*(K1^-1)*(X1^1)*(D1^-1)");
        strings[1] = strdup("X2. = V1*(K1^-1)*(X1^1)*(D1^-1) - V2*(K2^-1)*(X2^1)*(D2^-1)");
        strings[2] = strdup("Ds = 1 + (KS^-2)*(S^2)*(K12^-2)*(X2^2)");
        strings[3] = strdup("D1 = 1 + (K1^-1)*(X1^1)");
        strings[4] = strdup("D2 = 1 + (K2^-1)*(X2^1)");

        DSDesignSpace * ds;
        DSExpression ** expr = NULL;
        DSExpression * anExpression;
        DSUInteger * termSignature, caseNumber;
        DSCase *aCase;
        DSVariablePool * parameters;
        
        ds = DSDesignSpaceByParsingStrings(strings, NULL, 5);
        expr = DSDesignSpaceEquations(ds);

        printf("System's equations are:\n");
        for (i = 0; i < DSDesignSpaceNumberOfEquations(ds); i++) {
                DSExpressionPrint(expr[i]);
                DSExpressionFree(expr[i]);
        }
        
        DSSecureFree(expr);
        printf("Number of valid cases is: %i\n", DSDesignSpaceNumberOfValidCases(ds));
        printValidityOfCases(ds);


//        printf("Showing parameter values for case 8 \n");
//
//        caseNumber = 8;
//        termSignature = DSCaseSignatureForCaseNumber(caseNumber, ds->gma);
//        aCase = DSCaseWithTermsFromDesignSpace(ds, termSignature, DSDesignSpaceCasePrefix(ds));
//        parameters = DSCaseConsistentParameterAndStateSet(aCase);
//        printPool(parameters);

        return 0;
}



void printValidityOfCases(DSDesignSpace * ds){

    DSUInteger * termSignature, caseNumber, i;
    DSCase *aCase;
    DSVariablePool * parameters;


    printf("The total number of cases is: %u \n", DSDesignSpaceNumberOfCases(ds));
    for(i=0; i<DSDesignSpaceNumberOfCases(ds); i++){

    caseNumber = i + 1;
    termSignature = DSCaseSignatureForCaseNumber(caseNumber, ds->gma);
    aCase = DSCaseWithTermsFromDesignSpace(ds, termSignature, DSDesignSpaceCasePrefix(ds));

    printf("Case Nr. %u: The case is consistent?: %s. The conditions are valid?: %s. The case is valid?: %s \n",
    caseNumber, DSCaseIsConsistent(aCase)?"true":"false", DSCaseConditionsAreValid(aCase)?"true":"false",
    DSCaseIsValid(aCase, true)?"true":"false");

//    printf("The blow up signature for the case is: " );
//    blowUpSignature(i+1);
//    printf("\n");



    printf("Showing parameter values\n");
    parameters = DSCaseConsistentParameterAndStateSet(aCase);
//    parameters = DSCaseValidParameterAndStateSet(aCase);
    printPool(parameters);
    printf("-----------\n\n");


    }

}

void printPool(DSVariablePool *pool){

        for(int j = 0; j < pool->numberOfVariables; j++)
        {
        printf("Name %s ; Value %f  \n", pool->variables[j]->name, pool->variables[j]->value );
        }
}

//DSUInteger * blowUpSignature(const DSUInteger caseNumber, const DSUInteger numberOfDependentVariables){
//
//      DSUInteger num, *signature = NULL, i;
//
//      num = caseNumber - 1;
//      for (i = 0; i < 3; i++){
//            signature[] = num % 2 +1;
//            num /= 2;
//            printf("%u", signature);
//      }
//}
