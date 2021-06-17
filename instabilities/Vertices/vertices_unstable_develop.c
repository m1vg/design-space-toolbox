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
#include <glpk.h>
#include "vertices_unstable.h"


int main(int argc, const char ** argv) {
        int i;

        // simplest case
//        char * strings[2] = {"\0"};
//        strings[0] = strdup("X1. = Vo - Vm*(D^-1)*K^(-1)*X1");
//        strings[1] = strdup("D = 1 + (K^-1)*X1");

        // cyclical case
//        char * strings[3] = {"\0"};
//        strings[0] = strdup("X1. = a11 + 2*b31*X3 - b11*X1 - 2*b12*(X1^2)");
//        strings[1] = strdup("X2. = b12*(X1^2) - b23*X2 - b22*X2");
//        strings[2] = strdup("X3. = a31 + b23*X2 - b31*X3 - b33*X3");


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
        DSVariablePool *lowerBounds, *upperBounds;
        char  *xVariable, *yVariable;
        DSVertices * vertices = NULL;
        int casenumber = 3;
        
        ds = DSDesignSpaceByParsingStrings(strings, NULL, 5);
        expr = DSDesignSpaceEquations(ds);

        printf("System's equations are:\n");
        for (i = 0; i < DSDesignSpaceNumberOfEquations(ds); i++) {
                DSExpressionPrint(expr[i]);
                DSExpressionFree(expr[i]);
        }

        DSSecureFree(expr);
        printf("The total number of cases is: %u \n", DSDesignSpaceNumberOfCases(ds));
        printf("Number of valid cases is: %i\n", DSDesignSpaceNumberOfValidCases(ds));


        // let's create a case and then print vertices for that case.
        termSignature = DSCaseSignatureForCaseNumber(casenumber, ds->gma);
        printf("Term signature successfully constructed \n");

        aCase = DSCaseWithTermsFromDesignSpace(ds, termSignature, DSDesignSpaceCasePrefix(ds));
        printf("Case successfully constructed \n");

//        // construct boundary matrices
        DSCaseRecalculateBoundaryMatrices(aCase);
        printf("Boundary matrices successfully calculated \n");


        if (DSuCaseIsValid(aCase, true) == true){
                            DSUnstableCaseAddBoundaryMatrices(aCase);
                            DSSSystemSetIsUnstable((DSSSystem *)DSCaseSSystem(aCase), true);
        }

        printf("The Matrix U is: \n");
        DSMatrixPrint(aCase->U);

        printf("The Matrix Zeta is: \n");
        DSMatrixPrint(aCase->zeta);

//        printf("entering the loop \n");
//        loopValid(aCase);

        printf("Generting Case Boundaries \n");
        DSExpression ** caseboundaries = DSCaseBoundaries(aCase);
        DSExpressionPrint(caseboundaries[1]);

        printf("printing variable pool Xd_b of aCase: \n");
        printPool(aCase->ssys->Xd_b);

        DSSSystem *ssys = aCase->ssys;
        DSSSystem *ssys_copy = DSSSystemCopy(ssys);
        printf("The number of variables of Xd_b in the original syste is %u, and in the copy is %u \n",
                DSVariablePoolNumberOfVariables(ssys->Xd_b),
                DSVariablePoolNumberOfVariables(ssys_copy->Xd_b));




        printf("starting to process aCase \n");
        DSSSystem *collapsed = DSSSystemByRemovingAlgebraicConstraints(ssys);
        printf("aCase was sucessfully processed \n");
        printf("Number of variables in Xd_b in collapsed is %u \n", DSVariablePoolNumberOfVariables(collapsed->Xd_b));


        printf("testing function DSuSSystemSolution() \n");
        DSExpression ** expression = DSuSSystemLogarithmicSolution(aCase->ssys);
        printf("test passed \n");
        printf("printing 0: \n");

        DSExpressionPrint(expression[1]);
        printf("printing test passed \n");

//
//        DSUInteger numberOfBlowing = DSVariablePoolNumberOfVariables(ssys->Xd_b);
//        DSUInteger *indices = DSVariablePoolIndicesOfSubPool(ssys->Xd, ssys->Xd_b);










        // allocate x and yVariables and assign.
        xVariable = DSSecureCalloc(sizeof(char), 100);
        yVariable = DSSecureCalloc(sizeof(char), 100);
        sprintf(xVariable, "S");
        sprintf(yVariable, "V1");

//        xVariable = "V1";
//        yVariable = "S";

        // create DSVariablePool that contains boundaries.
        lowerBounds = DSVariablePoolAlloc();
        upperBounds = DSVariablePoolAlloc();

        // Create variable pool
        DSVariablePoolAddVariableWithName(lowerBounds, xVariable);
        DSVariablePoolAddVariableWithName(lowerBounds, yVariable);
        DSVariablePoolSetValueForVariableWithName(lowerBounds, xVariable, 1E-3);
        DSVariablePoolSetValueForVariableWithName(lowerBounds, yVariable, 1E-3);


        DSVariablePoolAddVariableWithName(upperBounds, xVariable);
        DSVariablePoolAddVariableWithName(upperBounds, yVariable);
        DSVariablePoolSetValueForVariableWithName(upperBounds, xVariable, 1E3);
        DSVariablePoolSetValueForVariableWithName(upperBounds, yVariable, 1E3);


        // call function to set all other parameter values.
        setRemainingParameterValues(lowerBounds, upperBounds);

//        if (DSCaseIsValidAtSlice(aCase, lowerBounds, upperBounds, false) == true) {
//                vertices = DSCaseVerticesFor2DSlice(aCase, lowerBounds, upperBounds, xVariable, yVariable);
//                printf("Showing %u Vertices for case %u \n", vertices->numberOfVertices, casenumber);
//                DSVerticesPrint(vertices);
//        } else{
//                printf("Case %u is not valid in slice \n", casenumber);
//
//        }

        return 0;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void printDictionary(DSDictionary *dictionary){

    DSUInteger i;
    printf("The total number of elements in the dictionary is: %u \n", dictionary->count);
    printf("Cases in the dictionary: \n");
    if (dictionary->count != 0){
            for (i=0; i<dictionary->count; i++){
                printf("%s \n", dictionary->names[i]);
            }
    }


    // print properties of a specific case.
    const char *case_number = strdup("4");
    printf("showing Xd_b pool for case 4 \n");
    DSUnstableCase *case_i = DSDictionaryValueForName(dictionary, case_number);
    printPool(case_i->Xd_b);
    printPool(case_i->knifeEdge);

}


void printPool(const DSVariablePool *pool){

        printf("showing %u elements for pool \n", pool->numberOfVariables);
        for(int j = 0; j < pool->numberOfVariables; j++)
        {
        printf("Name %s ; Value %f  \n", pool->variables[j]->name, pool->variables[j]->value );
        }
}

void setRemainingParameterValues(DSVariablePool *lowerBounds, DSVariablePool *upperBounds){

    char * variables[6] = {"K12",
                            "V2",
                            "KS",
                            "K2",
                            "Vs",
                            "K1"};
    int i;


    for (i=0; i<6; i++){


                DSVariablePoolAddVariableWithName(lowerBounds, variables[i]);
                DSVariablePoolSetValueForVariableWithName(lowerBounds, variables[i], 1);


                DSVariablePoolAddVariableWithName(upperBounds, variables[i]);
                DSVariablePoolSetValueForVariableWithName(upperBounds, variables[i], 1);

    }
}

void loopValid(DSCase *aCase){

    DSUnstableCase * unstableCase = DSSecureCalloc(sizeof(DSUnstableCase), 1);
    unstableCase->originalCase = DSCaseCopy(aCase);

    // perform main calculations
    // this function should populate Xd_b, Xd_e, pInverse and delta_unstable
    DSUnstableCaseIdentifyBlowingDependentVariables(unstableCase);

    // this function should set fields KnifeEdge, Cd_unstable, Ci_unstable and delta_unstable.
    //    DSUnstableCaseExpandConditionMatrices(unstableCase);

    DSUInteger *bSignature;
    DSUInteger numberOfCases, i;
    DSVariablePool *Xd_b = unstableCase->Xd_b;
    bool isValid;

    printf("showing Xd_b: \n");
    printPool(Xd_b);

    printf("showing Xd_e: \n");
    printPool(unstableCase->Xd_e);

    numberOfCases = pow(2, Xd_b->numberOfVariables);

    for (i=1; i <= numberOfCases; i++){
            // This is an array of length Xb->numberOfVariables containing values of 1 for blow-down and 2 for blow up.
            bSignature = DSBlowUpSignatureForCaseNumber(i, Xd_b);
            // construct case and return validity
            isValid = DSUnstableCaseConditionsAreValid(unstableCase, bSignature);
            DSSecureFree(bSignature);

            // If it has a strict solution, assign matrices and break.
            if ( isValid == true ){
                    DSUnstableCaseCreateBoundaryMatrices(unstableCase);
                    aCase->U = DSMatrixCopy(unstableCase->U);
                    aCase->zeta = DSMatrixCopy(unstableCase->Zeta);
                    printf("Subcase %u is valid \n ", i);


//                    break;
            }
    }
 }





