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
#include "instabilities.h"


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
        
        ds = DSDesignSpaceByParsingStrings(strings, NULL, 5);
        expr = DSDesignSpaceEquations(ds);

        printf("System's equations are:\n");
        for (i = 0; i < DSDesignSpaceNumberOfEquations(ds); i++) {
                DSExpressionPrint(expr[i]);
                DSExpressionFree(expr[i]);
        }

        printf("This is the file \n");
        
        DSSecureFree(expr);
        printf("The total number of cases is: %u \n", DSDesignSpaceNumberOfCases(ds));
        printf("Number of valid cases is: %i\n", DSDesignSpaceNumberOfValidCases(ds));
//        printValidityOfCases(ds);


        // The idea is to develop a function that identifies all blow-ups and saves them in the data structure
        // 'DSDesignSpace' as a dictionary.

        // testing routines:
        if (1==0){
        // print validity of cases
        printValidityOfCases(ds);

        // testing function DSBlowUpSignatureForCaseNumber.
        const DSUInteger * signature_1;
        signature_1 = DSBlowUpSignatureForCaseNumber(1, ds->Xd);
        printf("the signature for the case 1 is: %u%u%u%u%u \n", signature_1[0], signature_1[1], signature_1[2],
        signature_1[3], signature_1[4]);

        // Implement and Test Pseudo-inverse for case 3.
        printf("Entering pseudo inverse calculation \n");
        DSUInteger * termSignature3;
        const DSMatrix * pInverse;

        termSignature3 = DSCaseSignatureForCaseNumber(3, ds->gma);
        printf("Term signature successfully constructed \n");

        DSCase * case3 = DSCaseWithTermsFromDesignSpace(ds, termSignature3, DSDesignSpaceCasePrefix(ds) );
//      printf("DSSSystemAd(case3->ssys) is NULL: %s\n", DSSSystemAd(case3->ssys)==NULL?"true":"false");

        DSSSystem * collapsedSystem = NULL;
        collapsedSystem = DSSSystemByRemovingAlgebraicConstraints(case3->ssys);
        printf("The collapsed Ad matrix is: \n");
        DSMatrixPrint(DSSSystemAd(collapsedSystem));

        pInverse = DSMatrixPseudoInverse(  DSSSystemAd(collapsedSystem) );
        printf("The Pseudoinverse of the collapsed matrix Ad is: \n");
        DSMatrixPrint(pInverse);

        // testing DSCaseIsCyclical(const DSCase *aCase)
        printf("Case 3 is cyclical: %s \n", DSCaseIsCyclical(case3)?"true":"false");
        printf("The variable pool Xi is: \n");
        printPool(case3->Xi);


        // random testing
        DSUInteger *termSignature_i = DSCaseSignatureForCaseNumber(9, ds->gma);
        DSCase * case_i = DSCaseWithTermsFromDesignSpace(ds, termSignature_i, DSDesignSpaceCasePrefix(ds) );
        DSCaseIsCyclical(case_i);

        };

        // call main function to identify and analyze unstable cases and save them in the dictionary

        DSDesignSpaceCalculateUnstableCases(ds);
        printf("Showing features of field unstableCases: \n");
        printDictionary(ds->unstableCases);

        return 0;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void printDictionary(DSDictionary *dictionary)
{
    DSUnstableCase *uCase;
    DSUInteger i;
    printf("The total number of elements in the dictionary is: %u \n", dictionary->count);
    printf("Cases in the dictionary: \n");
    if (dictionary->count != 0){
            for (i=0; i<dictionary->count; i++){
                printf("%s \n", dictionary->names[i]);

//                uCase = DSDictionaryValueForName(dictionary, dictionary->names[i]);

//                printf("showing vertices for case %u: \n", uCase->originalCase->caseNumber );
//                DSGetVertices(uCase->originalCase);
            }
    }


    // print properties of a specific case.
    const char *case_number = strdup("7");
    printf("showing properties for case 7 \n");
    uCase = DSDictionaryValueForName(dictionary, case_number);
    printf("Matrix U is: \n");
    DSMatrixPrint(uCase->originalCase->U);
    DSGetVertices(uCase);
//    printPool(case_i->Xd_b);
//    printPool(case_i->knifeEdge);

}


void printValidityOfCases(DSDesignSpace * ds)
{

    DSUInteger * termSignature, caseNumber, i;
    DSCase *aCase;
    DSVariablePool * parameters;


    printf("The total number of cases is: %u \n", DSDesignSpaceNumberOfCases(ds));
    for(i=0; i<DSDesignSpaceNumberOfCases(ds); i++){

    caseNumber = i + 1;
    termSignature = DSCaseSignatureForCaseNumber(caseNumber, ds->gma);
    aCase = DSCaseWithTermsFromDesignSpace(ds, termSignature, DSDesignSpaceCasePrefix(ds));

    printf("Case Nr. %u: The case is consistent?: %s. The conditions are valid?: %s. The case is cyclical?: %s. The case is valid?: %s \n",
    caseNumber, DSCaseIsConsistent(aCase)?"true":"false", DSCaseConditionsAreValid(aCase)?"true":"false",
    DSCaseIsCyclical(aCase)?"true":"false",
    DSCaseIsValid(aCase, true)?"true":"false");

//    printf("The blow up signature for the case is: " );
//    blowUpSignature(i+1);
//    printf("\n");

    printf("Showing parameter values\n");
    parameters = DSCaseConsistentParameterAndStateSet(aCase);
//    parameters = DSCaseValidParameterAndStateSet(aCase);
//    printPool(parameters);
    printf("-----------\n\n");


    }

}

void printPool(const DSVariablePool *pool)
{

        printf("showing %u elements for pool \n", pool->numberOfVariables);
        for(int j = 0; j < pool->numberOfVariables; j++)
        {
        printf("Name %s ; Value %f  \n", pool->variables[j]->name, pool->variables[j]->value );
        }
}


// 10. Get and Print Vertices (24)
void DSGetVertices(DSUnstableCase *uCase)
{
    DSVertices *vertices = NULL;
    char *xVariable, *yVariable;
    DSVariablePool * parameters;
    DSVariablePool *lowerBounds, *upperBounds;


    // allocate x and yVariables and assign.
    xVariable = DSSecureCalloc(sizeof(char), 100);
    yVariable = DSSecureCalloc(sizeof(char), 100);
    sprintf(xVariable, "S");
    sprintf(yVariable, "V1");


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

    // change certain parameter values: K12 and V2
    DSVariablePoolSetValueForVariableWithName(lowerBounds, "K12", 10);
    DSVariablePoolSetValueForVariableWithName(upperBounds, "K12", 10);

    DSVariablePoolSetValueForVariableWithName(lowerBounds, "V2", 0.01);
    DSVariablePoolSetValueForVariableWithName(upperBounds, "V2", 0.01);


    // parameters for the simplest case:
    if (1==0){
//
//        // allocate x and yVariables and assign.
//    xVariable = DSSecureCalloc(sizeof(char), 100);
//    yVariable = DSSecureCalloc(sizeof(char), 100);
//    sprintf(xVariable, "Vo");
//    sprintf(yVariable, "Vm");
//
//
//    // create DSVariablePool that contains boundaries.
//    lowerBounds = DSVariablePoolAlloc();
//    upperBounds = DSVariablePoolAlloc();
//
//    // Create variable pool
//    DSVariablePoolAddVariableWithName(lowerBounds, xVariable);
//    DSVariablePoolAddVariableWithName(lowerBounds, yVariable);
//    DSVariablePoolSetValueForVariableWithName(lowerBounds, xVariable, 1E-3);
//    DSVariablePoolSetValueForVariableWithName(lowerBounds, yVariable, 1E-3);
//
//    DSVariablePoolAddVariableWithName(upperBounds, xVariable);
//    DSVariablePoolAddVariableWithName(upperBounds, yVariable);
//    DSVariablePoolSetValueForVariableWithName(upperBounds, xVariable, 1E3);
//    DSVariablePoolSetValueForVariableWithName(upperBounds, yVariable, 1E3);
//
//    // call function to set all other parameter values.
//    char * variables[1] = {"K"};
//    int i;
//    for (i=0; i<1; i++){
//                DSVariablePoolAddVariableWithName(lowerBounds, variables[i]);
//                DSVariablePoolSetValueForVariableWithName(lowerBounds, variables[i], 1);
//
//                DSVariablePoolAddVariableWithName(upperBounds, variables[i]);
//                DSVariablePoolSetValueForVariableWithName(upperBounds, variables[i], 1);
//    }
//
////    // change certain parameter values: K12 and V2
////    DSVariablePoolSetValueForVariableWithName(lowerBounds, "K12", 10);
////    DSVariablePoolSetValueForVariableWithName(upperBounds, "K12", 10);
////
////    DSVariablePoolSetValueForVariableWithName(lowerBounds, "V2", 0.01);
////    DSVariablePoolSetValueForVariableWithName(upperBounds, "V2", 0.01);
//
    }

//    // place matrices U and Zeta in the right place:
//    uCase->originalCase->U = uCase->U;
//    uCase->originalCase->zeta = uCase->Zeta;


    if (DSCaseIsValidAtSlice(uCase->originalCase, lowerBounds, upperBounds, true) == true) {
            vertices = DSCaseVerticesFor2DSlice(uCase->originalCase, lowerBounds, upperBounds, xVariable, yVariable);
            //    vertices = DSUnstableCaseVerticesFor2DSlice(uCase, lowerBounds, upperBounds, xVariable, yVariable);
            printf("The vertices are: \n ");
            DSVerticesPrint(vertices);
    }

//    printf("The pseudoinverse of the matrix is: \n");
//    DSMatrixPrint(uCase->pInverse);

//    printf("B is: \n");
//    DSMatrixPrint( DSSSystemB(uCase->originalCase->ssys));
//
//    printf("The independent pool is: \n");
//    printPool(uCase->originalCase->Xi);
//
//    printf("The matrix Cd_knife is: \n");
//    DSMatrixPrint(uCase->Cd_knife);
//
//    printf("The matrix Ci_knife is: \n");
//    DSMatrixPrint(uCase->Ci_knife);
//
//    printf("The matrix delta_knife is: \n");
//    DSMatrixPrint(uCase->delta_knife);
//
//    printf("The matrix U is: \n");
//    DSMatrixPrint(uCase->U);
////
//    printf("The matrix Zeta is: \n");
//    DSMatrixPrint(uCase->Zeta);

//    printf("The equality 0 matrix is: \n");
//    DSMatrixPrint( DSMatrixArrayMatrix(  uCase->Equality, 0));
//
//    printf("The equality 1 matrix is: \n");
//    DSMatrixPrint( DSMatrixArrayMatrix(  uCase->Equality, 1));

//    printf("The solution is: \n");
//    DSSSystem *collapsedSystem = DSSSystemByRemovingAlgebraicConstraints(uCase->originalCase->ssys);
//    DSMatrixPrint( DSMatrixByMultiplyingMatrix(uCase->pInverse, DSSSystemAi(collapsedSystem)    ));
//
    printf("Xd_b is: \n ");
    printPool(uCase->Xd_b);
//
//    printf("knife edge is: \n");
//    printPool(uCase->knifeEdge);

}
