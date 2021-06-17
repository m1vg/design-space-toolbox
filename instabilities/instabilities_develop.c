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
//        printf("Showing features of field unstableCases: \n");
//        printDictionary(ds->unstableCases);

        return 0;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void printDictionary(DSDictionary *dictionary)
{

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


// Functions to be implemented in the Toolbox:

// 1. Blowup signature (1)
DSUInteger * DSBlowUpSignatureForCaseNumber(const DSUInteger caseNumber, const DSVariablePool *freeDependentVariable)
{

      DSUInteger num, *signature = NULL, i, numberOfFreeVariables;

      numberOfFreeVariables = freeDependentVariable->numberOfVariables;
      signature = DSSecureMalloc(sizeof(DSUInteger)*numberOfFreeVariables);

      num = caseNumber - 1;
      for (i = 0; i < numberOfFreeVariables; i++){
            signature[i] = num % 2 +1;
            num /= 2;
      }

      return signature;
}


// 2.1 DSUnstableCaseGetSubSetPseudoInverse (2)
const DSMatrix *DSUnstableCaseGetSubSetPseudoInverse(const DSVariablePool *Xd_t, const DSVariablePool *Xd_e,
                                                     const DSVariablePool *Xd_b,
                                                     const DSMatrix *Ad)
{

    // The idea of this function is to control the way how the pseudo-inverse is being calculated and how it affects the solution
    // of the pool Xd_e. The idea is to calculate the pseudo-inverse from a subset of the matrix Ad that only contains
    // rows for Xd_e. Rows for Xd_b should be set to zero.

    DSUInteger *indices, nRows = Xd_b->numberOfVariables, i, j;
    DSMatrix *Ad_subset;
    const DSMatrix *pInverse = NULL;

    Ad_subset = DSMatrixCopy(Ad);

    // get indices of Xd_b
    indices = DSVariablePoolIndicesOfSubPool(Xd_t, Xd_b);

    // set rows contained in indices to zero.
    for (i=0; i<nRows; i++){
            for(j=0; j<DSMatrixColumns(Ad); j++){
                    DSMatrixSetDoubleValue(Ad_subset, indices[i], j, 0);
            }
    }

    pInverse = DSMatrixPseudoInverse(Ad_subset);
    DSSecureFree(indices);

    return pInverse;

}

// 2. Pseudoinverse (3)
const DSMatrix * DSMatrixPseudoInverse(const DSMatrix *matrix)
{

    DSMatrix *pMatrix = NULL;
    DSInteger n_row = DSMatrixRows(matrix);
	DSInteger n_col = DSMatrixColumns(matrix);
    DSUInteger i, j;
    DSMatrix *U, *V, *S, *invS, *trU, *pA, *VinvS;  //*trV,
    DSMatrixArray *array;

    if (matrix == NULL)
        goto bail;

    array = DSMatrixSVD(matrix);

    // unpack matrices
    S = DSMatrixArrayMatrix(array, 0);
    U = DSMatrixArrayMatrix(array, 1);
    V = DSMatrixArrayMatrix(array, 2);

    // calculate inverse of S, invS
    invS = DSMatrixCalloc(n_row, n_row);
    for (i = 0; i < n_row; i++){
        if (DSMatrixDoubleValue(S, 0, i) != 0){
            DSMatrixSetDoubleValue(invS, i, i, 1/DSMatrixDoubleValue(S, 0, i));
        }
    }

    // Transpose V
//    trV = DSMatrixTranspose(V);

    // Transpose U
    trU = DSMatrixTranspose(U);

    // Calculate pseudo inverse of matrix
    VinvS = DSMatrixByMultiplyingMatrix(V, invS);
    pMatrix = DSMatrixByMultiplyingMatrix(VinvS, trU);

    // free matrices
    DSMatrixFree(S);
	DSMatrixFree(U);
	DSMatrixFree(V);
	DSMatrixFree(invS);
//	DSMatrixFree(trV);
	DSMatrixFree(trU);
	DSMatrixFree(VinvS);

//	printf("The pseudoinverse is: \n");
//	DSMatrixPrint(pMatrix);

bail:
    return pMatrix;
}

// 3.14 DSUnstableCaseDetermineBlowingBehavior(lp, uCase, bSignature) (4)
void DSUnstableCaseDetermineBlowingBehavior(glp_prob *lp, DSUnstableCase *uCase, const DSUInteger *bSignature)
{

    const char *name;
    DSVariablePool *Xd_b = uCase->Xd_b;
    DSVariablePool *knifeEdge = DSVariablePoolAlloc();
    DSUInteger i, numberOfKnifes = Xd_b->numberOfVariables, k, j;
    DSUInteger numberOfXd = uCase->originalCase->Xd->numberOfVariables;
    double value;
    DSMatrix *gAi = DSMatrixArrayMatrix(uCase->Knife, 0);
    DSMatrix *gB = DSMatrixArrayMatrix(uCase->Knife, 1);
    DSMatrix *Cd_knife, *Ci_knife, *delta_knife;


    k = DSMatrixRows(uCase->originalCase->Ci) + uCase->Xd_b->numberOfVariables + uCase->Xd_e->numberOfVariables;

    // initialize matrices Ci_knife, delta_knife & Cd_knife.
    Ci_knife = DSMatrixCopy(gAi);
    Cd_knife = DSMatrixCalloc(numberOfKnifes, numberOfXd);
    delta_knife = DSMatrixCopy(gB);


    //loop over bSignature to assign blowing behavior.
    for (i=0; i<numberOfKnifes; i++){

        name = DSVariablePoolAllVariableNames(uCase->Xd_b)[i];
        DSVariablePoolAddVariableWithName(knifeEdge, name);
        value = glp_get_row_prim(lp, k+i+1) + DSMatrixDoubleValue(gB, i, 0);
        printf("value is the result of %f - %f \n", glp_get_row_prim(lp, k+i+1), DSMatrixDoubleValue(gB, i, 0) );
        DSVariablePoolSetValueForVariableWithName(knifeEdge, name, value);

        // set blowing behavior.
        if (bSignature[i]==1){
                DSVariablePoolSetValueForVariableWithName(Xd_b, name, 1E12);
        } else {
                DSVariablePoolSetValueForVariableWithName(Xd_b, name, 1E-12);
        }

        // adjust matrices Ci_knife and delta_knife.
        if (value < 0){

            for (j=0 ; j<DSMatrixColumns(gAi); j++){
                    DSMatrixSetDoubleValue(Ci_knife, i, j, DSMatrixDoubleValue(Ci_knife, i, j)*-1.0);
            }

            DSMatrixSetDoubleValue(delta_knife, i, 0, DSMatrixDoubleValue(delta_knife, i, 0)*-1.0);

        }
    }

    uCase->knifeEdge = knifeEdge;
    uCase->Cd_knife = Cd_knife;
    uCase->Ci_knife = Ci_knife;
    uCase->delta_knife = delta_knife;
}

// 3.13 DSUnstableCasePrintLinearProblem(lp, uCase, bSignature) (5)
void DSUnstableCasePrintLinearProblem(glp_prob *lp, const DSUnstableCase *uCase, const DSUInteger *bSignature)
{

    DSUInteger i;
    DSUInteger numberOfXd = uCase->originalCase->Xd->numberOfVariables;
    DSUInteger numberOfXi = uCase->originalCase->Xi->numberOfVariables;
    DSUInteger numberOfKnifes = uCase->Xd_b->numberOfVariables;
    DSUInteger k;
    const DSVariablePool *Xd = uCase->originalCase->Xd;
    const DSVariablePool *Xi = uCase->originalCase->Xi;
    DSVariablePool *Knifes = uCase->Xd_b;
    char *name;
    double value;

    printf("The values of the dependent variables for case %u are: \n", uCase->originalCase->caseNumber);

    for(i=0; i<numberOfXd; i++){

            name = DSVariableName(DSVariablePoolVariableAtIndex(Xd, i));
            value = glp_get_col_prim(lp, i+1) ;
            printf("%s = %f \n", name, value);
    }

    printf("The independent parameters are: \n");

    for(i=0; i<numberOfXi; i++){

            name = DSVariableName(DSVariablePoolVariableAtIndex(Xi, i));
            value = glp_get_col_prim(lp, numberOfXd+i+1) ;
            printf("%s = %f \n", name, value);
    }

    printf("The values of the knifes are: \n");

    k = DSMatrixRows(uCase->originalCase->Ci) + uCase->Xd_b->numberOfVariables + uCase->Xd_e->numberOfVariables;

    for(i=0; i<numberOfKnifes; i++){

            name = DSVariableName(DSVariablePoolVariableAtIndex(Knifes, i));
            value = glp_get_row_prim(lp, k+i+1) ;
            printf("%s = %f \n", name, value);
    }


}


// 3.12 dsUnstableCaseLinearProblemForMatrices (6)
glp_prob * dsUnstableCaseLinearProblemForMatrices(const DSMatrix *A_, const DSMatrix *B_, DSUnstableCase *uCase,
                                                  const DSUInteger * bSignature)
{
        glp_prob *linearProblem = NULL;
        int * ia = NULL, *ja = NULL;
        double *ar = NULL;
        DSUInteger i, numberOfXi, numberOfBoundaries;
        DSUInteger numberOfEqualities = uCase->Xd_e->numberOfVariables , numberOfKnifes = uCase->Xd_b->numberOfVariables;
        DSUInteger numberOfOriginalBounds = DSMatrixRows(uCase->originalCase->Cd) + numberOfKnifes;
        DSMatrixArray *lp_AB1 = NULL, *lp_AB2 = NULL;
        DSUInteger blowingIndex;
        char * name;
        const DSVariablePool *Xd = uCase->originalCase->Xd;
        DSMatrix *A, *B;

        glp_term_out(GLP_OFF);
        linearProblem = glp_create_prob();
        if (linearProblem == NULL) {
                DSError(M_DS_NULL ": Linear problem is NULL", A_DS_ERROR);
                goto bail;
        }
        if (numberOfKnifes == 0) {
                DSError(M_DS_NULL ": The case does not have any knife-edge constraints", A_DS_ERROR);
                goto bail;
        }

        //Then add Equality constraints if Xd_e is not zero.
        if (uCase->Xd_e->numberOfVariables !=0){
            lp_AB1 = dsUnstableCaseLinearProblemAddEqualityConstraints(A_, B_, uCase);
            // Then add Knife edge constraints
            lp_AB2 = dsUnstableCaseLinearProblemAddKnifeEdgeConditions(DSMatrixArrayMatrix(lp_AB1, 0),
                                                                       DSMatrixArrayMatrix(lp_AB1, 1), uCase);
         }else{
            // Then add Knife edge constraints
            lp_AB2 = dsUnstableCaseLinearProblemAddKnifeEdgeConditions(A_, B_, uCase);
         }

        A = DSMatrixArrayMatrix(lp_AB2, 0);
        B = DSMatrixArrayMatrix(lp_AB2, 1);

//        printf("Showing Matrix A: \n");
//        DSMatrixPrint(A);
//
//        printf("Showing Matrix B: \n");
//        DSMatrixPrint(B);

        numberOfXi = DSMatrixColumns(A);
        numberOfBoundaries = DSMatrixRows(A);

        ia = DSMatrixRowsForGLPK(A);
        ja = DSMatrixColumnsForGLPK(A);
        ar = DSMatrixDataForGLPK(A);

        glp_add_rows(linearProblem, numberOfBoundaries);
        glp_add_cols(linearProblem, numberOfXi);
        glp_set_obj_dir(linearProblem, GLP_MIN);
        glp_load_matrix(linearProblem, numberOfBoundaries*numberOfXi, ia, ja, ar);

        // set bounds for original boundaries and additional constrains coming from blow up/down assumption.
        for (i = 0; i < numberOfOriginalBounds; i++) {
                glp_set_row_bnds(linearProblem, i+1, GLP_UP, 0.0,
                                 DSMatrixDoubleValue(B, i, 0));
        }
        // set bounds for equality constraints
        for (i = numberOfOriginalBounds; i < numberOfOriginalBounds + numberOfEqualities; i++) {
                glp_set_row_bnds(linearProblem, i+1, GLP_FX, DSMatrixDoubleValue(B, i, 0),
                                 DSMatrixDoubleValue(B, i, 0));
        }
        // set bounds for knife-edges.
        for (i = numberOfOriginalBounds + numberOfEqualities;
             i < numberOfOriginalBounds + numberOfEqualities + numberOfKnifes; i++) {

                glp_set_row_bnds(linearProblem, i+1, GLP_FR, 0.0, 0.0);
        }

        // set bounds for parameter values:
        for (i = 0; i < numberOfXi; i++)
                glp_set_col_bnds(linearProblem, i+1, GLP_DB, -3.0, 3.0);

        //set value for blowing variables. Loop over number of Xd_b
        for(i=0; i<numberOfKnifes; i++){

                name = DSVariableName(DSVariablePoolVariableAtIndex(uCase->Xd_b, i));
                blowingIndex = DSVariablePoolIndexOfVariableWithName(Xd, name) ;
                if (bSignature[i] == 1){
                    glp_set_col_bnds(linearProblem, blowingIndex+1, GLP_FX, 12, 12);
                }else{
                    glp_set_col_bnds(linearProblem, blowingIndex+1, GLP_FX, -12, -12);
                }
        }

        if (ia != NULL)
                DSSecureFree(ia);
        if (ja != NULL)
                DSSecureFree(ja);
        if (ar != NULL)
                DSSecureFree(ar);
        if (uCase->Xd_e->numberOfVariables != 0){
                DSMatrixArrayFree(lp_AB1);
        }
        DSMatrixArrayFree(lp_AB2);
bail:
        return linearProblem;
}

// 3.11 DSUnstableCaseGetEqualityAndKnife(pInverse, gaussArray, uCase); (7)
void DSUnstableCaseGetEqualityAndKnife(const DSMatrix *pInverse, const DSMatrixArray *gaussArray, DSUnstableCase *uCase)
{

    //This function should populate fields uCase->Equality and uCase->knife.
    DSSSystem * o_ssys = uCase->originalCase->ssys, *collapsedSystem;
    DSMatrix * pInv_Ai, *pInv_Ai_all, * pInv_B, *pInv_B_all, *Knife = NULL, *pInv_Ad, *pInv_Ad_all;
    DSMatrixArray *array = NULL, *knife = NULL;
    DSUInteger i, numberOfEqualityConstraints = uCase->Xd_e->numberOfVariables;
    DSUInteger numberOfBlowing = uCase->Xd_b->numberOfVariables;
    DSUInteger *Xd_e_indices, *knife_row_indices, count=0;
    DSMatrix *gAi, *gAd, *gB, *gAi_all, *gAd_all, *gB_all;

    // merge algebraic constraints if necessary and get indices of Xd_e
    if( DSSSystemXd_a(o_ssys) != 0){

            collapsedSystem = DSSSystemByRemovingAlgebraicConstraints(o_ssys);

            if (numberOfEqualityConstraints != 0){
                    Xd_e_indices = DSVariablePoolIndicesOfSubPool(collapsedSystem->Xd, uCase->Xd_e);
                    // populate Equality first.
                    pInv_Ad_all = DSMatrixByMultiplyingMatrix(pInverse, DSSSystemAd(collapsedSystem));
                    pInv_Ai_all = DSMatrixByMultiplyingMatrix(pInverse, DSSSystemAi(collapsedSystem));
                    pInv_B_all  = DSMatrixByMultiplyingMatrix(pInverse, DSSSystemB(collapsedSystem));


                    // get a subset of pInv_Ai_all to construct pInv_Ai. Loop over Xd_e;
                    pInv_Ad = DSMatrixSubMatrixIncludingRows(pInv_Ad_all, numberOfEqualityConstraints, Xd_e_indices);
                    pInv_Ai = DSMatrixSubMatrixIncludingRows(pInv_Ai_all, numberOfEqualityConstraints, Xd_e_indices);
                    pInv_B  = DSMatrixSubMatrixIncludingRows(pInv_B_all, numberOfEqualityConstraints, Xd_e_indices);

                    array = DSMatrixArrayAlloc();
                    DSMatrixArrayAddMatrix(array, pInv_Ad);
                    DSMatrixArrayAddMatrix(array, pInv_Ai);
                    DSMatrixArrayAddMatrix(array, pInv_B);

                    // print matrices to debug.

//                    printf("The matrix pInv_Ad_all is: \n");
//                    DSMatrixPrint(pInv_Ad_all);
//
//                    printf("The matrix pInv_Ad is: \n");
//                    DSMatrixPrint(pInv_Ad);
//
//                    printf("The matrix pInv_Ai is: \n");
//                    DSMatrixPrint(pInv_Ai);
//
//                    printf("The matrix pInv_B is: \n");
//                    DSMatrixPrint(pInv_B);

                    //free variables
                    DSSecureFree(Xd_e_indices);
                    DSMatrixFree(pInv_Ai_all);
                    DSMatrixFree(pInv_B_all);
            }

            // populate Knife Edge
            // first, unpack gauss matrices
            gAd_all = DSMatrixArrayMatrix(gaussArray, 0);
            gAi_all = DSMatrixArrayMatrix(gaussArray, 1);
            gB_all  = DSMatrixArrayMatrix(gaussArray, 2);

            // then get row indices for which the matrix gAd has only zeros. and get only corresponding rows of gAdi
//            old method does not work: knife_row_indices = DSVariablePoolIndicesOfSubPool(collapsedSystem->Xd, uCase->Xd_b);
            knife_row_indices = DSSecureMalloc(numberOfBlowing* sizeof(DSUInteger));
            for (i=0; i<DSMatrixRows(gAd_all); i++){
                if (DSMatrixFirstNonZeroIndexAtRow(gAd_all,i) == 65000){
                        if (count < numberOfBlowing){
                                                knife_row_indices[count] = i;
                                                count ++;
                        }
                }
            }

            gAi = DSMatrixSubMatrixIncludingRows(gAi_all, numberOfBlowing, knife_row_indices);
            gB  = DSMatrixSubMatrixIncludingRows(gB_all, numberOfBlowing, knife_row_indices);

            // add matrices to the array.
            knife = DSMatrixArrayAlloc();
            DSMatrixArrayAddMatrix(knife, gAi);
            DSMatrixArrayAddMatrix(knife, gB);

            // free variables.
            DSSSystemFree(collapsedSystem);
            DSSecureFree(knife_row_indices);

    } else { // else, the system does not have algebraic constraints.
            if (numberOfEqualityConstraints != 0){

                    Xd_e_indices = DSVariablePoolIndicesOfSubPool(o_ssys->Xd, uCase->Xd_e);
                    // populate Equality first.
                    pInv_Ad_all = DSMatrixByMultiplyingMatrix(pInverse, DSSSystemAd(o_ssys));
                    pInv_Ai_all = DSMatrixByMultiplyingMatrix(pInverse, DSSSystemAi(o_ssys));
                    pInv_B_all  = DSMatrixByMultiplyingMatrix(pInverse, DSSSystemB(o_ssys));
                    // get a subset of pInv_Ai_all to construct pInv_Ai. Loop over Xd_e;
                    pInv_Ad = DSMatrixSubMatrixIncludingRows(pInv_Ad_all, numberOfEqualityConstraints, Xd_e_indices);
                    pInv_Ai = DSMatrixSubMatrixIncludingRows(pInv_Ai_all, numberOfEqualityConstraints, Xd_e_indices);
                    pInv_B  = DSMatrixSubMatrixIncludingRows(pInv_B_all, numberOfEqualityConstraints, Xd_e_indices);

                    array = DSMatrixArrayAlloc();
                    DSMatrixArrayAddMatrix(array, pInv_Ad);
                    DSMatrixArrayAddMatrix(array, pInv_Ai);
                    DSMatrixArrayAddMatrix(array, pInv_B);

                    //free variables
                    DSSecureFree(Xd_e_indices);
                    DSMatrixFree(pInv_Ad_all);
                    DSMatrixFree(pInv_Ai_all);
                    DSMatrixFree(pInv_B_all);
            }

            // populate Knife Edge
            // first, unpack gauss matrices
            gAd_all = DSMatrixArrayMatrix(gaussArray, 0);
            gAi_all = DSMatrixArrayMatrix(gaussArray, 1);
            gB_all  = DSMatrixArrayMatrix(gaussArray, 2);


            // then get row indices for which the matrix gAd has only zeros. and get only corresponding rows of gAdi
//            old method does not work: knife_row_indices = DSVariablePoolIndicesOfSubPool(collapsedSystem->Xd, uCase->Xd_b);
            knife_row_indices = DSSecureMalloc(numberOfBlowing* sizeof(DSUInteger));
            for (i=0; i<DSMatrixRows(gAd_all); i++){
                if (DSMatrixFirstNonZeroIndexAtRow(gAd_all,i) == 65000){
                        if (count < numberOfBlowing){
                                                knife_row_indices[count] = i;
                                                count ++;
                        }
                }
            }

            gAi = DSMatrixSubMatrixIncludingRows(gAi_all, numberOfBlowing, knife_row_indices);
            gB  = DSMatrixSubMatrixIncludingRows(gB_all, numberOfBlowing, knife_row_indices);

            // add matrices to the array.
            knife = DSMatrixArrayAlloc();
            DSMatrixArrayAddMatrix(knife, gAi);
            DSMatrixArrayAddMatrix(knife, gB);
            DSSecureFree(knife_row_indices);

    }

    uCase->Equality = array;
    uCase->Knife = knife;

}

// 3.-10 dsUnstableCaseLinearProblemAddKnifeEdgeConditions (8)
DSMatrixArray * dsUnstableCaseLinearProblemAddKnifeEdgeConditions(const DSMatrix *A_,
                                                            const DSMatrix *B_, const DSUnstableCase *uCase)
{

    // this function adds knife edge conditions. Note that the knife edge is a function of independent variables
    // and numerical coefficients. We will ignore numerical coefficients for now.

    // the merging will be in the form [[0], gAi, [0]] corresponding to matrix of zeros for dependent variables Xd
    // (including ) algebraic constraints, gauss Ai and matrix of zeros for slack.

    DSMatrix *zeros1, *zeros2, *A, *B, *temp1, *temp2;
    DSMatrixArray *array = NULL, *Knife = uCase->Knife;
    DSUInteger numberOfKnifes = uCase->Xd_b->numberOfVariables;
    DSUInteger numberOfXd = uCase->originalCase->Xd->numberOfVariables;

    if (numberOfKnifes == 0) {
                DSError(M_DS_NULL ": Case does not has any knife-edge constraints", A_DS_ERROR);
                goto bail;
    }

    zeros1 = DSMatrixCalloc(numberOfKnifes, numberOfXd);
    zeros2 = DSMatrixCalloc(numberOfKnifes, 1);

    temp1 = DSMatrixAppendMatrices(zeros1, DSMatrixArrayMatrix(Knife, 0), true);
    DSMatrixFree(zeros1);

    temp2 = DSMatrixAppendMatrices(temp1, zeros2, true);
    DSMatrixFree(temp1);

    A = DSMatrixAppendMatrices(A_, temp2, false);
    B = DSMatrixAppendMatrices(B_, zeros2, false);
    DSMatrixFree(zeros2);


    array = DSMatrixArrayAlloc();
    DSMatrixArrayAddMatrix(array, A);
    DSMatrixArrayAddMatrix(array, B);

bail:

    return array;
}

// 3.-9 dsUnstableCaseLinearProblemAddEqualityConstraints (9)
DSMatrixArray * dsUnstableCaseLinearProblemAddEqualityConstraints(const DSMatrix *A_, const DSMatrix *B_, const DSUnstableCase *uCase)
{

    DSUInteger numberOfXi, numberOfEqualities, numberOfAlgebraicConstraints, preExt;
    DSMatrix *Ae, *temp1, *temp2, *zeros1, *zeros2, *Be, *A, *B;
    DSMatrixArray *Equality = uCase->Equality, *array;

    numberOfAlgebraicConstraints = uCase->originalCase->Xd_a->numberOfVariables;
    numberOfEqualities = uCase->Xd_e->numberOfVariables;

    if (uCase == NULL) {
                DSError(M_DS_NULL ": uCase is NULL", A_DS_ERROR);
                goto bail;
    }

    if (numberOfEqualities == 0) {
                DSError(M_DS_NULL ": Case does not has any equality constraint", A_DS_ERROR);
                goto bail;
    }

    // debugging.
    if (uCase->originalCase->caseNumber == 3000){

        printf("Showing Equality Matrices for Case 3 \n");
        DSMatrixPrint(DSMatrixArrayMatrix(uCase->Equality, 0));
        DSMatrixPrint(DSMatrixArrayMatrix(uCase->Equality, 1));
        DSMatrixPrint(DSMatrixArrayMatrix(uCase->Equality, 2));

    }

    // Merge matrices [pInv * Ad] and [pInv * Ai], adding the correct amount of zeros to adjust for algebraic constraints
    // and slack variables. A = [ Equality[0]; [0]; Equality[1]; [0] ].

    if (numberOfAlgebraicConstraints != 0){
            zeros1 = DSMatrixCalloc(numberOfEqualities, numberOfAlgebraicConstraints);
            temp1 = DSMatrixAppendMatrices(DSMatrixArrayMatrix(Equality, 0), zeros1, true);
            DSMatrixFree(zeros1);
    }else{
            temp1 = DSMatrixArrayMatrix(Equality, 0);
    }
    zeros2 = DSMatrixCalloc(numberOfEqualities, 1);
    temp2 = DSMatrixAppendMatrices( DSMatrixArrayMatrix(Equality, 1), zeros2, true);

    Ae = DSMatrixAppendMatrices(temp1, temp2, true);
    A = DSMatrixAppendMatrices(A_, Ae, false);
    DSMatrixFree(Ae);

    Be = DSMatrixArrayMatrix(Equality, 2);
    B = DSMatrixAppendMatrices(B_, Be, false);

    array = DSMatrixArrayAlloc();
    DSMatrixArrayAddMatrix(array, A);
    DSMatrixArrayAddMatrix(array, B);

    if (numberOfAlgebraicConstraints != 0)
                DSMatrixFree(temp1);
    DSMatrixFree(zeros2);
    DSMatrixFree(temp2);

bail:
    return array;

}

// 3.-8 dsCaseLinearProblemForCaseValidity (10)
glp_prob * dsUnstableCaseLinearProblemForCaseValidity(const DSMatrix * U, const DSMatrix *zeta, DSUnstableCase *uCase,
                                                      const DSUInteger * bSignature )
{
        glp_prob *linearProblem = NULL;
        DSMatrix *slacks = NULL, * coefficients;
        DSUInteger numberOfXi, numberOfBoundaries;

        if (zeta == NULL) {
                DSError(M_DS_MAT_NULL, A_DS_ERROR);
                goto bail;
        }
        if (U == NULL)
                numberOfXi = 0;
        else
                numberOfXi = DSMatrixColumns(U);

        numberOfBoundaries = DSMatrixRows(zeta);
        if (numberOfXi > 0) {
                slacks = DSMatrixAlloc(numberOfBoundaries, 1);
                DSMatrixSetDoubleValueAll(slacks, 1.0);
                coefficients = DSMatrixAppendMatrices(U, slacks, true);
                DSMatrixMultiplyByScalar(coefficients, -1.0);
        } else {
                coefficients = DSMatrixAlloc(numberOfBoundaries, 1);
                DSMatrixSetDoubleValueAll(coefficients, -1.0);

        }

        linearProblem = dsUnstableCaseLinearProblemForMatrices(coefficients, zeta, uCase, bSignature);
//        printf("linear problem was built for case %u! \n", uCase->originalCase->caseNumber);
        glp_set_col_bnds(linearProblem, glp_get_num_cols(linearProblem), GLP_LO, -1.0, 0.0);
        glp_set_obj_coef(linearProblem, glp_get_num_cols(linearProblem), 1.0);

        DSMatrixFree(coefficients);
        if (slacks != NULL)
                DSMatrixFree(slacks);
bail:
        return linearProblem;
}

// 3.-7 dsUnstableCaseGetAdditionalConstraintMatrices(uCase, bSignature) (11)
void dsUnstableCaseGetAdditionalConstraintMatrices(DSUnstableCase *uCase, const DSUInteger *bSignature)
{

     if (uCase == NULL || uCase->originalCase == NULL) {
        DSError(M_DS_MAT_NULL, A_DS_ERROR);
     }

    DSMatrix *Cd_unstable, *Ci_unstable, *delta_unstable, *zeros, *temp;
    DSMatrix *Cd_unstable_all, *Ci_unstable_all, *delta_unstable_all;
    DSUInteger rows = uCase->Xd_b->numberOfVariables;
    DSUInteger i, n;
    DSSSystem *collapsedSystem, *o_ssys = uCase->originalCase->ssys;
    DSUInteger *Xd_b_indices;

    // merge algebraic constraints if necessary and get indices of Xd_b
    if( DSSSystemXd_a(o_ssys) != 0){
            collapsedSystem = DSSSystemByRemovingAlgebraicConstraints(o_ssys);
            Xd_b_indices = DSVariablePoolIndicesOfSubPool(collapsedSystem->Xd, uCase->Xd_b);

            temp = DSSSystemAd(collapsedSystem);
            Ci_unstable_all = DSSSystemAi(collapsedSystem);
            delta_unstable_all = DSSSystemB(collapsedSystem);

            // let's add a block of zeros with dimentions rows = DSMatrixRows(Cd_unstable_all) and columns = number of Xd_a
            zeros = DSMatrixCalloc(DSMatrixRows(temp), DSSSystemXd_a(o_ssys)->numberOfVariables);
            Cd_unstable_all = DSMatrixAppendMatrices(temp, zeros, true);

            DSMatrixFree(zeros);
            DSMatrixFree(temp);
            DSSSystemFree(collapsedSystem);
    }
    else {
            Xd_b_indices = DSVariablePoolIndicesOfSubPool(o_ssys->Xd, uCase->Xd_b);

            Cd_unstable_all = DSSSystemAd(o_ssys);
            Ci_unstable_all = DSSSystemAi(o_ssys);
            delta_unstable_all = DSSSystemB(o_ssys);
    }

    // set signs of rows of all matrices depending blowing up or down behavior. blowing up = 1
    // blowing down = 2. see  page 374 & 425
    for (i=0 ; i<rows; i++){
        if (bSignature[i] == 2){
            for (n=0; n<DSMatrixColumns(Cd_unstable_all); n++){
                DSMatrixSetDoubleValue(Cd_unstable_all, Xd_b_indices[i], n,
                                       DSMatrixDoubleValue(Cd_unstable_all, Xd_b_indices[i], n)*-1.0);
            }

            for (n=0; n<DSMatrixColumns(Ci_unstable_all); n++){
                DSMatrixSetDoubleValue(Ci_unstable_all, Xd_b_indices[i], n,
                                       DSMatrixDoubleValue(Ci_unstable_all, Xd_b_indices[i], n)*-1.0);
            }

        } else {

            for (n=0; n<DSMatrixColumns(delta_unstable_all); n++){
                DSMatrixSetDoubleValue(delta_unstable_all, Xd_b_indices[i], n,
                                       DSMatrixDoubleValue(delta_unstable_all, Xd_b_indices[i], n)*-1.0);
            }

        }
    }

    uCase->Cd_unstable = DSMatrixSubMatrixIncludingRows(Cd_unstable_all, rows, Xd_b_indices);
    uCase->Ci_unstable = DSMatrixSubMatrixIncludingRows(Ci_unstable_all, rows, Xd_b_indices);
    uCase->delta_unstable = DSMatrixSubMatrixIncludingRows(delta_unstable_all, rows, Xd_b_indices);

    DSMatrixFree(Cd_unstable_all);
    DSMatrixFree(Ci_unstable_all);
    DSMatrixFree(delta_unstable_all);
    DSSecureFree(Xd_b_indices);
}

// 3.-6 DSUnstableCaseConditionsAreValid(uCase, bSignature) constructs the case and returns its validity. (12)
const bool DSUnstableCaseConditionsAreValid(DSUnstableCase *uCase, const DSUInteger *bSignature)
{

        bool isValid = false;
        glp_prob *linearProblem = NULL;
        DSMatrix * U1, *U2, *U, * Zeta1, *Zeta;
        DSCase *aCase = uCase->originalCase;

        if (uCase == NULL || aCase == NULL ) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }

        // populate matrices Cd_unstable, Ci_unstable and delta_unstable
        dsUnstableCaseGetAdditionalConstraintMatrices(uCase, bSignature);

        // generate Matrices U and Z.
//        printf("Checking case %u . \n", uCase->originalCase->caseNumber);
//        printf("Showing Matrix Cd: \n");
//        DSMatrixPrint(DSCaseCd(aCase));
//
//        printf("Showing Matrix Ci: \n");
//        DSMatrixPrint(DSCaseCi(aCase));
//
//        printf("Showing Matrix Cd_unstable: \n");
//        DSMatrixPrint(uCase->Cd_unstable);
//
//        printf("Showing Matrix Ci_unstable: \n");
//        DSMatrixPrint(uCase->Ci_unstable);
//
//        printf("Showing Vector Delta: \n");
//        DSMatrixPrint(DSCaseDelta(aCase));
//
//        printf("Showing Vector Delta_unstable: \n");
//        DSMatrixPrint(uCase->delta_unstable);

        U1 = DSMatrixAppendMatrices(DSCaseCd(aCase), DSCaseCi(aCase), true);
        U2 = DSMatrixAppendMatrices(uCase->Cd_unstable, uCase->Ci_unstable, true);
        U = DSMatrixAppendMatrices(U1, U2, false);
//        DSMatrixMultiplyByScalar(U, -1.0);
        DSMatrixFree(U1);
        DSMatrixFree(U2);

        Zeta1 = DSCaseDelta(aCase);
        Zeta = DSMatrixAppendMatrices(Zeta1, uCase->delta_unstable, false);

//        printf("Showing Matrix U: \n");
//        DSMatrixPrint(U);
//
//        printf("Showing Matrix Zeta: \n");
//        DSMatrixPrint(Zeta);

        linearProblem = dsUnstableCaseLinearProblemForCaseValidity(U, Zeta, uCase, bSignature);
        DSMatrixFree(U);
        DSMatrixFree(Zeta);

        if (linearProblem != NULL) {
                glp_simplex(linearProblem, NULL);

                if (glp_get_obj_val(linearProblem) <= -1E-14 && glp_get_prim_stat(linearProblem) == GLP_FEAS) {
                        isValid = true;
                        DSUnstableCaseDetermineBlowingBehavior(linearProblem, uCase, bSignature);
//                        printf("Showing Parameter Values for Valid Case: \n");
//                        DSUnstableCasePrintLinearProblem(linearProblem, uCase, bSignature);
                }
                glp_delete_prob(linearProblem);
        }
//        printf("The problem is valid?: %s \n", isValid?"True":"False");
bail:
        return isValid;

}

// 3.-5 DSUnstableCaseExpandConditionMatrices (13)
void DSUnstableCaseExpandConditionMatrices(DSUnstableCase *uCase)
{

    DSUInteger *bSignature;
    DSUInteger numberOfCases, i, w;
    DSVariablePool *Xd_b = uCase->Xd_b;
    bool isValid = false;
    char * string = NULL;


    numberOfCases = pow(2, Xd_b->numberOfVariables);
    uCase->ValidCases = DSDictionaryAlloc();

    for (i=1; i <= numberOfCases; i++){

        // This is an array of length Xb->numberOfVariables containing values of 1 for blow-down and 2 for blow up.
        bSignature = DSBlowUpSignatureForCaseNumber(i, Xd_b);

        // construct case and return validity
        isValid = DSUnstableCaseConditionsAreValid(uCase, bSignature);
        DSSecureFree(bSignature);
        printf("Analyzing case %u. Subcase %u of %u is valid?: %s \n", uCase->originalCase->caseNumber, i ,
               numberOfCases, isValid?"True":"False");

        // If it has a solution, store the case number in the field uCase->validCases
        if ( isValid == true ) {
            string = DSSecureCalloc(sizeof(char), 100);
            sprintf(string, "%d", i);
            DSDictionaryAddValueWithName(uCase->ValidCases, string, (void*)1);

            DSUnstableCaseCreateBoundaryMatrices(uCase);

            // assign boundary matrices to original case.
            uCase->originalCase->U      = uCase->U;
            uCase->originalCase->zeta   = uCase->Zeta;

            DSGetVertices(uCase);

            DSSecureFree(string);
        }
    }



    // for de-bugin.
    // printf("The number of cases is %u \n", numberOfCases);

    // delete variables
}

// 3.-4 DSUnstableCaseGetXd (14)
void DSUnstableCaseGetXd_b(DSUnstableCase *uCase, DSMatrixArray *gaussArray)
{

    DSVariablePool *Xd_b, *knifeEdge, *Xd_e;
    DSUInteger i, index;
    const DSVariable *blow_variable;
    const DSVariable *equal_variable;
    const DSVariablePool *Xd;

    Xd_b = DSVariablePoolAlloc();
    Xd_e = DSVariablePoolAlloc();
    Xd = uCase->originalCase->Xd;

    // get the index of gAd where all elements of a column are zero. Use for that
    // DSMatrixFirstNonZeroIndexAtRow(const DSMatrix *matrix, const DSUInteger row). If Row contains only zeros, then
    // function returns 65000. Note that the matrix need to be transposed.

//    for debugging.
//    printf("showing matrix gAd \n");
//    DSMatrixPrint(DSMatrixArrayMatrix(gaussArray, 0));

    for (i=0; i<DSMatrixRows(DSMatrixArrayMatrix(gaussArray, 0)); i++){

        index = DSMatrixFirstNonZeroIndexAtRow( DSMatrixTranspose(DSMatrixArrayMatrix(gaussArray, 0)) , i );
        // if row contains zeros, add corresponding variable of Xd to to Xd_b, else it has a solution and belong to Xd_e
        if (index == 65000){
            blow_variable = DSVariablePoolAllVariables(Xd)[i];
            DSVariablePoolAddVariableWithName(Xd_b, DSVariableName(blow_variable));
        }else{
            equal_variable = DSVariablePoolAllVariables(Xd)[i];
            DSVariablePoolAddVariableWithName(Xd_e, DSVariableName(equal_variable));
        }
    }
    // assign Xd_b and Xd_e
    uCase->Xd_b = Xd_b;
    uCase->Xd_e = Xd_e;

// for debugging
//    printf("showing pool Xb: \n");
//    printPool(Xd_b);
}

// 3.-3 DSUnstableCaseIdentifyBlowingDependentVariables (15)
void DSUnstableCaseIdentifyBlowingDependentVariables(DSUnstableCase *uCase)
{

    DSVariablePool *Xb;
    DSMatrix *Ad, *Ai, *B, *gAd, *gAi, *gB;
    const DSMatrix *pInverse;
    DSSSystem *collapsedSystem, *o_ssys = uCase->originalCase->ssys;
    DSMatrixArray *gaussArray;

    // Calculation of the pseudo inverse of Ad.
    // check if system needs to be collapsed and then calculate pseudoinverse.
    if( DSSSystemXd_a(o_ssys) != 0)
        {
            collapsedSystem = DSSSystemByRemovingAlgebraicConstraints(o_ssys);

            // prepare matrices for gaussian elimination
            Ad = DSSSystemAd(collapsedSystem);
            Ai = DSSSystemAi(collapsedSystem);
            DSMatrixMultiplyByScalar(Ai, -1.0);
            B = DSSSystemB(collapsedSystem);
        }
        else {
            // prepare matrices for gaussian elimination
            Ad = DSSSystemAd(o_ssys);
            Ai = DSSSystemAi(o_ssys);
            DSMatrixMultiplyByScalar(Ai, -1.0);
            B = DSSSystemB(o_ssys);
        }

    // Gaussian Elimination
    gaussArray = DSMatrixGaussElimination(Ad, Ai, B);

    //Unpack Matrices
//    gAd = DSMatrixArrayMatrix(gaussArray, 0);
//    gAi = DSMatrixArrayMatrix(gaussArray, 1);
//    gB  = DSMatrixArrayMatrix(gaussArray, 2);

    // Identify Xb and assign.
    DSUnstableCaseGetXd_b(uCase, gaussArray);

    // Calculate pseudo-inverse of Ad by only considering a subset of the matrix Ad.
    pInverse = DSUnstableCaseGetSubSetPseudoInverse(uCase->originalCase->ssys->Xd_t, uCase->Xd_e,
                                                     uCase->Xd_b, Ad);
    uCase->pInverse = pInverse;

    // Get matrices Equality and KnifeEdge
    DSUnstableCaseGetEqualityAndKnife(pInverse, gaussArray, uCase);

    // delete variables.
    // consider deleting gaussArray.
    if( DSSSystemXd_a(o_ssys) != 0)
        DSSSystemFree(collapsedSystem);
    DSMatrixFree(Ad);
    DSMatrixFree(Ai);
    DSMatrixFree(B);
    DSMatrixArrayFree(gaussArray);
}

// 3.-2 (16)
DSUnstableCase * DSCalculateUnstableCase(DSDesignSpace *ds, DSCase *aCase)
{
        DSUnstableCase * unstableCase = NULL;
        if (ds == NULL) {
                DSError(M_DS_DESIGN_SPACE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }

        if (DSCaseIsValid(aCase, true) == true) {
                goto bail;
        }
        if (DSCaseConditionsAreValid(aCase) == false) {
                goto bail;
        }
        // add condition to rule out cyclical cases and/or cases with solution but not valid constraints.
        if (DSCaseIsCyclical(aCase) == true){
            goto bail;
        }

        unstableCase = DSSecureCalloc(sizeof(DSUnstableCase), 1);
        unstableCase->originalCase = aCase;

        // perform main calculations

        // this function should populate Xb, pInverse and delta_unstable
        DSUnstableCaseIdentifyBlowingDependentVariables(unstableCase);

        // this function should set fields KnifeEdge, Cd_unstable, Ci_unstable and delta_unstable.
        DSUnstableCaseExpandConditionMatrices(unstableCase);



        // fill fields of DSUnstableCase
//        unstableCase->Xb = ;
//        unstableCase->pInverse = ;

//        unstableCase->KnifeEdge = ;
//        unstableCase->Cd_unstable = ;
//        unstableCase->Ci_unstable = ;
//        unstableCase->delta_unstable = ;

bail:
        return unstableCase;
}

// 3.-1 (17)
void DSDesignSpaceCalculateUnstableCase(DSDesignSpace *ds, DSCase *aCase)
{

        DSUInteger caseNumber;
        char * string = NULL;
        DSUnstableCase * unstableCase = NULL;

        if (ds == NULL) {
                DSError(M_DS_DESIGN_SPACE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (DSCaseConditionsAreValid(aCase) == false) {
                goto bail;
        }
        string = DSSecureCalloc(sizeof(char), 100);
        caseNumber = DSCaseNumber(aCase);
        sprintf(string, "%d", caseNumber);
        if (DSDictionaryValueForName(DSDSUnstable(ds), string) == NULL) {
                unstableCase = DSCalculateUnstableCase(ds, aCase);
                if (unstableCase != NULL && unstableCase->ValidCases->count !=0){
//                        DSUnstableCaseCreateBoundaryMatrices(unstableCase);
                        DSDictionaryAddValueWithName(DSDSUnstable(ds), string, unstableCase);
//                        DSGetVertices(unstableCase);
//                        DSDictionaryAddValueWithName(DSDSUnstable(ds), string, (void*)1);
                }

        }
        if (string != NULL)
                DSSecureFree(string);
bail:
        return;

}

// 3. create main function DSDesignSpaceCalculateUnstableCases (18)
void DSDesignSpaceCalculateUnstableCases(DSDesignSpace *ds)
{

        DSUInteger i, caseNumber, numberOfCases, * termSignature;
        DSCase * aCase = NULL;
        DSDSUnstable(ds) = DSDictionaryAlloc();

        if (ds == NULL) {
                DSError(M_DS_DESIGN_SPACE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (DSDSGMA(ds) == NULL) {
                DSError(M_DS_GMA_NULL, A_DS_ERROR);
                goto bail;
        }
        if (DSGMASystemSignature(DSDSGMA(ds)) == NULL) {
                DSError(M_DS_WRONG ": GMA signature is NULL", A_DS_ERROR);
                goto bail;
        }
        numberOfCases = DSDesignSpaceNumberOfCases(ds);
        if (numberOfCases == 0) {
                DSError(M_DS_WRONG ": Number of cases to process must be more than 0", A_DS_ERROR);
                goto bail;
        }
        for (i = 0; i < numberOfCases; i++) {

                caseNumber = i+1;
//                if (caseNumber==6){

                if (caseNumber == 0)
                        continue;
                if (caseNumber > DSDesignSpaceNumberOfCases(ds)) {
                        DSError(M_DS_WRONG ": Case number out of bounds", A_DS_ERROR);
                        continue;
                }
                termSignature = DSCaseSignatureForCaseNumber(caseNumber, ds->gma);
                if (termSignature != NULL) {
                        aCase = DSCaseWithTermsFromDesignSpace(ds, termSignature, DSDesignSpaceCasePrefix(ds));
                        if (aCase != NULL) {
                                DSDesignSpaceCalculateUnstableCase(ds, aCase);
                                DSCaseFree(aCase);
                        }
                        DSSecureFree(termSignature);
                }

//                } // if(caseNumber==3)

        }

bail:
        return;


}

// 4. Test if case is cyclical// if matrices are consistent. (19)
const bool DSCaseIsCyclical(const DSCase *aCase)
{

        DSMatrix *A, *B, *Ai, *Ad;
        DSSSystem * collapsedSystem;
        DSUInteger i, j, rank_Ad, rank_Ai;
        bool isCyclical = false;
        DSMatrixArray *array;


        // get Matrices A and b corresponding to the equation A x = b.
        // These matrices will be subjected to gaussian elimination to determine if the system is cyclical
        // - in case that the case didn't pass the test DSCaseIsValid  -
        // note that A corresponds to the matrix Ad and b = (log(beta)-log(alpha)) - Ai.
        // if the case has algebraic constraints Xd_a != 0, they need to be merged first.

        if( DSSSystemXd_a(aCase->ssys) != 0)
        {
            collapsedSystem = DSSSystemByRemovingAlgebraicConstraints(aCase->ssys);
            Ad = DSSSystemAd(collapsedSystem);
            Ai = DSSSystemAi(collapsedSystem);
            DSMatrixMultiplyByScalar(Ai, -1.0);
            B = DSSSystemB(collapsedSystem);
            DSSSystemFree(collapsedSystem);
        }
        else {
            Ad = DSSSystemAd(aCase->ssys);
            Ai = DSSSystemAi(aCase->ssys);
            DSMatrixMultiplyByScalar(Ai, -1.0);
            B = DSSSystemB(aCase->ssys);
        }

//        printf("The matrix Ad is: \n");
//        DSMatrixPrint(Ad);
//
//        printf("The matrix b is: \n");
//        DSMatrixPrint(b);

        // call gaussian elimination on A (Ad) and b, eventually also on B!:
        array = DSMatrixGaussElimination(Ad, Ai, B);
        // calculate ranks for the diagonalized matrix Ad and b.
        rank_Ad = DSMatrixRank(DSMatrixArrayMatrix(array,0));
        if (DSMatrixRows(DSMatrixArrayMatrix(array,1)) < DSMatrixColumns(DSMatrixArrayMatrix(array,1))){
            rank_Ai = DSMatrixRank(DSMatrixTranspose(DSMatrixArrayMatrix(array,1)));
        }
        else{
            rank_Ai = DSMatrixRank(DSMatrixArrayMatrix(array,1));
        }

        // set value for isCyclical
        if(rank_Ad == rank_Ai){
            isCyclical = true;
        }
        // print matrices for de-bugging.
//        printf("The diagonalized left-hand side matrix A is: \n");
//        DSMatrixPrint(DSMatrixArrayMatrix(array,0));
//        printf("The rank of the left-hand side matrix A is %u: \n", DSMatrixRank(DSMatrixArrayMatrix(array,0)));
//
//        printf("The diagonialized right-hand side matrix b is: \n");
//        DSMatrixPrint(DSMatrixArrayMatrix(array,1));
//        printf("The rank of the left-hand side matrix b is %u: \n", DSMatrixRank(DSMatrixTranspose(DSMatrixArrayMatrix(array,1))));

        // delete variables
        DSMatrixFree(B);
        DSMatrixFree(Ai);
        DSMatrixFree(Ad);
        DSMatrixArrayFree(array);

        return isCyclical;
}

// 5. Gaussian elimination. (20)
DSMatrixArray * DSMatrixGaussElimination(const DSMatrix *Ad, const DSMatrix *Ai, const DSMatrix *B)
{
    DSUInteger h=0, k=0, i, j, m, n, jj, i_max, *swappingVector=NULL;
    double f, max_value;
    DSMatrix *column_k, *rhs, *lhs, *b, *rhs_b, *lhs_b, *b_b;
    m = DSMatrixRows(Ad);
    n = DSMatrixColumns(Ad);
    DSMatrixArray *array = NULL;

    // make copies of matrices:
    lhs = DSMatrixCopy(Ad);
    rhs = DSMatrixCopy(Ai);
    b = DSMatrixCopy(B);

    while( h < m && k <n ) {
        // get k-column
        column_k = DSMatrixSubMatrixIncludingColumnList(lhs, 1, k);
        // absolute values
        DSMatrixApplyFunction(column_k, fabs);
        // get the max value of that column
        max_value = maximumValue(column_k, false);
        // get the index of that value
        i_max = maximumValueIndex(column_k, max_value);
        if(max_value == 0.0){
            k = k +1;
            continue;
        } else{
            // swap rows of all three matrices if i_max and h indices differ:
            if (i_max != h){
                swappingVector = DSMatrixSwapRows(lhs, i_max, h);
                DSMatrixSwitchRows(rhs, i_max, h);
                DSMatrixSwitchRows(b, i_max, h);
            }
            for (i=h+1; i<m; i++ ){
                f = DSMatrixDoubleValue(lhs, i, k) / DSMatrixDoubleValue(lhs, h, k);
                // fill with zeros
                DSMatrixSetDoubleValue(lhs, i, k, 0);
                // for all elements in current row
                for(j=k+1; j<n; j++){
                    DSMatrixSetDoubleValue(lhs, i, j, (DSMatrixDoubleValue(lhs, i, j) - DSMatrixDoubleValue(lhs, h, j)*f));
                }
                // modify Ai accordingly
                for(jj=0; jj<DSMatrixColumns(rhs); jj++){
                    DSMatrixSetDoubleValue(rhs, i, jj, (DSMatrixDoubleValue(rhs, i, jj) - DSMatrixDoubleValue(rhs, h, jj)*f));
                }
                // modify B accordingly
                for(jj=0; jj<DSMatrixColumns(b); jj++){
                    DSMatrixSetDoubleValue(b, i, jj, (DSMatrixDoubleValue(b, i, jj) - DSMatrixDoubleValue(b, h, jj)*f));
                }

            }

        h++;
        k++;

        }
    }

    // swap matrices back and delete swappingVector and swapped matrices.
    if (swappingVector != NULL){

        lhs_b = DSMatrixSwapRowsBack(lhs, swappingVector);
        rhs_b = DSMatrixSwapRowsBack(rhs, swappingVector);
        b_b = DSMatrixSwapRowsBack(b, swappingVector);

        DSMatrixFree(lhs);
        DSMatrixFree(rhs);
        DSMatrixFree(b);
        DSSecureFree(swappingVector);

        array = DSMatrixArrayAlloc();
        DSMatrixArrayAddMatrix(array, lhs_b);
        DSMatrixArrayAddMatrix(array, rhs_b);
        DSMatrixArrayAddMatrix(array, b_b);

    }else{
        array = DSMatrixArrayAlloc();
        DSMatrixArrayAddMatrix(array, lhs);
        DSMatrixArrayAddMatrix(array, rhs);
        DSMatrixArrayAddMatrix(array, b);
    }
    return array;
}

// 7. Swaping function DSMatrixSwapRows (21)
DSUInteger * DSMatrixSwapRows(DSMatrix *matrix, DSUInteger rowA, DSUInteger rowB)
{
    // This functions swaps two given rows and writes the swapped matrix in *matrix. The function returns
    // a swappingVector as output that can be used to back swap the matrix.
    DSUInteger *swappingVector = DSSecureMalloc(sizeof(DSUInteger)*DSMatrixRows(matrix));
    DSUInteger i;
    DSMatrixSwitchRows(matrix, rowA, rowB);

    for(i=0; i<DSMatrixRows(matrix); i++){
        swappingVector[i] =i;
    }

    swappingVector[rowA] = rowB;
    swappingVector[rowB] = rowA;

    return swappingVector;
}

// 8. Swap back (22)
DSMatrix * DSMatrixSwapRowsBack(DSMatrix *matrix, DSUInteger *swappingVector)
{

    DSMatrix * back_matrix;
    back_matrix = DSMatrixSubMatrixIncludingRows(matrix, DSMatrixRows(matrix), swappingVector);

    return back_matrix;
}

// 9. create boundary matrices U and Z for unstable cases. (23)
void DSUnstableCaseCreateBoundaryMatrices(DSUnstableCase *uCase)
{

    // The idea of this function is to generate matrices U and zeta for the unstable case. Refer to
    // function dsCaseCreateBoundaryMatrices for analogous version for a stable case.

    if (uCase == NULL) {
                DSError(M_DS_CASE_NULL ": Pointer to unstable case is NUll", A_DS_ERROR);
                goto bail;
    }

    if (uCase->pInverse == NULL) {
                DSError(M_DS_CASE_NULL ": Pointer to pInverse is NUll", A_DS_ERROR);
                goto bail;
    }


    DSUInteger numberOfXi = 0, *Xd_t_indices, *Xd_b_indices;
    DSUInteger numberOfBlowing = uCase->Xd_b->numberOfVariables, i, j;
    DSMatrix * W = NULL, *Ai, *Zeta, *U, *Cd_xd_t, *B_xd_t, *temp, *delta;
    DSCase *aCase = uCase->originalCase;

    DSMatrix *Cd, *Ci;
    DSMatrix *B = DSSSystemB(DSCaseSSys(aCase));

    DSSSystem *o_ssys, *collapsedSSystem;
    const DSVariablePool *Xd_t = uCase->originalCase->ssys->Xd_t, *Xd = uCase->originalCase->ssys->Xd;
    const DSUInteger numberOfColumns = Xd_t->numberOfVariables;

    o_ssys = uCase->originalCase->ssys;
    numberOfXi = DSVariablePoolNumberOfVariables(DSCaseXi(aCase));

    // we first construct Cd, Ci and delta by merging different matrices.
    temp = DSMatrixAppendMatrices(DSCaseCd(aCase), uCase->Cd_unstable, false);
    Cd = DSMatrixAppendMatrices(temp, uCase->Cd_knife, false);
    DSMatrixFree(temp);

    temp = DSMatrixAppendMatrices(DSCaseCi(aCase), uCase->Ci_unstable, false);
    Ci = DSMatrixAppendMatrices(temp, uCase->Ci_knife, false);
    DSMatrixFree(temp);

    temp = DSMatrixAppendMatrices(DSCaseDelta(aCase), uCase->delta_unstable, false);
    delta = DSMatrixAppendMatrices(temp, uCase->delta_knife, false);
    DSMatrixFree(temp);

    // let's create the same matrices but without cd_knife.
//    Cd = DSMatrixAppendMatrices(DSCaseCd(aCase), uCase->Cd_unstable, false);
//    Ci = DSMatrixAppendMatrices(DSCaseCi(aCase), uCase->Ci_unstable, false);
//    delta = DSMatrixAppendMatrices(DSCaseDelta(aCase), uCase->delta_unstable, false);

    // let's create the same matrices but without original values.
//    Cd = DSMatrixAppendMatrices(uCase->Cd_unstable, uCase->Cd_knife, false);
//    Ci = DSMatrixAppendMatrices(uCase->Ci_unstable, uCase->Ci_knife, false);
//    delta = DSMatrixAppendMatrices(uCase->delta_unstable, uCase->delta_knife, false);



//    printf("Showing matrices Cd unstable  \n");
//    DSMatrixPrint(uCase->Cd_unstable);
//    printf("Showing matrices Ci unstable \n");
//    DSMatrixPrint(uCase->Ci_unstable);
//    printf("Showing matrices delta unstable \n");
//    DSMatrixPrint(uCase->delta_unstable);

    // Matrices Cd and B of the original system need to be cut to get rid of Xd_a,
    // this is because my pInverse only contains Xd_t.
    if( DSSSystemXd_a(o_ssys) != 0){

            collapsedSSystem = DSSSystemByRemovingAlgebraicConstraints(o_ssys);
            // Generate matrix Cd_xd_t as a subset of matrix Cd.
            Xd_t_indices = DSVariablePoolIndicesOfSubPool(Xd, Xd_t);
            Cd_xd_t = DSMatrixSubMatrixIncludingColumns(Cd, numberOfColumns, Xd_t_indices);

            // Generate vector B_xd_t as a subset of vector B
            B_xd_t = DSMatrixSubMatrixIncludingRows(B, numberOfColumns, Xd_t_indices);

            // now calculate W & Zeta
            W = DSMatrixByMultiplyingMatrix(Cd_xd_t, uCase->pInverse);
            Zeta = DSMatrixByMultiplyingMatrix(W, B_xd_t);
            DSMatrixAddByMatrix(Zeta, delta);
            if (numberOfXi != 0) {
                        Ai = DSSSystemAi(collapsedSSystem);
                        U = DSMatrixByMultiplyingMatrix(W, Ai);
                        if (DSCaseCi(aCase) != NULL)
                                DSMatrixSubstractByMatrix(U, Ci);
                        DSMatrixMultiplyByScalar(U, -1.0);
                        DSMatrixFree(Ai);
                }

            // now rows of Zeta are adjusted to consider for numerical values of blowing dependent variables
            // in the condition matrices.
            Xd_b_indices = DSVariablePoolIndicesOfSubPool(Xd_t, uCase->Xd_b);

            // loop over columns of Cd_xd_t that contain blowing variables
            for(i=0; i<numberOfBlowing; i++){
                // loop over rows of those columns.
                for (j=0; j< DSMatrixRows(Cd_xd_t); j++){
                        if( DSMatrixDoubleValue(Cd_xd_t, j, Xd_b_indices[i]) != 0.0 ){
                                DSMatrixSetDoubleValue(Zeta, j, 0,
                                     DSMatrixDoubleValue(Zeta, j, 0) + log10(uCase->Xd_b->variables[i]->value));
                        }
                }
            }

            DSSSystemFree(collapsedSSystem);
            DSMatrixFree(Cd_xd_t);
            DSMatrixFree(B_xd_t);
            DSMatrixFree(W);
            DSMatrixFree(B);
            DSSecureFree(Xd_t_indices);
            DSSecureFree(Xd_b_indices);

    } else {
            // do normal calculations.
            B = DSSSystemB(DSCaseSSys(aCase));
            numberOfXi = DSVariablePoolNumberOfVariables(DSCaseXi(aCase));
            W = DSMatrixByMultiplyingMatrix(Cd, uCase->pInverse);
            Zeta = DSMatrixByMultiplyingMatrix(W, B);
            DSMatrixAddByMatrix(Zeta, delta);
                if (numberOfXi != 0) {
                        Ai = DSSSystemAi(DSCaseSSys(aCase));
                        U = DSMatrixByMultiplyingMatrix(W, Ai);
                        if (DSCaseCi(aCase) != NULL)
                                DSMatrixSubstractByMatrix(U, Ci);
                        DSMatrixMultiplyByScalar(U, -1.0);
                        DSMatrixFree(Ai);
                }

            // now rows of Zeta are adjusted to consider for numerical values of blowing dependent variables
            // in the condition matrices.
            Xd_b_indices = DSVariablePoolIndicesOfSubPool(Xd_t, uCase->Xd_b);

            // loop over columns
            for(i=0; i<numberOfBlowing; i++){
                // loop over rows
                for (j=0; j< DSMatrixRows(Cd); j++){
                        if( DSMatrixDoubleValue(Cd, j, Xd_b_indices[i]) != 0.0 ){
                                DSMatrixSetDoubleValue(Zeta, j, 0,
                                     DSMatrixDoubleValue(Zeta, j, 0) + log10(uCase->Xd_b->variables[i]->value));
                        }
                }
            }

            DSMatrixFree(W);
            DSMatrixFree(B);
            DSSecureFree(Xd_b_indices);

    }

    uCase->U = U;
    uCase->Zeta = Zeta;


//    // for debugging:
//
//      printf("The pool of independent variables is: \n");
//      printPool(aCase->Xi);
//
//      printf("Showing Matrix Cd_xd_t: \n ");
//      DSMatrixPrint(Cd_xd_t);
//
//      printf("Showing matrix U: \n");
//      DSMatrixPrint(U);
////////
//      printf("Showing Matrix Zeta: \n");
//      DSMatrixPrint(Zeta);
//
//      printf("Showing the variable pool knifeEdge \n");
//      printPool(uCase->knifeEdge);

//     printf("The matrix Ci_knife is: \n");
//     DSMatrixPrint(uCase->Ci_knife);


bail:
    return;

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
    if (1==1){
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
            printf("The vertices for case %u are: \n ", uCase->originalCase->caseNumber);
            vertices = DSCaseVerticesFor2DSlice(uCase->originalCase, lowerBounds, upperBounds, xVariable, yVariable);
            //    vertices = DSUnstableCaseVerticesFor2DSlice(uCase, lowerBounds, upperBounds, xVariable, yVariable);
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
