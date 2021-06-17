//
//  DSUnstableCase.c
//  designspace
//
//  Created by Miguel Angel Valderrama Gomez on 12/4/18.
//

#include <stdio.h>
#include <string.h>
#include <glpk.h>
#include "DSMemoryManager.h"
#include "DSCase.h"
#include "DSUnstableCase.h"
#include "DSVariable.h"
#include "DSGMASystem.h"
#include "DSSSystem.h"
#include "DSDesignSpace.h"
#include "DSExpression.h"
#include "DSMatrix.h"
#include "DSMatrixArray.h"
#include "DSCyclicalCase.h"
#include "DSGMASystemParsingAux.h"


extern void DSuCaseFree(DSUnstableCase *uCase)
{
    if (uCase == NULL) {
        DSError(M_DS_CASE_NULL, A_DS_ERROR);
        goto bail;
    }
    
    if(uCase->originalCase != NULL)
        DSCaseFree(uCase->originalCase);
    
    if (uCase->Xd_b != NULL){
        DSVariablePoolSetReadWriteAdd((DSVariablePool *)uCase->Xd_b);
        DSVariablePoolFree(uCase->Xd_b);
    }
    
    if (uCase->Xd_e != NULL){
        DSVariablePoolSetReadWriteAdd((DSVariablePool *)uCase->Xd_e);
        DSVariablePoolFree(uCase->Xd_e);
    }
    
    if (uCase->knifeEdge != NULL){
        DSVariablePoolSetReadWriteAdd((DSVariablePool *)uCase->knifeEdge);
        DSVariablePoolFree(uCase->knifeEdge);
    }
    
    if (uCase->pInverse != NULL)
        DSMatrixFree((DSMatrix *)uCase->pInverse);
    if (uCase->Equality != NULL)
        DSMatrixArrayFree(uCase->Equality);
    if (uCase->Knife != NULL)
        DSMatrixArrayFree(uCase->Knife);
    
    if (uCase->Cd_unstable != NULL)
        DSMatrixFree(uCase->Cd_unstable);
    if (uCase->Ci_unstable != NULL)
        DSMatrixFree(uCase->Ci_unstable);
    if (uCase->delta_unstable != NULL)
        DSMatrixFree(uCase->delta_unstable);
    
    if (uCase->Cd_knife != NULL)
        DSMatrixFree(uCase->Cd_knife);
    if (uCase->Ci_knife != NULL)
        DSMatrixFree(uCase->Ci_knife);
    if (uCase->delta_knife != NULL)
        DSMatrixFree(uCase->delta_knife);
    if (uCase->ValidCases != NULL)
        DSDictionaryFreeWithFunction(uCase->ValidCases, DSCaseFree);
    
    if (uCase->U != NULL)
        DSMatrixFree(uCase->U);
    if (uCase->Zeta != NULL)
        DSMatrixFree(uCase->Zeta);
    
    DSSecureFree(uCase);
bail:
    return;
    
    
}

extern bool DSuCaseIsValid(DSCase *aCase, bool strict)
{
    bool isValid = false;
    DSUnstableCase * unstableCase = NULL;
    DSSSystem *reduced = NULL;

    if (aCase == NULL) {
        DSError(M_DS_CASE_NULL, A_DS_ERROR);
        goto bail;
    }
    if (DSCaseHasSolution(aCase) == true) {
        goto bail;
    }
    if (DSCaseCd(aCase) != NULL || DSCaseCi(aCase) != NULL){
        if (DSCaseConditionsAreValid(aCase) == false) {
            goto bail;
        }
    }
    reduced =DSSSystemByRemovingAlgebraicConstraints(aCase->ssys);
    if (reduced == NULL){
        printf("Skipping unstable %u \n", aCase->caseNumber);
        goto bail;
    } else{
        DSSSystemFree(reduced);
    }
    
    //  condition to rule out cyclical cases and/or cases with solution but not valid constraints.
    if (DSCaseIsCyclical(aCase) == true){
        goto bail;
    }
    if (aCase->Xd->numberOfVariables == 0){
        goto bail;
    }

    unstableCase = DSSecureCalloc(sizeof(DSUnstableCase), 1);
    unstableCase->originalCase = DSCaseCopy(aCase);
    
    // perform main calculations
    // this function should populate Xb, pInverse and delta_unstable
    DSUnstableCaseIdentifyBlowingDependentVariables(unstableCase);
    
    if (unstableCase->Xd_b->numberOfVariables == 0){
        goto bail;
    }

    DSUInteger *bSignature;
    DSUInteger numberOfCases = 1, i;
    DSVariablePool *Xd_b = unstableCase->Xd_b;
    numberOfCases = pow(2, Xd_b->numberOfVariables);
    
    for (i=1; i <= numberOfCases; i++){
        // This is an array of length Xb->numberOfVariables containing values of 1 for blow-down and 2 for blow up.
        bSignature = DSBlowUpSignatureForCaseNumber(i, Xd_b);
        // construct case and return validity
        isValid = DSUnstableCaseConditionsAreValid(unstableCase, bSignature);
        if (bSignature != NULL)
            DSSecureFree(bSignature);
        // If it has a strict solution, break.
        if ( isValid == true ){
            break;
        }
    }
    
    if (unstableCase != NULL)
        DSuCaseFree(unstableCase);
    
bail:
    return isValid;
}

void dsCaseCalculateUnstableCaseIdentifier(DSCase * aCase, DSCase * nCase, DSGMASystem *gma, const char * prefix, DSUInteger i)
{
    DSUInteger caseNumber = 0;
    char temp[1000] = {'\0'};
    if (aCase == NULL) {
        DSError(M_DS_CASE_NULL, A_DS_ERROR);
        goto bail;
    }
    caseNumber = aCase->caseNumber;
    sprintf(temp, "%i", caseNumber);
    if (prefix != NULL) {
        DSCaseId(nCase) = DSSecureCalloc(sizeof(char), strlen(prefix)+2+strlen(temp));
        sprintf(DSCaseId(nCase), "%s_%i.%i", prefix, caseNumber, i);
    } else {
        DSCaseId(nCase) = DSSecureCalloc(sizeof(char), strlen(temp)+1);
        sprintf(DSCaseId(nCase), "%i.%i", caseNumber, i);
    }
bail:
    return;
}

void dsUnstableCaseCloneaCase(DSCase *aCase, DSCase *nCase, const DSDesignSpace *ds, DSUInteger i){
    
//    nCase->Xd = aCase->Xd;
//    nCase->Xd_a = aCase->Xd_a;
//    nCase->Xi = aCase->Xi;
//    DSCaseSig(nCase) = DSCaseSig(aCase);
//    nCase->signature_3d = aCase->signature_3d;
//    DSCaseCd(nCase) = DSCaseCd(aCase);
//    DSCaseCi(nCase) = DSCaseCi(aCase);
//    DSCaseZeta(nCase) = DSCaseZeta(aCase);
//    DSCaseDelta(nCase) = DSCaseDelta(aCase);
//    DSCaseU(nCase) = DSCaseU(aCase);
//    nCase->ssys = DSSSystemCopy(aCase->ssys);
    dsCaseCalculateUnstableCaseIdentifier(aCase, nCase, ds->gma, ds->casePrefix, i);
    
}

extern DSUnstableCase * DSUnstableCaseAddBoundaryMatrices(DSCase *aCase, const DSDesignSpace *ds)
{
    DSUnstableCase *unstableCase = NULL;
    DSVariablePool *Xd_b = NULL;
    DSCase * nCase = NULL;
    DSUInteger *bSignature = NULL;
    DSUInteger numberOfCases, i;
    bool isValid = false, isValidSolution = false;
    char name[1000];
        
    unstableCase = DSSecureCalloc(sizeof(DSUnstableCase), 1);
    unstableCase->originalCase = DSCaseCopy(aCase);
    

    // allocate dictionary
    unstableCase->ValidCases = DSDictionaryAlloc();
    
    // perform main calculations
    // this function should populate Xd_b, Xd_e, pInverse and delta_unstable
    DSUnstableCaseIdentifyBlowingDependentVariables(unstableCase);
    Xd_b = unstableCase->Xd_b;
    
    
    numberOfCases = pow(2, Xd_b->numberOfVariables); // why not number of knife edges? comment by Miguel on 13 jan 2020
    
    for (i=1; i <= numberOfCases; i++){
                    
            // This is an array of length Xb->numberOfVariables containing values of 1 for blow-down and 2 for blow up.
            bSignature = DSBlowUpSignatureForCaseNumber(i, Xd_b);
                
            // construct case and return validity
            isValid = DSUnstableCaseConditionsAreValid(unstableCase, bSignature);
        
            if (bSignature != NULL)
                DSSecureFree(bSignature);

            // If it has a strict solution, assign matrices and save to the dictionary
            if (isValid == true){
                        
                    DSUnstableCaseCreateBoundaryMatrices_alt2(unstableCase);
                
                    //////////////////////////
                    // Allocate memory for new Case and link data structures to uCase->originalCase
//                    nCase = DSSecureCalloc(sizeof(DSCase), 1);
//                    nCase->freeVariables = false;
                    nCase = DSCaseCopy(unstableCase->originalCase);
                    dsUnstableCaseCloneaCase(unstableCase->originalCase, nCase, ds, i);
                    nCase->U = DSMatrixCopy(unstableCase->U);
                    nCase->zeta = DSMatrixCopy(unstableCase->Zeta);
                    nCase->ssys->Xd_b = DSVariablePoolCopy(unstableCase->Xd_b);
                    
                    // check if case with merged solution is valid strict.
                    isValidSolution = DSCaseIsValid(nCase, true);
                    if (isValidSolution == false){
                        DSCaseFree(nCase);
                        continue;
                    }

//                    if (DSCaseU(nCase) != NULL && DSCaseZeta(nCase) != NULL)
//                        DSCaseRemoveRedundantBoundaries(nCase);
                
                    DSSSystemSetIsUnstable((DSSSystem *)DSCaseSSystem(nCase), true);
                
                    sprintf(name, "%i", i);
                    // save nCase to unstableCase->ValidCases. Name is: Case_identifier_i
                    if (DSDictionaryValueForName(unstableCase->ValidCases, name) == NULL)
                            DSDictionaryAddValueWithName(unstableCase->ValidCases, name, nCase);
                    /////////////////////////
            }
    }
    //DSuCaseFree(unstableCase);
    return unstableCase;
}


extern void dsUnstableCaseGetAdditionalConstraintMatrices(DSUnstableCase *uCase, const DSUInteger *bSignature)
{
    
        if (uCase == NULL || uCase->originalCase == NULL) {
            DSError(M_DS_MAT_NULL, A_DS_ERROR);
        }

        DSMatrix *zeros, *temp;
        DSMatrix *Cd_unstable_all, *Ci_unstable_all, *delta_unstable_all;
        DSUInteger rows = uCase->Xd_b->numberOfVariables;
        DSUInteger i, n;
        DSSSystem *collapsedSystem, *o_ssys = uCase->originalCase->ssys;
        DSUInteger *Xd_b_indices;

        // merge algebraic constraints if necessary and get indices of Xd_b
        if( DSVariablePoolNumberOfVariables(DSSSystemXd_a(o_ssys)) != 0){

            collapsedSystem = DSSSystemByRemovingAlgebraicConstraints(o_ssys);
            Xd_b_indices = DSVariablePoolIndicesOfSubPool(collapsedSystem->Xd, uCase->Xd_b);
            
            temp = DSSSystemAd(collapsedSystem);
            Ci_unstable_all = DSSSystemAi(collapsedSystem);
            delta_unstable_all = DSSSystemB(collapsedSystem);
            
            // let's add a block of zeros with dimentions rows = DSMatrixRows(Cd_unstable_all) and columns = number of Xd_a
            zeros = DSMatrixCalloc(DSMatrixRows(temp), DSSSystemXd_a(o_ssys)->numberOfVariables);
            Cd_unstable_all = DSMatrixAppendMatrices(temp, zeros, true);
            if (zeros != NULL)
                DSMatrixFree(zeros);
            if (temp != NULL)
                DSMatrixFree(temp);
            if (collapsedSystem != NULL)
                DSSSystemFree(collapsedSystem);
        }else {
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
        if (Cd_unstable_all != NULL)
            DSMatrixFree(Cd_unstable_all);
        if (Ci_unstable_all != NULL)
            DSMatrixFree(Ci_unstable_all);
        if (delta_unstable_all != NULL)
            DSMatrixFree(delta_unstable_all);
        if (Xd_b_indices != NULL)
            DSSecureFree(Xd_b_indices);
}

extern void DSUnstableCaseExpandConditionMatrices(DSUnstableCase *uCase)
{
    
    DSUInteger *bSignature;
    DSUInteger numberOfCases, i;
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
                uCase->originalCase->U      = DSMatrixCopy(uCase->U) ;
                uCase->originalCase->zeta   = DSMatrixCopy(uCase->Zeta);
//                DSGetVertices(uCase);
                DSSecureFree(string);
        }
    }
}

extern void DSUnstableCaseGetKnife(DSUnstableCase *uCase, const DSMatrix *gAd)
{
        // this function should populate the field DSVariablePool *knifeEdge.
        // the gauss of the matrix Ad is used to identify the location of real knife edges. A knife edge is present if the row for gAd only contains zeros.
    
        DSVariablePool *knifeEdge = DSVariablePoolAlloc();
        DSUInteger i;
        const char *name;
    
        for (i = 0; i < DSMatrixRows(gAd); i++){
            if (DSMatrixFirstNonZeroIndexAtRow(gAd, i) == 65000){
                name = DSVariablePoolAllVariableNames(uCase->originalCase->ssys->Xd_t)[i];
                DSVariablePoolAddVariableWithName(knifeEdge, name);
            }
        }
        uCase->knifeEdge = knifeEdge;
        if (DSVariablePoolNumberOfVariables(uCase->knifeEdge) == 0.0)
            printf("Warning! knifeEdge is zero! \n");
}

void dsUnstableCaseGetXd_b_2(DSVariablePool *Xd_b, DSVariablePool *Xd_e_temp,
                             DSVariablePool *Xd_e,
                             DSVariablePool *knifeEdge_temp, DSUInteger n)
{
        // This case corresponds to the scenario where the number of free variables is larger than the number of zeros in the diagonal.
        // pick n last number of variables from from Xd_e_temp to be knifes and Xd_b. The rest are Xd_e
//        printf("the number of variables n is %u \n", n);
        DSUInteger i;
        const DSVariable *current_var;
    
        for (i = 0; i < DSVariablePoolNumberOfVariables(Xd_e_temp); i++){
            current_var = DSVariablePoolAllVariables(Xd_e_temp)[i];
            if (i > DSVariablePoolNumberOfVariables(Xd_e_temp) - n - 1){
                DSVariablePoolAddVariableWithName(Xd_b, DSVariableName(current_var));
                DSVariablePoolAddVariableWithName(knifeEdge_temp, DSVariableName(current_var));
            } else {
                DSVariablePoolAddVariableWithName(Xd_e, DSVariableName(current_var));
            }
        }
}


void dsUnstableCaseGetXd_b_3(DSVariablePool *knifeEdge, DSVariablePool *knifeEdge_temp, DSUInteger n)
{
        // This case corresponds to the scenario where the number of zeros in the diagonal is
        // larger than the number of free variables. Pick n knife_edges from knife_edges temp and free knifeEdge_temp.
//        printf("the number of variables n is %u \n", n);
        DSUInteger i;
        const DSVariable *current_var;
    
        for (i = 0; i < DSVariablePoolNumberOfVariables(knifeEdge_temp); i++){
            current_var = DSVariablePoolAllVariables(knifeEdge_temp)[i];
            if (i <= n-1){
                DSVariablePoolAddVariableWithName(knifeEdge, DSVariableName(current_var));
            }
        }
}


extern void DSUnstableCaseGetXd_b(DSUnstableCase *uCase, DSMatrix *Ad)
{
    
        DSVariablePool *Xd_b = DSVariablePoolAlloc();
        DSVariablePool *Xd_e_temp = DSVariablePoolAlloc(), *Xd_e = DSVariablePoolAlloc();
        DSVariablePool  *knifeEdge = DSVariablePoolAlloc(), *knifeEdge_temp  = DSVariablePoolAlloc();
        DSUInteger i, free_variable, n, case_id, zeros_diag;
        const DSVariable *blow_variable, *equal_variable;
        const DSVariablePool *Xd = uCase->originalCase->Xd;
        double value;
    
        // check the diagonal of the Ad matrix. If it is zero, variable belongs to the Xd_b and knifeEdge pool. Else it belongs to Xd_e
        for (i=0; i<DSMatrixRows(Ad); i++){
                    value = DSMatrixDoubleValue(Ad, i, i);
                    // if value is zero, add corresponding variable of Xd to Xd_b, else it has a solution and belong to Xd_e
                    if (value == 0.0){
                            blow_variable = DSVariablePoolAllVariables(Xd)[i];
                            DSVariablePoolAddVariableWithName(Xd_b, DSVariableName(blow_variable));
                            DSVariablePoolAddVariableWithName(knifeEdge_temp, DSVariableName(blow_variable));
                    }else{
                                equal_variable = DSVariablePoolAllVariables(Xd)[i];
                                DSVariablePoolAddVariableWithName(Xd_e_temp, DSVariableName(equal_variable));
                         }
        }
    
        // check if the number of zeros in the diagonal (number of variables in Xd_b) is smaller than the number of free variables
        free_variable = DSMatrixRows(Ad) - DSMatrixRank(Ad);
        zeros_diag = DSVariablePoolNumberOfVariables(Xd_b);
    
        case_id = (free_variable == zeros_diag) ? 1: ( (free_variable > zeros_diag) ? 2 : 3);
        n = free_variable - DSVariablePoolNumberOfVariables(Xd_b);
    
        switch (case_id) {
            case 1:
                uCase->knifeEdge = knifeEdge_temp;
                uCase->Xd_b = Xd_b;
                uCase->Xd_e = Xd_e_temp;
                if (knifeEdge != NULL)
                    DSVariablePoolFree(knifeEdge);
                if (Xd_e != NULL)
                    DSVariablePoolFree(Xd_e);
                
                break;
            
            case 2:
                dsUnstableCaseGetXd_b_2(Xd_b, Xd_e_temp, Xd_e, knifeEdge_temp, n);
                uCase->knifeEdge = knifeEdge_temp;
                uCase->Xd_b = Xd_b;
                uCase->Xd_e = Xd_e;
                if (knifeEdge != NULL)
                    DSVariablePoolFree(knifeEdge);
                if (Xd_e_temp != NULL)
                    DSVariablePoolFree(Xd_e_temp);
                break;
                
            case 3:
                n = DSVariablePoolNumberOfVariables(Xd_b) - free_variable;
                dsUnstableCaseGetXd_b_3(knifeEdge, knifeEdge_temp, n);
                uCase->knifeEdge = knifeEdge;
                uCase->Xd_b = Xd_b;
                uCase->Xd_e = Xd_e;
                if (knifeEdge_temp != NULL)
                    DSVariablePoolFree(knifeEdge_temp);
                if (Xd_e_temp != NULL)
                    DSVariablePoolFree(Xd_e_temp);
                break;
                
            default:
                printf("an error ocurred while setting Xd_b \n");
                break;
        }
    
        if(Xd_b->numberOfVariables ==0 ){
            printf("Warning. Anomalous case. See Case number %u \n", uCase->originalCase->caseNumber);
            printf("The matrix Ad is: \n");
            DSMatrixPrint(Ad);
        }

}

extern void DSUnstableCaseIdentifyBlowingDependentVariables(DSUnstableCase *uCase)
{
    
        DSMatrix *Ad, *Ai, *B;
        const DSMatrix *pInverse = NULL;
        DSSSystem *collapsedSystem = NULL, *o_ssys = uCase->originalCase->ssys;
        DSMatrixArray *gaussArray;
    
        // Calculation of the pseudo inverse of Ad.
        // check if system needs to be collapsed and then calculate pseudoinverse.
        if( DSVariablePoolNumberOfVariables(DSSSystemXd_a(o_ssys)) != 0){
            
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
    
        // Identify Xd_b and Xd_e and assign.
        DSUnstableCaseGetXd_b(uCase, Ad);
    
        if (uCase->Xd_b->numberOfVariables == 0)
            goto bail;
    
    
        // Calculate pseudo-inverse of Ad by only considering a subset of the matrix Ad.
        pInverse = DSUnstableCaseGetSubSetPseudoInverse(uCase->originalCase->ssys->Xd_t,
                                                        uCase->Xd_e,
                                                        uCase->Xd_b,
                                                        Ad);
        uCase->pInverse = pInverse;
    
        // Get matrices Equality and KnifeEdge
        DSUnstableCaseGetEqualityAndKnife(pInverse, gaussArray, uCase);
    
        // delete variables.
        if( DSVariablePoolNumberOfVariables(DSSSystemXd_a(o_ssys))!= 0 && collapsedSystem != NULL)
            DSSSystemFree(collapsedSystem);
        if (Ad != NULL)
            DSMatrixFree(Ad);
        if (Ai != NULL)
            DSMatrixFree(Ai);
        if (B != NULL)
            DSMatrixFree(B);
        if (gaussArray != NULL)
            DSMatrixArrayFree(gaussArray);
bail:
    return;
}

extern void DSUnstableCasePrintLinearProblem(glp_prob *lp, const DSUnstableCase *uCase, const DSUInteger *bSignature)
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

extern void DSUnstableCaseDetermineBlowingBehavior(glp_prob *lp, DSUnstableCase *uCase, const DSUInteger *bSignature)
{
    
        const char *name;
        DSVariablePool *Xd_b = uCase->Xd_b;
        DSVariablePool *knifeEdge = uCase->knifeEdge;
        DSUInteger i, numberOfBlowing = Xd_b->numberOfVariables, k, j;
        DSUInteger numberOfKnifes = uCase->knifeEdge->numberOfVariables;
        DSUInteger numberOfXd = uCase->originalCase->Xd->numberOfVariables;
        double value;
        DSMatrix *gAi = DSMatrixArrayMatrix(uCase->Knife, 0);
        DSMatrix *gB = DSMatrixArrayMatrix(uCase->Knife, 1);
        DSMatrix *Cd_knife, *Ci_knife, *delta_knife;
    
        if (uCase->originalCase->delta != NULL){
            k = DSMatrixRows(uCase->originalCase->delta) + uCase->Xd_b->numberOfVariables + uCase->Xd_e->numberOfVariables;
        } else {
            k = uCase->Xd_b->numberOfVariables + uCase->Xd_e->numberOfVariables;
        }
    
        // initialize matrices Ci_knife, delta_knife & Cd_knife.
        Ci_knife = DSMatrixCopy(gAi);
        Cd_knife = DSMatrixCalloc(numberOfKnifes, numberOfXd);
        delta_knife = DSMatrixCopy(gB);
    
    
        //loop over bSignature to assign blowing behavior.
        for (i=0; i<numberOfKnifes; i++){
            
            name = DSVariablePoolAllVariableNames(uCase->knifeEdge)[i];
            value = glp_get_row_prim(lp, k+i+1) + DSMatrixDoubleValue(gB, i, 0);
            //        printf("value is the result of %f - %f \n", glp_get_row_prim(lp, k+i+1), DSMatrixDoubleValue(gB, i, 0) );
            DSVariablePoolSetValueForVariableWithName(knifeEdge, name, value);
            
            // adjust matrices Ci_knife and delta_knife.
            if (value < 0){
                for (j=0 ; j<DSMatrixColumns(gAi); j++){
                    DSMatrixSetDoubleValue(Ci_knife, i, j, DSMatrixDoubleValue(Ci_knife, i, j)*-1.0);
                }
                DSMatrixSetDoubleValue(delta_knife, i, 0, DSMatrixDoubleValue(delta_knife, i, 0)*-1.0);
            }
        }
    
        for (i=0; i<numberOfBlowing; i++){
                name = DSVariablePoolAllVariableNames(uCase->Xd_b)[i];
                // set blowing behavior.
                if (bSignature[i] == 1){
                    DSVariablePoolSetValueForVariableWithName(Xd_b, name, 1.0E12);
                } else {
                    DSVariablePoolSetValueForVariableWithName(Xd_b, name, 1.0E-12);
                }
        }
        uCase->Cd_knife = Cd_knife;
        uCase->Ci_knife = Ci_knife;
        uCase->delta_knife = delta_knife;
        uCase->originalCase->ssys->Xd_b = DSVariablePoolCopy(Xd_b) ;
}

extern void DSUnstableCaseCreateBoundaryMatrices(DSUnstableCase *uCase)
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
    
    DSUInteger numberOfXi = 0, *Xd_t_indices, *Xd_b_indices, aa;
    DSUInteger numberOfBlowing = uCase->Xd_b->numberOfVariables, i, j;
    DSMatrix * W = NULL, *Ai, *Zeta, *U, *Cd_xd_t, *B_xd_t, *temp, *delta;
    DSCase *aCase = uCase->originalCase;
    double x;
    
    DSMatrix *Cd, *Ci;
    DSMatrix *B = DSSSystemB(DSCaseSSys(aCase));
    
    DSSSystem *o_ssys, *collapsedSSystem;
    const DSVariablePool *Xd_t = uCase->originalCase->ssys->Xd_t, *Xd = uCase->originalCase->ssys->Xd;
    const DSUInteger numberOfColumns = Xd_t->numberOfVariables;
    
    o_ssys = uCase->originalCase->ssys;
    numberOfXi = DSVariablePoolNumberOfVariables(DSCaseXi(aCase));
    
    // we first construct Cd, Ci and delta by merging different matrices.
    if ( DSCaseDelta(aCase) != NULL){
        temp = DSMatrixAppendMatrices(DSCaseCd(aCase), uCase->Cd_unstable, false);
        Cd = DSMatrixAppendMatrices(temp, uCase->Cd_knife, false);
        DSMatrixFree(temp);
        
        temp = DSMatrixAppendMatrices(DSCaseCi(aCase), uCase->Ci_unstable, false);
        Ci = DSMatrixAppendMatrices(temp, uCase->Ci_knife, false);
        DSMatrixFree(temp);
        
        temp = DSMatrixAppendMatrices(DSCaseDelta(aCase), uCase->delta_unstable, false);
        delta = DSMatrixAppendMatrices(temp, uCase->delta_knife, false);
        DSMatrixFree(temp);
    } else {
        Cd = DSMatrixAppendMatrices(uCase->Cd_unstable, uCase->Cd_knife, false);
        Ci = DSMatrixAppendMatrices(uCase->Ci_unstable, uCase->Ci_knife, false);
        delta = DSMatrixAppendMatrices(uCase->delta_unstable, uCase->delta_knife, false);
    }
//    printf("Reporting from within DSUnstableCaseCreateBoundaryMatrices. The matrix Cd, Ci and delta are: \n");
//    DSMatrixPrint(Cd);
//    printf("Ci: \n");
//    DSMatrixPrint(Ci);
//    printf("delta: \n");
//    DSMatrixPrint(delta);
    
    // Matrices Cd and B of the original system need to be cut to get rid of Xd_a,
    // this is because my pInverse only contains Xd_t.
    if( DSVariablePoolNumberOfVariables(DSSSystemXd_a(o_ssys)) != 0){
        
        collapsedSSystem = DSSSystemByRemovingAlgebraicConstraints(o_ssys);
        // Generate matrix Cd_xd_t as a subset of matrix Cd.
        Xd_t_indices = DSVariablePoolIndicesOfSubPool(Xd, Xd_t);
        Cd_xd_t = DSMatrixSubMatrixIncludingColumns(Cd, numberOfColumns, Xd_t_indices);
        
//        printf("-----------Reporting from DSUnstableCaseCreateBoundaryMatrices. Matrix Cd for case %i: \n ", uCase->originalCase->caseNumber);
//        DSMatrixPrint(Cd);
//        printf("Matrix Ad is: \n");
//        DSMatrixPrint(DSSSystemAd(o_ssys));
//        printf("The pseudoinverse of Ad is: \n");
//        DSMatrixPrint(DSMatrixPseudoInverse(DSSSystemAd(o_ssys)));
//        printf("Variable pool Xd is: \n");
//        
//        for(aa=0; aa<DSVariablePoolNumberOfVariables(Xd); aa++)
//            printf("%s \n",DSVariablePoolAllVariableNames(uCase->originalCase->Xd)[aa]);
        
        // Generate vector B_xd_t as a subset of vector B
        B_xd_t = DSMatrixSubMatrixIncludingRows(B, numberOfColumns, Xd_t_indices);
        
        // now calculate W & Zeta
        W = DSMatrixByMultiplyingMatrix(Cd_xd_t, uCase->pInverse);
        Zeta = DSMatrixByMultiplyingMatrix(W, B_xd_t);
        DSMatrixAddByMatrix(Zeta, delta);
        if (numberOfXi != 0) {
            Ai = DSSSystemAi(collapsedSSystem);
            U = DSMatrixByMultiplyingMatrix(W, Ai);
            if (Ci != NULL)
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
                x = DSMatrixDoubleValue(Cd_xd_t, j, Xd_b_indices[i]);
                if( x != 0.0 ){
                    DSMatrixSetDoubleValue(Zeta, j, 0,
                                           DSMatrixDoubleValue(Zeta, j, 0) + x*log10(uCase->Xd_b->variables[i]->value));
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
                    if (Ci != NULL)
                        DSMatrixSubstractByMatrix(U, Ci);
                    DSMatrixMultiplyByScalar(U, -1.0);
                    DSMatrixFree(Ai);
                }
        
                // now rows of Zeta are adjusted to consider for numerical values of blowing dependent variables
                // in the condition matrices.
                Xd_b_indices = DSVariablePoolIndicesOfSubPool(Xd_t, uCase->Xd_b);
//                printf("reporting from DSUnstableCaseCreateBoundaryMatrices .... \n");
//                printf("Xd_b indices are: %u %u \n", Xd_b_indices[0], Xd_b_indices[1] );
//                printf("The variable pool Xd_b is: \n");
//                DSVariablePoolPrint(uCase->Xd_b);
    
                // loop over columns
                for(i=0; i<numberOfBlowing; i++){
                    // loop over rows
                    for (j=0; j< DSMatrixRows(Cd); j++){
                        x = DSMatrixDoubleValue(Cd, j, Xd_b_indices[i]);
                        if( x != 0.0 ){
                            DSMatrixSetDoubleValue(Zeta, j, 0,
                                                   DSMatrixDoubleValue(Zeta, j, 0) + x*log10(uCase->Xd_b->variables[i]->value));
//                            printf("variable was assigned to blow in Z matrix \n");
                        }
                    }
                }
        
//                printf("The pseudoinverse of Ad is: \n");
//                DSMatrixPrint(uCase->pInverse);
//                printf("The matrix Ai is: \n");
//                DSMatrixPrint(Ai);
//                printf("The matrix B is: \n");
//                DSMatrixPrint(B);
        
        
                DSMatrixFree(W);
                DSMatrixFree(B);
                DSSecureFree(Xd_b_indices);
    }
    uCase->U = U;
    uCase->Zeta = Zeta;
    
bail:
    return;
    
}

static void dsUnstableCaseMatricesByPartitioningAuxiliaryVariables(const DSSSystem * ssystem,
                                                                   DSMatrix ** ADat, DSMatrix ** ADaa,
                                                                   DSMatrix ** AIa, DSMatrix ** Ba)
{
    DSUInteger * indices;
    const DSVariablePool * Xd_t, *Xd_a, *Xd;
    DSMatrix *Ad_a, *temp;
    if (ssystem == NULL) {
        DSError(M_DS_SSYS_NULL, A_DS_ERROR);
        goto bail;
    }
    if ((ADat && ADaa && AIa && Ba) == false) {
        DSError(M_DS_MAT_NULL, A_DS_ERROR);
        goto bail;
    }
    *ADat = NULL;
    *ADaa = NULL;
    *AIa = NULL;
    *Ba = NULL;
    if (DSSSystemXd_a(ssystem) == NULL) {
        DSError(M_DS_VAR_NULL ": Xd_a variable pool is NULL", A_DS_ERROR);
        goto bail;
    }
    if (DSVariablePoolNumberOfVariables(DSSSystemXd_a(ssystem)) > DSVariablePoolNumberOfVariables(DSSSystemXd(ssystem))) {
        DSError(M_DS_WRONG ": Number of algebraic variables exceeds number of total variables", A_DS_ERROR);
        goto bail;
    }
    if (DSVariablePoolNumberOfVariables(DSSSystemXd_a(ssystem)) == 0) {
        goto bail;
    }
    if (DSVariablePoolNumberOfVariables(DSSSystemXd_t(ssystem)) == 0) {
        DSError(M_DS_WRONG ": System does not have dynamic variables.", A_DS_WARN);
        goto bail;
    }
    Xd = DSSSystemXd(ssystem);
    Xd_t = DSSSystemXd_t(ssystem);
    Xd_a = DSSSystemXd_a(ssystem);
    indices = DSVariablePoolIndicesOfSubPool(Xd, Xd_a);
    temp = DSSSystemAd(ssystem);
    Ad_a = DSMatrixSubMatrixIncludingRows(temp, DSVariablePoolNumberOfVariables(Xd_a), indices);
    DSMatrixFree(temp);
    temp = DSSSystemAi(ssystem);
    *AIa = DSMatrixSubMatrixIncludingRows(temp, DSVariablePoolNumberOfVariables(Xd_a), indices);
    DSMatrixFree(temp);
    temp = DSSSystemB(ssystem);
    *Ba = DSMatrixSubMatrixIncludingRows(temp, DSVariablePoolNumberOfVariables(Xd_a), indices);
    DSMatrixFree(temp);
    *ADat = DSMatrixSubMatrixExcludingColumns(Ad_a, DSVariablePoolNumberOfVariables(Xd_a), indices);
    *ADaa = DSMatrixSubMatrixIncludingColumns(Ad_a, DSVariablePoolNumberOfVariables(Xd_a), indices);
    
    DSMatrixFree(Ad_a);
    DSSecureFree(indices);
    
    
bail:
    return;
}

static void dsUnstableCaseSolutionForAlgebraicConstraints(const DSSSystem * ssystem,
                                                          const DSMatrix * M_a,
                                                          DSMatrix ** Ad_at,
                                                          DSMatrix ** Ai_a,
                                                          DSMatrix ** B_a
                                                          ) {
        DSMatrix * LHS, * RHS;
        if (ssystem == NULL) {
            DSError(M_DS_SSYS_NULL, A_DS_ERROR);
            goto bail;
        }
        if ((M_a && Ad_at && Ai_a && B_a) == false) {
            DSError(M_DS_MAT_NULL, A_DS_ERROR);
            goto bail;
        }
        if (DSSSystemXd_a(ssystem) == NULL) {
            DSError(M_DS_VAR_NULL ": Xd_a variable pool is NULL", A_DS_ERROR);
            goto bail;
        }
        if (DSVariablePoolNumberOfVariables(DSSSystemXd_a(ssystem)) > DSVariablePoolNumberOfVariables(DSSSystemXd(ssystem))) {
            DSError(M_DS_WRONG ": Number of algebraic variables exceeds number of total variables", A_DS_ERROR);
            goto bail;
        }
        if (DSVariablePoolNumberOfVariables(DSSSystemXd_a(ssystem)) == 0) {
            goto bail;
        }
        if (DSVariablePoolNumberOfVariables(DSSSystemXd_t(ssystem)) == 0) {
            DSError(M_DS_WRONG ": System does not have dynamic variables.", A_DS_WARN);
            goto bail;
        }
    
        RHS = *Ad_at;
        LHS = DSMatrixByMultiplyingMatrix(M_a, RHS);
        *Ad_at = DSMatrixByMultiplyingScalar(LHS, -1.0);
        DSMatrixFree(LHS);
        DSMatrixFree(RHS);
    
        RHS = *Ai_a;
        LHS = DSMatrixByMultiplyingMatrix(M_a, RHS);
        *Ai_a = DSMatrixByMultiplyingScalar(LHS, -1.0);
        DSMatrixFree(LHS);
        DSMatrixFree(RHS);
    
        RHS = *B_a;
        LHS = DSMatrixByMultiplyingMatrix(M_a, RHS);
        *B_a = DSMatrixByMultiplyingScalar(LHS, -1.0);
        DSMatrixFree(LHS);
        DSMatrixFree(RHS);
    bail:
        return;
}

static void dsUnstableCaseConstraintsByRemovingAlgebraicConstraintsInternal(DSMatrix ** Cd,
                                                                            DSMatrix ** Ci,
                                                                            DSMatrix **delta,
                                                                            DSSSystem *ssys,
                                                                            const DSMatrix * MAd_at,
                                                                            const DSMatrix * MAi_a,
                                                                            const DSMatrix * MB_a)
{
        DSMatrix *Kd, *Ki, *Kd_t, *Kd_a, *a, *LHS, *RHS, *Cd_a, *temp;
        DSUInteger * indices, auxiliary_count;
        const DSVariablePool * Xd_t, *Xd_a, *Xd;

        if (MAd_at == NULL || MAi_a == NULL || MB_a == NULL) {
            DSError(M_DS_MAT_NULL, A_DS_ERROR);
            goto bail;
        }
    
        Xd = DSSSystemXd(ssys);
        Xd_t = DSSSystemXd_t(ssys);
        Xd_a = DSSSystemXd_a(ssys);
        indices = DSVariablePoolIndicesOfSubPool(Xd, Xd_a);
        auxiliary_count = DSVariablePoolNumberOfVariables(Xd_a);
    
        Kd = *Cd;
        Ki = *Ci;
        a = *delta;

        Kd_t = DSMatrixSubMatrixExcludingColumns(Kd, auxiliary_count, indices);
        Kd_a = DSMatrixSubMatrixIncludingColumns(Kd, auxiliary_count, indices);
        DSMatrixFree(Kd);
    
        LHS = Kd_t;
        RHS = DSMatrixByMultiplyingMatrix(Kd_a, MAd_at);
        Kd_t = DSMatrixByAddingMatrix(LHS, RHS);
        DSMatrixFree(LHS);
        DSMatrixFree(RHS);
    
        LHS = Ki;
        RHS = DSMatrixByMultiplyingMatrix(Kd_a, MAi_a);
        Ki = DSMatrixByAddingMatrix(LHS, RHS);
        DSMatrixFree(LHS);
        DSMatrixFree(RHS);
    
        LHS = a;
    
        // adjust sign of matrix Kd_a
        temp = DSMatrixByMultiplyingScalar(Kd_a, -1);
        DSMatrixFree(Kd_a);
        Kd_a = temp;
    
        // transforming matrix LHS to logarithmic coordinates is not neccesary because it is already in log
//        for (ii=0; ii<DSMatrixRows(LHS); ii++ )
//            DSMatrixSetDoubleValue(LHS, ii, 0, log10(DSMatrixDoubleValue(LHS, ii, 0)));
    
        RHS = DSMatrixByMultiplyingMatrix(Kd_a, MB_a);
        a = DSMatrixByAddingMatrix(LHS, RHS);
    
        // transform a back into cartesian coordinates is also not necessary!
//        for (ii=0; ii<DSMatrixRows(a); ii++ )
//            DSMatrixSetDoubleValue(a, ii, 0, pow(10, DSMatrixDoubleValue(a, ii, 0)));
    
        DSMatrixFree(LHS);
        DSMatrixFree(RHS);
    
        Cd_a = DSMatrixCalloc(DSMatrixRows(Kd_a), DSMatrixColumns(Kd_a));
        DSMatrixFree(Kd_a);
    
        *Cd = DSMatrixAppendMatrices(Kd_t, Cd_a, true) ;
        *Ci = Ki;
        *delta = a;
    
        DSMatrixFree(Kd_t);
        DSMatrixFree(Cd_a);
        DSSecureFree(indices);
    
bail:
    return;

}


static void dsUnstableCaseConditionMatricesByMergingAlgebraicConstraints(DSMatrix **Cd, DSMatrix **Ci, DSMatrix ** delta, DSCase *aCase)
{
    
        DSSSystem *ssys = aCase->ssys;
        DSMatrix * M_a, * Ad_at = NULL, * Ad_aa = NULL, * Ai_a = NULL, * B_a = NULL;
        if (aCase->ssys == NULL) {
            DSError(M_DS_SSYS_NULL, A_DS_ERROR);
            goto bail;
        }
        if (DSSSystemXd_a(ssys) == NULL) {
            DSError(M_DS_VAR_NULL ": Xd_a variable pool is NULL", A_DS_ERROR);
            goto bail;
        }
        if (DSVariablePoolNumberOfVariables(DSSSystemXd_a(ssys)) > DSVariablePoolNumberOfVariables(DSSSystemXd(ssys))) {
            DSError(M_DS_WRONG ": Number of algebraic variables exceeds number of total variables", A_DS_ERROR);
            goto bail;
        }
        if (DSVariablePoolNumberOfVariables(DSSSystemXd_t(ssys)) == 0) {
            DSError(M_DS_WRONG ": System does not have dynamic variables.", A_DS_WARN);
            goto bail;
        }
        dsUnstableCaseMatricesByPartitioningAuxiliaryVariables(ssys,
                                                               &Ad_at,&Ad_aa,&Ai_a,&B_a);
        if (Ad_aa == NULL) {
            goto bail;
        }
        M_a = DSMatrixInverse(Ad_aa);
        if (M_a == NULL) {
            goto bail;
        }
        dsUnstableCaseSolutionForAlgebraicConstraints(ssys,
                                                      M_a,
                                                      &Ad_at,
                                                      &Ai_a,
                                                      &B_a);
    
        DSMatrixFree(M_a);
        dsUnstableCaseConstraintsByRemovingAlgebraicConstraintsInternal(Cd, Ci, delta, ssys,
                                                                        Ad_at,
                                                                        Ai_a,
                                                                        B_a);
    
    bail:
        if (Ad_aa != NULL)
            DSMatrixFree(Ad_aa);
        if (Ad_at != NULL)
            DSMatrixFree(Ad_at);
        if (Ai_a != NULL)
            DSMatrixFree(Ai_a);
        if (B_a != NULL)
            DSMatrixFree(B_a);
}

extern void DSUnstableCaseCreateBoundaryMatrices_alt(DSUnstableCase *uCase)
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
    
        DSUInteger numberOfXi = 0, *Xd_b_indices = NULL;
        DSUInteger numberOfBlowing = uCase->Xd_b->numberOfVariables, i, j;
        DSMatrix * W = NULL, *Ai = NULL, *Zeta = NULL, *U = NULL, *temp = NULL, *delta = NULL, *Ad_cut = NULL;
        DSMatrix * Cd_xd_no_b = NULL;
        DSCase *aCase = uCase->originalCase;
        DSMatrix *Cd = NULL, *Ci = NULL;
        DSMatrix *B = NULL, *InvAd = NULL;
        DSSSystem *o_ssys;
        DSVariablePool *Xd_b = uCase->Xd_b;
        DSVariablePool *Xd = uCase->originalCase->ssys->Xd;
        double x;
        
        o_ssys = uCase->originalCase->ssys;
        numberOfXi = DSVariablePoolNumberOfVariables(DSCaseXi(aCase));
        
        // we first construct Cd, Ci and delta by merging different matrices.
        if ( DSCaseDelta(aCase) != NULL){
            temp = DSMatrixAppendMatrices(DSCaseCd(aCase), uCase->Cd_unstable, false);
            Cd = DSMatrixAppendMatrices(temp, uCase->Cd_knife, false);
            DSMatrixFree(temp);
            
            temp = DSMatrixAppendMatrices(DSCaseCi(aCase), uCase->Ci_unstable, false);
            Ci = DSMatrixAppendMatrices(temp, uCase->Ci_knife, false);
            DSMatrixFree(temp);
            
            temp = DSMatrixAppendMatrices(DSCaseDelta(aCase), uCase->delta_unstable, false);
            delta = DSMatrixAppendMatrices(temp, uCase->delta_knife, false);
            DSMatrixFree(temp);
        } else {
            Cd = DSMatrixAppendMatrices(uCase->Cd_unstable, uCase->Cd_knife, false);
            Ci = DSMatrixAppendMatrices(uCase->Ci_unstable, uCase->Ci_knife, false);
            delta = DSMatrixAppendMatrices(uCase->delta_unstable, uCase->delta_knife, false);
        }
    
        // If the system has algebraic constraints, process constraint matrices Cd, Ci and delta to merge those algebraic constraints.
        if (DSVariablePoolNumberOfVariables(DSSSystemXd_a(o_ssys)) != 0){
            dsUnstableCaseConditionMatricesByMergingAlgebraicConstraints(&Cd, &Ci, &delta, aCase);
        }
    
        // Matrices Cd and Ad of the original system need to be cut to get rid of blowing variables
        // Generate matrix Cd_xd_t as a subset of matrix Cd.
        Xd_b_indices= DSVariablePoolIndicesOfSubPool(Xd, Xd_b);
        Cd_xd_no_b = DSMatrixSubMatrixExcludingColumns(Cd, numberOfBlowing, Xd_b_indices);
    
        // Now cut Ad and get its inverse.
        Ad_cut = DSMatrixSubMatrixExcludingRowsAndColumns(DSSSystemAd(o_ssys), numberOfBlowing, numberOfBlowing,    Xd_b_indices, Xd_b_indices);
        InvAd = DSMatrixInverse(Ad_cut);
        if (Ad_cut != NULL)
            DSMatrixFree(Ad_cut);

        // now calculate W & Zeta
        B = DSMatrixSubMatrixExcludingRows(DSSSystemB(DSCaseSSys(aCase)), numberOfBlowing, Xd_b_indices);
        if (Cd_xd_no_b != NULL && InvAd != NULL)
            W = DSMatrixByMultiplyingMatrix(Cd_xd_no_b, InvAd);
        if (W != NULL && B != NULL)
            Zeta = DSMatrixByMultiplyingMatrix(W, B);

        if (Zeta != NULL)
            DSMatrixAddByMatrix(Zeta, delta);
        else
            Zeta = DSMatrixCopy(delta);

        if (numberOfXi != 0) {
            
                Ai = DSMatrixSubMatrixExcludingRows(DSSSystemAi(o_ssys), numberOfBlowing, Xd_b_indices);
                if (Ai != NULL && W != NULL)
                    U = DSMatrixByMultiplyingMatrix(W, Ai);
                if (Ci != NULL && U != NULL)
                    DSMatrixSubstractByMatrix(U, Ci);
                if (U == NULL){
                    U = DSMatrixCopy(Ci);
                    DSMatrixMultiplyByScalar(U, -1.0);
                }
                DSMatrixMultiplyByScalar(U, -1.0);
                if (Ai != NULL)
                    DSMatrixFree(Ai);
        }
    
    ///////////////////////////
        // now rows of Zeta are adjusted to consider for numerical values of blowing dependent variables
        // in the condition matrices.
    
        // loop over columns of Cd that contain blowing variables
        for(i=0; i<numberOfBlowing; i++){
            // loop over rows of those columns.
            for (j=0; j< DSMatrixRows(Cd); j++){
                x = DSMatrixDoubleValue(Cd, j, Xd_b_indices[i]);
                if( x != 0.0 ){
                    DSMatrixSetDoubleValue(Zeta, j, 0,
                                           DSMatrixDoubleValue(Zeta, j, 0) + x*log10(uCase->Xd_b->variables[i]->value));
                }
            }
        }

        if (delta != NULL)
            DSMatrixFree(delta);
        if (W != NULL)
            DSMatrixFree(W);
        if (B != NULL)
            DSMatrixFree(B);
        if (InvAd != NULL)
            DSMatrixFree(InvAd);
        if (Xd_b_indices != NULL)
            DSSecureFree(Xd_b_indices);
    
        uCase->U = U;
        uCase->Zeta = Zeta;
    
    bail:
        return;
    
}

static void dsUnstableCaseProcessMB_blow(DSMatrix *MB_blow){
    
    // The idea is to process the vector MB_blow to restrict values to [-12 and 12]
    DSUInteger i;
    for (i = 0; i<DSMatrixRows(MB_blow); i++){
        if (DSMatrixDoubleValue(MB_blow, i, 0) > 12.0)
            DSMatrixSetDoubleValue(MB_blow, i, 0, 12.0);
        if (DSMatrixDoubleValue(MB_blow, i, 0) < -12.0)
            DSMatrixSetDoubleValue(MB_blow, i, 0, -12.0);
    }
}


extern void DSUnstableCaseCreateBoundaryMatrices_alt2(DSUnstableCase *uCase)
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
    
    DSUInteger numberOfXi = 0, *Xd_b_indices = NULL;
    DSUInteger numberOfBlowing = uCase->Xd_b->numberOfVariables, i, n, col, row;
    DSMatrix * W = NULL, *Ai = NULL, *Zeta = NULL, *U = NULL, *temp = NULL, *delta = NULL, *w;
    DSMatrix *pInvAd, *I_pInvA, *Identity ;
    const DSMatrix *pInverse;
    DSCase *aCase = uCase->originalCase;
    DSMatrix *Cd = NULL, *Ci = NULL, *Ad, *Ad2;
    DSMatrix *B = NULL, *InvAd = NULL, * MB_blow = NULL;
    DSSSystem *o_ssys;
    DSVariablePool *Xd_b = uCase->Xd_b;
    DSVariablePool *Xd = uCase->originalCase->ssys->Xd;
    
    o_ssys = uCase->originalCase->ssys;
    numberOfXi = DSVariablePoolNumberOfVariables(DSCaseXi(aCase));
    
    // we first construct Cd, Ci and delta by merging different matrices.
    if ( DSCaseDelta(aCase) != NULL){
        temp = DSMatrixAppendMatrices(DSCaseCd(aCase), uCase->Cd_unstable, false);
        Cd = DSMatrixAppendMatrices(temp, uCase->Cd_knife, false);
        if (temp != NULL)
            DSMatrixFree(temp);
        
        temp = DSMatrixAppendMatrices(DSCaseCi(aCase), uCase->Ci_unstable, false);
        Ci = DSMatrixAppendMatrices(temp, uCase->Ci_knife, false);
        if (temp != NULL)
            DSMatrixFree(temp);
        
        temp = DSMatrixAppendMatrices(DSCaseDelta(aCase), uCase->delta_unstable, false);
        delta = DSMatrixAppendMatrices(temp, uCase->delta_knife, false);
        if (temp != NULL)
            DSMatrixFree(temp);
    } else {
        Cd = DSMatrixAppendMatrices(uCase->Cd_unstable, uCase->Cd_knife, false);
        Ci = DSMatrixAppendMatrices(uCase->Ci_unstable, uCase->Ci_knife, false);
        delta = DSMatrixAppendMatrices(uCase->delta_unstable, uCase->delta_knife, false);
    }
    
    // If the system has algebraic constraints, process constraint matrices Cd, Ci and delta to merge those algebraic constraints.
    if (DSVariablePoolNumberOfVariables(DSSSystemXd_a(o_ssys)) != 0){
        dsUnstableCaseConditionMatricesByMergingAlgebraicConstraints(&Cd, &Ci, &delta, aCase);
    }
    
    // Matrices Cd and Ad of the original system need to be cut to get rid of blowing variables
    // Generate matrix Cd_xd_t as a subset of matrix Cd.
    Xd_b_indices= DSVariablePoolIndicesOfSubPool(Xd, Xd_b);
    Ad = DSSSystemAd(o_ssys);
    Ad2 = DSSSystemAd(o_ssys);
    
    // Set both rows and colums of blowing variables to zero.
    for(n = 0; n<numberOfBlowing; n++){
        for(col = 0; col<DSMatrixColumns(Ad2); col++)
            DSMatrixSetDoubleValue(Ad2, Xd_b_indices[n], col, 0.0);
        for(row = 0; row<DSMatrixRows(Ad2); row++)
            DSMatrixSetDoubleValue(Ad2, row, Xd_b_indices[n], 0.0);
    }

    pInverse = DSMatrixPseudoInverse(Ad2);
    if (Ad2 != NULL)
        DSMatrixFree(Ad2);
    pInvAd = DSMatrixByMultiplyingMatrix(pInverse, Ad);
    Identity = DSMatrixIdentity(DSMatrixRows(pInvAd));
    I_pInvA = DSMatrixBySubstractingMatrix(Identity, pInvAd);
    
    if (Ad != NULL)
        DSMatrixFree(Ad);
    if (pInvAd != NULL)
        DSMatrixFree(pInvAd);
    if (Identity != NULL)
        DSMatrixFree(Identity);
    
    // now calculate W & Zeta
    B = DSSSystemB(DSCaseSSys(aCase));
    if (Cd != NULL && pInverse != NULL){
        W = DSMatrixByMultiplyingMatrix(Cd, pInverse);
        if (pInverse != NULL)
            DSMatrixFree((DSMatrix *)pInverse);
    }
    if (W != NULL && B != NULL)
        Zeta = DSMatrixByMultiplyingMatrix(W, B);
    if (Zeta != NULL)
        DSMatrixAddByMatrix(Zeta, delta);
    else
        Zeta = DSMatrixCopy(delta);
    
    if (numberOfXi != 0) {
        Ai = DSSSystemAi(o_ssys);
        if (Ai != NULL && W != NULL)
            U = DSMatrixByMultiplyingMatrix(W, Ai);
        if (Ci != NULL && U != NULL)
            DSMatrixSubstractByMatrix(U, Ci);
        if (U == NULL){
            U = DSMatrixCopy(Ci);
            DSMatrixMultiplyByScalar(U, -1.0);
        }
        DSMatrixMultiplyByScalar(U, -1.0);
        if (Ai != NULL)
            DSMatrixFree(Ai);
    }
    
    
    ///////////////////////////
    // now rows of Zeta are adjusted to consider for numerical values of blowing dependent variables
    // in the condition matrices.
    
    w = DSMatrixCalloc(DSMatrixRows(I_pInvA), 1);
    for (i = 0; i < numberOfBlowing; i++){
        DSMatrixSetDoubleValue(w, Xd_b_indices[i], 0,
                               log10(DSVariablePoolValueForVariableWithName(Xd_b,   DSVariableName(DSVariablePoolAllVariables(Xd_b)[i]))));
    }
    
    if (I_pInvA != NULL && Cd != NULL && w!= NULL){
        MB_blow = DSMatrixByMultiplyingMatrix(I_pInvA, w);
        temp = DSMatrixByMultiplyingMatrix(Cd, MB_blow);
        DSMatrixAddByMatrix(Zeta, temp);
        if (MB_blow != NULL)
            DSMatrixFree(MB_blow);
        if (temp != NULL)
            DSMatrixFree(temp);
    }
    
    
    if (w != NULL)
        DSMatrixFree(w);
    if (I_pInvA != NULL)
        DSMatrixFree(I_pInvA);
    if (delta != NULL)
        DSMatrixFree(delta);
    if (W != NULL)
        DSMatrixFree(W);
    if (B != NULL)
        DSMatrixFree(B);
    if (InvAd != NULL)
        DSMatrixFree(InvAd);
    if (Xd_b_indices != NULL)
        DSSecureFree(Xd_b_indices);
    if (Cd != NULL)
        DSMatrixFree(Cd);
    if (Ci != NULL)
        DSMatrixFree(Ci);
    
    uCase->U = U;
    uCase->Zeta = Zeta;
    
bail:
    return;
    
}


extern DSMatrixArray * dsUnstableCaseLinearProblemAddEqualityConstraints(const DSMatrix *A_, const DSMatrix *B_,
                                                                         const DSUnstableCase *uCase)
{

    DSUInteger numberOfEqualities, numberOfAlgebraicConstraints;
    DSMatrix *Ae, *temp1, *temp2, *zeros1, *zeros2, *Be, *A, *B;
    DSMatrixArray *Equality = uCase->Equality, *array=NULL;
    
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
//    if (uCase->originalCase->caseNumber == 3000){
    
//        printf("Showing Equality Matrices for Case 3 \n");
//        DSMatrixPrint(DSMatrixArrayMatrix(uCase->Equality, 0));
//        DSMatrixPrint(DSMatrixArrayMatrix(uCase->Equality, 1));
//        DSMatrixPrint(DSMatrixArrayMatrix(uCase->Equality, 2));
        
//    }
    
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
    
//    printf("--------Reporting from within function dsUnstableCaseLinearProblemAddEqualityConstraints \n ");
//    printf("Matrix A is: \n");
//    DSMatrixPrint(A);
//    printf("Matrix B is: \n");
//    DSMatrixPrint(B);
    
    if (numberOfAlgebraicConstraints != 0)
        DSMatrixFree(temp1);
    DSMatrixFree(zeros2);
    DSMatrixFree(temp2);
    
bail:
    return array;
    
}

extern DSMatrixArray * dsUnstableCaseLinearProblemAddKnifeEdgeConditions(const DSMatrix *A_,
                                                                         const DSMatrix *B_,
                                                                         const DSUnstableCase *uCase)
{
    
    // this function adds knife edge conditions. Note that the knife edge is a function of independent variables
    // and numerical coefficients. We will ignore numerical coefficients for now.
    
    // the merging will be in the form [[0], gAi, [0]] corresponding to matrix of zeros for dependent variables Xd
    // (including ) algebraic constraints, gauss Ai and matrix of zeros for slack.
    
    DSMatrix *zeros1, *zeros2, *A, *B, *temp1, *temp2;
    DSMatrixArray *array = NULL, *Knife = uCase->Knife;
    DSUInteger numberOfKnifes = uCase->knifeEdge->numberOfVariables;
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

extern const DSMatrix *DSUnstableCaseGetSubSetPseudoInverse(const DSVariablePool *Xd_t, const DSVariablePool *Xd_e,
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
            DSMatrixSetDoubleValue(Ad_subset, indices[i], j, 0.0);
        }
    }
    
    pInverse = DSMatrixPseudoInverse(Ad_subset);
    if (indices != NULL)
        DSSecureFree(indices);
    if(Ad_subset != NULL)
        DSMatrixFree(Ad_subset);
    
    return pInverse;
    
}

extern void DSUnstableCaseGetEqualityAndKnife(const DSMatrix *pInverse,
                                              const DSMatrixArray *gaussArray,
                                              DSUnstableCase *uCase)
{
    
        //This function should populate fields uCase->Equality and uCase->knife.
        DSSSystem * o_ssys = uCase->originalCase->ssys, *collapsedSystem;
        DSMatrix * pInv_Ai, *pInv_Ai_all, * pInv_B, *pInv_B_all, *pInv_Ad, *pInv_Ad_all;
        DSMatrixArray *array = NULL, *knife = NULL;
        DSUInteger numberOfEqualityConstraints = uCase->Xd_e->numberOfVariables;
        DSUInteger numberOfKnifes = uCase->knifeEdge->numberOfVariables;
        DSUInteger *Xd_e_indices, *knife_row_indices;
        DSMatrix *gAi, *gB, *gAi_all, *gB_all;
        DSMatrix *Ad, *Ai, *B;
    
        // merge algebraic constraints if necessary and get indices of Xd_e
        if( DSVariablePoolNumberOfVariables(DSSSystemXd_a(o_ssys)) != 0){
            
            collapsedSystem = DSSSystemByRemovingAlgebraicConstraints(o_ssys);
            
            if (numberOfEqualityConstraints != 0){
                Xd_e_indices = DSVariablePoolIndicesOfSubPool(collapsedSystem->Xd, uCase->Xd_e);
                
                // populate Equality first.
                Ad = DSSSystemAd(collapsedSystem);
                pInv_Ad_all = DSMatrixByMultiplyingMatrix(pInverse, Ad );
                if (Ad != NULL)
                    DSMatrixFree(Ad);
                
                Ai = DSSSystemAi(collapsedSystem);
                pInv_Ai_all = DSMatrixByMultiplyingMatrix(pInverse, Ai);
                if (Ai != NULL)
                    DSMatrixFree(Ai);

                B = DSSSystemB(collapsedSystem);
                pInv_B_all  = DSMatrixByMultiplyingMatrix(pInverse, B);
                if (B != NULL)
                    DSMatrixFree(B);
                
                // get a subset of pInv_Ai_all to construct pInv_Ai. Loop over Xd_e;
                pInv_Ad = DSMatrixSubMatrixIncludingRows(pInv_Ad_all, numberOfEqualityConstraints, Xd_e_indices);
                pInv_Ai = DSMatrixSubMatrixIncludingRows(pInv_Ai_all, numberOfEqualityConstraints, Xd_e_indices);
                pInv_B  = DSMatrixSubMatrixIncludingRows(pInv_B_all, numberOfEqualityConstraints, Xd_e_indices);
                
                array = DSMatrixArrayAlloc();
                DSMatrixArrayAddMatrix(array, pInv_Ad);
                DSMatrixArrayAddMatrix(array, pInv_Ai);
                DSMatrixArrayAddMatrix(array, pInv_B);
                
                //free variables
                if (Xd_e_indices != NULL)
                    DSSecureFree(Xd_e_indices);
                if (pInv_Ad_all != NULL)
                    DSMatrixFree(pInv_Ad_all);
                if (pInv_Ai_all != NULL)
                    DSMatrixFree(pInv_Ai_all);
                if (pInv_B_all != NULL)
                    DSMatrixFree(pInv_B_all);
            }
            
            // populate Knife Edge
            // first, unpack gauss matrices
            //gAd_all = DSMatrixArrayMatrix(gaussArray, 0);
            gAi_all = DSMatrixArrayMatrix(gaussArray, 1);
            gB_all  = DSMatrixArrayMatrix(gaussArray, 2);
            
            knife_row_indices = DSVariablePoolIndicesOfSubPool(collapsedSystem->Xd, uCase->knifeEdge);
        
            gAi = DSMatrixSubMatrixIncludingRows(gAi_all, numberOfKnifes, knife_row_indices);
            gB  = DSMatrixSubMatrixIncludingRows(gB_all, numberOfKnifes, knife_row_indices);
            
            // add matrices to the array.
            knife = DSMatrixArrayAlloc();
            DSMatrixArrayAddMatrix(knife, gAi);
            DSMatrixArrayAddMatrix(knife, gB);
            
            // free variables.
            if (collapsedSystem != NULL)
                DSSSystemFree(collapsedSystem);
            if (knife_row_indices != NULL)
                DSSecureFree(knife_row_indices);
            
        } else { // else, the system does not have algebraic constraints.
            if (numberOfEqualityConstraints != 0){
                
                Xd_e_indices = DSVariablePoolIndicesOfSubPool(o_ssys->Xd, uCase->Xd_e);
                // populate Equality first.
                
                Ad = DSSSystemAd(o_ssys);
                pInv_Ad_all = DSMatrixByMultiplyingMatrix( pInverse, Ad);
                if (Ad != NULL)
                    DSMatrixFree(Ad);
                
                Ai = DSSSystemAi(o_ssys);
                pInv_Ai_all = DSMatrixByMultiplyingMatrix(pInverse, Ai);
                if (Ai != NULL)
                    DSMatrixFree(Ai);
                
                B = DSSSystemB(o_ssys);
                pInv_B_all  = DSMatrixByMultiplyingMatrix(pInverse, B);
                if (B != NULL)
                    DSMatrixFree(B);
                
                // get a subset of pInv_Ai_all to construct pInv_Ai. Loop over Xd_e;
                pInv_Ad = DSMatrixSubMatrixIncludingRows(pInv_Ad_all, numberOfEqualityConstraints, Xd_e_indices);
                pInv_Ai = DSMatrixSubMatrixIncludingRows(pInv_Ai_all, numberOfEqualityConstraints, Xd_e_indices);
                pInv_B  = DSMatrixSubMatrixIncludingRows(pInv_B_all, numberOfEqualityConstraints, Xd_e_indices);
                
                array = DSMatrixArrayAlloc();
                DSMatrixArrayAddMatrix(array, pInv_Ad);
                DSMatrixArrayAddMatrix(array, pInv_Ai);
                DSMatrixArrayAddMatrix(array, pInv_B);
                
                //free variables
                if (Xd_e_indices != NULL)
                    DSSecureFree(Xd_e_indices);
                if (pInv_Ad_all != NULL)
                    DSMatrixFree(pInv_Ad_all);
                if (pInv_Ai_all != NULL)
                    DSMatrixFree(pInv_Ai_all);
                if (pInv_B_all != NULL)
                    DSMatrixFree(pInv_B_all);
            }
            
            // populate Knife Edge
            // first, unpack gauss matrices
            //gAd_all = DSMatrixArrayMatrix(gaussArray, 0);
            gAi_all = DSMatrixArrayMatrix(gaussArray, 1);
            gB_all  = DSMatrixArrayMatrix(gaussArray, 2);
            
    //        knife_row_indices = DSVariablePoolIndicesOfSubPool(o_ssys->Xd, uCase->Xd_b);
            knife_row_indices = DSVariablePoolIndicesOfSubPool(o_ssys->Xd, uCase->knifeEdge);
            gAi = DSMatrixSubMatrixIncludingRows(gAi_all, numberOfKnifes, knife_row_indices);
            gB  = DSMatrixSubMatrixIncludingRows(gB_all, numberOfKnifes, knife_row_indices);
            
            // add matrices to the array.
            knife = DSMatrixArrayAlloc();
            DSMatrixArrayAddMatrix(knife, gAi);
            DSMatrixArrayAddMatrix(knife, gB);
            DSSecureFree(knife_row_indices);
            
        }
    
        uCase->Equality = array;
        uCase->Knife = knife;
    
}


extern const DSCase ** DSUnstableCaseListAllSubcases(const DSCase * aCase, const DSDesignSpace *ds)
{
    
    const DSCase ** validSubCases = NULL;
    DSUnstableCase * uCase = NULL;
    DSUInteger i, numberValid = 0;
    if (ds == NULL) {
        DSError(M_DS_DESIGN_SPACE_NULL, A_DS_ERROR);
        goto bail;
    }
    
    uCase = DSDictionaryValueForName(ds->unstableCases, aCase->caseIdentifier);
    
    if (uCase == NULL){
        DSError(M_DS_NULL ": Unstable case not present in dictionary!", A_DS_ERROR);
        goto bail;
    }
    numberValid = uCase->ValidCases->count;
    if (numberValid == 0)
        goto bail;
    validSubCases = DSSecureCalloc(sizeof(DSCase *), numberValid);
    
    for (i = 0; i < numberValid; i++) {
        validSubCases[i] = DSCaseCopy(DSDictionaryValueForName(uCase->ValidCases, uCase->ValidCases->names[i]));
    }
    
bail:
    return validSubCases;
    
}

extern DSUInteger DSUnstableCaseSubcasesCount(const DSCase * aCase, const DSDesignSpace *ds)
{
    DSUInteger count = 0;
    DSUnstableCase *uCase = DSDictionaryValueForName(ds->unstableCases, aCase->caseIdentifier);
    DSDictionary *aDictionary = uCase->ValidCases;
    
    if (aDictionary == NULL) {
        DSError(M_DS_DICTIONARY_NULL, A_DS_ERROR);
        goto bail;
    }
    count = aDictionary->count;
bail:
    return  count;
}
