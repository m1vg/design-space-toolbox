/**
 * \file DSCase.c
 * \brief Implementation file with functions for dealing with cases in design space.
 *
 * \details 
 *
 * Copyright (C) 2011-2014 Jason Lomnitz.\n\n
 *
 * This file is part of the Design Space Toolbox V2 (C Library).
 *
 * The Design Space Toolbox V2 is free software: you can redistribute it and/or 
 * modify it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The Design Space Toolbox V2 is distributed in the hope that it will be 
 * useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with the Design Space Toolbox. If not, see 
 * <http://www.gnu.org/licenses/>.
 *
 * \author Jason Lomnitz.
 * \date 2011
 */

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

#if defined (__APPLE__) && defined (__MACH__)
#pragma mark - DSCase Global behavior
#endif

/* A + B > 0 => B > -A:S if S > 0*/

static char endian;

extern void DSCaseSetEndianness(char endianness)
{
        if (endianness != DS_CASE_NUMBER_BIG_ENDIAN && endianness != DS_CASE_NUMBER_SMALL_ENDIAN) {
                DSError(M_DS_WRONG ": Endianness must be big or small", A_DS_ERROR);
                goto bail;
        }
        endian = endianness;
bail:
        return;
}

extern char DSCaseEndianness(void)
{
        return endian;
}
#if defined (__APPLE__) && defined (__MACH__)
#pragma mark - Allocation, deallocation and initialization
#endif

static DSCase * DSCaseAlloc(void)
{
        DSCase * aCase = NULL;
        aCase = DSSecureCalloc(sizeof(DSCase), 1);
        aCase->freeVariables = false;
        return aCase;
}

static DSMassBalanceData * dsCaseDominantCaseCopyMassBalanceData(const DSCase *aCase){
    
    DSMassBalanceData * dominant_mass_balances = NULL;
    
    if (DSCaseDominantMassBalances(aCase) == NULL)
        goto bail;
    
    DSUInteger n, *signature=NULL, i, length = 1000, *signature_source = NULL;
    char **fin = NULL, **fout = NULL, **fin_source = NULL, **fout_source=NULL, *aux = NULL;
    char **fin_Xd = NULL, **fout_Xd = NULL, **fin_Xd_source = NULL, **fout_Xd_source=NULL;
    
    dominant_mass_balances = DSSecureCalloc(1, sizeof(DSMassBalanceData));
    fin_source = DSCaseDominantMassBalances(aCase)->fin;
    fout_source = DSCaseDominantMassBalances(aCase)->fout;
    fin_Xd_source = DSCaseDominantMassBalances(aCase)->fin_Xd;
    fout_Xd_source = DSCaseDominantMassBalances(aCase)->fout_Xd;
    signature_source = DSCaseDominantMassBalances(aCase)->signature;
    
    if (fin_source!=NULL)
        fin = DSSecureCalloc(DSCaseDominantMassBalances(aCase)->n, sizeof(char *));
    if (fout_source != NULL)
        fout = DSSecureCalloc(DSCaseDominantMassBalances(aCase)->n, sizeof(char *));
    if (fin_Xd_source!=NULL)
        fin_Xd = DSSecureCalloc(DSCaseDominantMassBalances(aCase)->n, sizeof(char *));
    if (fout_Xd_source != NULL)
        fout_Xd = DSSecureCalloc(DSCaseDominantMassBalances(aCase)->n, sizeof(char *));
    if (signature_source != NULL)
        signature = DSSecureCalloc(2*DSCaseDominantMassBalances(aCase)->n, sizeof(DSUInteger));
    
    for (i=0; i<2*DSCaseDominantMassBalances(aCase)->n; i++){
        
        if (signature_source != NULL)
            signature[i] = signature_source[i];

        if (i<DSCaseDominantMassBalances(aCase)->n){
            
            if (fin_source != NULL){
                aux = DSSecureCalloc(1, sizeof(char)*length);
                strcpy(aux, fin_source[i]);
                fin[i] = aux;
            }
            if (fout_source != NULL){
                aux = DSSecureCalloc(1, sizeof(char)*length);
                strcpy(aux, fout_source[i]);
                fout[i] = aux;
            }
            if (fin_Xd_source != NULL){
                aux = DSSecureCalloc(1, sizeof(char)*length);
                strcpy(aux, fin_Xd_source[i]);
                fin_Xd[i] = aux;
            }
            if (fout_Xd_source != NULL){
                aux = DSSecureCalloc(1, sizeof(char)*length);
                strcpy(aux, fout_Xd_source[i]);
                fout_Xd[i] = aux;
            }
            
        }
    }
    
    dominant_mass_balances->fin = fin;
    dominant_mass_balances->fout = fout;
    dominant_mass_balances->fin_Xd = fin_Xd;
    dominant_mass_balances->fout_Xd = fout_Xd;
    dominant_mass_balances->signature = signature;
    dominant_mass_balances->n = DSCaseDominantMassBalances(aCase)->n;
    
bail:
    return dominant_mass_balances;
}

extern DSCase * DSCaseCopy(const DSCase * aCase)
{
        DSCase * newCase = NULL;
        DSUInteger i, numberOfEquations;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        newCase = DSCaseAlloc();
        DSCaseSSys(newCase) = DSSSystemCopy(DSCaseSSystem(aCase));
        DSCaseNum(newCase) = DSCaseNum(aCase);
        numberOfEquations = DSCaseNumberOfEquations(aCase) + aCase->numberInheritedConservations;
        if (DSCaseSig(aCase) != NULL) {
                DSCaseSig(newCase) = DSSecureCalloc(sizeof(DSUInteger), numberOfEquations*2);
                for (i = 0; i < numberOfEquations*2; i++) {
                        DSCaseSig(newCase)[i] = DSCaseSig(aCase)[i];
                }
        }
        if (DSCase3Sig(aCase) != NULL) {
            DSCase3Sig(newCase) = DSSecureCalloc(sizeof(DSUInteger), numberOfEquations*3);
            for (i = 0; i < numberOfEquations*3; i++) {
                DSCase3Sig(newCase)[i] = DSCase3Sig(aCase)[i];
            }
        }
        if (DSCaseSigCons(aCase) != NULL){
            numberOfEquations = DSCaseNumberOfConservations(aCase) + DSCaseNumberOfEquations(aCase);
            DSCaseSigCons(newCase) = DSSecureCalloc(sizeof(DSUInteger), numberOfEquations*2);
            for (i = 0; i< numberOfEquations*2; i++)
                DSCaseSigCons(newCase)[i] = DSCaseSigCons(aCase)[i];
        }
        if (DSCaseCd(aCase) != NULL)
                DSCaseCd(newCase) = DSMatrixCopy(DSCaseCd(aCase));
        if (DSCaseCi(aCase) != NULL)
                DSCaseCi(newCase) = DSMatrixCopy(DSCaseCi(aCase));
        if (DSCaseZeta(aCase) != NULL)
                DSCaseZeta(newCase) = DSMatrixCopy(DSCaseZeta(aCase));
        if (DSCaseDelta(aCase) != NULL)
                DSCaseDelta(newCase) = DSMatrixCopy(DSCaseDelta(aCase));
        if (DSCaseU(aCase) != NULL)
                DSCaseU(newCase) = DSMatrixCopy(DSCaseU(aCase));
        if (DSCaseDominantMassBalances(aCase) != NULL)
            DSCaseDominantMassBalances(newCase) = dsCaseDominantCaseCopyMassBalanceData(aCase);
        newCase->Xd = DSSSystemXd(DSCaseSSys(newCase));
        newCase->Xi = DSSSystemXi(DSCaseSSys(newCase));
        newCase->Xd_a =DSSSystemXd_a(DSCaseSSys(newCase));
        newCase->numberInheritedConservations = aCase->numberInheritedConservations;
        DSCaseId(newCase) = strdup(DSCaseId(aCase));
bail:
        return newCase;
}

extern void DSCaseFree(DSCase * aCase)
{
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (DSCaseSSys(aCase) != NULL) {
                DSSSystemFree(DSCaseSSys(aCase));
        }
        if (aCase->freeVariables == true){
                DSVariablePoolSetReadWriteAdd((DSVariablePool *)aCase->Xi);
                DSVariablePoolSetReadWriteAdd((DSVariablePool *)aCase->Xd);
                DSVariablePoolSetReadWriteAdd((DSVariablePool *)aCase->Xd_a);
                if (aCase->Xi != NULL)
                    DSVariablePoolFree((DSVariablePool *)aCase->Xi);
                if (aCase->Xd != NULL)
                    DSVariablePoolFree((DSVariablePool *)aCase->Xd);
                if (aCase->Xd_a != NULL)
                    DSVariablePoolFree((DSVariablePool *)aCase->Xd_a);
        }
        if (DSCaseDominantMassBalances(aCase) != NULL)
            DSCaseMassBalanceDataFree(aCase);
        if (DSCaseSig(aCase) != NULL)
                DSSecureFree(DSCaseSig(aCase));
        if (DSCase3Sig(aCase) != NULL)
                DSSecureFree(DSCase3Sig(aCase));
        if (DSCaseSigCons(aCase) != NULL)
            DSSecureFree(DSCaseSigCons(aCase));
        if (DSCaseCd(aCase) != NULL)
                DSMatrixFree(DSCaseCd(aCase));
        if (DSCaseCi(aCase) != NULL)
                DSMatrixFree(DSCaseCi(aCase));
        if (DSCaseZeta(aCase) != NULL)
                DSMatrixFree(DSCaseZeta(aCase));
        if (DSCaseDelta(aCase) != NULL)
                DSMatrixFree(DSCaseDelta(aCase));
        if (DSCaseU(aCase) != NULL)
                DSMatrixFree(DSCaseU(aCase));
        if (DSCaseId(aCase) != NULL)
                DSSecureFree(DSCaseId(aCase));
        if ((DSCase *)aCase != NULL)
            DSSecureFree(aCase);
bail:
        return;
}

extern void DSCaseVolumeFree(DSCaseVolume *caseVolume){
    
    if (caseVolume == NULL) {
            DSError(M_DS_CASE_NULL, A_DS_ERROR);
            goto bail;
    }
    
    if (caseVolume->vertices != NULL)
        DSMatrixFree(caseVolume->vertices);
    DSSecureFree(caseVolume);
    
bail:
    return;
}

void DSCaseMassBalanceDataFree(DSCase *aCase){
    if (aCase == NULL || DSCaseDominantMassBalances(aCase) == NULL)
        goto bail;
    DSMassBalanceData *massbalances;
    DSUInteger i;
    massbalances = DSCaseDominantMassBalances(aCase);
    
    
    for (i=0; i<massbalances->n; i++){
        if (massbalances->fin[i] != NULL)
            DSSecureFree(massbalances->fin[i]);
        if (massbalances->fout[i] != NULL)
            DSSecureFree(massbalances->fout[i]);
        if (massbalances->fin_Xd[i] != NULL)
            DSSecureFree(massbalances->fin_Xd[i]);
        if (massbalances->fout_Xd[i] != NULL)
            DSSecureFree(massbalances->fout_Xd[i]);
    }
    
    if (massbalances->signature != NULL)
        DSSecureFree(massbalances->signature);
    if (massbalances->fin != NULL)
        DSSecureFree(massbalances->fin);
    if (massbalances->fout != NULL)
        DSSecureFree(massbalances->fout);
    if (massbalances->fin_Xd != NULL)
        DSSecureFree(massbalances->fin_Xd);
    if (massbalances->fout_Xd != NULL)
        DSSecureFree(massbalances->fout_Xd);
    
bail:
    return;
}


#if defined (__APPLE__) && defined (__MACH__)
#pragma mark - Factory functions
#endif

extern void DSCaseRemoveRedundantBoundaries(DSCase *aCase)
{
        DSMatrix *temp1, *temp2;
        temp1 = DSMatrixAppendMatrices(DSCaseU(aCase), DSCaseZeta(aCase), true);
        temp2 = DSMatrixWithUniqueRows(temp1);
        DSMatrixFree(temp1);
        if (temp2 == NULL)
                goto bail;
        DSMatrixFree(DSCaseU(aCase));
        DSMatrixFree(DSCaseZeta(aCase));
        DSCaseU(aCase) = DSMatrixSubMatrixExcludingColumnList(temp2, 1, DSMatrixColumns(temp2)-1);
        DSCaseZeta(aCase) = DSMatrixSubMatrixIncludingColumnList(temp2, 1, DSMatrixColumns(temp2)-1);
        DSMatrixFree(temp2);
bail:
        return;
}

extern void DSCaseAdjustStoichiometryOfCodominantCase(DSCase *aCase){
    
    if (aCase == NULL) {
            DSError(M_DS_CASE_NULL, A_DS_ERROR);
            goto bail;
    }
    
    if (DSSSystemAdjustCodominantStoichiometry(DSCaseSSystem(aCase)) == false)
        goto bail;
    
    DSSSystemAdjustStoichiometryOfCodominantCase(DSCaseSSystem(aCase));
    
bail:
    return;

}

static char ** dsCaseCodominantCaseCreateMassBalanceInequalities(const DSUInteger *signature_i,
                                                                 const DSUInteger *mass_balance_signature,
                                                                 const char ** fin,
                                                                 const char **fout,
                                                                 DSUInteger n_mass_balances,
                                                                 DSUInteger n_inequalities){
    

    char ** inequalities = NULL, *non_dominant_string = NULL, *dominant_string = NULL, *aux=NULL ;
    DSUInteger i, j,  number_of_terms, counter=0;
    DSExpression *fin_i = NULL, *fout_i = NULL, *fin_i_rhs = NULL, *fout_i_rhs = NULL;
    
    if (signature_i == NULL || fin == NULL || fout == NULL || n_inequalities == 0)
        goto bail;
    
    inequalities = DSSecureCalloc(n_inequalities, sizeof(char *));
    
    for (i=0; i<n_mass_balances*2; i++ ){
            if (mass_balance_signature[i] <= 1)
                    continue;
            if (i%2 == 0){
                    fin_i = DSExpressionByParsingString(fin[i/2]);
                    fin_i_rhs = DSExpressionEquationRHSExpression(fin_i);
                    dominant_string = DSExpressionBranchAtIndexAsString(fin_i_rhs, signature_i[i]);
                    if (fin_i != NULL){
                        DSExpressionFree(fin_i);
                        fin_i = NULL;
                    }
            }else{
                    fout_i = DSExpressionByParsingString(fout[i/2]);
                    fout_i_rhs = DSExpressionEquationRHSExpression(fout_i);
                    dominant_string = DSExpressionBranchAtIndexAsString(fout_i_rhs, signature_i[i]);
                    if (fout_i != NULL){
                        DSExpressionFree(fout_i);
                        fout_i = NULL;
                    }
            }
            for (j=0; j<mass_balance_signature[i]; j++){
                    aux = DSSecureMalloc(sizeof(char)*1000);
                    if (j+1 == signature_i[i])
                        continue;
                    
                    if (i%2 == 0){
                        non_dominant_string = DSExpressionBranchAtIndexAsString(fin_i_rhs, j+1);
                    }else{
                        non_dominant_string = DSExpressionBranchAtIndexAsString(fout_i_rhs, j+1);
                    }
                    sprintf(aux, "%s > %s", dominant_string, non_dominant_string);
                    inequalities[counter] = aux;
                    counter += 1;
            }
        
            if (fout_i_rhs != NULL){
                DSExpressionFree(fout_i_rhs);
                fout_i_rhs = NULL;
            }
            if (fin_i_rhs != NULL){
                DSExpressionFree(fin_i_rhs);
                fin_i_rhs = NULL;
            }
            if (dominant_string!=NULL){
                DSSecureFree(dominant_string);
                dominant_string = NULL;
            }
            if (non_dominant_string!=NULL){
                DSSecureFree(non_dominant_string);
                non_dominant_string = NULL;
            }

    }

    
bail:
    return inequalities;
}

static DSExpression * dsExpressionProcessExpressionByReplacingDependentVariables(DSExpression * aux_rhs,
                                                                                 DSCase * aCase_i){
    
    DSExpression * processed_expression = NULL, *target = NULL, *substitute = NULL, * aux = NULL;
    DSExpression ** solution = NULL;
    DSUInteger i;
    bool contain_dependent_variables = true;
    DSVariablePool *Xi=NULL, *Xd=NULL;
    char **names, *aux_string=NULL;
    
    if (aux_rhs == NULL)
        goto bail;
    if (aCase_i == NULL){
        DSError(M_DS_CASE_NULL, A_DS_ERROR);
        goto bail;
    }
        
    solution = DSCaseSolution(aCase_i);
    processed_expression = DSExpressionCopy(aux_rhs);
    Xd = DSCaseXd(aCase_i);
    
    while (contain_dependent_variables == true){
            for (i=0; i<DSSSystemNumberOfEquations(DSCaseSSystem(aCase_i)); i++){
                
//                    printf("expression before round of substitutions of dependent variables [%u]: %s \n", i, DSExpressionAsString(processed_expression));
                
                    target = DSExpressionEquationLHSExpression(solution[i]);
                    substitute = DSExpressionEquationRHSExpression(solution[i]);
                    aux = DSExpressionByReplacingSubExpression(processed_expression, target, substitute);
                    aux_string = DSExpressionAsString(aux);
                
                    if (processed_expression != NULL)
                        DSExpressionFree(processed_expression);
                    if (target != NULL)
                        DSExpressionFree(target);
                    if (substitute != NULL)
                        DSExpressionFree(substitute);
                    if (aux != NULL)
                        DSExpressionFree(aux);
                    processed_expression = DSExpressionByParsingString(aux_string);
                    if (aux_string != NULL)
                        DSSecureFree(aux_string);
                
//                printf("expression after round of substitutions of dependent variables [%u]: %s \n", i, DSExpressionAsString(processed_expression));
                
            }
            contain_dependent_variables = false;
//            printf("Expression before being simplified %s \n", DSExpressionAsString(processed_expression));
            aux = DSExpressionSimplifyExpressionAsExpression(processed_expression);
            if (processed_expression != NULL)
                DSExpressionFree(processed_expression);
//            printf("Expression after being simplified %s \n", DSExpressionAsString(aux));
            processed_expression = aux;
            Xi = DSExpressionVariablesInExpression(processed_expression);
            names = DSVariablePoolAllVariableNames(Xi);
            for (i=0; i<DSVariablePoolNumberOfVariables(Xi); i++){
//                    printf("checking if %s is a dependent variable! \n", names[i]);
                if (DSVariablePoolHasVariableWithName(Xd, names[i]) == true){
                        contain_dependent_variables = true;
//                        printf("yes, it is! \n");
                }
            }
            
            if (Xi != NULL){
                DSVariablePoolFree(Xi);
                Xi = NULL;
            }
            if (names != NULL){
                DSSecureFree(names);
                names = NULL;
            }
    }
    

    
    
bail:
    return processed_expression;
}

static bool dsCaseCodominatCaseMassBalancesAreValid(DSCase * aCase_i,
                                                    const DSUInteger * signature_i,
                                                    const DSUInteger * mass_balance_signature,
                                                    DSUInteger n_mass_balances,
                                                    const char **fin,
                                                    const char **fout){
    if (aCase_i == NULL){
        DSError(M_DS_CASE_NULL, A_DS_ERROR);
        goto bail;
    }
    
    if (signature_i == NULL || fin == NULL || fout == NULL)
        goto bail;
    
    bool massBalanceisValid = false;
    DSExpression *aux = NULL, *aux_rhs = NULL, * aux_rhs_dom = NULL, * ratio_simplified = NULL;
    DSExpression *fin_dom = NULL, *fout_dom = NULL;
    DSUInteger i;
    char * ratio = NULL, **fin_dom_char=NULL, **fout_dom_char=NULL, **fin_dom_Xd=NULL, **fout_dom_Xd=NULL;
    DSVariablePool * Xi = NULL;
    DSMassBalanceData *dominant_exchange_fluxes = NULL;
    
    dominant_exchange_fluxes = DSSecureCalloc(sizeof(DSMassBalanceData), 1);
    dominant_exchange_fluxes->n = n_mass_balances;
    fin_dom_char = DSSecureCalloc(n_mass_balances, sizeof(char *));
    fout_dom_char = DSSecureCalloc(n_mass_balances, sizeof(char *));
    fin_dom_Xd = DSSecureCalloc(n_mass_balances, sizeof(char *));
    fout_dom_Xd = DSSecureCalloc(n_mass_balances, sizeof(char *));
    
    
    // We first need to identify the dominant fin and dominant fout from the signatures and the fin and fout expressions
//    printf("About to test %u metabolic blocks: \n", n_mass_balances);
    for (i=0; i<2*n_mass_balances; i++){
        
            if (i%2 == 0){
//                printf("The inputs are: %s \n", fin[i/2]);
                aux = DSExpressionByParsingString(fin[i/2]);
                aux_rhs = DSExpressionEquationRHSExpression(aux);
                if (mass_balance_signature[i] == 1)
                    aux_rhs_dom = aux_rhs;
                else
                    aux_rhs_dom = DSExpressionBranchAtIndex(aux_rhs, signature_i[i]);
//                printf("input without replacing dependent variables: %s \n", DSExpressionAsString(aux_rhs_dom));
                fin_dom = dsExpressionProcessExpressionByReplacingDependentVariables(aux_rhs_dom, aCase_i);
                fin_dom_Xd[i/2] = DSExpressionAsString(aux_rhs_dom);
            
            }else{
                aux = DSExpressionByParsingString(fout[i/2]);
                aux_rhs = DSExpressionEquationRHSExpression(aux);
                if (mass_balance_signature[i] == 1)
                    aux_rhs_dom = aux_rhs;
                else
                    aux_rhs_dom = DSExpressionBranchAtIndex(aux_rhs, signature_i[i]);
//                printf("output without replacing dependent variables: %s \n", DSExpressionAsString(aux_rhs_dom));
                fout_dom = dsExpressionProcessExpressionByReplacingDependentVariables(aux_rhs_dom, aCase_i);
                ratio = DSSecureMalloc(sizeof(char)*1000);
                sprintf(ratio, "%s/(%s)", DSExpressionAsString(fin_dom) , DSExpressionAsString(fout_dom));
//                printf("input = %s; output = %s \n", DSExpressionAsString(fin_dom), DSExpressionAsString(fout_dom));
                ratio_simplified = DSExpressionByParsingStringSimplifyExpressionAsExpression(ratio);
//                printf("The simplified ratio is %s \n", DSExpressionAsString(ratio_simplified));
                Xi = DSExpressionVariablesInExpression(ratio_simplified);
                if (DSVariablePoolNumberOfVariables(Xi) == 0)
                        massBalanceisValid = true;
                else
                        massBalanceisValid = false;
                
                fin_dom_char[i/2] = DSExpressionAsString(fin_dom);
                fout_dom_char[i/2] = DSExpressionAsString(fout_dom);
                fout_dom_Xd[i/2] = DSExpressionAsString(aux_rhs_dom);
            
                if (ratio != NULL)
                    DSSecureFree(ratio);
                if (ratio_simplified != NULL)
                    DSExpressionFree(ratio_simplified);
                if (Xi != NULL)
                    DSVariablePoolFree(Xi);
                if (fin_dom != NULL)
                    DSExpressionFree(fin_dom);
                if (fout_dom != NULL)
                    DSExpressionFree(fout_dom);
            }
            
            if (aux !=NULL)
                DSExpressionFree(aux);
            if (aux_rhs != NULL)
                DSExpressionFree(aux_rhs);
    }
    
    dominant_exchange_fluxes->fin = fin_dom_char;
    dominant_exchange_fluxes->fout = fout_dom_char;
    dominant_exchange_fluxes->fin_Xd = fin_dom_Xd;
    dominant_exchange_fluxes->fout_Xd = fout_dom_Xd;
    DSCaseDominantMassBalances(aCase_i) = dominant_exchange_fluxes;
    
bail:
    return massBalanceisValid;
}

static bool dominant_input_output_are_connected(const DSDesignSpace *ds,
                                                DSCase *newCase,
                                                const char *fin,
                                                const char *fout,
                                                const char **rxns,
                                                DSUInteger * fin_index,
                                                DSUInteger * fout_index,
                                                const DSMatrix *S){
    
    
    DSUInteger i;
    bool check_internal = false;

    
    // we first find the indices of fin and fout in the columns of S.
    for (i=0; i<DSMatrixColumns(S); i++){
        if(DSStringIsEqualToStringIgnoringOrder(fin, rxns[i]) == true)
            *fin_index = i;
        if(DSStringIsEqualToStringIgnoringOrder(fout, rxns[i]) == true)
            *fout_index = i;
    }
    
    
    // We then check the stoichiometric matrix to see if internal mass balances should be considered for this metabolic block
    for (i=0; i<DSMatrixRows(S); i++){
        if (DSMatrixDoubleValue(S, i, *fin_index) != 0.0){
            if (DSMatrixDoubleValue(S, i, *fout_index) == 0.0){
                check_internal = true;
                break;
            }
        }
    }
    
    return check_internal;
}

static DSMatrix * reduce_S_rxns_Xd(DSMatrix *S,
                                   DSUInteger fin_index,
                                   DSUInteger fout_index,
                                   char ** rxns,
                                   char **rxns_reduced,
                                   DSVariablePool *Xd,
                                   DSVariablePool *Xd_reduced){
    
    
    DSMatrix *S_reduced = NULL;
    DSUInteger i, n, numberOfRows = 0, numberOfColumns = 2, *rows = DSSecureMalloc(sizeof(DSUInteger)*1);
    DSUInteger *columns = DSSecureMalloc(sizeof(DSUInteger)*2);
    
    if (fin_index > DSMatrixColumns(S))
        goto bail;
    
    if (fout_index > DSMatrixColumns(S))
        goto bail;
    
    columns[0] = fin_index;
    columns[1] = fout_index;
    
    // we first delete the row(s) containing the pool fed by fin_index.
    for (i=0; i<DSMatrixRows(S); i++){
        // capture indices of rows to delete
        if (DSMatrixDoubleValue(S, i, fin_index) != 0){

            if (numberOfRows > 0)
                DSSecureRealloc(rows, sizeof(DSUInteger)*(numberOfRows + 1));
            rows[numberOfRows] = i;
            numberOfRows++;
                
        // capture dependent elements that will survive.
        } else{
            DSVariablePoolAddVariableWithName(Xd_reduced,
                                              DSVariableName(DSVariablePoolAllVariables(Xd)[i]));
        }
    }
    
    printf("@reduce_S_rxns_Xd. The elements of rows are: \n");
    for (i=0; i<numberOfRows; i++)
        printf("i=%u: row[i] = %u \n", i, rows[i]);
    
    S_reduced = DSMatrixSubMatrixExcludingRowsAndColumns(S, numberOfRows, numberOfColumns, rows, columns);
    
    
    for (i=0, n=0; i<DSMatrixColumns(S); i++){
        if (i == fin_index || i == fout_index)
            continue;
        rxns_reduced[n] = strdup(rxns[i]);
        n++;
    }
    
bail:
    return S_reduced;
    
}

static DSUInteger * determine_mass_balances(const DSMatrix * S_reduced, const char ** rxns_reduced,
                                            const DSVariablePool *Xd_reduced, const char *fin, const char *fout){
    
    DSUInteger * mass_balance_signature = NULL;
    
    
bail:
    return mass_balance_signature;
}

static DSUInteger * determine_dominant_input_output(const DSCase * newCase, const DSUInteger * mass_balance_signature){
    
    DSUInteger * dominant_signature = NULL;
    
    
bail:
    return dominant_signature;
}

static void adjust_internal_stoichiometry(DSCase * newCase, DSMatrix * alpha, DSUInteger * dominant_signature,
                                          const char * fin, const char * fout){
    
    
    return;
}

static void DSCaseCodominantCaseAdjustInternalStoichiometry_recursive(DSMatrix * alpha,
                                                                      const DSDesignSpace *ds,
                                                                      DSCase *newCase,
                                                                      const char *fin,
                                                                      const char *fout,
                                                                      char **rxns,
                                                                      DSVariablePool *Xd_i,
                                                                      DSMatrix *S_i,
                                                                      double factor){
    
    bool check_internal = false;
    DSUInteger fin_index, fout_index, *dominant_signature = NULL, *mass_balance_signature = NULL;
    DSMatrix *S_reduced = NULL;
    char **allvariables, *fin_i = NULL, *fout_i = NULL;
    
    
    // 0. Check number of rows in S. If 1 or less, return.
    if (S_i == NULL)
        goto bail;
    if (DSMatrixRows(S_i) <= 1)
        goto bail;
    
    // 1. Adjust internal stoichiometries if dominant input and output do not share a common pool.
    check_internal = dominant_input_output_are_connected(ds, newCase, fin, fout, rxns, &fin_index, &fout_index, S_i);
    if (check_internal == false)
        goto bail;
    
    char **rxns_reduced = DSSecureMalloc(sizeof(char *)* (DSMatrixColumns(S_i) - 2));
    DSVariablePool *Xd_reduced = DSVariablePoolAlloc();
      
    //2. Adjust the stoichiometric matrix: remove pool(s) (rows) affected by fin as well as columns for fin_index and fout_index. rxns_reduced is obtained from rxns. We will also probably have to modify the pools.
    S_reduced = reduce_S_rxns_Xd(S_i, fin_index, fout_index, rxns,  rxns_reduced, Xd_i, Xd_reduced);
    
    fin_i = DSSecureMalloc(sizeof(char )* 1000);
    fout_i = DSSecureMalloc(sizeof(char )* 1000);
    //3. Determine metabolic blocks and their corresponding mass balances.
    mass_balance_signature = determine_mass_balances(S_reduced, rxns_reduced, Xd_reduced, fin_i, fout_i);
    
    if (mass_balance_signature == NULL)
        goto bail;
    
    //4. determine dominant input and output fluxes using linear programming
    dominant_signature = determine_dominant_input_output(newCase, mass_balance_signature);
    
    if (dominant_signature == NULL)
        goto bail;
    
    //5. adjust s-system equations as well as boundaries for both dominant input and dominant output
    adjust_internal_stoichiometry(newCase, alpha, dominant_signature, fin_i, fout_i);
    
    //6. call this same function.
    if (DSMatrixRows(S_reduced) > 1 && (factor - 1.0) > 0.0)
        DSCaseCodominantCaseAdjustInternalStoichiometry_recursive(alpha, ds, newCase, fin_i, fout_i, rxns_reduced,
                                                                  Xd_reduced, S_reduced, factor - 1.0);
    
    
bail:
    if (S_reduced != NULL)
        DSMatrixFree(S_reduced);
    if (rxns_reduced != NULL)
        DSSecureFree(rxns_reduced);
    if (Xd_reduced != NULL)
        DSVariablePoolFree(Xd_reduced);
    if (mass_balance_signature != NULL)
        DSSecureFree(mass_balance_signature);
    if (dominant_signature != NULL)
        DSSecureFree(dominant_signature);
    if (fin_i != NULL)
        DSSecureFree(fin_i);
    if(fout_i != NULL)
        DSSecureFree(fout_i);

    return;
}


static DSMatrix * s_for_metabolic_block(DSDesignSpace *ds, DSMatrix *S, DSVariablePool *Xd_i, DSUInteger block_i){
    
    DSMatrix * S_i =NULL;
    DSUInteger numberOfRows = 0, *rows = DSSecureCalloc(sizeof(DSUInteger), 1), i, value;
    DSVariablePool *X_metabolic_blocks = ds->massBalances->metabolicBlocks;
    char *name;
    
    
    // We will only include rows for certain metabolic species.
    for (i=0; i<DSVariablePoolNumberOfVariables(X_metabolic_blocks); i++){
        name = DSVariableName(DSVariablePoolAllVariables(X_metabolic_blocks)[i]);
        value = (DSUInteger )DSVariablePoolValueForVariableWithName(X_metabolic_blocks, name);
        if (value != block_i)
            continue;
        DSVariablePoolAddVariableWithName(Xd_i, name);
        
        if (numberOfRows > 0)
            DSSecureRealloc(rows, sizeof(DSUInteger)*(numberOfRows+1));
        rows[numberOfRows] = i;
        numberOfRows++;
        
    }
    
    printf("@s_for_metabolic_block. The elements of rows (%u) are: \n", numberOfRows);
    for (i=0; i<numberOfRows; i++)
        printf("i=%u: row[i] = %u \n", i, rows[i]);
    
    printf("The metabolic pool Xd_i is: \n");
    DSVariablePoolPrint(X_metabolic_blocks);
    
    S_i = DSMatrixSubMatrixIncludingRows(S, numberOfRows, rows);
    
bail:
    return S_i;
}

static void DSCaseCodominantCaseAdjustInternalStoichiometryUsingDominantMassBalance(DSMatrix * alpha,
                                                                                    DSCase *newCase,
                                                                                    const DSUInteger * signature_i,
                                                                                    DSUInteger *numberZeroBoundariesPerBlock,
                                                                                    const DSDesignSpace *ds,
                                                                                    const double *factors){
    
    
    
    // The elements that I need here are the stoichiometric matrix, the definition of metabolic blocks, dominant intput and output fluxes, as well as the vector of reaction fluxes.
    
    if (alpha == NULL || newCase == NULL || ds == NULL){
        DSError(M_DS_CASE_NULL, A_DS_ERROR);
        goto bail;
    }
    
//    DSMatrix *S = ds->massBalances->S;
//
//    DSUInteger n_mass_balances = ds->massBalances->n;
//    char *fin, *fout, **rxns = ds->massBalances->rxns;
//    DSUInteger i;
//    double factor;
//    DSVariablePool *Xd = DSCaseXd(newCase), *Xd_i;
//    DSMatrix *S_i = NULL;
//
//    // for each mass balances (metabolic block), define if internal stoichiometries need to be adjusted. If so, do this in a recursive way by modifying the stoichiometric matrix.
//
//    for (i=0; i<n_mass_balances; i++){
//        Xd_i = DSVariablePoolAlloc();
//        S_i = s_for_metabolic_block(ds, S, Xd_i, i);
//        fin = newCase->dominant_mass_balances->fin_Xd[i];
//        fout = newCase->dominant_mass_balances->fout_Xd[i];
//        factor = (double)numberZeroBoundariesPerBlock[i];
//        DSCaseCodominantCaseAdjustInternalStoichiometry_recursive(alpha, ds, newCase, fin, fout, rxns, Xd_i, S, factor);
//
//        if (S_i != NULL)
//            DSMatrixFree(S_i);
//        DSVariablePoolFree(Xd_i);
//    }
    
        
        
    // delete after debugging. Manual adjustment of delta matrices for k7 network. original stoichometry
//
//    char aux1[100], aux2[100], aux3[100], aux4[100];
//    double previous_value;
//    DSUInteger index;
//    DSMatrix *delta;
//
//    sprintf(aux1, "%s", "7330_3");
//    sprintf(aux2, "%s", "6");
//    sprintf(aux3, "%s", "9");
//    sprintf(aux4, "%s", "10");
    
//    if (strcmp(aux1, newCase->caseIdentifier) == 0 ){
        
//        index = 2;
//        previous_value = DSMatrixDoubleValue(alpha, index, 0);
//        DSMatrixSetDoubleValue(alpha, index, 0, previous_value * 1.0);
        
//        index = 2;
//        DSMatrixSetDoubleValue(DSSSystemBeta(DSCaseSSystem(newCase)), index, 0,  3.0);
//
//        index = 0;
//        previous_value = DSMatrixDoubleValue(alpha, index, 0);
//        DSMatrixSetDoubleValue(alpha, index, 0, previous_value * 3.0);

//
//        delta = DSCaseDelta(newCase);
//        previous_value = DSMatrixDoubleValue(delta, index, 0);
//        DSMatrixSetDoubleValue(delta, index, 0, log10(pow(10, previous_value) * 3.0 ));
        
//    }
    
    
    
//        if (strcmp(aux1, newCase->caseIdentifier) == 0 ){
//
//            previous_value = DSMatrixDoubleValue(alpha, 0, 0);
//            DSMatrixSetDoubleValue(alpha, 0, 0, previous_value * 0.5);
//
//
//            delta = DSCaseDelta(newCase);
//            previous_value = DSMatrixDoubleValue(delta, 0, 0);
//            DSMatrixSetDoubleValue(delta, 0, 0, log10(pow(10, previous_value) * 2.0 ));
//
//            previous_value = DSMatrixDoubleValue(delta, 2, 0);
//            DSMatrixSetDoubleValue(delta, 2, 0, log10(pow(10, previous_value) * 2.0 ));
//
//
//            if (strcmp(aux1, newCase->caseIdentifier) == 0 ||
//                strcmp(aux3, newCase->caseIdentifier) == 0){
//
//                previous_value = DSMatrixDoubleValue(alpha, 1, 0);
//                DSMatrixSetDoubleValue(alpha, 1, 0, previous_value * 2.0);
//
//            }
//
//
//        }else{
//            delta = DSCaseDelta(newCase);
//            previous_value = DSMatrixDoubleValue(delta, 2, 0);
//            DSMatrixSetDoubleValue(delta, 2, 0, log10(pow(10, previous_value) * 0.5 ));
//
//            previous_value = DSMatrixDoubleValue(delta, 3, 0);
//            DSMatrixSetDoubleValue(delta, 3, 0, log10(pow(10, previous_value) * 0.5 ));
//
//            previous_value = DSMatrixDoubleValue(delta, 4, 0);
//            DSMatrixSetDoubleValue(delta, 4, 0, log10(pow(10, previous_value) * 0.5 ));
//
//            previous_value = DSMatrixDoubleValue(delta, 5, 0);
//            DSMatrixSetDoubleValue(delta, 5, 0, log10(pow(10, previous_value) * 0.5 ));
//        }
    /// end delete **
    

bail:
    return;
}

static void DSCaseCodominantCaseAdjustExternalStoichiometryUsingDominantMassBalance(DSMatrix * alpha,
                                                                                    DSCase *newCase,
                                                                                    const DSUInteger * signature_i,
                                                                                    const DSUInteger * mass_balance_signature,
                                                                                    const char ** fin,
                                                                                    const char ** fout,
                                                                                    DSUInteger n_mass_balances,
                                                                                    DSUInteger *numberZeroBoundariesPerBlock,
                                                                                    const DSDesignSpace *ds,
                                                                                    const double *factors
                                                                                    ){
    
    DSUInteger i, eq, nr_equations, term_i, index_i = 0;
    bool adjusted;
    DSExpression * fin_i = NULL, *fin_i_rhs = NULL, *dominant_input = NULL, ** equations = NULL, *equation_i=NULL, *term_i_expression;
    double previous_value, factor;
    
    equations = DSCaseEquations(newCase);
    nr_equations = DSCaseNumberOfEquations(newCase);
    
    for (i=0; i<n_mass_balances*2; i++ ){
            if (i%2 != 0)
                continue;
            adjusted = false;
            fin_i = DSExpressionByParsingString(fin[i/2]);
            fin_i_rhs = DSExpressionEquationRHSExpression(fin_i);
            if (mass_balance_signature[i] == 1)
                dominant_input = fin_i_rhs;
            else
                dominant_input = DSExpressionBranchAtIndex(fin_i_rhs, signature_i[i]);
            for (eq=0; eq<nr_equations; eq++ ){
                    equation_i = DSExpressionEquationRHSExpression(equations[eq]);
                    for(term_i=0; term_i<DSExpressionNumberOfTerms(equation_i); term_i++){
                            term_i_expression = DSExpressionBranchAtIndex(equation_i, term_i);
                            if (DSExpressionIsEqualToExpressionIgnoringOrder(dominant_input, term_i_expression) == true){
                                    index_i = eq;
                                    adjusted = true;
                                    break;
                            }
                    }
                    if (equation_i != NULL)
                        DSExpressionFree(equation_i);
                    if (adjusted == true)
                        break;
            }
            if (adjusted == false)
               continue;
            if (fin_i != NULL){
                DSExpressionFree(fin_i);
                fin_i = NULL;
            }

            // we use the index to modify the matrix
            previous_value = DSMatrixDoubleValue(DSSSystemAlpha(DSCaseSSystem(newCase)), index_i, 0);
            factor = numberZeroBoundariesPerBlock[i/2] + 1.0;
            DSMatrixSetDoubleValue(alpha, index_i, 0, previous_value/factor);
    
            // we will try a new approach that uses the information contained in the dominant input.
            if (DSDesignSpaceShouldAdjustCodominantBoundaries(ds) == true)
                dsCyclicalDesignSpaceCodominantCaseAdjustConditionMatrices_DominantInput(newCase,
                                                                                         dominant_input,
                                                                                         1.0/factor);
            if (fin_i_rhs != NULL){
                DSExpressionFree(fin_i_rhs);
                fin_i_rhs = NULL;
            }
            
    }
    

    
    if (equations != NULL)
        DSSecureFree(equations);
    
}

extern void DSCaseCodominantCaseAdjustStoichiometryUsingDominantMassBalance(DSCase * newCase,
                                                                            const DSUInteger * signature_i,
                                                                            DSUInteger *numberZeroBoundariesPerBlock,
                                                                            const DSDesignSpace *ds,
                                                                            const double *factors){
    
    if (newCase == NULL){
        DSError(M_DS_CASE_NULL, A_DS_ERROR);
        goto bail;
    }
    
    if (ds->massBalances == NULL)
        goto bail;
    
    const char ** fin = ds->massBalances->fin;
    const char ** fout = ds->massBalances->fout;
    const DSUInteger * mass_balance_signature = ds->massBalances->signature;
    DSUInteger n_mass_balances = ds->massBalances->n;
    
    if (signature_i == NULL || fin == NULL || fout == NULL)
        goto bail;

    DSMatrix *alpha = NULL;
    bool adjust_internal_stoichiometry = false;
    
    // we first take care of the extenal stoichiometry, which is contained in the vector numberZeroBoundariesPerBlock.
    alpha = DSMatrixCopy(DSSSystemAlpha(DSCaseSSystem(newCase)));
    DSCaseCodominantCaseAdjustExternalStoichiometryUsingDominantMassBalance(alpha, newCase, signature_i,
                                                                            mass_balance_signature, fin, fout,
                                                                            n_mass_balances,
                                                                            numberZeroBoundariesPerBlock, ds,
                                                                            factors);
    
    // we now take care of the internal stoichiometry. We will postpone the development of this function for a later time point.
//    DSCaseCodominantCaseAdjustInternalStoichiometryUsingDominantMassBalance(alpha,
//                                                                            newCase,
//                                                                            signature_i,
//                                                                            numberZeroBoundariesPerBlock,
//                                                                            ds,
//                                                                            factors);
    
    newCase->ssys->alpha_adjusted = alpha;
    if (DSDesignSpaceShouldAdjustCodominantBoundaries(ds) == true){
        if (DSSSystemAlpha(DSCaseSSystem(newCase)) != NULL)
            DSMatrixFree(DSSSystemAlpha(DSCaseSSystem(newCase)));
        newCase->ssys->alpha = DSMatrixCopy(newCase->ssys->alpha_adjusted);
    }else{
        DSSSystemSetAdjustCodominantStoichiometry(DSCaseSSystem(newCase), true);
    }
    

bail:
    return;
}



extern bool DSCaseCodominantCaseFulfillMassBalances(const DSDesignSpace *ds,
                                                    DSCase * newCase,
                                                    DSUInteger numberZeroBoundaries,
                                                    DSUInteger *zeroBoundaries,
                                                    DSUInteger *numberZeroBoundariesPerBlock,
                                                    const double *factors){
    
    if (ds == NULL){
        DSError(M_DS_DESIGN_SPACE_NULL, A_DS_ERROR);
        goto bail;
    }
    if (newCase == NULL)
        goto bail;
    
    DSUInteger i, ii, hh, n_subcases = 1, *mass_balance_signature = NULL, n_mass_balances, kk;
    DSUInteger num, * signature_i = NULL, numberOfConstraints = 0;
    char ** inequalities = NULL, **fin, **fout;
    DSCase *aCase_i = NULL;
    bool isValid = false, massBalanceisValid = false, fulfill_mass_balances = false;
    char *id[100];
    
    n_mass_balances = ds->massBalances->n;
    mass_balance_signature = ds->massBalances->signature;
    fin = ds->massBalances->fin;
    fout = ds->massBalances->fout;
    
    
    // Determine the number of subcases to be analyzed and the number of dominance constraints
    for (i=0; i<n_mass_balances*2; i++){
        n_subcases *= mass_balance_signature[i];
        numberOfConstraints += mass_balance_signature[i] - 1;
    }
//    printf("Testing mass balances for case %s \n", newCase->caseIdentifier);
    
    // loop over each subcase. Construct additional inequalities. If the case is valid, check if mass balance conditions are fulfilled. If so, break the loop.
    for (i=0; i<n_subcases; i++){
            
            // we first create signatures for each subcase.
            num = i;
            signature_i = DSSecureMalloc(sizeof(DSUInteger)*2*n_mass_balances);
            for (ii = 0; ii < 2*n_mass_balances; ii++) {
                    signature_i[ii] = num % mass_balance_signature[ii] +1;
                    num /= mass_balance_signature[ii];
            }
        
            // We initialice a subcase
            aCase_i = DSCaseCopy(newCase);
        
            // Ignore codominant boundaries.
            for (hh = 0; hh < numberZeroBoundaries; hh++)
                    DSMatrixSetDoubleValue(DSCaseDelta(aCase_i), zeroBoundaries[hh], 0, log10(2.0));
            DSCaseRecalculateBoundaryMatrices(aCase_i);
                
            // now, we need to create inequalities for each signature.
            inequalities = dsCaseCodominantCaseCreateMassBalanceInequalities(signature_i,
                                                                             mass_balance_signature,
                                                                             fin,
                                                                             fout,
                                                                             n_mass_balances,
                                                                             numberOfConstraints);
            if (inequalities != NULL)
                DSCaseAddConstraints(aCase_i, inequalities, numberOfConstraints);
        
            // Test the validity of the case
            if (DSCaseIsValid(aCase_i, true) == false){
                DSSecureFree(aCase_i);
                aCase_i = NULL;
                continue;
            }
            
//            printf("Mass balances are valid for case %s, subcase %u. About to test validity of the ratio \n", aCase_i->caseIdentifier, i);
        
            // Test the validity of the mass balance fin/fout
            massBalanceisValid = dsCaseCodominatCaseMassBalancesAreValid(aCase_i,
                                                                         signature_i,
                                                                         mass_balance_signature,
                                                                         n_mass_balances,
                                                                         fin,
                                                                         fout);
            if (massBalanceisValid == false){
                DSSecureFree(aCase_i);
                aCase_i = NULL;
                continue;
            }
                
            fulfill_mass_balances = true;
            DSCaseDominantMassBalances(newCase) = DSCaseDominantMassBalances(aCase_i);
            DSCaseDominantMassBalances(aCase_i) = NULL;
        
            DSCaseCodominantCaseAdjustStoichiometryUsingDominantMassBalance(newCase,
                                                                            signature_i,
                                                                            numberZeroBoundariesPerBlock,
                                                                            ds,
                                                                            factors);
            
            if (signature_i != NULL)
                DSSecureFree(signature_i);
            if (inequalities != NULL)
                DSSecureFree(inequalities);
            if (aCase_i != NULL)
                DSCaseFree(aCase_i);
    }
    
    
bail:
    return fulfill_mass_balances;
}


extern void DSCaseRemoveZeroBoundaries(DSCase *aCase){
        DSUInteger * zeroRows = NULL;
        DSUInteger i, j, numberOfZeroRows = 0, maxSize = 1000;
        DSMatrix *temp1, *temp2;
        if (DSCaseU(aCase) == NULL || DSCaseZeta(aCase) == NULL) {
                goto bail;
        }
        temp1 = DSMatrixAppendMatrices(DSCaseU(aCase), DSCaseZeta(aCase), true);
        for (i = 0; i < DSMatrixRows(temp1); i++) {
                for (j = 0; j < DSMatrixColumns(temp1); j++) {
                        if (DSMatrixDoubleValue(temp1, i, j) != 0.0)
                                break;
                }
                if (j == DSMatrixColumns(temp1)) {
                        numberOfZeroRows++;
                        if (numberOfZeroRows == 1) {
                                maxSize = 1000;
                                zeroRows = DSSecureMalloc(sizeof(DSUInteger)*maxSize);
                        } else if (numberOfZeroRows == maxSize) {
                                maxSize += 1000;
                                zeroRows = DSSecureRealloc(&zeroRows, sizeof(DSUInteger)*maxSize);
                        }
                        zeroRows[numberOfZeroRows-1] = i;
                        break;
                }
        }
        if (numberOfZeroRows == 0) {
                DSMatrixFree(temp1);
                goto bail;
        }
        temp2 = DSMatrixSubMatrixExcludingRows(temp1, numberOfZeroRows, zeroRows);
        DSMatrixFree(temp1);
        if (temp2 == NULL)
                goto bail;
        DSMatrixFree(DSCaseU(aCase));
        DSMatrixFree(DSCaseZeta(aCase));
        DSCaseU(aCase) = DSMatrixSubMatrixExcludingColumnList(temp2, 1, DSMatrixColumns(temp2)-1);
        DSCaseZeta(aCase) = DSMatrixSubMatrixIncludingColumnList(temp2, 1, DSMatrixColumns(temp2)-1);
        DSMatrixFree(temp2);
bail:
        if (zeroRows != NULL)
                DSSecureFree(zeroRows);
        return;
}

static void dsCaseCreateBoundaryMatrices(DSCase *aCase)
{
        DSUInteger numberOfXi = 0;
        DSMatrix * W = NULL, *B, *Ai;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (DSSSystemHasSolution(DSCaseSSys(aCase)) == false) {
                goto bail;
        }
        if (DSCaseCd(aCase) == NULL) {
                goto bail;
        }
        B = DSSSystemB(DSCaseSSys(aCase));
        numberOfXi =DSVariablePoolNumberOfVariables(DSCaseXi(aCase));
        W = DSMatrixByMultiplyingMatrix(DSCaseCd(aCase), DSSSystemM(DSCaseSSys(aCase)));
        DSCaseZeta(aCase) = DSMatrixByMultiplyingMatrix(W, B);
        DSMatrixAddByMatrix(DSCaseZeta(aCase), DSCaseDelta(aCase));
        if (numberOfXi != 0) {
                Ai = DSSSystemAi(DSCaseSSys(aCase));
                DSCaseU(aCase) = DSMatrixByMultiplyingMatrix(W, Ai);
                if (DSCaseCi(aCase) != NULL)
                        DSMatrixSubstractByMatrix(DSCaseU(aCase), DSCaseCi(aCase));
                DSMatrixMultiplyByScalar(DSCaseU(aCase), -1.0);
//                DSCaseRemoveRedundantBoundaries(aCase);
                DSMatrixFree(Ai);
        }
        DSMatrixFree(W);
        DSMatrixFree(B);
bail:
        return;
}

extern void DSCaseRecalculateBoundaryMatrices(DSCase *aCase)
{
        DSMatrix * U, *zeta;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (DSSSystemHasSolution(DSCaseSSys(aCase)) == false) {
                goto bail;
        }
        if (DSCaseCd(aCase) == NULL) {
                goto bail;
        }
        U = DSCaseU(aCase);
        zeta = DSCaseZeta(aCase);
        if (U != NULL) {
                DSMatrixFree(U);
                DSMatrixFree(zeta);
                DSCaseU(aCase) = NULL;
                DSCaseZeta(aCase) = NULL;
                dsCaseCreateBoundaryMatrices(aCase);
        }
bail:
        return;
}


static void dsCaseCreateConditionMatrices(DSCase *aCase, const DSGMASystem * gma)
{
        DSUInteger i, j, k, l, numberOfConditions, numberOfEquations;
        DSUInteger numberOfXi, numberOfXd;
        const DSUInteger *termArray;
        double value;
        const DSMatrix * (*a)(const DSGMASystem * gma);
        const DSMatrixArray * (*kd)(const DSGMASystem * gma);
        const DSMatrixArray * (*ki)(const DSGMASystem * gma);
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (gma == NULL) {
                DSError(M_DS_GMA_NULL, A_DS_ERROR);
                goto bail;
        }
        numberOfEquations = DSVariablePoolNumberOfVariables(DSCaseXd(aCase));
        numberOfXd = numberOfEquations;
        numberOfXi = DSVariablePoolNumberOfVariables(DSCaseXi(aCase));
        numberOfConditions = 0;
        termArray = DSCaseSig(aCase);
        for (i = 0; i < 2*numberOfEquations; i++)
                numberOfConditions += DSGMASystemSignature(gma)[i]-1;
        if (numberOfConditions == 0) {
                goto bail;
        }
        DSCaseCd(aCase) = DSMatrixCalloc(numberOfConditions, numberOfXd);
        DSCaseCi(aCase) = DSMatrixCalloc(numberOfConditions, numberOfXi);
        DSCaseDelta(aCase) = DSMatrixCalloc(numberOfConditions, 1);
        for (i = 0, l = 0; i < 2*numberOfEquations; i++) {
                if (i % 2 == 0) {
                        a =  DSGMASystemAlpha;
                        kd = DSGMASystemGd;
                        ki = DSGMASystemGi;
                } else {
                        a =  DSGMASystemBeta;
                        kd = DSGMASystemHd;
                        ki = DSGMASystemHi;
                }
                for (j = 0; j < DSGMASystemSignature(gma)[i]; j++) {
                        if (j == termArray[i]-1)
                                continue;
                        value = log10(DSMatrixDoubleValue(a(gma), i/2, termArray[i]-1)
                                 /DSMatrixDoubleValue(a(gma), i/2, j));
                        DSMatrixSetDoubleValue(DSCaseDelta(aCase), l, 0, value);
                        for (k = 0; k < numberOfXd; k++) {
                                value = DSMatrixArrayDoubleWithIndices(kd(gma), i/2, termArray[i]-1, k);
                                value -= DSMatrixArrayDoubleWithIndices(kd(gma), i/2, j, k);
                                DSMatrixSetDoubleValue(DSCaseCd(aCase), l, k, value);
                        }
                        for (k = 0; k < numberOfXi; k++) {
                                value = DSMatrixArrayDoubleWithIndices(ki(gma), i/2, termArray[i]-1, k);
                                value -= DSMatrixArrayDoubleWithIndices(ki(gma), i/2, j, k);
                                DSMatrixSetDoubleValue(DSCaseCi(aCase), l, k, value);
                        }
                        l++;
                }
        }
//        printf("Reporting from within funciton dsCaseCreateConditionMatrices! \n");
//        printf("Matrix Delta is: \n");
//        DSMatrixPrint(DSCaseDelta(aCase));
//
//        printf("Matrix Cd is: \n");
//        DSMatrixPrint(DSCaseCd(aCase));
//
//        printf("Matrix Ci is: \n");
//        DSMatrixPrint(DSCaseCi(aCase));
    
bail:
        return;
}

static void dsCaseCalculateCaseNumber(DSCase * aCase, const DSGMASystem * gma, const char endianness)
{
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (gma == NULL) {
                DSError(M_DS_GMA_NULL, A_DS_ERROR);
                goto bail;
        }
        DSCaseNum(aCase) = DSCaseNumberForSignature(DSCaseSig(aCase), gma);
bail:
        return;
}

static void dsCaseCalculateCaseIdentifier(DSCase * aCase, const DSGMASystem * gma, const char endianness, const char * prefix)
{
        DSUInteger caseNumber = 0;
        char temp[1000] = {'\0'};
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (gma == NULL) {
                DSError(M_DS_GMA_NULL, A_DS_ERROR);
                goto bail;
        }
        caseNumber = DSCaseNumberForSignature(DSCaseSig(aCase), gma);
        sprintf(temp, "%i", caseNumber);
        if (prefix != NULL) {
                DSCaseId(aCase) = DSSecureCalloc(sizeof(char), strlen(prefix)+2+strlen(temp));
                sprintf(DSCaseId(aCase), "%s_%i", prefix, caseNumber);
        } else {
                DSCaseId(aCase) = DSSecureCalloc(sizeof(char), strlen(temp)+1);
                sprintf(DSCaseId(aCase), "%i", caseNumber);
        }
bail:
        return;
}

extern DSCase * DSCaseWithTermsFromGMA(const DSGMASystem * gma, const DSUInteger * termArray, const char * prefix)
{
        DSCase *aCase = NULL;
        DSUInteger i, term1, term2, numberOfEquations;
        if (gma == NULL) {
                DSError(M_DS_GMA_NULL ": Template GMA to make S-System is NULL", A_DS_ERROR);
                goto bail;
        }
        if (termArray == NULL) {
                DSError(M_DS_NULL ": Array of dominant terms is NULL", A_DS_ERROR);
                goto bail;
        }
        aCase = DSCaseAlloc();
        DSCaseSSys(aCase) = DSSSystemWithTermsFromGMA(gma, termArray);
        aCase->Xi = DSSSystemXi(DSCaseSSys(aCase));
        aCase->Xd = DSSSystemXd(DSCaseSSys(aCase));
        numberOfEquations = DSGMASystemNumberOfEquations(gma);
        DSCaseSig(aCase) = DSSecureMalloc(sizeof(DSUInteger)*(2*numberOfEquations));
        for (i = 0; i < 2*numberOfEquations; i+=2) {
                term1 = termArray[i];
                term2 = termArray[i+1];
                DSCaseSig(aCase)[i] = term1;
                DSCaseSig(aCase)[i+1] = term2;
                if (term1 > DSGMASystemSignature(gma)[i] || term2 > DSGMASystemSignature(gma)[i+1])
                        break;
                if (term1 <= 0 || term2 <= 0)
                        break;
        }
        if (i == 2*numberOfEquations) {
                dsCaseCreateConditionMatrices(aCase, gma);
                dsCaseCreateBoundaryMatrices(aCase);
                dsCaseCalculateCaseNumber(aCase, gma, endian);
                dsCaseCalculateCaseIdentifier(aCase, gma, endian, NULL);
        } else {
                DSCaseFree(aCase);
                aCase = NULL;
        }
bail:
        return aCase;
}

static void dsCaseAppendDesignSpaceConditions(DSCase * aCase, const DSDesignSpace * ds)
{
        DSMatrix * temp;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (ds == NULL) {
                DSError(M_DS_DESIGN_SPACE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (ds->Cd == NULL) 
                goto bail;
        if (DSCaseCd(aCase) == NULL) {
                DSCaseCd(aCase) = DSMatrixCopy(ds->Cd);
                DSCaseDelta(aCase) = DSMatrixCopy(ds->delta);
                if (DSVariablePoolNumberOfVariables(DSCaseXi(aCase)) > 0)
                        DSCaseCi(aCase) = DSMatrixCopy(ds->Ci);
        } else {
                temp = DSMatrixAppendMatrices(DSCaseCd(aCase), ds->Cd, false);
                DSMatrixFree(DSCaseCd(aCase));
                DSCaseCd(aCase) = temp;
                temp = DSMatrixAppendMatrices(DSCaseDelta(aCase), ds->delta, false);
                DSMatrixFree(DSCaseDelta(aCase));
                DSCaseDelta(aCase) = temp;
                if (DSVariablePoolNumberOfVariables(DSCaseXi(aCase)) > 0) {
                        temp = DSMatrixAppendMatrices(DSCaseCi(aCase), ds->Ci, false);
                        DSMatrixFree(DSCaseCi(aCase));
                        DSCaseCi(aCase) = temp;
                }

        }
bail:
        return;
}

static void dsSSystemSolveEquations(DSSSystem *ssys)
{
    DSMatrix *M, *Ad;
    if (ssys == NULL) {
        DSError(M_DS_NULL ": S-System being modified is NULL", A_DS_ERROR);
        goto bail;
    }
    DSSSystemSetIsSingular(ssys, true);
    Ad = DSMatrixBySubstractingMatrix(ssys->Gd, ssys->Hd);
    M = DSMatrixInverse(Ad);
    if (M != NULL) {
        DSSSystemSetIsSingular(ssys, false);
        ssys->M = M;
    }
    DSMatrixFree(Ad);
bail:
    return;
}

static void dsSSystemReshapeMatrices(DSSSystem *ssys, DSUInteger numberOfAssociatedVariables , DSUInteger numberOfConservations,
                                     DSUInteger *variables_to_delete, DSUInteger *conservation_indices, DSDictionary *swap)


{
    
        DSMatrix *alpha_temp = NULL, *beta_temp = NULL, *Gd_temp = NULL;
        DSMatrix *Gi_temp = NULL, *Hd_temp = NULL, *Hi_temp = NULL;
        DSUInteger i, col1, col2;
    
        // Delete Rows contained in variables_to_delete of all matrices! For matrices Gd and Hd also delete columns of Xcn
        if (ssys->alpha != NULL){
            alpha_temp = DSMatrixSubMatrixExcludingRows(ssys->alpha, numberOfConservations, variables_to_delete);
            DSMatrixFree(ssys->alpha);
            ssys->alpha = alpha_temp;
        }
    
        if (ssys->beta != NULL){
            beta_temp = DSMatrixSubMatrixExcludingRows(ssys->beta, numberOfConservations, variables_to_delete);
            DSMatrixFree(ssys->beta);
            ssys->beta = beta_temp;
        }
    
        // Matrices Gd need to be treated in a special way!
        if (ssys->Gd != NULL){
            // We first switch colums corresponding to realationships Xc <-> X1. These relationships are contained in swap dictionary
            for (i=0; i<DSDictionaryCount(swap); i++){
                col1 = DSVariablePoolIndexOfVariableWithName(ssys->Xd, DSDictionaryNames(swap)[i]);
                col2 = DSVariablePoolIndexOfVariableWithName(ssys->Xd, DSDictionaryValueForName(swap, DSDictionaryNames(swap)[i]));
                DSMatrixSwitchColumns(ssys->Gd, col1, col2);
            }
            Gd_temp = DSMatrixSubMatrixExcludingRowsAndColumns(ssys->Gd, numberOfAssociatedVariables , numberOfAssociatedVariables, variables_to_delete, variables_to_delete);
            DSMatrixFree(ssys->Gd);
            ssys->Gd = Gd_temp;
        }
    
        if (ssys->Gi != NULL){
            Gi_temp = DSMatrixSubMatrixExcludingRows(ssys->Gi, numberOfConservations, variables_to_delete);
            DSMatrixFree(ssys->Gi);
            ssys->Gi = Gi_temp;
        }
    
        // Matrices Hd need to be treated in a special way!
        if (ssys->Hd != NULL){
            // We first switch colums corresponding to realationships Xc <-> X1. These relationships are contained in swap dictionary
            for (i=0; i<DSDictionaryCount(swap); i++){
                col1 = DSVariablePoolIndexOfVariableWithName(ssys->Xd, DSDictionaryNames(swap)[i]);
                col2 = DSVariablePoolIndexOfVariableWithName(ssys->Xd, DSDictionaryValueForName(swap, DSDictionaryNames(swap)[i]));
                DSMatrixSwitchColumns(ssys->Hd, col1, col2);
            }
            Hd_temp = DSMatrixSubMatrixExcludingRowsAndColumns(ssys->Hd, numberOfAssociatedVariables , numberOfAssociatedVariables, variables_to_delete, variables_to_delete);
            DSMatrixFree(ssys->Hd);
            ssys->Hd = Hd_temp;
        }
    
        if (ssys->Hi != NULL){
            Hi_temp = DSMatrixSubMatrixExcludingRows(ssys->Hi, numberOfConservations, variables_to_delete);
            DSMatrixFree(ssys->Hi);
            ssys->Hi = Hi_temp;
        }
    
    
}

static void dsSSystemReshapeVariablePools(DSSSystem * ssys, const DSDesignSpace *ds,
                                          const DSDictionary * variables_to_delete_dic,
                                          const DSDictionary * conservation_variables_dic)
{
    
    // Now re-shape variable pools. For each Xci, do:
    // Eliminate Xci from Xd
    // Change identity of associated Xd from Xd_t to Xd_a. Add that Xd to Xd_a_c
    
    // The idea would be to loop over ssys->Xd and then create new pools depending on the existence of that certain variable on either dictionary.
    
    // todo: check memory management of Xd and corresponding sub-pools!. SSystem should delete these pools!
    
    DSVariablePool *Xd_gma = ds->gma->Xd;
    DSVariablePool *Xd_a_gma = ds->gma->Xd_a;
    DSVariablePool *Xd_new = NULL;
    DSVariablePool *Xd_t_new = NULL;
    DSVariablePool *Xd_a_new = NULL;
    DSVariablePool *Xd_a_c = NULL;
    DSUInteger i, count=0;
    const char ** variableNames = NULL;
    const char *associatedVariable = NULL;
    
    variableNames = DSVariablePoolAllVariableNames(Xd_gma);
    Xd_new = DSVariablePoolAlloc();
    Xd_t_new = DSVariablePoolAlloc();
    Xd_a_new = DSVariablePoolAlloc();
    Xd_a_c = DSVariablePoolAlloc();
    
    
    //    printf("dictionary variables_to_delete_dic contains following variables: \n");
    //    DSDictionaryPrint(variables_to_delete_dic);
    //    printf("Showing %u names saved in variables_to_delete_dic : \n ", DSDictionaryCount(variables_to_delete_dic));
    //    for(i=0; i<DSDictionaryCount(variables_to_delete_dic); i++){
    //        printf(" %s \n", DSDictionaryNames(variables_to_delete_dic)[i]);
    //    }
    
    
    
    //    printf("dictionary conservation_variables_dic contains following variables: \n");
    //    DSDictionaryPrint(conservation_variables_dic);
    //    printf("Showing %u names saved in conservation_variables_dic : \n ", DSDictionaryCount(conservation_variables_dic));
    //    for(i=0; i<DSDictionaryCount(conservation_variables_dic); i++){
    //        printf(" %s \n", DSDictionaryNames(conservation_variables_dic)[i]);
    //    }
    
    
    
    // loop over Xd_gma and populate Xd_new
    for (i=0; i<DSVariablePoolNumberOfVariables(Xd_gma); i++){
        
        //        printf("Analyzing variable with name %s .. \n ",variableNames[i] );
        
        if( DSDictionaryValueForName(variables_to_delete_dic, variableNames[i]) != NULL ){
            //            printf("Skipping variable %s \n", variableNames[i]);
            continue;
        }
        
        if( DSDictionaryValueForName(conservation_variables_dic, variableNames[i]) != NULL ){
            // if variable is a conservation variable, add associated variable to the dictionary Xc --> X1
            associatedVariable = DSDictionaryValueForName(conservation_variables_dic, variableNames[i]);
            
            if (DSVariablePoolHasVariableWithName(Xd_new, associatedVariable) == false){
                DSVariablePoolAddVariableWithName(Xd_new, associatedVariable);
                
            }
            continue;
        }
        //        printf("About to add variable %s without any check \n", variableNames[i] );
        DSVariablePoolAddVariableWithName(Xd_new, variableNames[i]);
        if (DSVariablePoolHasVariableWithName(Xd_a_gma, variableNames[i]) == false){
            DSVariablePoolSetValueForVariableWithName(Xd_new, variableNames[i], count);
            count++;
        }
        
    }
    if (variableNames != NULL)
        DSSecureFree(variableNames);
    
    // First assign true Xd_a from Xd_a_gma to Xd_a_new.
    // Loop over Xd_a_gma and assign to Xd_a_new only if variable is not present in dictionary "conservation_variables_dic"
    variableNames = DSVariablePoolAllVariableNames(Xd_a_gma);
    for (i=0; i<DSVariablePoolNumberOfVariables(Xd_a_gma); i++ ){
        if (DSDictionaryValueForName(conservation_variables_dic, variableNames[i]) == NULL)
            DSVariablePoolAddVariableWithName(Xd_a_new, variableNames[i]);
    }
    if (variableNames != NULL)
        DSSecureFree(variableNames);
    
    // now loop over the newly created Xd_new and assign each variable to either Xd_t_new or Xd_a_new.
    // Assign variable to Xd_a_new if variable is present in dictionary "variables_to_delete_dic", else to Xd_t only if variable is not present in Xd_a_new already
    variableNames = DSVariablePoolAllVariableNames(Xd_new);
    for (i=0; i<DSVariablePoolNumberOfVariables(Xd_new); i++){
        
        if (DSDictionaryValueForName(variables_to_delete_dic, variableNames[i]) != NULL){
            DSVariablePoolAddVariableWithName(Xd_a_new, variableNames[i]);
            DSVariablePoolAddVariableWithName(Xd_a_c, variableNames[i]);
        } else {
            if (DSVariablePoolHasVariableWithName(Xd_t_new, variableNames[i]) == false && DSVariablePoolHasVariableWithName(Xd_a_gma, variableNames[i]) == false )
                DSVariablePoolAddVariableWithName(Xd_t_new, variableNames[i]);
        }
    }
    if (variableNames != NULL)
        DSSecureFree(variableNames);
    
    // now let's print results!
    //    printf("Xd_gma is: \n");
    //    DSVariablePoolPrint(Xd_gma);
    //    printf("Xd_new is: \n");
    //    DSVariablePoolPrint(Xd_new);
    //
    //    printf("Xd_t_gma is: \n");
    //    DSVariablePoolPrint(ds->gma->Xd_t);
    //    printf("Xd_t_new is: \n");
    //    DSVariablePoolPrint(Xd_t_new);
    //
    //    printf("Xd_a_gma is: \n");
    //    DSVariablePoolPrint(ds->gma->Xd_a);
    //    printf("Xd_a_new is: \n");
    //    DSVariablePoolPrint(Xd_a_new);
    
    // and now let's assign the newly created pools to the ssystem
    ssys->Xd = Xd_new;
    ssys->Xd_t = Xd_t_new;
    ssys->Xd_a = Xd_a_new;
    ssys->Xd_a_c = Xd_a_c;
    DSSSystemSetShouldFreeXd(ssys, true);
    
    
}

static void dsSSystemReshapeVariablePools_future(DSSSystem * ssys, const DSDesignSpace *ds,
                                          const DSDictionary * variables_to_delete_dic,
                                          const DSDictionary * conservation_variables_dic,
                                          char * conservedKey)
{
    
        // Now re-shape variable pools. For each Xci, do:
        // Eliminate Xci from Xd
        // Change identity of associated Xd from Xd_t to Xd_a. Add that Xd to Xd_a_c
        // The idea would be to loop over ssys->Xd and then create new pools depending on the existence of that certain variable on either dictionary.
    
        DSVariablePool *Xd_gma = ds->gma->Xd;
        DSVariablePool *Xd_a_gma = ds->gma->Xd_a;
        DSVariablePool *Xd_new = NULL;
        DSVariablePool *Xd_t_new = NULL;
        DSVariablePool *Xd_a_new = NULL;
        DSVariablePool *Xd_a_c = NULL;
        DSUInteger i, count=0;
        const char ** variableNames = NULL;
        const char *associatedVariable = NULL;
    
        if (DSDictionaryValueForName(ds->Xd_dic, conservedKey) != NULL ||
            DSDictionaryValueForName(ds->Xd_t_dic, conservedKey) != NULL ||
            DSDictionaryValueForName(ds->Xd_a_dic, conservedKey) != NULL ||
            DSDictionaryValueForName(ds->Xd_a_c_dic, conservedKey) != NULL)
            goto bail;
    
        variableNames = DSVariablePoolAllVariableNames(Xd_gma);
        Xd_new = DSVariablePoolAlloc();
        Xd_t_new = DSVariablePoolAlloc();
        Xd_a_new = DSVariablePoolAlloc();
        Xd_a_c = DSVariablePoolAlloc();
    
        // loop over Xd_gma and populate Xd_new
        for (i=0; i<DSVariablePoolNumberOfVariables(Xd_gma); i++){
                if( DSDictionaryValueForName(variables_to_delete_dic, variableNames[i]) != NULL ){
                    continue;
                }
            
                if( DSDictionaryValueForName(conservation_variables_dic, variableNames[i]) != NULL ){
                        // if variable is a conservation variable, add associated variable to the dictionary Xc --> X1
                    associatedVariable = DSDictionaryValueForName(conservation_variables_dic, variableNames[i]);
                    if (DSVariablePoolHasVariableWithName(Xd_new, associatedVariable) == false){
                                DSVariablePoolAddVariableWithName(Xd_new, associatedVariable);

                    }
                    continue;
                }
                DSVariablePoolAddVariableWithName(Xd_new, variableNames[i]);
                if (DSVariablePoolHasVariableWithName(Xd_a_gma, variableNames[i]) == false){
                    DSVariablePoolSetValueForVariableWithName(Xd_new, variableNames[i], count);
                    count++;
                }
            
        }
        if (variableNames != NULL)
            DSSecureFree(variableNames);
    
        // First assign true Xd_a from Xd_a_gma to Xd_a_new.
        // Loop over Xd_a_gma and assign to Xd_a_new only if variable is not present in dictionary "conservation_variables_dic"
        variableNames = DSVariablePoolAllVariableNames(Xd_a_gma);
        for (i=0; i<DSVariablePoolNumberOfVariables(Xd_a_gma); i++ ){
                if (DSDictionaryValueForName(conservation_variables_dic, variableNames[i]) == NULL)
                    DSVariablePoolAddVariableWithName(Xd_a_new, variableNames[i]);
        }
        if (variableNames != NULL)
            DSSecureFree(variableNames);
    
        // now loop over the newly created Xd_new and assign each variable to either Xd_t_new or Xd_a_new.
        // Assign variable to Xd_a_new if variable is present in dictionary "variables_to_delete_dic", else to Xd_t only if variable is not present in Xd_a_new already
        variableNames = DSVariablePoolAllVariableNames(Xd_new);
        for (i=0; i<DSVariablePoolNumberOfVariables(Xd_new); i++){
            
                if (DSDictionaryValueForName(variables_to_delete_dic, variableNames[i]) != NULL){
                    DSVariablePoolAddVariableWithName(Xd_a_new, variableNames[i]);
                    DSVariablePoolAddVariableWithName(Xd_a_c, variableNames[i]);
                } else {
                        if (DSVariablePoolHasVariableWithName(Xd_t_new, variableNames[i]) == false && DSVariablePoolHasVariableWithName(Xd_a_gma, variableNames[i]) == false )
                                DSVariablePoolAddVariableWithName(Xd_t_new, variableNames[i]);
                }
        }
        if (variableNames != NULL)
            DSSecureFree(variableNames);
    
//    if (DSDictionaryValueForName(ds->Xd_dic, conservedKey) == NULL ||
//        DSDictionaryValueForName(ds->Xd_t_dic, conservedKey) == NULL ||
//        DSDictionaryValueForName(ds->Xd_a_dic, conservedKey) == NULL ||
//        DSDictionaryValueForName(ds->Xd_a_c_dic, conservedKey) == NULL){
    
                DSDictionaryAddValueWithName(ds->Xd_dic, conservedKey, Xd_new);
                DSDictionaryAddValueWithName(ds->Xd_t_dic, conservedKey, Xd_t_new);
                DSDictionaryAddValueWithName(ds->Xd_a_dic, conservedKey, Xd_a_new);
                DSDictionaryAddValueWithName(ds->Xd_a_c_dic, conservedKey, Xd_a_c);
//    } else{
//                DSVariablePoolFree(Xd_new);
//                DSVariablePoolFree(Xd_t_new);
//                DSVariablePoolFree(Xd_a_new);
//                DSVariablePoolFree(Xd_a_c);
//    }
    
bail:
    
    ssys->Xd = DSDictionaryValueForName(ds->Xd_dic, conservedKey);
    ssys->Xd_t = DSDictionaryValueForName(ds->Xd_t_dic, conservedKey);
    ssys->Xd_a = DSDictionaryValueForName(ds->Xd_a_dic, conservedKey);
    ssys->Xd_a_c = DSDictionaryValueForName(ds->Xd_a_c_dic, conservedKey);

}

static void dSCaseReshapeConditionMatrices(DSCase *aCase, DSUInteger numberOfAssociatedVariables,
                                           DSUInteger *variables_to_delete, DSDictionary *swap)
{
    // this function is used to reshape one of the condition matrices aCase->Cd. The modification consist on first swapping colums associated in Xc <-> X1 and then deleting colums associated with  associated variables contained in *variables_to_delete.
    
        DSUInteger i, col1, col2;
        DSSSystem *ssys = aCase->ssys;
        DSMatrix * Cd_new;
    
        //swap
        for (i=0; i<DSDictionaryCount(swap); i++){
            col1 = DSVariablePoolIndexOfVariableWithName(ssys->Xd, DSDictionaryNames(swap)[i]);
            col2 = DSVariablePoolIndexOfVariableWithName(ssys->Xd, DSDictionaryValueForName(swap, DSDictionaryNames(swap)[i]));
            DSMatrixSwitchColumns(aCase->Cd, col1, col2);
        }
    
        //delete columns
        Cd_new = DSMatrixSubMatrixExcludingColumns(aCase->Cd, numberOfAssociatedVariables, variables_to_delete);
        DSMatrixFree(aCase->Cd);
        aCase->Cd = Cd_new;
    
}

static void dSCaseAdjustCaseSignature(DSCase *aCase,
                                      DSUInteger numberOfAssociatedVariables, DSUInteger *variables_to_delete)
{
    // First copy aCase->signature into aCase->conserved_sig and then modify it
    // Set 0 for all indices contained in variables_to_delete.

    
        DSUInteger ii, index;
        DSUInteger numberOfEquations;
    
        if (DSCaseSig(aCase) == NULL)
        goto bail;
    
        if (DSCase3Sig(aCase) != NULL)
        {
            goto bail;
        } else {
                numberOfEquations = DSCaseNumberOfConservations(aCase) + DSCaseNumberOfEquations(aCase);
                DSCaseSigCons(aCase) = DSSecureCalloc(sizeof(DSUInteger), numberOfEquations*2);
                for (ii = 0; ii< numberOfEquations*2; ii++)
                    DSCaseSigCons(aCase)[ii] = DSCaseSig(aCase)[ii];
            
                for (ii=0; ii<numberOfAssociatedVariables; ii++){
                    index = variables_to_delete[ii];
                    DSCaseSigCons(aCase)[2*index] = 0;
                    DSCaseSigCons(aCase)[2*index + 1] = 0;
                }
        }
    
    bail:
        return;
    
}


static void dsSSystemIntegrateConservedRelationships(DSCase *aCase, const DSDesignSpace *ds)
{
        // This function should help re-shape the matrices of the S-system so that auxiliary variables
        // can be sucessfully merged. Matrices that need to be re-shaped are Alpha, Beta, Gd, Hd, Gi, Hi.
        // The number of conservation relationships is contained in ds->numberOfConservations. The conservation
        // relationships correspond to the last n equations, being n = ds->numberOfConservations.
        // Conservation relationships have been saved till this point as Xd_a variables.
        // The idea of this function is to: a) reshape matrices and b) reshape variable pools.
    

        DSSSystem *ssys = aCase->ssys;
        DSUInteger *conservation_indices = NULL, *variables_to_delete = NULL;
        DSUInteger numberOfConservations, i, conservation_indx, associated_indx, numberOfAssociatedVariables = 0;
        DSVariablePool * Xd = ssys->Xd;
        DSMatrix *Hd = ssys->Hd;
        char conservedVar[100], *associatedVar, conservedKey[1000], buff[10];
        DSDictionary *variables_to_delete_dic = NULL, *conservation_variables_dic = NULL;
    
        if (DSDesignSpaceConserved(ds) == false){
            aCase->ssys->numberOfConservations = 0;
            goto bail;
        }
    
        if (aCase->ssys == NULL)
            goto bail;
    
        aCase->ssys->numberOfConservations = ds->numberOfConservations;
        DSSSystemSetIsConserved(aCase->ssys, true);
    
        variables_to_delete_dic = DSDictionaryAlloc();
        conservation_variables_dic = DSDictionaryAlloc();
        numberOfConservations = ds->numberOfConservations;
    
        // find indices for conservation variables Xc1, Xc2,.., Xcn
        // find indices of rows of matrices that need to be deleted. Check for non-zero values in Matrix Hd
        conservation_indices = DSSecureMalloc(sizeof(DSUInteger)*numberOfConservations);
        variables_to_delete = DSSecureMalloc(sizeof(DSUInteger)*numberOfConservations);
        for (i=0; i<numberOfConservations; i++){
                sprintf(conservedVar, "Xc%u", i+1);
                conservation_indx = DSVariablePoolIndexOfVariableWithName(Xd, conservedVar);
                associated_indx = DSMatrixFirstNonZeroIndexAtRow(Hd, conservation_indx);
                if (associated_indx == 65000 ){
                    char error[1000];
                    sprintf(error, "Parsing Error in Conservation Eqn. # %u"": It must be a f() of variables with a differential equation", i+1);
                    DSError(error, A_DS_ERROR);
                    goto bail;
                }
                if (i==0)
                    sprintf(conservedKey, "%u%u", DSCaseSignature(aCase)[2*conservation_indx],
                                                  DSCaseSignature(aCase)[2*conservation_indx + 1]);
                else{
                    sprintf(buff, "%u%u", DSCaseSignature(aCase)[2*conservation_indx],
                            DSCaseSignature(aCase)[2*conservation_indx + 1]);
                    strcat(conservedKey, buff);
                }
                associatedVar = DSVariablePoolVariableAtIndex(Xd, associated_indx)->name;
                DSDictionaryAddValueWithName(conservation_variables_dic, conservedVar, associatedVar);
                if (DSDictionaryValueForName(variables_to_delete_dic, associatedVar) == NULL){
                        DSDictionaryAddValueWithName(variables_to_delete_dic, associatedVar, (void*)1);
                        variables_to_delete[numberOfAssociatedVariables] = associated_indx;
                        numberOfAssociatedVariables++;
                }
                conservation_indices[i] = conservation_indx;
        }
        if (numberOfConservations != numberOfAssociatedVariables ){
            DSSSystemSetIsConserved(aCase->ssys, false);
            aCase->ssys->numberOfConservations = 0;
            goto bail;
        }
    
        dsSSystemReshapeMatrices(ssys, numberOfAssociatedVariables ,
                                 numberOfConservations, variables_to_delete, conservation_indices,
                                 conservation_variables_dic);
        dSCaseReshapeConditionMatrices(aCase, numberOfAssociatedVariables, variables_to_delete,     conservation_variables_dic);
//        dsSSystemReshapeVariablePools_future(ssys, ds, variables_to_delete_dic, conservation_variables_dic, conservedKey);
            dsSSystemReshapeVariablePools(ssys, ds, variables_to_delete_dic, conservation_variables_dic);
    
        if (DSDesignSpaceCyclical(ds) == false)
            dSCaseAdjustCaseSignature(aCase, numberOfAssociatedVariables, variables_to_delete);
        dsSSystemSolveEquations(ssys);
    
        aCase->Xi = DSSSystemXi(DSCaseSSys(aCase));
        aCase->Xd = DSSSystemXd(DSCaseSSys(aCase));
        aCase->Xd_a = DSSSystemXd_a(DSCaseSSys(aCase));
    
    bail:
        if (variables_to_delete_dic != NULL)
            DSDictionaryFree(variables_to_delete_dic);
        if (conservation_variables_dic != NULL)
            DSDictionaryFree(conservation_variables_dic);
        if (variables_to_delete != NULL)
            DSSecureFree(variables_to_delete);
        if (conservation_indices != NULL)
            DSSecureFree(conservation_indices);

        return;
}

static void dsSubCaseSignatureProcessLocationThreeDigitMatrices(DSuIntegerMatrix *location,
                                                                DSuIntegerMatrix *three_digit,
                                                                DSUInteger *parentSignature,
                                                                DSUInteger numberCycles,
                                                                DSUInteger numberInheritedConservations,
                                                                DSUInteger numberOfEqutions,
                                                                bool print)
{
    
    // This function should modify integer matrices *location and *three_digit. The number of modifications is dictated by the number of inherited conservations and by the location of the conservations. This location can be extracted from the parentSignature. A triplet of zeros symbolize the location of a conservation.
    
    DSUInteger * conservations_location = NULL;
    DSUInteger i, ii, count = 0, increase_location, increase_three_digit, value_location, value_three_digit;
    
    conservations_location = DSSecureCalloc(sizeof(DSUInteger), numberInheritedConservations);
    
    // Generate conservations_location containing indices for the location of the conservations.
    for (i = 0; i < numberOfEqutions; i++){
        if (parentSignature[3*i] == 0 && parentSignature[3*i + 1] == 0 && parentSignature[3*i + 2] == 0){
            conservations_location[count] = i;
            count++;
        }
    }
    
    // loop over location and three_digit and modify if conservations are located before cycles.
    for (i = 0; i < numberCycles; i++){
            increase_location = 0;
            increase_three_digit = 0;
        
            for (ii = 0; ii<numberInheritedConservations; ii++){
                if(conservations_location[ii] <= DSuIntegerMatrixValue(location, i, 0))
                    increase_location++;
                if(conservations_location[ii] <= DSuIntegerMatrixValue(three_digit, i, 0))
                    increase_three_digit++;
            }
        
            // loc: add the number of conservations before loc.
            // three_digit: add the number of conservations before loc to the first column.
            value_location = DSuIntegerMatrixValue(location, i, 0) + increase_location;
            value_three_digit = DSuIntegerMatrixValue(three_digit, i, 0) + increase_three_digit;
            DSuIntegerMatrixSetValue(location, i, 0, value_location);
            DSuIntegerMatrixSetValue(three_digit, i, 0, value_three_digit);
    }
    
    if (conservations_location != NULL)
        DSSecureFree(conservations_location);
    
    return;
}

static void DSSubCaseAllMainAndSecondaryVariables(const DSDesignSpace * ds,
                                                  DSUInteger  *all_main_and_secondary_variables,
                                                  DSUInteger *count){
    
    DSUInteger numberCycles, numberAllSecondaryVariables=0, i, jj, *aux =NULL;
    DSDesignSpace *previous = NULL;
    
    if (ds == NULL || ds->extensionData == NULL)
        goto bail;
    
    numberCycles = ds->extensionData->numberCycles;
    previous = ds->extensionData->parent_ds;
    
    for (i=0; i<numberCycles; i++)
        numberAllSecondaryVariables += ds->extensionData->numberSecondaryVariables[i];
    
    if (*count == 0){
        all_main_and_secondary_variables = (DSUInteger *)DSSecureRealloc(all_main_and_secondary_variables,
                                                                         sizeof(DSUInteger)*(numberCycles + numberAllSecondaryVariables));
    }else{
        all_main_and_secondary_variables = (DSUInteger *)DSSecureRealloc(all_main_and_secondary_variables,
                                                                         (*count + numberCycles + numberAllSecondaryVariables)*sizeof(DSUInteger));
    }
        
    for (i = 0; i<numberCycles; i++){
        all_main_and_secondary_variables[*count] = ds->extensionData->mainCycleVariables[i];
        *count = *count + 1;
        for (jj=0; jj<ds->extensionData->numberSecondaryVariables[i]; jj++){
            all_main_and_secondary_variables[*count] = ds->extensionData->allSecondaryVariables[i][jj];
            *count = *count + 1;
        }
    }
    
    if (previous != NULL)
        DSSubCaseAllMainAndSecondaryVariables(previous, all_main_and_secondary_variables, count);
    
bail:
    return;
}



extern void DSSubCaseGenerate3dSignature(const DSDesignSpace *ds, DSCase *aCase, DSuIntegerMatrix *three_digit, DSuIntegerMatrix *location){
    
        DSUInteger numberOfEquations, numberCycles;
        DSUInteger *all_main_and_secondary_variables, vector_count = 0;
        DSUInteger ii, loc, jj;
        bool is_main_secondary = false;
        numberOfEquations = DSGMASystemNumberOfEquations(DSDesignSpaceGMASystem(ds)) + ds->numberInheritedConservations;

        if (DSDesignSpaceCyclical(ds) == true && ds->parent3DigitsSignature != NULL && three_digit != NULL) {
            numberCycles = ds->extensionData->numberCycles;
            // we first copy the 3d signature from the parent design space.
            aCase->signature_3d = DSSecureCalloc(sizeof(DSUInteger), numberOfEquations*3);
            for (ii = 0; ii< numberOfEquations*3; ii++){
                aCase->signature_3d[ii] = ds->parent3DigitsSignature[ii];
            }
                        
            if(ds->numberInheritedConservations != 0)
                dsSubCaseSignatureProcessLocationThreeDigitMatrices(location, three_digit,
                                                                    aCase->signature_3d, numberCycles,
                                                                    ds->numberInheritedConservations,
                                                                    numberOfEquations, false);
            
            // then we use information located in *three_digit and *location to modify that signature.
            for (ii = 0; ii < numberCycles; ii++ ){
                loc = DSuIntegerMatrixValue(location, ii, 0);
                aCase->signature_3d[loc*3] = DSuIntegerMatrixValue(three_digit, ii, 0) + 1;
                aCase->signature_3d[loc*3 + 1] = DSuIntegerMatrixValue(three_digit, ii, 1) + 1;
                aCase->signature_3d[loc*3 + 2] = DSuIntegerMatrixValue(three_digit, ii, 2) + 1;
            }
            
            // finally, we copy the internal signature from the case into the 3d signature for variables that are not main cyclical or secondary and have more than 1 element in the respective position
            all_main_and_secondary_variables = DSSecureCalloc(sizeof(DSUInteger), 1);
            DSSubCaseAllMainAndSecondaryVariables(ds, all_main_and_secondary_variables, &vector_count);
            for (ii = 0; ii<numberOfEquations; ii++){
                is_main_secondary = false;
                for (jj=0; jj<vector_count; jj++){
                    if (all_main_and_secondary_variables[jj] == ii){
                        is_main_secondary = true;
                    }
                }
                if (is_main_secondary == true)
                    continue;
                if (DSGMASystemSignature(ds->gma)[ii] == 1)
                    continue;
                aCase->signature_3d[ii*3 + 1] = DSCaseSignature(aCase)[ii*2];
                aCase->signature_3d[ii*3 + 2] = DSCaseSignature(aCase)[ii*2 + 1];
            }
            
            
            DSuIntegerMatrixFree(location);
            DSuIntegerMatrixFree(three_digit);
            if (all_main_and_secondary_variables != NULL)
                DSSecureFree(all_main_and_secondary_variables);
        } else {
            if(three_digit != NULL)
                DSuIntegerMatrixFree(three_digit);
            if(location != NULL)
                DSuIntegerMatrixFree(location);
        }
}

extern DSCase * DSCaseWithTermsFromDesignSpace(const DSDesignSpace * ds,
                                               const DSUInteger * termArray,
                                               const char * prefix)
{
        DSCase *aCase = NULL;
        DSUInteger i, term1, term2, numberOfEquations;
        DSUnstableCase * uCase = NULL;
        char name[1000];
        DSUInteger numberCycles;
        DSuIntegerMatrix *three_digit = NULL, *location = NULL;
//        char aux[100];
    
        if (ds == NULL) {
                DSError(M_DS_NULL ": Template GMA to make S-System is NULL", A_DS_ERROR);
                goto bail;
        }
        if (termArray == NULL) {
                DSError(M_DS_DESIGN_SPACE_NULL ": Array of dominant terms is NULL", A_DS_ERROR);
                goto bail;
        }
        aCase = DSCaseAlloc();
        numberOfEquations = DSGMASystemNumberOfEquations(DSDesignSpaceGMASystem(ds));
        if (DSDesignSpaceCyclical(ds) == true) {
            numberCycles = ds->extensionData->numberCycles;
            three_digit = DSuIntegerMatrixCalloc(numberCycles, 3);
            location = DSuIntegerMatrixCalloc(numberCycles, 1);
//            if (ds->casePrefix != NULL)
//                sprintf(aux, "%s", ds->casePrefix);
//            printf("****4 Calling DSSSystemWithTermsFromGMACyclical for case %u (%s) \n", DSCaseNumberForSignature(termArray, ds->gma), aux);
            DSCaseSSys(aCase) = DSSSystemWithTermsFromGMACyclical(ds, termArray, three_digit, location);
            aCase->numberInheritedConservations = ds->numberInheritedConservations;
        }else{
            DSCaseSSys(aCase) = DSSSystemWithTermsFromGMA(DSDesignSpaceGMASystem(ds), termArray);
            aCase->numberInheritedConservations = 0;
        }
        aCase->Xi = DSSSystemXi(DSCaseSSys(aCase));
        aCase->Xd = DSSSystemXd(DSCaseSSys(aCase));
        aCase->Xd_a = DSSSystemXd_a(DSCaseSSys(aCase));
        DSCaseSig(aCase) = DSSecureMalloc(sizeof(DSUInteger)*(2*numberOfEquations));
        for (i = 0; i < 2*numberOfEquations; i+=2) {
            term1 = termArray[i];
            term2 = termArray[i+1];
            DSCaseSig(aCase)[i] = term1;
            DSCaseSig(aCase)[i+1] = term2;
            if (term1 > DSGMASystemSignature(DSDesignSpaceGMASystem(ds))[i] || term2 > DSGMASystemSignature(DSDesignSpaceGMASystem(ds))[i+1])
                break;
            if (term1 <= 0 || term2 <= 0)
                break;
        }
//        printf("****4 signatures created \n");
        DSSubCaseGenerate3dSignature(ds, aCase, three_digit, location);
        if (i == 2*numberOfEquations) {
                dsCaseCreateConditionMatrices(aCase, DSDesignSpaceGMASystem(ds));
                dsCaseAppendDesignSpaceConditions(aCase, ds);
                dsCaseCalculateCaseNumber(aCase, DSDesignSpaceGMASystem(ds), endian);
                dsSSystemIntegrateConservedRelationships(aCase, ds);
                dsCaseCreateBoundaryMatrices(aCase);
//                dsCaseCalculateCaseNumber(aCase, DSDesignSpaceGMASystem(ds), endian);
                dsCaseCalculateCaseIdentifier(aCase, DSDesignSpaceGMASystem(ds),
                                              endian, DSDesignSpaceCasePrefix(ds));
//                printf("****4 case identifier created \n");
                sprintf(name, "%i", aCase->caseNumber);
                if (DSSSystemHasSolution(DSCaseSSys(aCase)) == false &&
                    DSDesignSpaceUnstable(ds) == true &&
                    DSDictionaryValueForName(ds->unstableCases, name) == NULL) {
                    
                        if (DSuCaseIsValid(aCase, true) == true){
                                
                                uCase = DSUnstableCaseAddBoundaryMatrices(aCase, ds);
                                // Now save uCase in the dictionary DSDesignSpace->unstableCases
                                if (DSDictionaryValueForName(ds->unstableCases, name) == NULL){
                                        DSDictionaryAddValueWithName(ds->unstableCases, name, uCase);
                                }else{
                                    if (uCase != NULL )
                                        DSuCaseFree(uCase);
                                }
                            
                                DSSSystemSetIsUnstable((DSSSystem *)DSCaseSSystem(aCase), true);
                        }
                    
                }
//            printf("****4 matrices and other elements created \n");
        } else {
                if (aCase != NULL)
                    DSCaseFree(aCase);
                aCase = NULL;
        }
bail:
//    printf("****4 Finished function DSSSystemWithTermsFromGMACyclical for case %u (%s) \n", DSCaseNumberForSignature(termArray, ds->gma), aux);
//    if (aCase != NULL){
//        printf("The equations for this case are: \n");
//        DSCasePrintEquations(aCase);
//    }
        return aCase;
}


#if defined (__APPLE__) && defined (__MACH__)
#pragma mark - Getter functions
#endif

extern const bool DSCaseHasSolution(const DSCase *aCase)
{
        bool hasSolution = false;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        hasSolution = (DSCaseU(aCase) != NULL);
bail:
        return hasSolution;
}

extern const DSUInteger DSCaseNumberOfEquations(const DSCase *aCase)
{
        DSUInteger numberOfEquations = 0;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        numberOfEquations = DSSSystemNumberOfEquations(DSCaseSSys(aCase));
bail:
        return numberOfEquations;
}

extern const DSUInteger DSCaseNumberOfConservations(const DSCase *aCase)
{
    DSUInteger numberOfConservations = 0;
    if (aCase == NULL) {
        DSError(M_DS_CASE_NULL, A_DS_ERROR);
        goto bail;
    }
    numberOfConservations = DSSSystemNumberOfConservations(DSCaseSSys(aCase));
bail:
    return numberOfConservations;
}

extern const DSUInteger DSCaseNumberOfInheritedConservations(const DSCase *aCase)
{
    DSUInteger numberOfInheritedConservations = 0;
    if (aCase == NULL) {
        DSError(M_DS_CASE_NULL, A_DS_ERROR);
        goto bail;
    }
    numberOfInheritedConservations = aCase->numberInheritedConservations;
bail:
    return numberOfInheritedConservations;
}


extern DSExpression ** DSCaseEquations(const DSCase *aCase)
{
        DSExpression **equations = NULL;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (DSCaseSSys(aCase) == NULL) {
                DSError(M_DS_SSYS_NULL, A_DS_ERROR);
                goto bail;
        }
        equations = DSSSystemEquations(DSCaseSSys(aCase));
bail:
        return equations;
}

extern DSExpression ** DSCaseSolution(const DSCase *aCase)
{
        DSExpression **solution = NULL;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (DSCaseSSys(aCase) == NULL) {
                DSError(M_DS_SSYS_NULL, A_DS_ERROR);
                goto bail;
        }
        solution = DSSSystemSolution(DSCaseSSys(aCase));
bail:
        return solution;
}

extern DSExpression ** DSCaseLogarithmicSolution(const DSCase * aCase)
{
        DSExpression **solution = NULL;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (DSCaseSSys(aCase) == NULL) {
                DSError(M_DS_SSYS_NULL, A_DS_ERROR);
                goto bail;
        }
        solution = DSSSystemLogarithmicSolution(DSCaseSSys(aCase));
bail:
        return solution;
}

extern const DSUInteger DSCaseNumberOfConditions(const DSCase *aCase)
{
        DSUInteger numberOfConditions = 0;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (DSCaseCd(aCase) == NULL) {
                goto bail;
        }
        numberOfConditions = DSMatrixRows(DSCaseCd(aCase));
bail:
        return numberOfConditions;
}

extern const DSUInteger DSCaseNumberOfBoundaries(const DSCase *aCase)
{
        DSUInteger numberOfConditions = 0;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (DSCaseU(aCase) == NULL) {
                goto bail;
        }
        numberOfConditions = DSMatrixRows(DSCaseU(aCase));
bail:
        return numberOfConditions;
}

static void dsCaseConditionToString(const DSCase *aCase, 
                                    const DSUInteger condition, 
                                    char ** string, 
                                    DSUInteger *length, const bool inLog)
{
        DSUInteger i, numberOfXd, numberOfXi;
        char tempString[100];
        const char *name;
        double value;
        if (aCase == NULL) {
                DSError(M_DS_SSYS_NULL, A_DS_ERROR);
                goto bail;
        }
        numberOfXd = DSVariablePoolNumberOfVariables(DSCaseXd(aCase));
        if (condition >= DSMatrixRows(DSCaseCd(aCase))) {
                DSError("Equation does not exist: Check number of equations", A_DS_ERROR);
                goto bail;
        }
        if (string == NULL) {
                DSError(M_DS_NULL ": Pointer to string is NULL", A_DS_ERROR);
                goto bail;
        }
        if (length == NULL) {
                DSError(M_DS_NULL ": Pointer to length is NULL", A_DS_ERROR);
                goto bail;
        }
        if (*string == NULL) {
                DSError(M_DS_NULL ": String should be initialized", A_DS_ERROR);
                goto bail;
        }
        numberOfXi = DSVariablePoolNumberOfVariables(DSCaseXi(aCase));
        if (inLog == true) 
                sprintf(tempString, "%lf", DSMatrixDoubleValue(DSCaseDelta(aCase), condition, 0));
        else
                sprintf(tempString, "10^%lf", DSMatrixDoubleValue(DSCaseDelta(aCase), condition, 0));
        if (*length-strlen(*string) < 100) {
                *length += 1000;
                *string = DSSecureRealloc(*string, sizeof(char)**length);
        }
        strncat(*string, tempString, *length-strlen(*string));
        for (i = 0; i < numberOfXi+numberOfXd; i++) {
                if (*length-strlen(*string) < 100) {
                        *length += 1000;
                        *string = DSSecureRealloc(*string, sizeof(char)**length);
                }
                if (i < numberOfXi) {
                        name = DSVariableName(DSVariablePoolAllVariables(DSCaseXi(aCase))[i]);
                        value = DSMatrixDoubleValue(DSCaseCi(aCase), condition, i);
                } else {
                        name = DSVariableName(DSVariablePoolAllVariables(DSCaseXd(aCase))[i-numberOfXi]);
                        value = DSMatrixDoubleValue(DSCaseCd(aCase), condition, i-numberOfXi);
                }
                if (value == 0.0)
                        continue;
                if (inLog == true)
                        sprintf(tempString, "+%lf*log(%s)", value, name);
                else if (value == 1.0)
                        sprintf(tempString, "*%s", name);
                else
                        sprintf(tempString, "*%s^%lf", name, value);
                strncat(*string, tempString, *length-strlen(*string));
        }
bail:
        return;
}

extern DSExpression ** DSCaseConditions(const DSCase *aCase)
{
        DSUInteger i, numberOfConditions, length;
        DSExpression ** conditions = NULL;
        char *tempString, * equationString;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
    
        if (DSCaseDelta(aCase) == NULL || DSCaseCi(aCase) == NULL || DSCaseCd(aCase) == NULL)
            goto bail;

        numberOfConditions = DSMatrixRows(DSCaseCd(aCase));
        if (numberOfConditions == 0) {
                DSError("Case being accessed has no conditions", A_DS_ERROR);
                goto bail;
        }
        conditions = DSSecureCalloc(sizeof(DSExpression *), numberOfConditions);
        length = 1000;
        tempString = DSSecureCalloc(sizeof(char), length);
        for (i = 0; i < numberOfConditions; i++) {
                tempString[0] = '\0';
                dsCaseConditionToString(aCase, i, &tempString, &length, false);
                equationString = DSSecureCalloc(sizeof(char),
                                                strlen(tempString)+6);
                equationString = strcpy(equationString, tempString);
                equationString = strcat(equationString, " > 1");
                conditions[i] = DSExpressionByParsingString(equationString);
                DSSecureFree(equationString);
        }
        DSSecureFree(tempString);
bail:
        return conditions;
}

extern DSExpression ** DSCaseLogarithmicConditions(const DSCase *aCase)
{
        DSUInteger i, numberOfConditions, length;
        DSExpression ** conditions = NULL;
        char *tempString, * equationString;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (DSCaseCd(aCase) == NULL) {
                goto bail;
        }
        numberOfConditions = DSMatrixRows(DSCaseCd(aCase));
        if (numberOfConditions == 0) {
                DSError("Case being accessed has no conditions", A_DS_ERROR);
                goto bail;
        }
        conditions = DSSecureCalloc(sizeof(DSExpression *), numberOfConditions);
        length = 1000;
        tempString = DSSecureCalloc(sizeof(char), length);
        for (i = 0; i < numberOfConditions; i++) {
                tempString[0] = '\0';
                dsCaseConditionToString(aCase, i, &tempString, &length, true);
                equationString = DSSecureCalloc(sizeof(char),
                                                strlen(tempString)+6);
                equationString = strcpy(equationString, tempString);
                equationString = strcat(equationString, " > 0");
                conditions[i] = DSExpressionByParsingString(equationString);
                DSSecureFree(equationString);
        }
        DSSecureFree(tempString);
bail:
        return conditions;
}

static void dsCaseBoundaryToString(const DSCase *aCase, 
                                   const DSUInteger boundary, 
                                   char ** string, 
                                   DSUInteger *length, const bool inLog)
{
        DSUInteger i, numberOfXi;
        char tempString[100];
        const char *name;
        double value;
        if (aCase == NULL) {
                DSError(M_DS_SSYS_NULL, A_DS_ERROR);
                goto bail;
        }
        if (DSCaseU(aCase) == NULL) {
                goto bail;
        }
        if (boundary >= DSMatrixRows(DSCaseU(aCase))) {
                DSError("Equation does not exist: Check number of equations", A_DS_ERROR);
                goto bail;
        }
        if (string == NULL) {
                DSError(M_DS_NULL ": Pointer to string is NULL", A_DS_ERROR);
                goto bail;
        }
        if (length == NULL) {
                DSError(M_DS_NULL ": Pointer to length is NULL", A_DS_ERROR);
                goto bail;
        }
        if (*string == NULL) {
                DSError(M_DS_NULL ": String should be initialized", A_DS_ERROR);
                goto bail;
        }
        numberOfXi = DSVariablePoolNumberOfVariables(DSCaseXi(aCase));
        if (inLog == true) 
                sprintf(tempString, "%lf", DSMatrixDoubleValue(DSCaseZeta(aCase), boundary, 0));
        else
                sprintf(tempString, "10^%lf", DSMatrixDoubleValue(DSCaseZeta(aCase), boundary, 0));
        if (*length-strlen(*string) < 100) {
                *length += 1000;
                *string = DSSecureRealloc(*string, sizeof(char)**length);
        }
        strncat(*string, tempString, *length-strlen(*string));
        for (i = 0; i < numberOfXi; i++) {
                if (*length-strlen(*string) < 100) {
                        *length += 1000;
                        *string = DSSecureRealloc(*string, sizeof(char)**length);
                }
                name = DSVariableName(DSVariablePoolAllVariables(DSCaseXi(aCase))[i]);
                value = DSMatrixDoubleValue(DSCaseU(aCase), boundary, i);
                if (value == 0.0)
                        continue;
                if (inLog == true)
                        sprintf(tempString, "+%lf*log(%s)", value, name);
                else if (value == 1.0)
                        sprintf(tempString, "*%s", name);
                else
                        sprintf(tempString, "*%s^%lf", name, value);
                strncat(*string, tempString, *length-strlen(*string));
        }
bail:
        return;
}

extern DSExpression ** DSCaseBoundaries(const DSCase *aCase)
{
        DSUInteger i, numberOfConditions, length;
        DSExpression ** boundaries = NULL;
        char *tempString, * equationString;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (DSCaseU(aCase) == NULL) {
                goto bail;
        }
        numberOfConditions = DSMatrixRows(DSCaseU(aCase));
        if (numberOfConditions == 0) {
                DSError("Case being accessed has no conditions", A_DS_ERROR);
                goto bail;
        }
        boundaries = DSSecureCalloc(sizeof(DSExpression *), numberOfConditions);
        length = 1000;
        tempString = DSSecureCalloc(sizeof(char), length);
        for (i = 0; i < numberOfConditions; i++) {
                tempString[0] = '\0';
                dsCaseBoundaryToString(aCase, i, &tempString, &length, false);
                equationString = DSSecureCalloc(sizeof(char),
                                                strlen(tempString)+6);
                equationString = strcpy(equationString, tempString);
                equationString = strcat(equationString, " > 1");
                boundaries[i] = DSExpressionByParsingString(equationString);
                DSSecureFree(equationString);
        }
        DSSecureFree(tempString);
bail:
        return boundaries;
}

extern DSExpression ** DSCaseLogarithmicBoundaries(const DSCase *aCase)
{
        DSUInteger i, numberOfConditions, length;
        DSExpression ** boundaries = NULL;
        char *tempString, *equationString;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (DSCaseU(aCase) == NULL) {
                goto bail;
        }
        numberOfConditions = DSMatrixRows(DSCaseU(aCase));
        if (numberOfConditions == 0) {
                goto bail;
        }
        boundaries = DSSecureCalloc(sizeof(DSExpression *), numberOfConditions);
        length = 1000;
        tempString = DSSecureCalloc(sizeof(char), length);
        for (i = 0; i < numberOfConditions; i++) {
                tempString[0] = '\0';
                dsCaseBoundaryToString(aCase, i, &tempString, &length, true);
                equationString = DSSecureCalloc(sizeof(char),
                                                strlen(tempString)+6);
                equationString = strcpy(equationString, tempString);
                equationString = strcat(equationString, " > 0");
                boundaries[i] = DSExpressionByParsingString(equationString);
                DSSecureFree(equationString);
        }
        DSSecureFree(tempString);
bail:
        return boundaries;
}

extern DSUInteger DSCaseNumber(const DSCase * aCase)
{
        DSUInteger caseNumber = 0;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        caseNumber = DSCaseNum(aCase);
bail:
        return caseNumber;
}

extern const char * DSCaseIdentifier(const DSCase * aCase)
{
        const char * caseIdentifier = NULL;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        caseIdentifier = DSCaseId(aCase);
bail:
        return caseIdentifier;
}

extern const DSUInteger * DSCaseSignature(const DSCase * aCase)
{
        const DSUInteger *signature = NULL;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        signature = DSCaseSig(aCase);
bail:
        return signature;
}


extern const DSSSystem *DSCaseSSystem(const DSCase * aCase)
{
        DSSSystem * ssys = NULL;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        ssys = DSCaseSSys(aCase);
bail:
        return ssys;
}

extern double DSCaseLogarithmicGain(const DSCase *aCase, const char *XdName, const char *XiName)
{
        double logGain = INFINITY;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        logGain = DSSSystemLogarithmicGain(DSCaseSSys(aCase), XdName, XiName);
bail:
        return logGain;
}


extern DSStack * DSCaseVertexEquationsFor2DSlice(const DSCase *aCase, const DSVariablePool * lowerBounds, const DSVariablePool *upperBounds, const char * xVariable, const char *yVariable, const bool log_out)
{
        DSStack * equations = NULL;
        DSVertices *vertices = NULL;
        DSUInteger yIndex, xIndex;
        DSMatrix *U, *Zeta, *temp, * vars;
        DSMatrix * solution, *Us, *Up, *Um;
        DSUInteger i, j, k, count, * activeRows;
        const char * variables[2];
        DSVariablePool * Xd, * Xi;
        char * name, * string, * new;
        DSExpression ** expressions, *rhs;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        vertices = DSCaseVerticesFor2DSlice(aCase, lowerBounds, upperBounds, xVariable, yVariable);
        if (vertices == NULL) {
                goto bail;
        }
        yIndex = DSVariablePoolIndexOfVariableWithName(DSCaseXi(aCase), yVariable);
        xIndex = DSVariablePoolIndexOfVariableWithName(DSCaseXi(aCase), xVariable);
        temp = DSMatrixCalloc(4, DSMatrixColumns(DSCaseU(aCase)));
        DSMatrixSetDoubleValue(temp, 0, xIndex, 1.0);
        DSMatrixSetDoubleValue(temp, 1, xIndex, -1.0);
        DSMatrixSetDoubleValue(temp, 2, yIndex, 1.0);
        DSMatrixSetDoubleValue(temp, 3, yIndex, -1.0);
        U = DSMatrixAppendMatrices(DSCaseU(aCase), temp, false);
        DSMatrixFree(temp);
        temp = DSMatrixCalloc(4, 1);
        DSMatrixSetDoubleValue(temp, 0, 0, -log10(DSVariableValue(DSVariablePoolVariableWithName(lowerBounds, xVariable))));
        DSMatrixSetDoubleValue(temp, 1, 0, log10(DSVariableValue(DSVariablePoolVariableWithName(upperBounds, xVariable))));
        DSMatrixSetDoubleValue(temp, 2, 0, -log10(DSVariableValue(DSVariablePoolVariableWithName(lowerBounds, yVariable))));
        DSMatrixSetDoubleValue(temp, 3, 0, log10(DSVariableValue(DSVariablePoolVariableWithName(upperBounds, yVariable))));
        Zeta = DSMatrixAppendMatrices(DSCaseZeta(aCase), temp, false);
        DSMatrixFree(temp);
        equations = DSStackAlloc();
        vars = DSVariablePoolValuesAsVector(lowerBounds, false);
        DSMatrixApplyFunction(vars, log10);
        Xd = DSVariablePoolAlloc();
        Xi = DSVariablePoolAlloc();
        count = 2;
        activeRows = DSSecureCalloc(sizeof(DSUInteger), count);
        for (i = 0; i < DSVariablePoolNumberOfVariables(DSCaseXi(aCase)); i++) {
                if (i == xIndex || i == yIndex)
                        continue;
                name = DSVariableName(DSVariablePoolVariableAtIndex(DSCaseXi(aCase), i));
                DSVariablePoolAddVariableWithName(Xi, name);
        }
        if (xIndex < yIndex) {
                variables[0] = xVariable;
                variables[1] = yVariable;
        } else {
                variables[0] = yVariable;
                variables[1] = xVariable;
        }
        for (i = 0; i < vertices->numberOfVertices; i++) {
                DSMatrixSetDoubleValue(vars, xIndex, 0, vertices->vertices[i][0]);
                DSMatrixSetDoubleValue(vars, yIndex, 0, vertices->vertices[i][1]);
                solution = DSMatrixByMultiplyingMatrix(U, vars);
                DSMatrixAddByMatrix(solution, Zeta);
                k = 0;
                for (j = 0; j < DSMatrixRows(solution); j++) {
                        if (DSMatrixDoubleValue(solution, j, 0) < 1e-5) {
                                if (k == count) {
                                        count++;
                                        activeRows = DSSecureRealloc(activeRows, sizeof(DSUInteger)*count);
                                }
                                activeRows[k++] = j;
                        }
                }
                DSMatrixFree(solution);
                temp = DSMatrixSubMatrixIncludingRows(U, k, activeRows);
                if (k > 2) {
                        Us = DSMatrixSubMatrixIncludingColumnList(temp, 2, xIndex, yIndex);
                        for (j = 0; j < DSMatrixColumns(Us); j++) {
                                for (k = 0; k < DSMatrixRows(Us); k++) {
                                        if (DSMatrixDoubleValue(Us, k, j) != 0) {
                                                if (j == 1) {
                                                        if (activeRows[j-1] == k)
                                                                continue;
                                                }
                                                activeRows[j] = k;
                                        }
                                }
                        }
                        Us = DSMatrixSubMatrixIncludingRows(temp, 2, activeRows);
                        DSMatrixFree(temp);
                        temp = Us;
                }
                Us = DSMatrixSubMatrixIncludingColumnList(temp, 2, xIndex, yIndex);
                Up = DSMatrixSubMatrixExcludingColumnList(temp, 2, xIndex, yIndex);
                DSMatrixFree(temp);
                Um = DSMatrixInverse(Us);
                DSMatrixFree(Us);
                Us = DSMatrixByMultiplyingMatrix(Um, Up);
                DSMatrixFree(Up);
                DSMatrixMultiplyByScalar(Us, -1.);
                temp = DSMatrixSubMatrixIncludingRows(Zeta, 2, activeRows);
                Up = DSMatrixByMultiplyingMatrix(Um, temp);
                DSMatrixFree(temp);
                temp = Up;
                DSMatrixMultiplyByScalar(temp, -1.f);
                expressions = DSSecureCalloc(sizeof(DSExpression *), 2);
                for (k = 0; k < 2; k++) {
                        DSMatrixSetDoubleValue(temp, k, 0, pow(10, DSMatrixDoubleValue(temp, k, 0)));
                        if (log_out == false) {
                                rhs = DSExpressionFromPowerlawInMatrixForm(k, NULL, Xd, Us, Xi, Up);
                                string = DSExpressionAsString(rhs);
                                new = DSSecureCalloc(sizeof(char), strlen(string)+1000);
                                sprintf(new, "%s = %s", variables[k], string);
                                DSSecureFree(string);
                                DSExpressionFree(rhs);
                                expressions[k] = DSExpressionByParsingString(new);
                                DSSecureFree(new);
                        } else {
                                rhs = DSExpressionFromLogPowerlawInMatrixForm(k, NULL, Xd, Us, Xi, Up);
                                string = DSExpressionAsString(rhs);
                                new = DSSecureCalloc(sizeof(char), strlen(string)+1000);
                                sprintf(new, "log(%s) = %s", variables[k], string);
                                DSSecureFree(string);
                                DSExpressionFree(rhs);
                                expressions[k] = DSExpressionByParsingString(new);
                                DSSecureFree(new);
                        }
                }
                DSStackPush(equations, expressions);
                DSMatrixFree(temp);
                DSMatrixFree(Us);
        }
        DSVariablePoolFree(Xd);
        DSVariablePoolFree(Xi);
        DSMatrixFree(vars);
        DSMatrixFree(Zeta);
        DSMatrixFree(U);
bail:
        return equations;
}

#if defined (__APPLE__) && defined (__MACH__)
#pragma mark - Utility functions
#endif

extern DSMatrix * DSCaseDoubleValueBoundariesAtPoint(const DSCase * aCase, const DSVariablePool * point)
{
        DSMatrix * values = NULL;
        DSMatrix * U, *Zeta, *Xi;
        U = DSCaseU(aCase);
        Zeta = DSCaseZeta(aCase);
        Xi = DSVariablePoolValuesAsVector(point, false);
        values = DSMatrixByMultiplyingMatrix(U, Xi);
        DSMatrixAddByMatrix(values, Zeta);
        return values;
}

extern DSMatrix * DSCaseDoubleValueBoundariesAtPointSortXi(const DSCase * aCase,
                                                           const DSVariablePool * point)
{
    DSMatrix * values = NULL;
    DSVariablePool *Xi0 = NULL;
    DSUInteger i;
    const DSSSystem *ssys = DSCaseSSystem(aCase);
    const char *name;
    
    if (DSVariablePoolNumberOfVariables(ssys->Xi) != 0) {
        Xi0 = DSVariablePoolAlloc();
        for (i=0; i < DSVariablePoolNumberOfVariables(ssys->Xi); i++) {
            name =  DSVariableName(DSVariablePoolAllVariables(ssys->Xi)[i]);
            DSVariablePoolAddVariableWithName(Xi0, name);
            if (DSVariablePoolHasVariableWithName(point, name) == false) {
                DSVariablePoolFree(Xi0);
                goto bail;
            }
            DSVariablePoolSetValueForVariableWithName(Xi0,
                                                      name,
                                                      log10(DSVariableValue(DSVariablePoolVariableWithName(point, name))));
        }
        values = DSCaseDoubleValueBoundariesAtPoint(aCase, Xi0);
        DSVariablePoolFree(Xi0);
    } else {
        values = DSMatrixCopy(DSCaseZeta(aCase));
    }
    
bail:
    return values;
}


extern const DSVariablePool * DSCaseXd(const DSCase * aCase)
{
        const DSVariablePool * Xd = NULL;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        Xd = aCase->Xd;
bail:
        return Xd;
}
extern const DSVariablePool * DSCaseXd_a(const DSCase * aCase)
{
        const DSVariablePool * Xd_a = NULL;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        Xd_a = aCase->Xd_a;
bail:
        return Xd_a;
}

extern const DSVariablePool * DSCaseXi(const DSCase * aCase)
{
        const DSVariablePool * Xi = NULL;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        Xi = aCase->Xi;
bail:
        return Xi;
}

extern DSUInteger DSCaseNumberOfMassBalances(const DSCase *aCase){
    DSUInteger n = 0;
    
    if (aCase == NULL) {
            DSError(M_DS_CASE_NULL, A_DS_ERROR);
            goto bail;
    }
    if (DSCaseDominantMassBalances(aCase) == NULL)
        goto bail;
    
    n = DSCaseDominantMassBalances(aCase)->n;
    
bail:
    return n;
}

extern const char * DSCaseDominantFinAtIndex(const DSCase *aCase, DSUInteger i){
    char *fin = NULL;
    
    if (aCase == NULL) {
            DSError(M_DS_CASE_NULL, A_DS_ERROR);
            goto bail;
    }
    if (DSCaseDominantMassBalances(aCase) == NULL || i > DSCaseDominantMassBalances(aCase)->n )
        goto bail;
    
    fin = DSCaseDominantMassBalances(aCase)->fin[i];

bail:
    return fin;
    
}

extern const char * DSCaseDominantFoutAtIndex(const DSCase *aCase, DSUInteger i){
    char * fout = NULL;
    
    if (aCase == NULL) {
            DSError(M_DS_CASE_NULL, A_DS_ERROR);
            goto bail;
    }
    if (DSCaseDominantMassBalances(aCase) == NULL || i > DSCaseDominantMassBalances(aCase)->n)
        goto bail;
    
    fout = DSCaseDominantMassBalances(aCase)->fout[i];
    
bail:
    return fout;
    
}


#if defined (__APPLE__) && defined (__MACH__)
#pragma mark Additional Constraints
#endif

static void dsCaseAddBoundariesFromConditions(DSCase *aCase, const DSMatrix * Cd, const DSMatrix * Ci, const DSMatrix * delta)
{
        DSUInteger numberOfXi = 0;
        DSMatrix * W = NULL, *Zeta, *U = NULL, *B, *Ai, *temp;
        bool hasSSys;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
//        if (DSSSystemHasSolution(DSCaseSSys(aCase)) == false) {
//                goto bail;
//        }
        if (Cd == NULL) {
                goto bail;
        }
        numberOfXi = DSVariablePoolNumberOfVariables(DSCaseXi(aCase));
        hasSSys = DSCaseSSys(aCase) != NULL;
        if (hasSSys == false) {
                W = DSMatrixCalloc(DSMatrixRows(Cd), DSMatrixColumns(Cd));
                Zeta = DSMatrixCopy(delta);
                U = DSMatrixCopy(Ci);
        } else if (DSSSystemHasSolution(DSCaseSSys(aCase)) == false) {
                W = DSMatrixCalloc(DSMatrixRows(Cd), DSMatrixColumns(Cd));
                Zeta = DSMatrixCopy(delta);
                U = DSMatrixCopy(Ci);
        } else {
                
                W = DSMatrixByMultiplyingMatrix(Cd, DSSSystemM(DSCaseSSys(aCase)));
                B = DSSSystemB(DSCaseSSys(aCase));
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
                DSMatrixFree(B);
        }
        temp = DSCaseZeta(aCase);
        DSCaseZeta(aCase) = DSMatrixAppendMatrices(temp, Zeta, false);
        DSMatrixFree(temp);
        DSMatrixFree(Zeta);
        temp = DSCaseU(aCase);
        DSCaseU(aCase) = DSMatrixAppendMatrices(temp, U, false);
        DSMatrixFree(temp);
        if (U != NULL)
                DSMatrixFree(U);
        DSMatrixFree(W);
//        DSCaseRemoveRedundantBoundaries(aCase);
bail:
        return;
}

static void dsCaseAddConditions(DSCase *aCase, const DSMatrix * Cd, const DSMatrix * Ci, const DSMatrix * delta)
{
        DSMatrix *temp = NULL;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (Cd == NULL) {
                DSError(M_DS_MAT_NULL ": Cd is NULL", A_DS_ERROR);
                goto bail;
        }
        if (delta == NULL) {
                DSError(M_DS_MAT_NULL ": Delta is NULL", A_DS_ERROR);
                goto bail;
        }
        if (Ci == NULL && DSVariablePoolNumberOfVariables(DSCaseXi(aCase)) != 0) {
                DSError(M_DS_MAT_NULL ": Ci is NULL", A_DS_ERROR);
                goto bail;
        }
        if (DSMatrixColumns(Cd) != DSVariablePoolNumberOfVariables(DSCaseXd(aCase))) {
                DSError(M_DS_WRONG ": Number of dep. variables must match number of columns of Cd", A_DS_ERROR);
                goto bail;
        }
        if (Ci != NULL) {
                if (DSMatrixColumns(Ci) != DSVariablePoolNumberOfVariables(DSCaseXi(aCase))) {
                        DSError(M_DS_WRONG ": Number of indep. variables must match number of columns of Ci", A_DS_ERROR);
                        goto bail;
                }
                if (DSMatrixRows(Cd) != DSMatrixRows(Ci)) {
                        DSError(M_DS_WRONG ": Rows of Ci must match rows of Cd", A_DS_ERROR);
                        goto bail;
                }
        }
        if (DSMatrixRows(Cd) != DSMatrixRows(delta)) {
                DSError(M_DS_WRONG ": Rows of Cd must match rows of delta", A_DS_ERROR);
                goto bail;
        }
        if (DSCaseCd(aCase) == NULL) {
                DSCaseCd(aCase) = DSMatrixCopy(Cd);
                DSCaseDelta(aCase) = DSMatrixCopy(delta);
                if (Ci != NULL)
                        DSCaseCi(aCase) = DSMatrixCopy(Ci);
        } else {
                temp = DSMatrixAppendMatrices(DSCaseCd(aCase), Cd, false);
                DSMatrixFree(DSCaseCd(aCase));
                DSCaseCd(aCase) = temp;
                temp = DSMatrixAppendMatrices(DSCaseDelta(aCase), delta, false);
                DSMatrixFree(DSCaseDelta(aCase));
                DSCaseDelta(aCase) = temp;
                if (Ci != NULL) {
                        temp = DSMatrixAppendMatrices(DSCaseCi(aCase), Ci, false);
                        DSMatrixFree(DSCaseCi(aCase));
                        DSCaseCi(aCase) = temp;
                }
        }
bail:
        return;
}

static void dsCaseConstraintsProcessExponentBasePairs(const DSCase *aCase, gma_parseraux_t *current, DSInteger sign,
                                                             DSUInteger index, DSMatrix * Cd, DSMatrix * Ci, DSMatrix *delta)
{
        DSUInteger j, varIndex;
        const char *varName;
        double currentValue;
        if (current == NULL) {
                goto bail;
        }
        if (sign == AUX_SIGN_NEGATIVE) {
                sign = -1;
        } else {
                sign = 1;
        }
        for (j = 0; j < DSGMAParserAuxNumberOfBases(current); j++) {
                if (DSGMAParserAuxBaseAtIndexIsVariable(current, j) == false) {
                        currentValue = DSMatrixDoubleValue(delta, index, 0);
                        currentValue += sign * log10(DSGMAParseAuxsConstantBaseAtIndex(current, j));
                        DSMatrixSetDoubleValue(delta,
                                               index, 0,
                                               currentValue);
                        continue;
                }
                varName = DSGMAParserAuxVariableAtIndex(current, j);
                if (DSVariablePoolHasVariableWithName(DSCaseXd(aCase), varName) == true) {
                        varIndex = DSVariablePoolIndexOfVariableWithName(DSCaseXd(aCase), varName);
                        currentValue = DSMatrixDoubleValue(Cd, index, varIndex);
                        currentValue += sign * DSGMAParserAuxExponentAtIndex(current, j);
                        DSMatrixSetDoubleValue(Cd, index, varIndex, currentValue);
                } else if (DSVariablePoolHasVariableWithName(DSCaseXi(aCase), varName) == true) {
                        varIndex = DSVariablePoolIndexOfVariableWithName(DSCaseXi(aCase), varName);
                        currentValue = DSMatrixDoubleValue(Ci, index, varIndex);
                        currentValue += sign * DSGMAParserAuxExponentAtIndex(current, j);
                        DSMatrixSetDoubleValue(Ci, index, varIndex, currentValue);
                }
        }
bail:
        return;
}

static void dsCaseConstraintsCreateSystemMatrices(DSCase *aCase, DSUInteger numberOfConstraints, gma_parseraux_t **aux)
{
        gma_parseraux_t *current;
        DSUInteger i;
        DSMatrix * Cd, *Ci, *delta;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL ": Case being modified is NULL", A_DS_ERROR);
                goto bail;
        }
        if (aux == NULL) {
                DSError(M_DS_NULL ": Parser auxiliary data is NULL", A_DS_ERROR);
                goto bail;
        }
        if (DSCaseXd(aCase) == NULL || DSCaseXi(aCase) == NULL) {
                DSError(M_DS_WRONG ": GMA data is incomplete: Need Xi and Xd", A_DS_ERROR);
                goto bail;
        }
        Cd = DSMatrixCalloc(numberOfConstraints, DSVariablePoolNumberOfVariables(DSCaseXd(aCase)));
        Ci = DSMatrixCalloc(numberOfConstraints, DSVariablePoolNumberOfVariables(DSCaseXi(aCase)));
        delta = DSMatrixCalloc(numberOfConstraints, 1);
        for (i = 0; i < numberOfConstraints; i++) {
                current = aux[i];
                dsCaseConstraintsProcessExponentBasePairs(aCase, current, current->sign, i, Cd, Ci, delta);
                current = DSGMAParserAuxNextNode(current);
                dsCaseConstraintsProcessExponentBasePairs(aCase, current, current->sign, i, Cd, Ci, delta);
        }
        dsCaseAddConditions(aCase, Cd, Ci, delta);
        dsCaseAddBoundariesFromConditions(aCase, Cd, Ci, delta);
        DSMatrixFree(Cd);
        DSMatrixFree(Ci);
        DSMatrixFree(delta);
bail:
        return;
}

extern void DSCaseAddConstraints(DSCase * aCase, const char ** strings, DSUInteger numberOfConstraints)
{
        DSUInteger i;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        gma_parseraux_t **aux = NULL;
        aux = (gma_parseraux_t **)DSDesignSpaceTermListForAllStrings(strings, numberOfConstraints);
        if (aux == NULL) {
                goto bail;
        }
        dsCaseConstraintsCreateSystemMatrices(aCase, numberOfConstraints, aux);
        for (i=0; i < numberOfConstraints; i++) {
                if (aux[i] != NULL)
                        DSGMAParserAuxFree(aux[i]);
        }
        DSSecureFree(aux);
bail:
        return;
}


#if defined (__APPLE__) && defined (__MACH__)
#pragma mark Case signature and case number
#endif

extern DSUInteger * DSCaseSignatureForCaseNumber(const DSUInteger caseNumber, const DSGMASystem * gma)
{
        DSUInteger *signature = NULL;
        DSUInteger num;
        DSInteger i;
        DSUInteger numberOfEquations;
        if (gma == NULL) {
                DSError(M_DS_GMA_NULL, A_DS_ERROR);
                goto bail;
        }
        if (caseNumber == 0) {
                DSError(M_DS_WRONG ": Case number is 0", A_DS_ERROR);
                goto bail;
        }
        if (caseNumber > DSGMASystemNumberOfCases(gma)) {
            printf("The case number you are trying to access is %u and the number of cases is %u \n", caseNumber, DSGMASystemNumberOfCases(gma) );
                DSError(M_DS_WRONG ": Case number is out of bounds", A_DS_ERROR);
                goto bail;
        }
        numberOfEquations = DSGMASystemNumberOfEquations(gma);
        signature = DSSecureMalloc(sizeof(DSUInteger)*2*numberOfEquations);
        num = caseNumber-1;
        switch (endian) {
                case DS_CASE_NUMBER_SMALL_ENDIAN:
                        for (i = 0; i < 2*numberOfEquations; i++) {
                                signature[i] = num % DSGMASystemSignature(gma)[i] +1;
                                num /= DSGMASystemSignature(gma)[i];
                        }
                        break;
                case DS_CASE_NUMBER_BIG_ENDIAN:
                default:
                        for (i = 2*numberOfEquations-1; i >= 0; i--) {
                                signature[i] = num % DSGMASystemSignature(gma)[i] +1;
                                num /= DSGMASystemSignature(gma)[i];
                        }
                        break;
        }
bail:
        return signature;
}

extern const DSUIntegerVector * DSCaseGetSignatureNeighbors(const DSCase *aCase, const DSDesignSpace *ds)
{
    
    DSUInteger *neighbors_array = NULL, *case_sig = NULL, *neighbors_array_aux = NULL, *neighbor_sig = NULL;
    const DSUInteger * system_sig = NULL;
    DSUInteger max_nr_neighbors = 0, neighbor = 0, i, j, caseNr;
    DSCase *aNeighbor = NULL;
    const DSGMASystem *gma = DSDesignSpaceGMASystem(ds);
    DSUIntegerVector *neighbors_vector = NULL;
    
    caseNr = DSCaseNumber(aCase);
    
    /* Calculate the maximal number of neighbors  */
    system_sig = DSGMASystemSignature(gma);
    case_sig = DSCaseSignatureForCaseNumber(caseNr, gma);
    for(i=0; i<DSVariablePoolNumberOfVariables(DSGMASystemXd(gma))*2; i++){
        max_nr_neighbors += system_sig[i] - 1;
    }
        
    /* Allocate auxiliary vector to store integers for neighbors */
    neighbors_array_aux = DSSecureMalloc(sizeof(DSUInteger)*max_nr_neighbors);
    
    /* Generate neighbors*/
    for(i=0; i<DSVariablePoolNumberOfVariables(DSGMASystemXd(gma))*2; i++){
        
            for (j=1; j<=system_sig[i]; j++){
                    
                    if (case_sig[i] == j){
                        continue;
                    }
                    else{
                        neighbor_sig = DSCaseSignatureForCaseNumber(caseNr, gma);
                        neighbor_sig[i] = j;
                        
                        /* Create Neighbor case and see if it is valid. If valid, save and increase neighbor counter  */
                        aNeighbor = DSDesignSpaceCaseWithCaseSignature(ds, neighbor_sig);
                        if (DSCaseIsValid(aNeighbor, true) == true){
                            neighbors_array_aux[neighbor] = DSCaseNumber(aNeighbor);
                            neighbor++;
                        }
                        if (neighbor_sig != NULL)
                            DSSecureFree(neighbor_sig);
                        if (aNeighbor != NULL)
                            DSCaseFree(aNeighbor);
                    }
            }
    }
    
    if (neighbor == 0)
        goto bail;
    
    /* Allocate vector to store neighbors and delete variables */
    neighbors_array = DSSecureMalloc(sizeof(DSUInteger)*neighbor);
    neighbors_vector = DSSecureMalloc(sizeof(DSUIntegerVector));
    for (i=0; i<neighbor; i++){
        neighbors_array[i] = neighbors_array_aux[i];
    }
    neighbors_vector->vector = neighbors_array;
    neighbors_vector->dimension = neighbor;
    
bail:
    if (neighbors_array_aux != NULL)
        DSSecureFree(neighbors_array_aux);
    
    return neighbors_vector;
}

extern DSUInteger DSUIntegerVectorDimension(const DSUIntegerVector *aVector)
{
    DSUInteger dimension = 0;
    
    if (aVector == NULL)
        goto bail;
    
    dimension = aVector->dimension;
    
bail:
    return dimension;
}

extern DSUInteger DSUIntegerVectorValueAtIndex(const DSUIntegerVector *aVector, const DSUInteger index)
{
    DSUInteger value = 0;
    
    if (aVector == NULL || aVector == NULL)
        goto bail;
    
    value = aVector->vector[index];
    
bail:
    return value;
}


extern const DSUInteger DSCaseNumberForSignature(const DSUInteger * signature, const DSGMASystem * gma)
{
        DSInteger i;
        DSUInteger temp, numberOfEquations,  caseNumber = 0;
        if (signature == NULL) {
                DSError(M_DS_NULL ": Case Signature is NULL", A_DS_ERROR);
                goto bail;
        }
        if (gma == NULL) {
                DSError(M_DS_GMA_NULL ": Template GMA to make S-System is NULL", A_DS_ERROR);
                goto bail;
        }
        numberOfEquations = DSVariablePoolNumberOfVariables(DSGMASystemXd(gma));
        caseNumber = 1;
        temp = 1;
        switch (endian) {
                case DS_CASE_NUMBER_SMALL_ENDIAN:
                        for (i = 0; i < 2*numberOfEquations; i++) {
                                caseNumber += (signature[i]-1)*temp;
                                temp *= DSGMASystemSignature(gma)[i];
                        }
                        break;
                case DS_CASE_NUMBER_BIG_ENDIAN:
                default:
                        for (i = 2*numberOfEquations-1; i >= 0; i--) {
                                caseNumber += (signature[i]-1)*temp;
                                temp *= DSGMASystemSignature(gma)[i];
                        }
                        break;
        }
bail:
        return caseNumber;
}

extern char * DSCaseSignatureToString(const DSCase *aCase)
{
        char temp[100];
        char * string = NULL;
        char space[50];
        DSUInteger i, numberOfEquations;
        numberOfEquations = DSCaseNumberOfConservations(aCase) + DSCaseNumberOfEquations(aCase) + DSCaseNumberOfInheritedConservations(aCase);
    
//        if (strcmp(aCase->caseIdentifier, "6917_25_5") == 0 )
//        printf("Number of equations for case %u is %u. Number of equations is %u, number of conservations is %u, number of inherited conservations is %u \n", aCase->caseNumber, numberOfEquations, DSCaseNumberOfEquations(aCase), DSCaseNumberOfConservations(aCase),aCase->numberInheritedConservations );
    
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        string = DSSecureCalloc(sizeof(char), 5*numberOfEquations);
        strcpy(space, "  ");
        if (DSCaseSigCons(aCase) != NULL){
                    for (i = 0; i < 2*numberOfEquations; i++) {
                        if (DSCaseSigCons(aCase)[i] >= 10)
                            sprintf(temp, "(%i)", DSCaseSigCons(aCase)[i]);
                        else
                            sprintf(temp, "%i", DSCaseSigCons(aCase)[i]);
                        strncat(string, temp, 100-strlen(string));
                        if ((i+1)%2 == 0)
                            strncat(string, space, 100-strlen(space));
                    }
                    goto bail;
        }else if (DSCase3Sig(aCase) == NULL){
                    for (i = 0; i < 2*numberOfEquations; i++) {
                            if (DSCaseSig(aCase)[i] >= 10)
                                    sprintf(temp, "(%i)", DSCaseSig(aCase)[i]);
                            else
                                    sprintf(temp, "%i", DSCaseSig(aCase)[i]);
                            strncat(string, temp, 100-strlen(string));
                            if ((i+1)%2 == 0)
                                strncat(string, space, 100-strlen(space));
                    }
        } else {
                    for (i = 0; i < 3*numberOfEquations; i++) {
                        if (DSCase3Sig(aCase)[i] >= 10){
                                if (DSCase3Sig(aCase)[i] != 0){
                                    sprintf(temp, "(%i)", DSCase3Sig(aCase)[i]);
                                } else
                                    continue;
                        } else {
                            if (DSCase3Sig(aCase)[i] != 0){
                                    sprintf(temp, "%i", DSCase3Sig(aCase)[i]);
                            }else if ((i+2)%3 == 0 || (i+1)%3 == 0){
                                sprintf(temp, "%i", DSCase3Sig(aCase)[i]);
                            } else
                                continue;
                        }
                        strncat(string, temp, 100-strlen(string));
                        if ((i+1)%3 == 0)
                            strncat(string, space, 100-strlen(space));
                    }
        }
bail:
        return string;
}

#if defined (__APPLE__) && defined (__MACH__)
#pragma mark Printing functions
#endif

extern void DSCasePrint(const DSCase *aCase)
{
        DSUInteger i;
        int (*print)(const char *, ...);
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (DSPrintf == NULL)
                print = printf;
        else
                print = DSPrintf;
        print("\t  Case: %i\n\t   Sig: ",
              DSCaseNum(aCase));
        for (i = 0; i < DSCaseNumberOfEquations(aCase); i++) {
                if (DSCaseSig(aCase)[2*i] >= 10)
                        print("(");
                print("%i", DSCaseSig(aCase)[2*i]);
                if (DSCaseSig(aCase)[2*i] >= 10)
                        print(")");
                if (DSCaseSig(aCase)[2*i+1] >= 10)
                        print("(");
                print("%i", DSCaseSig(aCase)[2*i+1]);
                if (DSCaseSig(aCase)[2*i+1] >= 10)
                        print(")");
        }
        print("\n");
        DSSSystemPrint(DSCaseSSys(aCase));
bail:
        return;
}

extern void DSCasePrintEquations(const DSCase *aCase)
{
        DSUInteger i;
        DSExpression ** equations = NULL;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        equations = DSCaseEquations(aCase);
        if (equations != NULL) {
                for (i= 0; i < DSCaseNumberOfEquations(aCase); i++) {
                        DSExpressionPrint(equations[i]);
                        DSExpressionFree(equations[i]);
                }
                DSSecureFree(equations);
        }
bail:
        return;
}

extern void DSCasePrintSolution(const DSCase *aCase)
{
        DSUInteger i;
        DSExpression ** solution = NULL;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        solution = DSCaseSolution(aCase);
        if (solution != NULL) {
                for (i= 0; i < DSCaseNumberOfEquations(aCase); i++) {
                        DSExpressionPrint(solution[i]);
                        DSExpressionFree(solution[i]);
                }
                DSSecureFree(solution);
        }
bail:
        return;
}

extern void DSCasePrintLogarithmicSolution(const DSCase *aCase)
{
        DSUInteger i;
        DSExpression ** solution = NULL;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        solution = DSCaseLogarithmicSolution(aCase);
        if (solution != NULL) {
                for (i= 0; i < DSCaseNumberOfEquations(aCase); i++) {
                        DSExpressionPrint(solution[i]);
                        DSExpressionFree(solution[i]);
                }
                DSSecureFree(solution);
        }
bail:
        return;
}

extern void DSCasePrintConditions(const DSCase *aCase)
{
        DSUInteger i;
        DSExpression ** conditions = NULL;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        conditions = DSCaseConditions(aCase);
        if (conditions != NULL) {
                for (i= 0; i < DSCaseNumberOfConditions(aCase); i++) {
                        DSExpressionPrint(conditions[i]);
                        DSExpressionFree(conditions[i]);
                }
                DSSecureFree(conditions);
        }
bail:
        return;
}

extern void DSCasePrintLogarithmicConditions(const DSCase *aCase)
{
        DSUInteger i;
        DSExpression ** conditions = NULL;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        conditions = DSCaseLogarithmicConditions(aCase);
        if (conditions != NULL) {
                for (i= 0; i < DSCaseNumberOfConditions(aCase); i++) {
                        DSExpressionPrint(conditions[i]);
                        DSExpressionFree(conditions[i]);
                }
                DSSecureFree(conditions);
        }
bail:
        return;
}

extern void DSCasePrintBoundaries(const DSCase *aCase)
{
        DSUInteger i;
        DSExpression ** boundaries = NULL;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        boundaries = DSCaseBoundaries(aCase);
        if (boundaries != NULL) {
                for (i= 0; i < DSCaseNumberOfBoundaries(aCase); i++) {
                        printf("0 < ");
                        DSExpressionPrint(boundaries[i]);
                        DSExpressionFree(boundaries[i]);
                }
                DSSecureFree(boundaries);
        }
bail:
        return;
}

extern void DSCasePrintLogarithmicBoundaries(const DSCase *aCase)
{
        DSUInteger i;
        DSExpression ** boundaries = NULL;
        int (*print)(const char *, ...);
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        print = DSPrintf;
        if (print == NULL) 
                print = printf;
        boundaries = DSCaseLogarithmicBoundaries(aCase);
        if (boundaries != NULL) {
                for (i= 0; i < DSCaseNumberOfBoundaries(aCase); i++) {
                        print("0 < ");
                        DSExpressionPrint(boundaries[i]);
                        DSExpressionFree(boundaries[i]);
                }
                DSSecureFree(boundaries);
        }
bail:
        return;
}

#if defined(__APPLE__) && defined (__MACH__)
#pragma mark - Data Serialization
#endif


extern DSCaseMessage * DSCaseEncode(const DSCase * aCase)
{
        DSCaseMessage * message = NULL;
        DSUInteger i, numberOfEquations;
        numberOfEquations = DSCaseNumberOfEquations(aCase) + DSCaseNumberOfConservations(aCase) + DSCaseNumberOfInheritedConservations(aCase);
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        message = DSSecureMalloc(sizeof(DSCaseMessage));
        dscase_message__init(message);
        message->ssystem = DSSSystemEncode(DSCaseSSystem(aCase));
        message->casenumber = DSCaseNumber(aCase);
        message->cd = DSMatrixEncode(DSCaseCd(aCase));
        message->ci = DSMatrixEncode(DSCaseCi(aCase));
        message->n_signature = numberOfEquations*2;
        message->signature = DSSecureMalloc(sizeof(DSUInteger)*message->n_signature);
        message->delta = DSMatrixEncode(DSCaseDelta(aCase));
        for (i = 0; i < numberOfEquations*2; i++) {
                message->signature[i] = DSCaseSignature(aCase)[i];
        }
        if (aCase->signature_3d != NULL){
                message->n_signature_3d = numberOfEquations*3;
                message->signature_3d = DSSecureMalloc(sizeof(DSUInteger)*message->n_signature_3d);
                for (i = 0; i< message->n_signature_3d; i++)
                    message->signature_3d[i] = DSCase3Sig(aCase)[i];
        }
        if (DSSSystemHasSolution(DSCaseSSystem(aCase))) {
                message->u = DSMatrixEncode(DSCaseU(aCase));
                message->zeta = DSMatrixEncode(DSCaseZeta(aCase));
        } else {
                message->u = NULL;
                message->zeta = NULL;
        }
        message->caseidentifier = strdup(DSCaseId(aCase));
        if(aCase->conserved_sig != NULL){
            message->n_conserved_sig = numberOfEquations*2;
            message->conserved_sig = DSSecureMalloc(sizeof(DSUInteger)*message->n_conserved_sig);
            for (i = 0; i < message->n_conserved_sig; i++) {
                message->conserved_sig[i] = DSCaseSigCons(aCase)[i];
            }
        }
        message->has_numberinheritedconservations = true;
        message->numberinheritedconservations = DSCaseNumberOfInheritedConservations(aCase);
    
bail:
        return message;
}

extern DSCase * DSCaseFromCaseMessage(const DSCaseMessage * message)
{
        DSCase * aCase = NULL;
        DSUInteger i;
        if (message == NULL) {
                printf("message is NULL\n");
                goto bail;
        }
        aCase = DSCaseAlloc();
        aCase->caseNumber = message->casenumber;
        aCase->Cd = DSMatrixFromMatrixMessage(message->cd);
        aCase->Ci = DSMatrixFromMatrixMessage(message->ci);
        aCase->ssys = DSSSystemFromSSystemMessage(message->ssystem);
        aCase->Xd = DSSSystemXd(DSCaseSSystem(aCase));
        aCase->Xd_a= DSSSystemXd_a(DSCaseSSystem(aCase));
        aCase->Xi = DSSSystemXi(DSCaseSSystem(aCase));
        aCase->delta = DSMatrixFromMatrixMessage(message->delta);
        if (message->u != NULL) {
                aCase->U = DSMatrixFromMatrixMessage(message->u);
                aCase->zeta = DSMatrixFromMatrixMessage(message->zeta);
        } else {
                aCase->U = NULL;
                aCase->zeta = NULL;
        }
        aCase->signature = DSSecureMalloc(sizeof(DSUInteger)*message->n_signature);
        for (i = 0; i < message->n_signature; i++) {
                aCase->signature[i] = message->signature[i];
        }
        if (message->signature_3d != NULL){
                aCase->signature_3d = DSSecureMalloc(sizeof(DSUInteger)*message->n_signature_3d);
                for (i = 0; i < message->n_signature_3d; i++) {
                        aCase->signature_3d[i] = message->signature_3d[i];
                }
        }
        if (message->conserved_sig != NULL){
            aCase->conserved_sig = DSSecureMalloc(sizeof(DSUInteger)*message->n_conserved_sig);
            for (i = 0; i < message->n_conserved_sig; i++) {
                aCase->conserved_sig[i] = message->conserved_sig[i];
            }
        }
        DSCaseId(aCase) = strdup(message->caseidentifier);
        aCase->numberInheritedConservations = message->numberinheritedconservations;
    
bail:
        return aCase;
}

extern DSCase * DSCaseDecode(size_t length, const void * buffer)
{
        DSCase * aCase = NULL;
        DSCaseMessage * message;
        message = dscase_message__unpack(NULL, length, buffer);
        aCase = DSCaseFromCaseMessage(message);
        dscase_message__free_unpacked(message, NULL);
//bail:
        return aCase;
}

extern DSDesignSpace * DSCaseEigenSubspaces(const DSCase * aCase)
{
        DSDesignSpace * subspaceSystem = NULL;
        DSSSystem * reducedSSystem;
        DSExpression *rhs, * p, *n;
        DSUInteger i, j, numberOfEquations;
        DSVariablePool * Xi;
        DSMatrix * Ci, * Cd;
        DSUInteger * columns;
        char * temp, *name;
        char ** equations;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        reducedSSystem = DSSSystemByRemovingAlgebraicConstraints(DSCaseSSystem(aCase));
        numberOfEquations = DSSSystemNumberOfEquations(reducedSSystem);
        numberOfEquations += 2*numberOfEquations;
        equations = DSSecureCalloc(sizeof(char*), numberOfEquations);
        columns = DSSecureMalloc(sizeof(DSUInteger)*DSSSystemNumberOfEquations(reducedSSystem));
        Xi = DSVariablePoolCopy(DSCaseXi(aCase));
        DSVariablePoolSetReadWriteAdd(Xi);
        j = DSSSystemNumberOfEquations(reducedSSystem);
        for (i = 0; i < DSSSystemNumberOfEquations(reducedSSystem); i++) {
                p = DSExpressionFromPowerlawInMatrixForm(i,
                                                         DSSSystemGd(reducedSSystem),
                                                         DSSSystemXd(reducedSSystem),
                                                         DSSSystemGi(reducedSSystem),
                                                         DSSSystemXi(reducedSSystem),
                                                         DSSSystemAlpha(reducedSSystem));
                n = DSExpressionFromPowerlawInMatrixForm(i,
                                                         DSSSystemHd(reducedSSystem),
                                                         DSSSystemXd(reducedSSystem),
                                                         DSSSystemHi(reducedSSystem),
                                                         DSSSystemXi(reducedSSystem),
                                                         DSSSystemBeta(reducedSSystem));
                name = DSVariableName(DSVariablePoolVariableAtIndex(DSSSystemXd_t(reducedSSystem), i));
                temp = DSExpressionAsString(p);
                equations[j] = DSSecureCalloc(sizeof(char), strlen(temp) + strlen(name)+10);
                sprintf(equations[j++], "$e_%s_p = %s", name, temp);
                DSSecureFree(temp);

                temp = DSExpressionAsString(n);
                equations[j] = DSSecureCalloc(sizeof(char), strlen(temp) + strlen(name)+10);
                sprintf(equations[j++], "$e_%s_n = %s", name, temp);
                DSSecureFree(temp);
                rhs = DSExpressionAddExpressions(p, n);
                temp = DSExpressionAsString(rhs);

                equations[i] = DSSecureCalloc(sizeof(char), strlen(temp) + strlen(name)+10);
                sprintf(equations[i], "$e_%s = %s", name, temp);
                DSSecureFree(temp);
                DSExpressionFree(rhs);
                DSVariablePoolAddVariableWithName(Xi, name);
                columns[i] = DSVariablePoolIndexOfVariableWithName(DSCaseXd(aCase), name);
        }
        subspaceSystem = DSDesignSpaceByParsingStringsWithXi(equations, NULL, Xi, numberOfEquations);
        DSVariablePoolFree(Xi);
        Cd = DSMatrixSubMatrixIncludingColumns(DSCaseCd(aCase), DSSSystemNumberOfEquations(reducedSSystem), columns);
        Ci = DSMatrixAppendMatrices(DSCaseCi(aCase), Cd, true);
        DSMatrixFree(Cd);
        Cd = DSMatrixCalloc(DSMatrixRows(Ci), numberOfEquations);
        DSDesignSpaceAddConditions(subspaceSystem, Cd, Ci, DSCaseDelta(aCase));
        DSSecureFree(columns);
        for (i = 0; i < numberOfEquations; i++) {
                DSSecureFree(equations[i]);
        }
        DSSecureFree(equations);
        DSSSystemFree(reducedSSystem);
        DSMatrixFree(Cd);
        DSMatrixFree(Ci);
bail:
        return subspaceSystem;
}






