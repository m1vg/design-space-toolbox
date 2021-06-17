/**
 * \file DSCaseLinearProgramming.c
 * \brief Implementation file with functions for linear programming
 *        operations dealing with cases in design space.
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

/**
 * \todo Find/write a parallelizable linear programming package
 */

#include <stdio.h>
#include <string.h>
#include <glpk.h>

#include "qhull_ra.h"
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
#include "DSVertices.h"
#include "DSNVertexEnumeration.h"
#include "DSGMASystemParsingAux.h"
#include "DSExpressionTokenizer.h"
#include "DSCaseOptimizationFunctionGrammar.h"
#if defined (__APPLE__) && defined (__MACH__)
#pragma mark Linear programming functions
#endif


static glp_prob * dsCaseLinearProblemForMatrices(const DSMatrix *A, const DSMatrix *B)
{
        glp_prob *linearProblem = NULL;
        int * ia = NULL, *ja = NULL;
        double *ar = NULL;
        DSUInteger i, numberOfXi, numberOfBoundaries;
        
        glp_term_out(GLP_OFF);
        linearProblem = glp_create_prob();
        if (linearProblem == NULL) {
                DSError(M_DS_NULL ": Linear problem is NULL", A_DS_ERROR);
                goto bail;
        }
        numberOfXi = DSMatrixColumns(A);
        numberOfBoundaries = DSMatrixRows(A);
        
        ia = DSMatrixRowsForGLPK(A);
        ja = DSMatrixColumnsForGLPK(A);
        ar = DSMatrixDataForGLPK(A);
        
        glp_add_rows(linearProblem, numberOfBoundaries);
        glp_add_cols(linearProblem, numberOfXi);
        
        glp_set_obj_dir(linearProblem, GLP_MIN);
        glp_load_matrix(linearProblem, numberOfBoundaries*numberOfXi,
                        ia, ja, ar);
        for (i = 0; i < numberOfBoundaries; i++) {
                glp_set_row_bnds(linearProblem, i+1, GLP_UP, 0.0,
                                 DSMatrixDoubleValue(B, i, 0));
        }
        for (i = 0; i < numberOfXi; i++){
                glp_set_col_bnds(linearProblem, i+1, GLP_FR, 0., 0.); // original jason
//                glp_set_col_bnds(linearProblem, i+1, GLP_DB, -20.0, 20.0); // edits miguel
        }
        
        if (ia != NULL)
                DSSecureFree(ia);
        if (ja != NULL)
                DSSecureFree(ja);
        if (ar != NULL)
                DSSecureFree(ar);
bail:
        return linearProblem;
}

static glp_prob * dsCaseLinearProblemForCaseValidity(const DSMatrix * U, const DSMatrix *zeta)
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
        
        linearProblem = dsCaseLinearProblemForMatrices(coefficients, zeta);
        glp_set_col_bnds(linearProblem, glp_get_num_cols(linearProblem), GLP_LO, -1.0, 0.0);
        glp_set_obj_coef(linearProblem, glp_get_num_cols(linearProblem), 1.0);
        
        
        DSMatrixFree(coefficients);
        if (slacks != NULL)
                DSMatrixFree(slacks);
bail:
        return linearProblem;
}

extern const bool DSCaseConditionsAreValid(const DSCase *aCase)
{
        bool isValid = false;
        glp_prob *linearProblem = NULL;
        DSMatrix * U, * Zeta;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        U = DSMatrixAppendMatrices(DSCaseCd(aCase), DSCaseCi(aCase), true);
        Zeta = DSCaseDelta(aCase);
        linearProblem = dsCaseLinearProblemForCaseValidity(U, Zeta);
        DSMatrixFree(U);
        if (linearProblem != NULL) {
                glp_simplex(linearProblem, NULL);
                if (glp_get_obj_val(linearProblem) <= -1E-14 && glp_get_prim_stat(linearProblem) == GLP_FEAS) {
                        isValid = true;
                }
                glp_delete_prob(linearProblem);
        }
bail:
        return isValid;
}

extern const bool DSCaseHasSharedBoundaries(const DSCase * aCase1, const DSCase * aCase2, const bool intersecting){
    
        bool has_shared = false, has_shared_intersecting = false, has_shared_non_intersecting = false;
        glp_prob *linearProblem = NULL;
        const DSCase **pointer_cases = NULL;
        DSCase * aCase = NULL;
    
        if (aCase1 == NULL || aCase2 == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (DSCaseHasSolution(aCase1) == false || DSCaseHasSolution(aCase2) == false ) {
                goto bail;
        }
    
        pointer_cases = DSSecureMalloc(sizeof( DSCase *)*2);
        pointer_cases[0] = aCase1;
        pointer_cases[1] = aCase2;
        aCase = DSPseudoCaseFromIntersectionOfCases(2, pointer_cases);
    
        linearProblem = dsCaseLinearProblemForCaseValidity(DSCaseU(aCase), DSCaseZeta(aCase));
    
        if (linearProblem != NULL) {
                glp_simplex(linearProblem, NULL);
                if (glp_get_obj_val(linearProblem) <= -1E-14 && glp_get_prim_stat(linearProblem) == GLP_FEAS) {
                        has_shared_intersecting = true;
                }
                if ((glp_get_obj_val(linearProblem) >= -1E-14 || glp_get_obj_val(linearProblem) <= 1E-14 ) && glp_get_prim_stat(linearProblem) == GLP_FEAS) {
                        has_shared_non_intersecting = true;
                }
                glp_delete_prob(linearProblem);
        }
        if (intersecting == true)
            has_shared = has_shared_intersecting;
        else
            has_shared = has_shared_non_intersecting;
    
        if (pointer_cases != NULL)
            DSSecureFree(pointer_cases);
        if (aCase != NULL)
            DSSecureFree(aCase);
bail:
    return has_shared;
}

extern DSMatrix * DSCaseSharedBoundaries(const DSCase *aCase1, const DSCase *aCase2, const bool intersecting){
    
    bool has_shared = false;
    DSMatrix *aux1 = NULL, *aux2 = NULL, *aux3 = NULL, *aux3_t = NULL;
    DSMatrix *p1 = NULL, * p2 = NULL, * SharedBoundaries_aux = NULL, *SharedBoundaries = NULL;
    DSUInteger bound_1, bound_2, rank, counter = 0, * rows = NULL, i, max_val;
    
    
    //0. Initialice matrix SharedBoundaries_aux
    max_val = DSCaseNumberOfBoundaries(aCase1) * DSCaseNumberOfBoundaries(aCase2);

    SharedBoundaries_aux = DSMatrixAlloc(max_val, 2);
    
    //1. First, let's create a intersection from these two cases and see if they have the pontential to have shared boundaries
    has_shared = DSCaseHasSharedBoundaries(aCase1, aCase2, intersecting);
    
    // 2. If these two cases have the potential for shared boundaries, calculate pairwise rank of boundaries.
    if (has_shared == true){
        
        for (bound_1 = 0; bound_1 < DSCaseNumberOfBoundaries(aCase1); bound_1++){
            for (bound_2 = 0; bound_2 < DSCaseNumberOfBoundaries(aCase2); bound_2++){
                
                
                p1 = DSMatrixSubMatrixIncludingRowList(aCase1->U, 1, bound_1);
                p2 = DSMatrixSubMatrixIncludingRowList(aCase1->zeta, 1, bound_1);
                aux1 = DSMatrixAppendMatrices(p1, p2, true);
                if (p1 != NULL)
                    DSMatrixFree(p1);
                if (p2 != NULL)
                    DSMatrixFree(p2);
                
                p1 = DSMatrixSubMatrixIncludingRowList(aCase2->U, 1, bound_2);
                p2 = DSMatrixSubMatrixIncludingRowList(aCase2->zeta, 1, bound_2);
                aux2 = DSMatrixAppendMatrices(p1, p2, true);
                if (p1 != NULL)
                    DSMatrixFree(p1);
                if (p2 != NULL)
                    DSMatrixFree(p2);
                
                aux3 = DSMatrixAppendMatrices(aux1, aux2, false);
                aux3_t = DSMatrixTranspose(aux3);
                
                rank = DSMatrixRank(aux3_t);
            
                if (rank == 1){
                    DSMatrixSetDoubleValue(SharedBoundaries_aux, counter, 0, bound_1);
                    DSMatrixSetDoubleValue(SharedBoundaries_aux, counter, 1, bound_2);
                    counter ++;
                }
                
                if (aux1 != NULL)
                    DSMatrixFree(aux1);
                if (aux2 != NULL)
                    DSMatrixFree(aux2);
                if (aux3 != NULL)
                    DSMatrixFree(aux3);
                if (aux3_t != NULL)
                    DSMatrixFree(aux3_t);
            }
        }
    }
    
    // 3. Initialize the rows vector, cut matrix SharedBoundaries_aux according to rows.
    if (counter != 0){
        rows = DSSecureMalloc(sizeof(DSUInteger)*counter);
        for (i=0; i<counter; i++){
            rows[i] = i;
        }
        SharedBoundaries = DSMatrixSubMatrixIncludingRows(SharedBoundaries_aux, counter, rows);
    }
    
    if (SharedBoundaries_aux != NULL)
        DSMatrixFree(SharedBoundaries_aux);
    if (rows != NULL)
        DSSecureFree(rows);

    return SharedBoundaries;
}

extern const long int DSCaseSharedBoundariesNumberOfVertices(const DSCase *aCase1,
                                                             const DSCase *aCase2,
                                                             const DSVariablePool *lowerBounds,
                                                             const DSVariablePool *upperBounds,
                                                             const long int maxVertices,
                                                             const bool limitVertices){
    
    DSUInteger n_variables = 0, k;
    const DSCase **pointer_cases = NULL;
    DSCase * aCase = NULL;
    DSVariablePool *lowerBounds_int = NULL, *upperBounds_int = NULL;
    const DSVariable ** all_variables = NULL ;
    long int numberOfVertices = 0;
    
    //* The variable pools lowerBounds and upperBounds should assign values to the variables contained in the Xi pool (DSCaseXi(aCase))
    
    if (aCase1 == NULL || aCase2 == NULL){
        DSError(M_DS_CASE_NULL, A_DS_ERROR);
        goto bail;
    }

    if (lowerBounds == NULL || upperBounds == NULL){
        DSError(M_DS_NULL":Variable Pool is Null",A_DS_ERROR);
        goto bail;
    }
    

    pointer_cases = DSSecureMalloc(sizeof( DSCase *)*2);
    pointer_cases[0] = aCase1;
    pointer_cases[1] = aCase2;
    aCase = DSPseudoCaseFromIntersectionOfCases(2, pointer_cases);
        
    if (aCase == NULL) {
            DSError(M_DS_CASE_NULL, A_DS_ERROR);
            goto bail;
    }
    
    n_variables = DSVariablePoolNumberOfVariables(DSCaseXi(aCase));
    
    if (n_variables != DSVariablePoolNumberOfVariables(lowerBounds)
        || DSVariablePoolNumberOfVariables(lowerBounds) != DSVariablePoolNumberOfVariables(upperBounds)){
        DSError(M_DS_WRONG":Inconsistent number of variables in pool",A_DS_ERROR);
        goto bail;
    }
    
    //* Let's just re-assign the values from lowerBounds and upperBounds just to make sure they are set in the correct order.
    all_variables = DSVariablePoolAllVariables(DSCaseXi(aCase));
    lowerBounds_int = DSVariablePoolCopy(DSCaseXi(aCase));
    upperBounds_int = DSVariablePoolCopy(DSCaseXi(aCase));
    for (k=0; k<n_variables; k++){
        
        DSVariablePoolSetValueForVariableWithName(lowerBounds_int,
                                                  DSVariableName(all_variables[k]),
                                                  DSVariablePoolValueForVariableWithName(lowerBounds, DSVariableName(all_variables[k])));
        
        DSVariablePoolSetValueForVariableWithName(upperBounds_int,
                                                  DSVariableName(all_variables[k]),
                                                  DSVariablePoolValueForVariableWithName(upperBounds, DSVariableName(all_variables[k])));
        
    }
    
    DSCaseRemoveRedundantBoundaries(aCase);
    numberOfVertices = DSCaseNDVertexEnumerationNumberOfVertices(aCase,
                                                                lowerBounds_int,
                                                                upperBounds_int,
                                                                maxVertices,
                                                                limitVertices);
    
    if (pointer_cases != NULL)
        DSSecureFree(pointer_cases);
    if (aCase != NULL)
        DSCaseFree(aCase);
    if (lowerBounds_int != NULL)
        DSVariablePoolFree(lowerBounds_int);
    if (upperBounds_int != NULL)
        DSVariablePoolFree(upperBounds_int);

bail:
    return numberOfVertices;
    
}

extern const bool DSCaseIsValid(const DSCase *aCase, const bool strict)
{
        bool isValid = false;
        glp_prob *linearProblem = NULL;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (DSCaseHasSolution(aCase) == false) {
                goto bail;
        }
        linearProblem = dsCaseLinearProblemForCaseValidity(DSCaseU(aCase), DSCaseZeta(aCase));
        if (linearProblem != NULL) {
                glp_simplex(linearProblem, NULL);
                if (strict == true) {
                        if (glp_get_obj_val(linearProblem) <= -1E-14 && glp_get_prim_stat(linearProblem) == GLP_FEAS) {
                                isValid = true;
                        }
                } else {
                        if (glp_get_obj_val(linearProblem) <= 0.0f && glp_get_prim_stat(linearProblem) == GLP_FEAS) {
                                isValid = true;
                        }
                }
                glp_delete_prob(linearProblem);
        }
bail:
        return isValid;
}

extern const bool DSCaseSharedBoundariesIsValid(const DSCase *aCase)
{
        bool isValid = false;
        glp_prob *linearProblem = NULL;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (DSCaseHasSolution(aCase) == false) {
                goto bail;
        }
        linearProblem = dsCaseLinearProblemForCaseValidity(DSCaseU(aCase), DSCaseZeta(aCase));
        if (linearProblem != NULL) {
                glp_simplex(linearProblem, NULL);
                if (glp_get_obj_val(linearProblem) <= 1E-14 && glp_get_prim_stat(linearProblem) == GLP_FEAS) {
                        isValid = true;
                }
                glp_delete_prob(linearProblem);
        }
bail:
        return isValid;
}

extern const bool DSCasesSharedBoundariesIsValid(const DSCase *aCase1, const DSCase *aCase2){
    
    bool isValid = false;
    const DSCase **pointer_cases = NULL;
    DSCase * aCase = NULL;


    
    if (aCase1 == NULL || aCase2 == NULL){
        DSError(M_DS_CASE_NULL, A_DS_ERROR);
        goto bail;
    }
    
    pointer_cases = DSSecureMalloc(sizeof( DSCase *)*2);
    pointer_cases[0] = aCase1;
    pointer_cases[1] = aCase2;
    aCase = DSPseudoCaseFromIntersectionOfCases(2, pointer_cases);
    
    isValid = DSCaseSharedBoundariesIsValid(aCase);
    
    if (pointer_cases != NULL)
        DSSecureFree(pointer_cases);
    if (aCase != NULL)
        DSCaseFree(aCase);
    
bail:
    return isValid;
}

__deprecated extern const bool DSCaseIsValidInStateSpace(const DSCase *aCase) {
        return DSCaseIsConsistent(aCase);
}

extern const bool DSCaseIsConsistent(const DSCase *aCase)
{
        bool isValid = false;
        glp_prob *linearProblem = NULL;
        DSMatrix * C;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        C  = DSMatrixAppendMatrices(DSCaseCd(aCase), DSCaseCi(aCase), true);
        linearProblem = dsCaseLinearProblemForCaseValidity(C, DSCaseDelta(aCase));
        if (linearProblem != NULL) {
                glp_simplex(linearProblem, NULL);
                if (glp_get_obj_val(linearProblem) <= -1E-14 && glp_get_prim_stat(linearProblem) == GLP_FEAS) {
                        isValid = true;
                }
                glp_delete_prob(linearProblem);
        }
bail:
        return isValid;
}

//extern const bool DSCaseIsValidInStateSpaceAtSlice(const DSCase *aCase, const DSVariablePool * lower, const DSVariablePool * upper)
//{
//        bool isValid = false;
//        glp_prob *linearProblem = NULL;
//        DSMatrix * C;
//        if (aCase == NULL) {
//                DSError(M_DS_CASE_NULL, A_DS_ERROR);
//                goto bail;
//        }
//        C  = DSMatrixAppendMatrices(DSCaseCd(aCase), DSCaseCi(aCase), true);
//        linearProblem = dsCaseLinearProblemForCaseValidity(C, DSCaseDelta(aCase));
//        if (linearProblem != NULL) {
//                glp_simplex(linearProblem, NULL);
//                if (glp_get_obj_val(linearProblem) <= -1E-14 && glp_get_prim_stat(linearProblem) == GLP_FEAS) {
//                        isValid = true;
//                }
//                glp_delete_prob(linearProblem);
//        }
//bail:
//        return isValid;
//}

extern const bool DSCaseIsValidAtPoint(const DSCase *aCase, const DSVariablePool * variablesToFix)
{
        bool isValid = false;
        DSUInteger i, numberToRemove, indexOfVariable;
        DSMatrix *result, *Xi;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (DSCaseHasSolution(aCase) == false) {
                goto bail;
        }
        if (variablesToFix == NULL) {
                DSError(M_DS_VAR_NULL ": Variable pool with variables to fix is NULL", A_DS_ERROR);
                goto bail;
        }
        if (DSVariablePoolNumberOfVariables(variablesToFix) != DSVariablePoolNumberOfVariables(DSCaseXi(aCase))) {
                DSError(M_DS_WRONG ": Number of variables the same as the number Xi", A_DS_ERROR);
                goto bail;
        }
        if (DSVariablePoolNumberOfVariables(variablesToFix) == 0) {
                DSError(M_DS_WRONG ": Case has no independent variables", A_DS_WARN);
                isValid = DSCaseIsValid(aCase, false);
                goto bail;
        }
        numberToRemove = DSVariablePoolNumberOfVariables(variablesToFix);
        Xi = DSMatrixAlloc(DSVariablePoolNumberOfVariables(DSCaseXi(aCase)), 1);
        for (i = 0; i < numberToRemove; i++) {
                indexOfVariable = DSVariablePoolIndexOfVariableWithName(DSCaseXi(aCase), 
                                                                        DSVariableName(DSVariablePoolAllVariables(variablesToFix)[i]));
                if (indexOfVariable >= DSMatrixRows(Xi)) {
                        DSMatrixFree(Xi);
                        goto bail;
                }
                DSMatrixSetDoubleValue(Xi, indexOfVariable, 0, log10(DSVariableValue(DSVariablePoolAllVariables(variablesToFix)[i])));
        }
        result = DSMatrixByMultiplyingMatrix(DSCaseU(aCase), Xi);
        DSMatrixAddByMatrix(result, DSCaseZeta(aCase));
        
        
        for (i = 0; i < DSMatrixRows(result); i++)
                if (DSMatrixDoubleValue(result, i, 0) < 0)
                        break;
        if (i == DSMatrixRows(result))
                isValid = true;
        DSMatrixFree(result);
        DSMatrixFree(Xi);
bail:
        return isValid;
}

__deprecated extern const bool DSCaseIsValidInStateSpaceAtPoint(const DSCase *aCase, const DSVariablePool * Xd_p, const DSVariablePool * Xi_p) {
        return DSCaseIsConsistentAtPoint(aCase, Xd_p, Xi_p);
}

extern const bool DSCaseIsConsistentAtPoint(const DSCase *aCase, const DSVariablePool * Xd_p, const DSVariablePool * Xi_p)
{
        bool isValid = false;
        DSUInteger i, numberOfXi, numberOfXd, indexOfVariable;
        DSMatrix *result, *CdYd, *CiYi, *Yi, *Yd;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (Xd_p == NULL) {
                DSError(M_DS_VAR_NULL ": Variable pool with values for dependent variable is NULL", A_DS_ERROR);
                goto bail;
        }
        if (Xi_p == NULL) {
                DSError(M_DS_VAR_NULL ": Variable pool with values for independent variable is NULL", A_DS_ERROR);
                goto bail;
        }
        if (DSVariablePoolNumberOfVariables(Xd_p) != DSVariablePoolNumberOfVariables(DSCaseXd(aCase))) {
                DSError(M_DS_WRONG ": Inconsistent number of dependent variables", A_DS_ERROR);
                goto bail;
        }
        if (DSVariablePoolNumberOfVariables(Xi_p) != DSVariablePoolNumberOfVariables(DSCaseXi(aCase))) {
                DSError(M_DS_WRONG ": Inconsistent number of independent variables", A_DS_ERROR);
                goto bail;
        }
        if (DSVariablePoolNumberOfVariables(Xi_p) == 0) {
                DSError(M_DS_WRONG ": Case has no independent variables", A_DS_WARN);
                isValid = DSCaseIsValid(aCase, false);
                goto bail;
        }
        numberOfXi = DSVariablePoolNumberOfVariables(Xi_p);
        numberOfXd = DSVariablePoolNumberOfVariables(Xd_p);
        Yd = DSMatrixAlloc(numberOfXd, 1);
        Yi = DSMatrixAlloc(numberOfXi, 1);
        for (i = 0; i < numberOfXd; i++) {
                indexOfVariable = DSVariablePoolIndexOfVariableWithName(DSCaseXd(aCase),
                                                                        DSVariableName(DSVariablePoolAllVariables(Xd_p)[i]));
                if (indexOfVariable >= DSMatrixRows(Yd)) {
                        DSMatrixFree(Yi);
                        DSMatrixFree(Yd);
                        goto bail;
                }
                DSMatrixSetDoubleValue(Yd, indexOfVariable, 0, log10(DSVariableValue(DSVariablePoolAllVariables(Xd_p)[i])));
        }
        for (i = 0; i < numberOfXi; i++) {
                indexOfVariable = DSVariablePoolIndexOfVariableWithName(DSCaseXi(aCase),
                                                                        DSVariableName(DSVariablePoolAllVariables(Xi_p)[i]));
                if (indexOfVariable >= DSMatrixRows(Yi)) {
                        DSMatrixFree(Yi);
                        DSMatrixFree(Yd);
                        goto bail;
                }
                DSMatrixSetDoubleValue(Yi, indexOfVariable, 0, log10(DSVariableValue(DSVariablePoolAllVariables(Xi_p)[i])));
        }
        CdYd = DSMatrixByMultiplyingMatrix(DSCaseCd(aCase), Yd);
        CiYi = DSMatrixByMultiplyingMatrix(DSCaseCi(aCase), Yi);
        result = DSMatrixByAddingMatrix(CdYd, CiYi);
        DSMatrixAddByMatrix(result, DSCaseDelta(aCase));
        for (i = 0; i < DSMatrixRows(result); i++)
                if (DSMatrixDoubleValue(result, i, 0) < 0)
                        break;
        if (i == DSMatrixRows(result))
                isValid = true;
        DSMatrixFree(result);
        DSMatrixFree(Yi);
        DSMatrixFree(Yd);
        DSMatrixFree(CiYi);
        DSMatrixFree(CdYd);
bail:
        return isValid;
}

extern DSVariablePool * DSCaseConsistentParameterAndStateSet(const DSCase *aCase)
{
        DSVariablePool * Xi = NULL;
        glp_prob *linearProblem = NULL;
        DSUInteger i;
        DSMatrix * C;
        char * name;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        C = DSMatrixAppendMatrices(DSCaseCd(aCase), DSCaseCi(aCase), true);
        linearProblem = dsCaseLinearProblemForCaseValidity(C, DSCaseDelta(aCase));
        DSMatrixFree(C);
        if (linearProblem != NULL) {
                glp_simplex(linearProblem, NULL);
                if (glp_get_obj_val(linearProblem) <= -1E-14 && glp_get_prim_stat(linearProblem) == GLP_FEAS) {
                        Xi = DSVariablePoolAlloc();
                        for (i = 0; i < DSVariablePoolNumberOfVariables(DSCaseXd(aCase)); i++) {
                                name = DSVariableName(DSVariablePoolVariableAtIndex(DSCaseXd(aCase), i));
                                DSVariablePoolAddVariableWithName(Xi, name);
                                DSVariablePoolSetValueForVariableWithName(Xi, name, pow(10, glp_get_col_prim(linearProblem, i+1)));
                        }
                        for (i = 0; i < DSVariablePoolNumberOfVariables(DSCaseXi(aCase)); i++) {
                                name = DSVariableName(DSVariablePoolVariableAtIndex(DSCaseXi(aCase), i));
                                DSVariablePoolAddVariableWithName(Xi, name);
                                DSVariablePoolSetValueForVariableWithName(Xi, name, pow(10, glp_get_col_prim(linearProblem, DSVariablePoolNumberOfVariables(DSCaseXd(aCase))+i+1)));
                        }
                } else {
                        printf("invalid.\n");
                }
                glp_delete_prob(linearProblem);
        }
bail:
        return Xi;
}

extern DSVariablePool * DSCaseValidParameterAndStateSet(const DSCase *aCase)
{
        DSVariablePool * Xi = NULL;
        glp_prob *linearProblem = NULL;
        DSUInteger i;
        DSMatrix * C, *temp1, *temp2;
        char * name;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        
        temp1 = DSMatrixCalloc(DSMatrixRows(DSCaseU(aCase)), DSMatrixColumns(DSCaseCd(aCase)));
        temp2 = DSMatrixAppendMatrices(temp1, DSCaseU(aCase), true);
        DSMatrixFree(temp1);
        temp1 = DSMatrixAppendMatrices(DSCaseCd(aCase), DSCaseCi(aCase), true);
        C = DSMatrixAppendMatrices(temp1, temp2, false);
        DSMatrixFree(temp1);
        DSMatrixFree(temp2);
        temp1 = DSMatrixAppendMatrices(DSCaseDelta(aCase), DSCaseZeta(aCase), false);
        linearProblem = dsCaseLinearProblemForCaseValidity(C, temp1);
        DSMatrixFree(C);
        DSMatrixFree(temp1);
        if (linearProblem != NULL) {
                glp_simplex(linearProblem, NULL);
                if (glp_get_obj_val(linearProblem) <= -1E-14 && glp_get_prim_stat(linearProblem) == GLP_FEAS) {
                        Xi = DSVariablePoolAlloc();
                        for (i = 0; i < DSVariablePoolNumberOfVariables(DSCaseXd(aCase)); i++) {
                                name = DSVariableName(DSVariablePoolVariableAtIndex(DSCaseXd(aCase), i));
                                DSVariablePoolAddVariableWithName(Xi, name);
                                DSVariablePoolSetValueForVariableWithName(Xi, name, pow(10, glp_get_col_prim(linearProblem, i+1)));
                        }
                        for (i = 0; i < DSVariablePoolNumberOfVariables(DSCaseXi(aCase)); i++) {
                                name = DSVariableName(DSVariablePoolVariableAtIndex(DSCaseXi(aCase), i));
                                DSVariablePoolAddVariableWithName(Xi, name);
                                DSVariablePoolSetValueForVariableWithName(Xi, name, pow(10, glp_get_col_prim(linearProblem, DSVariablePoolNumberOfVariables(DSCaseXd(aCase))+i+1)));
                        }
                } else {
                        printf("invalid.\n");
                }
                glp_delete_prob(linearProblem);
        }
bail:
        return Xi;
}

extern DSVariablePool * DSCaseValidParameterSet(const DSCase *aCase)
{
        DSVariablePool * Xi = NULL;
        glp_prob *linearProblem = NULL;
        DSUInteger i;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (DSCaseIsValid(aCase, false) == false)
                goto bail;
        linearProblem = dsCaseLinearProblemForCaseValidity(DSCaseU(aCase), DSCaseZeta(aCase));
        if (linearProblem != NULL) {
                glp_simplex(linearProblem, NULL);
                if (glp_get_obj_val(linearProblem) <= -1E-14 && glp_get_prim_stat(linearProblem) == GLP_FEAS) {
                        Xi = DSVariablePoolCopy(DSCaseXi(aCase));
                        DSVariablePoolSetReadWriteAdd(Xi);
                        for (i = 0; i < DSVariablePoolNumberOfVariables(Xi); i++) {
                                DSVariableSetValue(DSVariablePoolAllVariables(Xi)[i], pow(10, glp_get_col_prim(linearProblem, i+1)));
                        }
                }
                glp_delete_prob(linearProblem);
        }
bail:
        return Xi;
}

extern DSVariablePool * DSCaseSharedBoundariesValidParameterSet(const DSCase *aCase)
{
        DSVariablePool * Xi = NULL;
        glp_prob *linearProblem = NULL;
        DSUInteger i;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (DSCaseSharedBoundariesIsValid(aCase) == false){
                goto bail;
        }
        linearProblem = dsCaseLinearProblemForCaseValidity(DSCaseU(aCase), DSCaseZeta(aCase));
        if (linearProblem != NULL) {
                glp_simplex(linearProblem, NULL);
                if (glp_get_obj_val(linearProblem) <= 1E-14 && glp_get_prim_stat(linearProblem) == GLP_FEAS) {
                        Xi = DSVariablePoolCopy(DSCaseXi(aCase));
                        DSVariablePoolSetReadWriteAdd(Xi);
                        for (i = 0; i < DSVariablePoolNumberOfVariables(Xi); i++) {
                                DSVariableSetValue(DSVariablePoolAllVariables(Xi)[i], pow(10, glp_get_col_prim(linearProblem, i+1)));
                        }
                }
                glp_delete_prob(linearProblem);
        }
        
bail:
        return Xi;
}

extern DSVariablePool * DSCaseValidParameterSetByOptimizingFunction(const DSCase *aCase, const char * function, const bool minimize)
{
        DSVariablePool * Xi = NULL;
        glp_prob *linearProblem = NULL;
        DSUInteger i;
        DSMatrixArray * objective;
        DSMatrix * Od, *U, *Zeta;
        DSMatrix * delta;
        DSExpression * expression;
        char * processedFunction;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (DSCaseIsValid(aCase, false) == false)
                goto bail;
        expression = DSExpressionByParsingString(function);
        if (expression == NULL) {
                DSError(M_DS_NULL ": Could not parse function string", A_DS_ERROR);
                goto bail;
        }
        processedFunction = DSExpressionAsString(expression);
        DSExpressionFree(expression);
        objective = DSCaseParseOptimizationFunction(aCase, processedFunction);
        DSSecureFree(processedFunction);
        if (objective == NULL) {
                goto bail;
        }
        U = DSMatrixCopy(DSCaseU(aCase));
        Zeta = DSMatrixCopy(DSCaseZeta(aCase));
        Od = DSMatrixArrayMatrix(objective, 0);
        delta = DSMatrixArrayMatrix(objective, 1);
        DSMatrixMultiplyByScalar(U, -1.0);
        linearProblem = dsCaseLinearProblemForMatrices(U, Zeta);
        DSMatrixFree(U);
        DSMatrixFree(Zeta);
        if (linearProblem == NULL) {
                DSMatrixArrayFree(objective);
                goto bail;
        }
        if (minimize == false) {
                glp_set_obj_dir(linearProblem, GLP_MAX);
        }
        for (i = 0; i < DSMatrixColumns(Od); i++) {
                glp_set_obj_coef(linearProblem, i+1, DSMatrixDoubleValue(Od, 0, i));
                // Limits on optimization bounded between 1e-20 and 1e20
                glp_set_col_bnds(linearProblem, i+1, GLP_DB, -20, 20);
        }
        glp_set_obj_coef(linearProblem, 0, DSMatrixDoubleValue(delta, 0, 0));
        glp_simplex(linearProblem, NULL);
        if (glp_get_status(linearProblem) != GLP_OPT) {
                glp_delete_prob(linearProblem);
                DSMatrixArrayFree(objective);
                goto bail;
        }
        Xi = DSVariablePoolCopy(DSCaseXi(aCase));
        DSVariablePoolSetReadWriteAdd(Xi);
        for (i = 0; i < DSVariablePoolNumberOfVariables(Xi); i++) {
                DSVariableSetValue(DSVariablePoolAllVariables(Xi)[i], pow(10, glp_get_col_prim(linearProblem, i+1)));
        }
        glp_delete_prob(linearProblem);
        DSMatrixArrayFree(objective);
bail:
        return Xi;
}

static DSUInteger dsCaseNumberOfFreeDependentVariablesForBounds(const DSCase *aCase, const DSVariablePool * lowerBounds, const DSVariablePool *upperBounds)
{
        DSUInteger i, variableIndex, freeVariables = 0;
        const DSVariable * lowVariable, *highVariable;
        double low, high;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (lowerBounds == NULL && upperBounds == NULL) {
                DSError(M_DS_VAR_NULL ": Variable pool with variables to fix is NULL", A_DS_ERROR);
                goto bail;
        }
        if (DSVariablePoolNumberOfVariables(lowerBounds) != DSVariablePoolNumberOfVariables(upperBounds)) {
                DSError(M_DS_WRONG ": Number of variables to bound must match", A_DS_ERROR);
                goto bail;
        }
        for (i = 0; i < DSVariablePoolNumberOfVariables(lowerBounds); i++) {
                
                lowVariable = DSVariablePoolAllVariables(lowerBounds)[i];
                highVariable = DSVariablePoolVariableWithName(upperBounds, DSVariableName(lowVariable));
                
                if (lowVariable == NULL || highVariable == NULL) {
                        DSError(M_DS_WRONG ": Variables to bound are not consistent", A_DS_WARN);
                        continue;
                }
                
                if (DSVariablePoolHasVariableWithName(DSCaseXd(aCase), DSVariableName(lowVariable)) == false) {
                        continue;
                }
                variableIndex = DSVariablePoolIndexOfVariableWithName(DSCaseXd(aCase),
                                                                      DSVariableName(lowVariable));
                low = DSVariableValue(lowVariable);
                high = DSVariableValue(highVariable);
                
                if (low > high) {
                        DSError(M_DS_WRONG ": Variable bounds are not consistent", A_DS_WARN);
                        continue;
                }
                
                if (variableIndex >= DSVariablePoolNumberOfVariables(DSCaseXd(aCase)))
                        continue;
                if (low == high)
                        continue;
                freeVariables++;
        }
bail:
        return freeVariables;
}

static DSUInteger dsCaseNumberOfFreeVariablesForBounds(const DSCase *aCase, const DSVariablePool * lowerBounds, const DSVariablePool *upperBounds)
{
        DSUInteger i, variableIndex, freeVariables = 0;
        const DSVariable * lowVariable, *highVariable;
        double low, high;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (lowerBounds == NULL && upperBounds == NULL) {
                DSError(M_DS_VAR_NULL ": Variable pool with variables to fix is NULL", A_DS_ERROR);
                goto bail;
        }
        if (DSVariablePoolNumberOfVariables(lowerBounds) != DSVariablePoolNumberOfVariables(upperBounds)) {
                DSError(M_DS_WRONG ": Number of variables to bound must match", A_DS_ERROR);
                goto bail;
        }
        for (i = 0; i < DSVariablePoolNumberOfVariables(lowerBounds); i++) {
                
                lowVariable = DSVariablePoolAllVariables(lowerBounds)[i];
                highVariable = DSVariablePoolVariableWithName(upperBounds, DSVariableName(lowVariable));
                
                if (lowVariable == NULL || highVariable == NULL) {
                        DSError(M_DS_WRONG ": Variables to bound are not consistent", A_DS_WARN);
                        continue;
                }
                
                if (DSVariablePoolHasVariableWithName(DSCaseXi(aCase), DSVariableName(lowVariable)) == false) {
                        continue;
                }
                variableIndex = DSVariablePoolIndexOfVariableWithName(DSCaseXi(aCase),
                                                                      DSVariableName(lowVariable));
                low = DSVariableValue(lowVariable);
                high = DSVariableValue(highVariable);
                
                if (low > high) {
                        DSError(M_DS_WRONG ": Variable bounds are not consistent", A_DS_WARN);
                        continue;
                }
                
                if (variableIndex >= DSVariablePoolNumberOfVariables(DSCaseXi(aCase)))
                        continue;
                if (low == high) 
                        continue;
                freeVariables++;
        }
bail:
        return freeVariables;
}

static DSUInteger dsCaseSetDependentAndIndependentVariableBoundsLinearProblem(const DSCase *aCase, glp_prob *linearProblem,  const DSVariablePool * lowerBounds, const DSVariablePool *upperBounds)
{
        DSUInteger i, variableIndex, freeVariables = 0;
        const DSVariable * lowVariable, *highVariable;
        double low, high;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (lowerBounds == NULL && upperBounds == NULL) {
                DSError(M_DS_VAR_NULL ": Variable pool with variables to fix is NULL", A_DS_ERROR);
                goto bail;
        }
        if (DSVariablePoolNumberOfVariables(lowerBounds) != DSVariablePoolNumberOfVariables(upperBounds)) {
                DSError(M_DS_WRONG ": Number of variables to bound must match", A_DS_ERROR);
                goto bail;
        }
        if (linearProblem == NULL) {
                DSError(M_DS_NULL ": Linear problem is NULL", A_DS_WARN);
                goto bail;
        }
        for (i = 0; i < glp_get_num_cols(linearProblem); i++) {
                glp_set_col_bnds(linearProblem, i+1, GLP_FR, 0.0, 0.0);
        }
        for (i = 0; i < DSVariablePoolNumberOfVariables(lowerBounds); i++) {
                
                lowVariable = DSVariablePoolAllVariables(lowerBounds)[i];
                highVariable = DSVariablePoolVariableWithName(upperBounds, DSVariableName(lowVariable));
                if (lowVariable == NULL || highVariable == NULL) {
                        DSError(M_DS_WRONG ": Variables to bound are not consistent", A_DS_WARN);
                        freeVariables = 0;
                        break;
                }
                if (DSVariablePoolHasVariableWithName(DSCaseXd(aCase), DSVariableName(lowVariable)) == true) {
                        variableIndex = DSVariablePoolIndexOfVariableWithName(DSCaseXd(aCase),
                                                                              DSVariableName(lowVariable));
                } else {
                        variableIndex = DSVariablePoolIndexOfVariableWithName(DSCaseXi(aCase),
                                                                              DSVariableName(lowVariable)) + DSVariablePoolNumberOfVariables(DSCaseXd(aCase));
                }
                
                low = DSVariableValue(lowVariable);
                high = DSVariableValue(highVariable);
                
                if (low > high) {
                        DSError(M_DS_WRONG ": Variable bounds are not consistent", A_DS_WARN);
                        freeVariables = 0;
                        break;
                }
                
                if (variableIndex >= DSVariablePoolNumberOfVariables(DSCaseXi(aCase)) + DSVariablePoolNumberOfVariables(DSCaseXd(aCase))) {
                        freeVariables = 0;
                        break;
                }
                if (low == -INFINITY && high == INFINITY)
                        glp_set_col_bnds(linearProblem, variableIndex+1, GLP_FR, 0.0, 0.0);
                else if (low == -INFINITY)
                        glp_set_col_bnds(linearProblem, variableIndex+1, GLP_UP, 0.0, log10(high));
                else if (high == INFINITY)
                        glp_set_col_bnds(linearProblem, variableIndex+1, GLP_LO, log10(low), 0.0);
                else if (low == high)
                        glp_set_col_bnds(linearProblem, variableIndex+1, GLP_FX, log10(low), 0.0);
                else
                        glp_set_col_bnds(linearProblem, variableIndex+1, GLP_DB, log10(low), log10(high));
                if (glp_get_col_type(linearProblem, variableIndex+1) != GLP_FX)
                        freeVariables++;
                
        }
bail:
        return freeVariables;
}

static DSUInteger dsCaseSetVariableBoundsLinearProblem(const DSCase *aCase, glp_prob *linearProblem,  const DSVariablePool * lowerBounds, const DSVariablePool *upperBounds)
{
        DSUInteger i, variableIndex, freeVariables = 0;
        const DSVariable * lowVariable, *highVariable;
        double low, high;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (lowerBounds == NULL && upperBounds == NULL) {
                DSError(M_DS_VAR_NULL ": Variable pool with variables to fix is NULL", A_DS_ERROR);
                goto bail;
        }
        if (DSVariablePoolNumberOfVariables(lowerBounds) != DSVariablePoolNumberOfVariables(upperBounds)) {
                DSError(M_DS_WRONG ": Number of variables to bound must match", A_DS_ERROR);
                goto bail;
        }
        if (linearProblem == NULL) {
                DSError(M_DS_NULL ": Linear problem is NULL", A_DS_WARN);
                goto bail;
        }
        for (i = 0; i < glp_get_num_cols(linearProblem); i++) {
                glp_set_col_bnds(linearProblem, i+1, GLP_FR, 0.0, 0.0);
        }
        for (i = 0; i < DSVariablePoolNumberOfVariables(lowerBounds); i++) {
                
                lowVariable = DSVariablePoolAllVariables(lowerBounds)[i];
                highVariable = DSVariablePoolVariableWithName(upperBounds, DSVariableName(lowVariable));
                if (lowVariable == NULL || highVariable == NULL) {
                        DSError(M_DS_WRONG ": Variables to bound are not consistent", A_DS_WARN);
                        freeVariables = 0;
                        break;
                }
                
                variableIndex = DSVariablePoolIndexOfVariableWithName(DSCaseXi(aCase), 
                                                                      DSVariableName(lowVariable));
                
                low = DSVariableValue(lowVariable);
                high = DSVariableValue(highVariable);
                
                if (low > high) {
                        DSError(M_DS_WRONG ": Variable bounds are not consistent", A_DS_WARN);
                        freeVariables = 0;
                        break;
                }
                
                if (variableIndex >= DSVariablePoolNumberOfVariables(DSCaseXi(aCase))) {
                        freeVariables = 0;
                        break;
                }
                if (low == -INFINITY && high == INFINITY)
                        glp_set_col_bnds(linearProblem, variableIndex+1, GLP_FR, 0.0, 0.0);
                else if (low == -INFINITY)
                        glp_set_col_bnds(linearProblem, variableIndex+1, GLP_UP, 0.0, log10(high));
                else if (high == INFINITY)
                        glp_set_col_bnds(linearProblem, variableIndex+1, GLP_LO, log10(low), 0.0);
                else if (low == high)
                        glp_set_col_bnds(linearProblem, variableIndex+1, GLP_FX, log10(low), 0.0);
                else
                        glp_set_col_bnds(linearProblem, variableIndex+1, GLP_DB, log10(low), log10(high));
                if (glp_get_col_type(linearProblem, variableIndex+1) != GLP_FX)
                        freeVariables++;
                
        }
bail:
        return freeVariables;
}

extern DSVariablePool * DSCaseValidParameterSetAtSlice(const DSCase *aCase, const DSVariablePool * lowerBounds, const DSVariablePool *upperBounds)
{
        bool isValid = false;
        glp_prob *linearProblem = NULL;
        DSVariablePool * Xi = NULL;
        DSUInteger i;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (DSCaseHasSolution(aCase) == false) {
                goto bail;
        }
        if (lowerBounds == NULL || upperBounds == NULL) {
                DSError(M_DS_VAR_NULL ": Variable pool with variables to fix is NULL", A_DS_ERROR);
                goto bail;
        }
        if (DSVariablePoolNumberOfVariables(lowerBounds) != DSVariablePoolNumberOfVariables(upperBounds)) {
                DSError(M_DS_WRONG ": Number of variables to bound must match", A_DS_ERROR);
                goto bail;
        }
        linearProblem = dsCaseLinearProblemForCaseValidity(DSCaseU(aCase), DSCaseZeta(aCase));
        if (linearProblem == NULL) {
                DSError(M_DS_NULL ": Linear problem was not created", A_DS_WARN);
                goto bail;
        }
        if (dsCaseSetVariableBoundsLinearProblem(aCase, linearProblem, lowerBounds, upperBounds) <= DSVariablePoolNumberOfVariables(DSCaseXi(aCase))) {
                glp_simplex(linearProblem, NULL);
                if (glp_get_obj_val(linearProblem) <= -1E-14 && glp_get_prim_stat(linearProblem) == GLP_FEAS)
                        isValid = true;
        }
        if (isValid == true) {
                Xi = DSVariablePoolCopy(DSCaseXi(aCase));
                DSVariablePoolSetReadWrite(Xi);
                for (i = 0; i < DSVariablePoolNumberOfVariables(Xi); i++) {
                        DSVariableSetValue(DSVariablePoolAllVariables(Xi)[i],
                                           pow(10, glp_get_col_prim(linearProblem, i+1)));
                }
        }
        glp_delete_prob(linearProblem);
bail:
        return Xi;
}

extern DSVariablePool * DSCaseValidParameterSetAtSliceByOptimizingFunction(const DSCase *aCase,
                                                                           const DSVariablePool * lowerBounds, const DSVariablePool *upperBounds,
                                                                           const char * function, const bool minimize)
{
        glp_prob *linearProblem = NULL;
        DSVariablePool * Xi = NULL;
        DSUInteger i;
        DSMatrixArray * objective;
        DSMatrix * Oi, *delta, *U, *Zeta;
        DSExpression * expression;
        char * processedFunction;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (DSCaseHasSolution(aCase) == false) {
                goto bail;
        }
        if (lowerBounds == NULL || upperBounds == NULL) {
                DSError(M_DS_VAR_NULL ": Variable pool with variables to fix is NULL", A_DS_ERROR);
                goto bail;
        }
        if (DSVariablePoolNumberOfVariables(lowerBounds) != DSVariablePoolNumberOfVariables(upperBounds)) {
                DSError(M_DS_WRONG ": Number of variables to bound must match", A_DS_ERROR);
                goto bail;
        }
        expression = DSExpressionByParsingString(function);
        if (expression == NULL) {
                DSError(M_DS_NULL ": Could not parse function string", A_DS_ERROR);
                goto bail;
        }
        processedFunction = DSExpressionAsString(expression);
        DSExpressionFree(expression);
        objective = DSCaseParseOptimizationFunction(aCase, processedFunction);
        DSSecureFree(processedFunction);
        if (objective == NULL) {
                goto bail;
        }
        Oi = DSMatrixArrayMatrix(objective, 0);
        delta = DSMatrixArrayMatrix(objective, 1);
        U = DSMatrixCopy(DSCaseU(aCase));
        Zeta = DSMatrixCopy(DSCaseZeta(aCase));
        DSMatrixMultiplyByScalar(U, -1.0);
        linearProblem = dsCaseLinearProblemForMatrices(U, Zeta);
        DSMatrixFree(U);
        DSMatrixFree(Zeta);        if (linearProblem == NULL) {
                DSError(M_DS_NULL ": Linear problem is null", A_DS_WARN);
                DSMatrixArrayFree(objective);
                goto bail;
        }
        if (dsCaseSetVariableBoundsLinearProblem(aCase, linearProblem, lowerBounds, upperBounds) > DSVariablePoolNumberOfVariables(DSCaseXi(aCase))) {
                glp_delete_prob(linearProblem);
                DSMatrixArrayFree(objective);
                goto bail;
        }
        for (i = 0; i < DSMatrixColumns(Oi); i++) {
                glp_set_obj_coef(linearProblem, i+1, DSMatrixDoubleValue(Oi, 0, i));
        }
        if (minimize == false) {
                glp_set_obj_dir(linearProblem, GLP_MAX);
        }
        glp_set_obj_coef(linearProblem, 0, DSMatrixDoubleValue(delta, 0, 0));
        if (linearProblem != NULL) {
                glp_simplex(linearProblem, NULL);
                if (glp_get_status(linearProblem) != GLP_OPT) {
                        glp_delete_prob(linearProblem);
                        DSMatrixArrayFree(objective);
                        goto bail;
                }
                Xi = DSVariablePoolCopy(DSCaseXi(aCase));
                DSVariablePoolSetReadWriteAdd(Xi);
                for (i = 0; i < DSVariablePoolNumberOfVariables(Xi); i++) {
                        DSVariableSetValue(DSVariablePoolAllVariables(Xi)[i], pow(10, glp_get_col_prim(linearProblem, i+1)));
                }
        }
        Xi = DSVariablePoolCopy(DSCaseXi(aCase));
        DSVariablePoolSetReadWrite(Xi);
        for (i = 0; i < DSVariablePoolNumberOfVariables(Xi); i++) {
                DSVariableSetValue(DSVariablePoolAllVariables(Xi)[i],
                                   pow(10, glp_get_col_prim(linearProblem, i+1)));
        }
        DSMatrixArrayFree(objective);
        glp_delete_prob(linearProblem);
bail:
        return Xi;
}

extern const bool DSCaseIsValidAtSlice(const DSCase *aCase, const DSVariablePool * lowerBounds, const DSVariablePool *upperBounds, const bool strict)
{
        bool isValid = false;
        glp_prob *linearProblem = NULL;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (DSCaseHasSolution(aCase) == false) {
                goto bail;
        }
        if (lowerBounds == NULL && upperBounds == NULL) {
                DSError(M_DS_VAR_NULL ": Variable pool with variables to fix is NULL", A_DS_ERROR);
                goto bail;
        }
        if (DSVariablePoolNumberOfVariables(lowerBounds) != DSVariablePoolNumberOfVariables(upperBounds)) {
                DSError(M_DS_WRONG ": Number of variables to bound must match", A_DS_ERROR);
                goto bail;
        }
        if (lowerBounds == upperBounds) {
                isValid = DSCaseIsValidAtPoint(aCase, lowerBounds);
                goto bail;
        }
        if (dsCaseNumberOfFreeVariablesForBounds(aCase, lowerBounds, upperBounds) == 0) {
                isValid = DSCaseIsValidAtPoint(aCase, lowerBounds);
                goto bail;
        }
        linearProblem = dsCaseLinearProblemForCaseValidity(DSCaseU(aCase), DSCaseZeta(aCase));
        if (linearProblem == NULL) {
                DSError(M_DS_NULL ": Linear problem was not created", A_DS_WARN);
                goto bail;
        }
        if (dsCaseSetVariableBoundsLinearProblem(aCase, linearProblem, lowerBounds, upperBounds) <= DSVariablePoolNumberOfVariables(DSCaseXi(aCase))) {
                glp_simplex(linearProblem, NULL);
                if (strict == true) {
                        if (glp_get_obj_val(linearProblem) <= -1E-14 && glp_get_prim_stat(linearProblem) == GLP_FEAS)
                                isValid = true;
                } else {
                        if (glp_get_obj_val(linearProblem) <= 0.0f && glp_get_prim_stat(linearProblem) == GLP_FEAS)
                                isValid = true;
                }
        }
        
        glp_delete_prob(linearProblem);
    
bail:
        return isValid;
}

extern const bool DSCaseIsConsistentAtSlice(const DSCase *aCase, const DSVariablePool * lowerBounds, const DSVariablePool *upperBounds, const bool strict)
{
        bool isValid = false;
        glp_prob *linearProblem = NULL;
        
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (DSCaseHasSolution(aCase) == false) {
                goto bail;
        }
        if (lowerBounds == NULL && upperBounds == NULL) {
                DSError(M_DS_VAR_NULL ": Variable pool with variables to fix is NULL", A_DS_ERROR);
                goto bail;
        }
        if (DSVariablePoolNumberOfVariables(lowerBounds) != DSVariablePoolNumberOfVariables(upperBounds)) {
                DSError(M_DS_WRONG ": Number of variables to bound must match", A_DS_ERROR);
                goto bail;
        }
        if (lowerBounds == upperBounds) {
//                isValid = DSCaseIsConsistentAtPoint(aCase, lowerBounds);
                goto bail;
        }
        if ((dsCaseNumberOfFreeVariablesForBounds(aCase, lowerBounds, upperBounds)+dsCaseNumberOfFreeDependentVariablesForBounds(aCase, lowerBounds, upperBounds)) == 0) {
//                isValid = DSCaseIsConsistentAtPoint(aCase, lowerBounds);
                goto bail;
        }
        DSMatrix * C  = DSMatrixAppendMatrices(DSCaseCd(aCase), DSCaseCi(aCase), true);
        linearProblem = dsCaseLinearProblemForCaseValidity(C, DSCaseDelta(aCase));
        if (linearProblem == NULL) {
                DSError(M_DS_NULL ": Linear problem was not created", A_DS_WARN);
                goto bail;
        }
        if (dsCaseSetDependentAndIndependentVariableBoundsLinearProblem(aCase, linearProblem, lowerBounds, upperBounds) <= DSVariablePoolNumberOfVariables(DSCaseXi(aCase)) + DSVariablePoolNumberOfVariables(DSCaseXd(aCase))) {
                glp_simplex(linearProblem, NULL);
                if (strict == true) {
                        if (glp_get_obj_val(linearProblem) <= -1E-14 && glp_get_prim_stat(linearProblem) == GLP_FEAS)
                                isValid = true;
                } else {
                        if (glp_get_obj_val(linearProblem) <= 0.0f && glp_get_prim_stat(linearProblem) == GLP_FEAS)
                                isValid = true;
                }
        }
        
        glp_delete_prob(linearProblem);
bail:
        return isValid;
}

static DSUInteger nchoosek(DSUInteger n, DSUInteger k)
{
        double denominator, numerator;
        DSUInteger i, lowerBound, upperBound, value = 0;
        if (n == 0 || k == 0)
                goto bail;
        lowerBound = (k > n-k) ? k : (n-k);
        upperBound = (k < n-k) ? k : (n-k);
        denominator = 1;
        for (i = 2; i <= upperBound; i++)
                denominator *= (double)i;
        numerator = lowerBound+1;
        for (i = lowerBound+2; i <= n; i++)
                numerator *= (double)i;
        value = (DSUInteger)(numerator/denominator);
bail:
        return value;
}

static DSVertices * dsCaseCalculateBoundingRange(const DSCase * aCase, glp_prob * linearProblem, const DSUInteger index)
{
        DSVertices *vertices = NULL;
        double minVal, maxVal, val[1] = {INFINITY};
        vertices = DSVerticesAlloc(1);
        glp_set_obj_coef(linearProblem, index+1, 1.0);
        glp_simplex(linearProblem, NULL);
        maxVal = glp_get_obj_val(linearProblem);
        if (glp_get_prim_stat(linearProblem) == GLP_FEAS) {
                val[0] = maxVal;
                DSVerticesAddVertex(vertices, val);
        }
        glp_set_obj_coef(linearProblem, index+1, -1.0);
        glp_simplex(linearProblem, NULL);
        if (glp_get_prim_stat(linearProblem) == GLP_FEAS) {
                minVal = -glp_get_obj_val(linearProblem);
                if (minVal != maxVal) {
                        val[0] = minVal;
                        DSVerticesAddVertex(vertices, val);
                }
        }
bail:
        return vertices;
}

extern DSVertices * DSCaseBoundingRangeForVariableWithConstraints(const DSCase *aCase, const char * variable, DSVariablePool * lowerBounds, DSVariablePool * upperBounds)
{
        DSVertices *vertices = NULL;
        DSUInteger index;
        DSMatrix *A = NULL, *Zeta = NULL, *temp;
        glp_prob * linearProblem = NULL;
        
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        
        index = DSVariablePoolIndexOfVariableWithName(DSCaseXi(aCase), variable);
        
        if (index >= DSVariablePoolNumberOfVariables(DSCaseXi(aCase))) {
                DSError(M_DS_WRONG ": Case does not have variable", A_DS_ERROR);
                goto bail;
        }
        
        temp = DSMatrixCalloc(2, DSVariablePoolNumberOfVariables(DSCaseXi(aCase)));
        DSMatrixSetDoubleValue(temp, 0, index, 1.0);
        DSMatrixSetDoubleValue(temp, 1, index, -1.0);
        A = DSMatrixAppendMatrices(DSCaseU(aCase), temp, false);
        DSMatrixFree(temp);
        temp = DSMatrixCalloc(2, 1);
        DSMatrixSetDoubleValue(temp, 0, 0, 15.0f);
        DSMatrixSetDoubleValue(temp, 1, 0, 15.0f);
        Zeta = DSMatrixAppendMatrices(DSCaseZeta(aCase), temp, false);
        DSMatrixFree(temp);
        DSMatrixMultiplyByScalar(A, -1.0);
        linearProblem = dsCaseLinearProblemForMatrices(A, Zeta);
        
        if (linearProblem == NULL) {
                DSError(M_DS_NULL ": Linear problem is NULL", A_DS_ERROR);
                goto bail;
        }
        
        if (dsCaseSetVariableBoundsLinearProblem(aCase, linearProblem, lowerBounds, upperBounds) == 0) {
                DSError(M_DS_WRONG ": Needs at least one free variables", A_DS_ERROR);
                glp_delete_prob(linearProblem);
                goto bail;
        }
        
        if (glp_get_col_type(linearProblem, DSVariablePoolIndexOfVariableWithName(DSCaseXi(aCase), variable)+1) == GLP_FX) {
                DSError(M_DS_WRONG ": variable is fixed", A_DS_ERROR);
                glp_delete_prob(linearProblem);
                goto bail;
        }
        
        vertices = dsCaseCalculateBoundingRange(aCase, linearProblem, index);
        glp_delete_prob(linearProblem);
bail:
        if (A != NULL)
                DSMatrixFree(A);
        if (Zeta != NULL)
                DSMatrixFree(Zeta);
        return vertices;
}

extern DSVertices * DSCaseBoundingRangeForVariable(const DSCase *aCase, const char * variable)
{
        DSVertices *vertices = NULL;
        DSUInteger i, index;
        DSMatrix *A = NULL, *Zeta = NULL, *temp;
        glp_prob * linearProblem = NULL;
        
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        
        index = DSVariablePoolIndexOfVariableWithName(DSCaseXi(aCase), variable);
        
        if (index >= DSVariablePoolNumberOfVariables(DSCaseXi(aCase))) {
                DSError(M_DS_WRONG ": Case does not have variable", A_DS_ERROR);
                goto bail;
        }
        
        temp = DSMatrixCalloc(2, DSVariablePoolNumberOfVariables(DSCaseXi(aCase)));
        DSMatrixSetDoubleValue(temp, 0, index, 1.0);
        DSMatrixSetDoubleValue(temp, 1, index, -1.0);
        A = DSMatrixAppendMatrices(DSCaseU(aCase), temp, false);
        DSMatrixFree(temp);
        temp = DSMatrixCalloc(2, 1);
        DSMatrixSetDoubleValue(temp, 0, 0, 15.0f);
        DSMatrixSetDoubleValue(temp, 1, 0, 15.0f);
        Zeta = DSMatrixAppendMatrices(DSCaseZeta(aCase), temp, false);
        DSMatrixFree(temp);
        DSMatrixMultiplyByScalar(A, -1.0);
        linearProblem = dsCaseLinearProblemForMatrices(A, Zeta);
        
        if (linearProblem == NULL) {
                DSError(M_DS_NULL ": Linear problem is NULL", A_DS_ERROR);
                goto bail;
        }
        
        for (i = 0; i < DSVariablePoolNumberOfVariables(DSCaseXd(aCase)); i++)
                glp_set_col_bnds(linearProblem, i+1, GLP_FR, 0.0, 0.0);
        
        vertices = dsCaseCalculateBoundingRange(aCase, linearProblem, index);
        glp_delete_prob(linearProblem);
bail:
        if (A != NULL)
                DSMatrixFree(A);
        if (Zeta != NULL)
                DSMatrixFree(Zeta);
        return vertices;
}

static DSVertices * dsCaseCalculate1DVertices(const DSCase * aCase, glp_prob * linearProblem, const DSMatrix * A, const DSMatrix *Zeta, const DSUInteger xIndex, const DSVariablePool * lower, const DSVariablePool * upper)
{
        DSVertices *vertices = NULL;
        double minVal = 0, maxVal = 0, val[1] = {INFINITY};
        vertices = DSVerticesAlloc(1);
        glp_set_obj_coef(linearProblem, xIndex+1, 1.0);
        glp_simplex(linearProblem, NULL);
        if (glp_get_prim_stat(linearProblem) == GLP_FEAS) {
                maxVal = glp_get_obj_val(linearProblem);
                val[0] = maxVal;
                DSVerticesAddVertex(vertices, val);
        }
        glp_set_obj_coef(linearProblem, xIndex+1, -1.0);
        glp_simplex(linearProblem, NULL);
        if (glp_get_prim_stat(linearProblem) == GLP_FEAS) {
                minVal = -glp_get_obj_val(linearProblem);
                if (minVal != maxVal) {
                        val[0] = minVal;
                        DSVerticesAddVertex(vertices, val);
                }
        }
bail:
        return vertices;
}


extern DSVertices * DSCaseVerticesFor1DSlice(const DSCase *aCase, const DSVariablePool * lowerBounds, const DSVariablePool *upperBounds, const char * xVariable)
{
        DSVertices *vertices = NULL;
        DSUInteger xIndex;
        DSMatrix *A, *Zeta, *temp;
        glp_prob * linearProblem = NULL;
        
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        
        if (dsCaseNumberOfFreeVariablesForBounds(aCase, lowerBounds, upperBounds) != 1) {
                DSError(M_DS_WRONG ": Must have only one free variables", A_DS_ERROR);
                goto bail;
        }
        
        xIndex = DSVariablePoolIndexOfVariableWithName(DSCaseXi(aCase), xVariable);
        
        if (xIndex >= DSVariablePoolNumberOfVariables(DSCaseXi(aCase))) {
                DSError(M_DS_WRONG ": Case does not have X variable", A_DS_ERROR);
                goto bail;
        }
        
        temp = DSMatrixCalloc(2, DSVariablePoolNumberOfVariables(DSCaseXi(aCase)));
        DSMatrixSetDoubleValue(temp, 0, xIndex, 1.0);
        DSMatrixSetDoubleValue(temp, 1, xIndex, -1.0);
        A = DSMatrixAppendMatrices(DSCaseU(aCase), temp, false);
        DSMatrixFree(temp);
        temp = DSMatrixCalloc(2, 1);
        DSMatrixSetDoubleValue(temp, 0, 0, -log10(DSVariableValue(DSVariablePoolVariableWithName(lowerBounds, xVariable))));
        DSMatrixSetDoubleValue(temp, 1, 0, log10(DSVariableValue(DSVariablePoolVariableWithName(upperBounds, xVariable))));
        Zeta = DSMatrixAppendMatrices(DSCaseZeta(aCase), temp, false);
        DSMatrixFree(temp);
        DSMatrixMultiplyByScalar(A, -1.0);
        linearProblem = dsCaseLinearProblemForMatrices(A, Zeta);
        
        if (linearProblem == NULL) {
                DSError(M_DS_NULL ": Linear problem is NULL", A_DS_ERROR);
                DSMatrixFree(A);
                DSMatrixFree(Zeta);
                goto bail;
        }
        
        if (dsCaseSetVariableBoundsLinearProblem(aCase, linearProblem, lowerBounds, upperBounds) != 1) {
                DSError(M_DS_WRONG ": Need one free variables", A_DS_ERROR);
                DSMatrixFree(A);
                DSMatrixFree(Zeta);
                glp_delete_prob(linearProblem);
                goto bail;
        }
        
        if (glp_get_col_type(linearProblem, DSVariablePoolIndexOfVariableWithName(DSCaseXi(aCase), xVariable)+1) != GLP_DB) {
                DSError(M_DS_WRONG ": X Variable is not double bound", A_DS_ERROR);
                DSMatrixFree(A);
                DSMatrixFree(Zeta);
                glp_delete_prob(linearProblem);
                goto bail;
        }
        
        vertices = dsCaseCalculate1DVertices(aCase, linearProblem, A, Zeta, xIndex, lowerBounds, upperBounds);
        DSMatrixFree(A);
        DSMatrixFree(Zeta);
        glp_delete_prob(linearProblem);
bail:
        return vertices;
}


static DSVertices * dsCaseCalculate3DVertices(const DSCase * aCase, glp_prob * linearProblem, const DSMatrix * A, const DSMatrix *Zeta, const DSUInteger xIndex, const DSUInteger yIndex, const DSUInteger zIndex)
{
        DSVertices *vertices = NULL;
        DSUInteger numberOfCombinations, i, firstIndex, secondIndex, thirdIndex;
        DSUInteger numberOfBoundaries;
        double bnd, xVal, yVal, zVal, vals[3];
        
        numberOfBoundaries = DSMatrixRows(A);
        numberOfCombinations = 0;
        
        vertices = DSVerticesAlloc(3);
        for (firstIndex = 0; firstIndex < numberOfBoundaries; firstIndex++) {
                for (secondIndex = firstIndex+1; secondIndex < numberOfBoundaries; secondIndex++) {
                        for (thirdIndex = secondIndex+1; thirdIndex < numberOfBoundaries; thirdIndex++) {

                                for (i = 0; i < numberOfBoundaries; i++) {
                                        bnd = glp_get_row_ub(linearProblem, i+1);
                                        glp_set_row_bnds(linearProblem, i+1, GLP_UP, 0.0, bnd);
                                }
                                bnd = glp_get_row_ub(linearProblem, firstIndex+1);
                                glp_set_row_bnds(linearProblem, firstIndex+1, GLP_FX, bnd, bnd);
                                bnd = glp_get_row_ub(linearProblem, secondIndex+1);
                                glp_set_row_bnds(linearProblem, secondIndex+1, GLP_FX, bnd, bnd);
                                bnd = glp_get_row_ub(linearProblem, thirdIndex+1);
                                glp_set_row_bnds(linearProblem, thirdIndex+1, GLP_FX, bnd, bnd);
                                for (i = 0; i < DSVariablePoolNumberOfVariables(DSCaseXi(aCase)); i++)
                                        glp_set_obj_coef(linearProblem, i+1, 0.0);
                                glp_set_obj_coef(linearProblem, xIndex+1, 1.0);
                                glp_simplex(linearProblem, NULL);
                                if (glp_get_prim_stat(linearProblem) != GLP_FEAS)
                                        continue;
                                xVal = glp_get_obj_val(linearProblem);
                                glp_set_obj_coef(linearProblem, xIndex+1, 0.0);
                                glp_set_obj_coef(linearProblem, yIndex+1, 1.0);
                                glp_simplex(linearProblem, NULL);
                                if (glp_get_prim_stat(linearProblem) != GLP_FEAS)
                                        continue;
                                yVal = glp_get_obj_val(linearProblem);
                                glp_set_obj_coef(linearProblem, yIndex+1, 0.0);
                                glp_set_obj_coef(linearProblem, zIndex+1, 1.0);
                                glp_simplex(linearProblem, NULL);
                                if (glp_get_prim_stat(linearProblem) != GLP_FEAS)
                                        continue;
                                zVal = glp_get_obj_val(linearProblem);
                                vals[0] = xVal;
                                vals[1] = yVal;
                                vals[2] = zVal;
                                DSVerticesAddVertex(vertices, vals);
                                numberOfCombinations++;
                        }
                }
        }
        return vertices;
}


static DSVertices * dsCaseCalculate2DVertices(const DSCase * aCase, glp_prob * linearProblem, const DSMatrix * A, const DSMatrix *Zeta, const DSUInteger xIndex, const DSUInteger yIndex)
{
        DSVertices *vertices = NULL;
        DSUInteger numberOfCombinations, i, j, firstIndex, secondIndex;
        DSUInteger numberOfBoundaries, activeIndex;
        double bnd, xVal, yVal, vals[2];
        
        numberOfBoundaries = DSMatrixRows(A);
        numberOfCombinations = nchoosek(numberOfBoundaries, 2);
        
        vertices = DSVerticesAlloc(2);
        for (i = 0; i < numberOfCombinations; i++) {
                
                for (j = 0; j < numberOfBoundaries; j++) {
                        bnd = glp_get_row_ub(linearProblem, j+1);
                        glp_set_row_bnds(linearProblem, j+1, GLP_UP, 0.0, bnd);
                }
                
                secondIndex = 0;
                for (j = 0, firstIndex = 1; j <= i; j += numberOfBoundaries-firstIndex, firstIndex++)
                        secondIndex += numberOfBoundaries-firstIndex;
                secondIndex -= i;
                secondIndex = numberOfBoundaries-secondIndex;
                firstIndex -= 2;
                
                bnd = glp_get_row_ub(linearProblem, firstIndex+1);
                glp_set_row_bnds(linearProblem, firstIndex+1, GLP_FX, bnd, bnd);
                bnd = glp_get_row_ub(linearProblem, secondIndex+1);
                glp_set_row_bnds(linearProblem, secondIndex+1, GLP_FX, bnd, bnd);
                
                for (j = 0; j < DSVariablePoolNumberOfVariables(DSCaseXi(aCase)); j++)
                        glp_set_obj_coef(linearProblem, j+1, 0.0);
                
                glp_set_obj_coef(linearProblem, xIndex+1, 1.0);
                
                if (fabs(DSMatrixDoubleValue(A, firstIndex, yIndex)) >= 1E-14)
                        activeIndex = firstIndex;
                else if (fabs(DSMatrixDoubleValue(A, secondIndex, yIndex)) >= 1E-14)
                        activeIndex = secondIndex;
                else
                        continue;
                
                glp_simplex(linearProblem, NULL);
                if (glp_get_prim_stat(linearProblem) != GLP_FEAS)
                        continue;
                xVal = glp_get_obj_val(linearProblem);
                yVal = -DSMatrixDoubleValue(Zeta, activeIndex, 0);
                for (j = 0; j < DSVariablePoolNumberOfVariables(DSCaseXi(aCase)); j++) {
                        if (j == yIndex)
                                continue;
                        if (j == xIndex) {
                                yVal += DSMatrixDoubleValue(A, activeIndex, j) * xVal;
                        } else {
                                yVal += DSMatrixDoubleValue(A, activeIndex, j) * glp_get_col_ub(linearProblem, j+1);
                        }
                }
                yVal /= -DSMatrixDoubleValue(A, activeIndex, yIndex);
                vals[0] = xVal;
                vals[1] = yVal;
                DSVerticesAddVertex(vertices, vals);
        }
        DSVerticesOrder2DVertices(vertices);
        return vertices;
}

extern DSMatrixArray * DSCaseFacesFor3DSliceAndConnectivity(const DSCase *aCase, const DSVariablePool * lowerBounds, const DSVariablePool *upperBounds, const char * xVariable, const char *yVariable, const char *zVariable)
{
        DSMatrixArray * faces = NULL;
        DSVertices * vertices;
        DSUInteger xIndex, yIndex, zIndex;
        vertices = DSCaseVerticesFor3DSlice(aCase, lowerBounds, upperBounds, xVariable, yVariable, zVariable);
        if (vertices == NULL) {
                goto exit;
        }
        yIndex = DSVariablePoolIndexOfVariableWithName(DSCaseXi(aCase), yVariable);
        xIndex = DSVariablePoolIndexOfVariableWithName(DSCaseXi(aCase), xVariable);
        zIndex = DSVariablePoolIndexOfVariableWithName(DSCaseXi(aCase), zVariable);
        faces = DSVertices3DFaces(vertices, aCase, lowerBounds, upperBounds, xIndex, yIndex, zIndex);
        DSVerticesFree(vertices);
exit:
        return faces;
}


extern DSMatrixArray * DSCaseVerticesFor3DSliceAndConnectivity(const DSCase *aCase, const DSVariablePool * lowerBounds, const DSVariablePool *upperBounds, const char * xVariable, const char *yVariable, const char *zVariable)
{
        DSMatrixArray * verticesAndConnectivity = NULL;
        DSVertices * vertices;
        DSMatrix * connectivity;
        DSUInteger xIndex, yIndex, zIndex;
        vertices = DSCaseVerticesFor3DSlice(aCase, lowerBounds, upperBounds, xVariable, yVariable, zVariable);
        if (vertices == NULL) {
                goto exit;
        }
        yIndex = DSVariablePoolIndexOfVariableWithName(DSCaseXi(aCase), yVariable);
        xIndex = DSVariablePoolIndexOfVariableWithName(DSCaseXi(aCase), xVariable);
        zIndex = DSVariablePoolIndexOfVariableWithName(DSCaseXi(aCase), zVariable);
        connectivity = DSVertices3DConnectivityMatrix(vertices, aCase, lowerBounds, upperBounds, xIndex, yIndex, zIndex);
        verticesAndConnectivity = DSMatrixArrayAlloc();
        DSMatrixArrayAddMatrix(verticesAndConnectivity, DSVerticesToMatrix(vertices));
        DSMatrixArrayAddMatrix(verticesAndConnectivity, connectivity);
        DSVerticesFree(vertices);
exit:
        return verticesAndConnectivity;
}
extern DSVertices * DSCaseVerticesFor3DSlice(const DSCase *aCase, const DSVariablePool * lowerBounds, const DSVariablePool *upperBounds, const char * xVariable, const char *yVariable, const char *zVariable)
{
        DSVertices *vertices = NULL;
        DSUInteger yIndex, xIndex, zIndex;
        DSMatrix *A, *Zeta, *temp;
        glp_prob * linearProblem = NULL;
        
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        
        if (dsCaseNumberOfFreeVariablesForBounds(aCase, lowerBounds, upperBounds) != 3) {
                DSError(M_DS_WRONG ": Must have only three free variables", A_DS_ERROR);
                goto bail;
        }
        
        yIndex = DSVariablePoolIndexOfVariableWithName(DSCaseXi(aCase), yVariable);
        xIndex = DSVariablePoolIndexOfVariableWithName(DSCaseXi(aCase), xVariable);
        zIndex = DSVariablePoolIndexOfVariableWithName(DSCaseXi(aCase), zVariable);
        
        if (xIndex >= DSVariablePoolNumberOfVariables(DSCaseXi(aCase))) {
                DSError(M_DS_WRONG ": Case does not have X variable", A_DS_ERROR);
                goto bail;
        }
        if (yIndex >= DSVariablePoolNumberOfVariables(DSCaseXi(aCase))) {
                DSError(M_DS_WRONG ": Case does not have Y variable", A_DS_ERROR);
                goto bail;
        }
        if (zIndex >= DSVariablePoolNumberOfVariables(DSCaseXi(aCase))) {
                DSError(M_DS_WRONG ": Case does not have Z variable", A_DS_ERROR);
                goto bail;
        }
        temp = DSMatrixCalloc(6, DSVariablePoolNumberOfVariables(DSCaseXi(aCase)));
        DSMatrixSetDoubleValue(temp, 0, xIndex, 1.0);
        DSMatrixSetDoubleValue(temp, 1, xIndex, -1.0);
        DSMatrixSetDoubleValue(temp, 2, yIndex, 1.0);
        DSMatrixSetDoubleValue(temp, 3, yIndex, -1.0);
        DSMatrixSetDoubleValue(temp, 4, zIndex, 1.0);
        DSMatrixSetDoubleValue(temp, 5, zIndex, -1.0);
        
        A = DSMatrixAppendMatrices(DSCaseU(aCase), temp, false);
        DSMatrixFree(temp);
        temp = DSMatrixCalloc(6, 1);
        DSMatrixSetDoubleValue(temp, 0, 0, -log10(DSVariableValue(DSVariablePoolVariableWithName(lowerBounds, xVariable))));
        DSMatrixSetDoubleValue(temp, 1, 0, log10(DSVariableValue(DSVariablePoolVariableWithName(upperBounds, xVariable))));
        DSMatrixSetDoubleValue(temp, 2, 0, -log10(DSVariableValue(DSVariablePoolVariableWithName(lowerBounds, yVariable))));
        DSMatrixSetDoubleValue(temp, 3, 0, log10(DSVariableValue(DSVariablePoolVariableWithName(upperBounds, yVariable))));
        DSMatrixSetDoubleValue(temp, 4, 0, -log10(DSVariableValue(DSVariablePoolVariableWithName(lowerBounds, zVariable))));
        DSMatrixSetDoubleValue(temp, 5, 0, log10(DSVariableValue(DSVariablePoolVariableWithName(upperBounds, zVariable))));
        Zeta = DSMatrixAppendMatrices(DSCaseZeta(aCase), temp, false);
        DSMatrixFree(temp);
        DSMatrixMultiplyByScalar(A, -1.0);
        linearProblem = dsCaseLinearProblemForMatrices(A, Zeta);
        
        if (linearProblem == NULL) {
                DSError(M_DS_NULL ": Linear problem is NULL", A_DS_ERROR);
                DSMatrixFree(A);
                DSMatrixFree(Zeta);
                goto bail;
        }
        
        if (dsCaseSetVariableBoundsLinearProblem(aCase, linearProblem, lowerBounds, upperBounds) != 3) {
                DSError(M_DS_WRONG ": Need three free variables", A_DS_ERROR);
                DSMatrixFree(A);
                DSMatrixFree(Zeta);
                glp_delete_prob(linearProblem);
                goto bail;
        }
        
        if (glp_get_col_type(linearProblem, DSVariablePoolIndexOfVariableWithName(DSCaseXi(aCase), xVariable)+1) == GLP_FX) {
                DSError(M_DS_WRONG ": X Variable is fixed", A_DS_ERROR);
                DSMatrixFree(A);
                DSMatrixFree(Zeta);
                glp_delete_prob(linearProblem);
                goto bail;
        }
        if (glp_get_col_type(linearProblem, DSVariablePoolIndexOfVariableWithName(DSCaseXi(aCase), yVariable)+1) == GLP_FX) {
                DSError(M_DS_WRONG ": Y Variable is fixed", A_DS_ERROR);
                DSMatrixFree(A);
                DSMatrixFree(Zeta);
                glp_delete_prob(linearProblem);
                goto bail;
        }
        if (glp_get_col_type(linearProblem, DSVariablePoolIndexOfVariableWithName(DSCaseXi(aCase), zVariable)+1) == GLP_FX) {
                DSError(M_DS_WRONG ": Z Variable is fixed", A_DS_ERROR);
                DSMatrixFree(A);
                DSMatrixFree(Zeta);
                glp_delete_prob(linearProblem);
                goto bail;
        }
        vertices = dsCaseCalculate3DVertices(aCase, linearProblem, A, Zeta, xIndex, yIndex, zIndex);
        DSMatrixFree(A);
        DSMatrixFree(Zeta);
        glp_delete_prob(linearProblem);
bail:
        return vertices;
}

extern DSVertices * DSCaseVerticesFor2DSlice(const DSCase *aCase, const DSVariablePool * lowerBounds, const DSVariablePool *upperBounds, const char * xVariable, const char *yVariable)
{
        DSVertices *vertices = NULL;
        DSUInteger yIndex, xIndex;
        DSMatrix *A, *Zeta, *temp;
        glp_prob * linearProblem = NULL;
        
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        
        if (dsCaseNumberOfFreeVariablesForBounds(aCase, lowerBounds, upperBounds) != 2) {
                DSError(M_DS_WRONG ": Must have only two free variables", A_DS_ERROR);
                goto bail;
        }
        
        yIndex = DSVariablePoolIndexOfVariableWithName(DSCaseXi(aCase), yVariable);
        xIndex = DSVariablePoolIndexOfVariableWithName(DSCaseXi(aCase), xVariable);
        
        if (xIndex >= DSVariablePoolNumberOfVariables(DSCaseXi(aCase))) {
                DSError(M_DS_WRONG ": Case does not have X variable", A_DS_ERROR);
                goto bail;
        }
        if (yIndex >= DSVariablePoolNumberOfVariables(DSCaseXi(aCase))) {
                DSError(M_DS_WRONG ": Case does not have Y variable", A_DS_ERROR);
                goto bail;
        }
        
        temp = DSMatrixCalloc(4, DSVariablePoolNumberOfVariables(DSCaseXi(aCase)));
        DSMatrixSetDoubleValue(temp, 0, xIndex, 1.0);
        DSMatrixSetDoubleValue(temp, 1, xIndex, -1.0);
        DSMatrixSetDoubleValue(temp, 2, yIndex, 1.0);
        DSMatrixSetDoubleValue(temp, 3, yIndex, -1.0);
        A = DSMatrixAppendMatrices(DSCaseU(aCase), temp, false);
        DSMatrixFree(temp);
        temp = DSMatrixCalloc(4, 1);
        DSMatrixSetDoubleValue(temp, 0, 0, -log10(DSVariableValue(DSVariablePoolVariableWithName(lowerBounds, xVariable))));
        DSMatrixSetDoubleValue(temp, 1, 0, log10(DSVariableValue(DSVariablePoolVariableWithName(upperBounds, xVariable))));
        DSMatrixSetDoubleValue(temp, 2, 0, -log10(DSVariableValue(DSVariablePoolVariableWithName(lowerBounds, yVariable))));
        DSMatrixSetDoubleValue(temp, 3, 0, log10(DSVariableValue(DSVariablePoolVariableWithName(upperBounds, yVariable))));
        Zeta = DSMatrixAppendMatrices(DSCaseZeta(aCase), temp, false);
        DSMatrixFree(temp);
        DSMatrixMultiplyByScalar(A, -1.0);
        linearProblem = dsCaseLinearProblemForMatrices(A, Zeta);
        
        if (linearProblem == NULL) {
                DSError(M_DS_NULL ": Linear problem is NULL", A_DS_ERROR);
                DSMatrixFree(A);
                DSMatrixFree(Zeta);
                goto bail;
        }
        
        if (dsCaseSetVariableBoundsLinearProblem(aCase, linearProblem, lowerBounds, upperBounds) != 2) {
                DSError(M_DS_WRONG ": Need two free variables", A_DS_ERROR);
                DSMatrixFree(A);
                DSMatrixFree(Zeta);
                glp_delete_prob(linearProblem);
                goto bail;
        }
        
        if (glp_get_col_type(linearProblem, DSVariablePoolIndexOfVariableWithName(DSCaseXi(aCase), xVariable)+1) == GLP_FX) {
                DSError(M_DS_WRONG ": X Variable is fixed", A_DS_ERROR);
                DSMatrixFree(A);
                DSMatrixFree(Zeta);
                glp_delete_prob(linearProblem);
                goto bail;
        }
        if (glp_get_col_type(linearProblem, DSVariablePoolIndexOfVariableWithName(DSCaseXi(aCase), yVariable)+1) == GLP_FX) {
                DSError(M_DS_WRONG ": Y Variable is fixed", A_DS_ERROR);
                DSMatrixFree(A);
                DSMatrixFree(Zeta);
                glp_delete_prob(linearProblem);
                goto bail;
        }
        
        vertices = dsCaseCalculate2DVertices(aCase, linearProblem, A, Zeta, xIndex, yIndex);
        DSMatrixFree(A);
        DSMatrixFree(Zeta);
        glp_delete_prob(linearProblem);
bail:
        return vertices;
}

extern DSMatrixArray * DSCaseVerticesForNDSlice(const DSCase *aCase, const DSVariablePool * lowerBounds, const DSVariablePool *upperBounds)
{
        DSUInteger i, j;
        DSVertices * vertices = NULL;
        DSMatrix * vertexMatrix;
        DSMatrix * cobasisMatrix;
        DSMatrixArray * matrixArray, *vertexAndConnectivity = NULL;
        double * coordinates;
        matrixArray = DSCaseNDVertexEnumeration(aCase, lowerBounds, upperBounds);
        if (matrixArray == NULL) {
                goto exit;
        }
        if (DSMatrixArrayNumberOfMatrices(matrixArray) < 2) {
                goto exit;
        }
        vertexMatrix = DSMatrixArrayMatrix(matrixArray, 0);
        cobasisMatrix = DSMatrixArrayMatrix(matrixArray, 1);
        if (vertexMatrix == NULL) {
                goto exit;
        }
        if (cobasisMatrix == NULL) {
                goto exit;
        }
        vertices = DSVerticesAlloc(DSVariablePoolNumberOfVariables(lowerBounds));
        coordinates = DSSecureCalloc(sizeof(double), DSVariablePoolNumberOfVariables(lowerBounds));
        for (i = 0; i < DSMatrixRows(vertexMatrix); i++) {
                for (j = 0; j < DSVariablePoolNumberOfVariables(lowerBounds); j++) {
                        coordinates[j] = DSMatrixDoubleValue(vertexMatrix, i, j);
                }
                DSVerticesAddVertex(vertices, coordinates);
        }
        vertexAndConnectivity = DSMatrixArrayAlloc();
        DSMatrixArrayAddMatrix(vertexAndConnectivity, DSMatrixCopy(vertexMatrix));
        DSMatrixArrayAddMatrix(vertexAndConnectivity, DSVerticesConnectivityMatrix(vertices, aCase, lowerBounds, upperBounds));
        DSSecureFree(coordinates);
        DSMatrixArrayFree(matrixArray);
        DSVerticesFree(vertices);
exit:
        return vertexAndConnectivity;
}

long int DSCaseVerticesNumberOfVerticesForNDSlice(const DSCase *aCase,
                                                  const DSVariablePool * lowerBounds,
                                                  const DSVariablePool *upperBounds,
                                                  const long int maxVertices,
                                                  const bool limitVertices){
        
        long int numberVertices;
    
        numberVertices = DSCaseNDVertexEnumerationNumberOfVertices(aCase,
                                                                   lowerBounds,
                                                                   upperBounds,
                                                                   maxVertices,
                                                                   limitVertices);
    
        return numberVertices;
}

extern DSVertices * DSCaseVerticesForSlice(const DSCase *aCase, const DSVariablePool * lowerBounds, const DSVariablePool *upperBounds, const DSUInteger numberOfVariables, const char **variables)
{
        DSVertices *vertices = NULL;
        DSUInteger i, numberOfFreeVariables = 0;
        
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (lowerBounds == NULL && upperBounds == NULL) {
                DSError(M_DS_VAR_NULL, A_DS_ERROR);
                goto bail;
        }
        if (variables == NULL) {
                DSError(M_DS_NULL ": String with variable names is NULL", A_DS_ERROR);
                goto bail;
        }
        numberOfFreeVariables = dsCaseNumberOfFreeVariablesForBounds(aCase, lowerBounds, upperBounds);
        if (numberOfFreeVariables != numberOfVariables) {
                DSError(M_DS_WRONG ": Number of free variables does not match number of variables", A_DS_ERROR);
                goto bail;
        }
        for (i = 0; i < numberOfFreeVariables; i++) {
                if (variables[i] == NULL) {
                        DSError(M_DS_NULL ": String with variable is NULL", A_DS_ERROR);
                        goto bail;
                }
                if (strlen(variables[i]) == 0) {
                        DSError(M_DS_WRONG ": String with variable is empty", A_DS_ERROR);
                        goto bail;
                }
                if (DSVariablePoolHasVariableWithName(DSCaseXi(aCase), variables[i]) == false) {
                        DSError(M_DS_WRONG ": Case does not have variable for slice", A_DS_ERROR);
                        goto bail;
                }
        }
       
        if (numberOfFreeVariables == 2) {
                vertices = DSCaseVerticesFor2DSlice(aCase, lowerBounds, upperBounds, variables[0], variables[1]);
        } else {
                DSError(M_DS_NOT_IMPL ": N-dimensional vertex enumeration not implemented", A_DS_WARN);
        }
bail:
        return vertices;
}

static gma_parseraux_t * dsCaseParseStringToTermList(const char * string)
{
        void *parser = NULL;
        struct expression_token *tokens, *current;
        gma_parseraux_t *root = NULL, *parser_aux;
        if (string == NULL) {
                DSError(M_DS_WRONG ": String to parse is NULL", A_DS_ERROR);
                goto bail;
        }
        if (strlen(string) == 0) {
                DSError(M_DS_WRONG ": String to parse is empty", A_DS_WARN);
                goto bail;
        }
        tokens = DSExpressionTokenizeString(string);
        if (tokens == NULL) {
                DSError(M_DS_PARSE ": Token stream is NULL", A_DS_ERROR);
                goto bail;
        }
        parser = DSCaseOptimizationFunctionParserAlloc(DSSecureMalloc);//DSGMASystemParserAlloc(DSSecureMalloc);
        root = DSGMAParserAuxAlloc();
        parser_aux = root;
        current = tokens;
        while (current != NULL) {
                if (DSExpressionTokenType(current) == DS_EXPRESSION_TOKEN_START) {
                        current = DSExpressionTokenNext(current);
                        continue;
                }
                DSCaseOptimizationFunctionParser(parser,
                                              DSExpressionTokenType(current),
                                              current,
                                              ((void**)&parser_aux));
                current = DSExpressionTokenNext(current);
        }
        DSCaseOptimizationFunctionParser(parser,
                                      0,
                                      NULL,
                                      ((void **)&parser_aux));
        DSCaseOptimizationFunctionParserFree(parser, DSSecureFree);
        DSExpressionTokenFree(tokens);
        if (DSGMAParserAuxParsingFailed(root) == true) {
                DSGMAParserAuxFree(root);
                root = NULL;
        }
bail:
        return root;
}

//extern void * DSDesignSpaceTermListForAllStrings(char * const * const strings, const DSUInteger numberOfEquations)
//{
//        DSUInteger i;
//        gma_parseraux_t **aux = NULL;
//        DSExpression *expr;
//        char *aString;
//        bool failed = false;
//        aux = DSSecureCalloc(sizeof(gma_parseraux_t *), numberOfEquations);
//        for (i = 0; i < numberOfEquations; i++) {
//                if (strings[i] == NULL) {
//                        DSError(M_DS_WRONG ": String to parse is NULL", A_DS_ERROR);
//                        failed = true;
//                        break;
//                }
//                if (strlen(strings[i]) == 0) {
//                        DSError(M_DS_WRONG ": String to parse is empty", A_DS_ERROR);
//                        failed = true;
//                        break;
//                }
//                expr = DSExpressionByParsingString(strings[i]);
//                if (expr != NULL) {
//                        aString = DSExpressionAsString(expr);
//                        aux[i] = dsDesignSpaceParseStringToTermList(aString);
//                        DSSecureFree(aString);
//                        DSExpressionFree(expr);
//                }
//                if (aux[i] == NULL) {
//                        DSError(M_DS_PARSE ": Expression not in GMA format", A_DS_ERROR);
//                        failed = true;
//                        break;
//                }
//        }
//        if (failed == true) {
//                for (i = 0; i < numberOfEquations; i++)
//                        if (aux[i] != NULL)
//                                DSGMAParserAuxFree(aux[i]);
//                DSSecureFree(aux);
//                aux = NULL;
//        }
//bail:
//        return aux;
//}
//

static void dsCaseOptimizationFunctionProcessExponentBasePairs(const DSCase *aCase, gma_parseraux_t *aux,
                                                               DSMatrix * Od, DSMatrix * Oi, DSMatrix *delta)
{
        DSUInteger j, varIndex;
        const char *varName;
        double currentValue;
        if (aux == NULL) {
                goto bail;
        }
        for (j = 0; j < DSGMAParserAuxNumberOfBases(aux); j++) {
                if (DSGMAParserAuxBaseAtIndexIsVariable(aux, j) == false) {
                        currentValue = DSMatrixDoubleValue(delta, index, 0);
                        currentValue += log10(DSGMAParseAuxsConstantBaseAtIndex(aux, j));
                        DSMatrixSetDoubleValue(delta,
                                               index, 0,
                                               currentValue);
                        continue;
                }
                varName = DSGMAParserAuxVariableAtIndex(aux, j);
                if (DSVariablePoolHasVariableWithName(DSCaseXd(aCase), varName) == true) {
                        varIndex = DSVariablePoolIndexOfVariableWithName(DSCaseXd(aCase), varName);
                        currentValue = DSMatrixDoubleValue(Od, 0, varIndex);
                        currentValue += DSGMAParserAuxExponentAtIndex(aux, j);
                        DSMatrixSetDoubleValue(Od, 0, varIndex, currentValue);
                } else if (DSVariablePoolHasVariableWithName(DSCaseXi(aCase), varName) == true) {
                        varIndex = DSVariablePoolIndexOfVariableWithName(DSCaseXi(aCase), varName);
                        currentValue = DSMatrixDoubleValue(Oi, 0, varIndex);
                        currentValue += DSGMAParserAuxExponentAtIndex(aux, j);
                        DSMatrixSetDoubleValue(Oi, 0, varIndex, currentValue);
                }
        }
bail:
        return;
}

static DSMatrixArray * dsCaseOptimiztionFunctionCreateMatrix(const DSCase *aCase, gma_parseraux_t *aux, bool hasXd)
{
        DSMatrixArray * optimizationMatrices = NULL;
        DSMatrix * Od, *Oi, *delta;
        const DSSSystem * ssystem;
        const DSVariablePool * Xd, *Xi;
        DSMatrix * Ai, * B, * MAi, * Mb;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (aux == NULL) {
                DSError(M_DS_NULL ": Parser auxiliary data is NULL", A_DS_ERROR);
                goto bail;
        }
        Xd = DSCaseXd(aCase);
        Xi = DSCaseXi(aCase);
        if (Xd == NULL || Xi == NULL) {
                DSError(M_DS_WRONG ": Need Xi and Xd", A_DS_ERROR);
                goto bail;
        }
        ssystem = DSCaseSSystem(aCase);
        Od = DSMatrixCalloc(1, DSVariablePoolNumberOfVariables(Xd));
        Oi = DSMatrixCalloc(1, DSVariablePoolNumberOfVariables(Xi));
        delta = DSMatrixCalloc(1, 1);
        dsCaseOptimizationFunctionProcessExponentBasePairs(aCase, aux, Od, Oi, delta);
        if (hasXd == true) {
                if (DSSSystemHasSolution(ssystem) == false) {
                        DSMatrixFree(Od);
                        DSMatrixFree(Oi);
                        DSMatrixFree(delta);
                        goto bail;
                }
                Ai = DSSSystemAi(ssystem);
                B = DSSSystemB(ssystem);
                MAi = DSMatrixByMultiplyingMatrix(DSSSystemM(ssystem), Ai);
                Mb = DSMatrixByMultiplyingMatrix(DSSSystemM(ssystem), B);
                DSMatrixFree(Ai);
                DSMatrixFree(B);
                Ai = MAi;
                MAi = DSMatrixByMultiplyingMatrix(Od, MAi);
                DSMatrixFree(Ai);
                B = Mb;
                Mb = DSMatrixByMultiplyingMatrix(Od, Mb);
                DSMatrixFree(B);
                DSMatrixSubstractByMatrix(Oi, MAi);
                DSMatrixAddByMatrix(delta, Mb);
                DSMatrixFree(MAi);
                DSMatrixFree(Mb);
        }
        DSMatrixFree(Od);
        optimizationMatrices = DSMatrixArrayAlloc();
        DSMatrixArrayAddMatrix(optimizationMatrices, Oi);
        DSMatrixArrayAddMatrix(optimizationMatrices, delta);
bail:
        return optimizationMatrices;
}

extern DSMatrixArray * DSCaseParseOptimizationFunction(const DSCase * aCase, const char * string)
{
        DSMatrixArray * O = NULL;
        DSVariablePool * eqVars = NULL;
        DSExpression * expr;
        bool hasXd = false;
        DSUInteger i;
        char * name;
        if (aCase == NULL) {
                DSError(M_DS_DESIGN_SPACE_NULL, A_DS_ERROR);
                goto bail;
        }
        gma_parseraux_t *aux = NULL;
        aux = dsCaseParseStringToTermList(string);
        if (aux == NULL) {
                goto bail;
        }
        expr = DSExpressionByParsingString(string);
        eqVars = DSExpressionVariablesInExpression(expr);
        for (i = 0; i < DSVariablePoolNumberOfVariables(eqVars); i++) {
                name = DSVariableName(DSVariablePoolVariableAtIndex(eqVars, i));
                if (DSVariablePoolHasVariableWithName(DSCaseXd(aCase), name)) {
                        hasXd = true;
                        break;
                }
        }
        DSExpressionFree(expr);
        DSVariablePoolFree(eqVars);
        O = dsCaseOptimiztionFunctionCreateMatrix(aCase, aux, hasXd);
        DSGMAParserAuxFree(aux);

bail:
        return O;
}

extern DSCaseVolume * DSCaseVolume_lrs(const DSCase *aCase,
                                       const DSVariablePool *lowerBounds,
                                       const DSVariablePool *upperBounds,
                                       const long int maxNumberVertices,
                                       const bool limitVertices,
                                       const bool return_vertices_matrix){
    
    /* The calculation of the volume using the lrs library involves two steps. First, we calculate the vertices of the polytope using the function DSCaseVerticesNDSlice(). In a second step we go back to the H-Representation again using the lrs library. */

    DSCaseVolume *Volume = NULL;
    DSMatrixArray *matrixArray = NULL;
    DSMatrix * Vertices = NULL;
    
    if (aCase == NULL) {
            DSError(M_DS_CASE_NULL, A_DS_ERROR);
            goto bail;
    }
    Volume = DSSecureMalloc(sizeof(DSCaseVolume));
    matrixArray = DSCaseNDVertexEnumerationVertices(aCase,
                                                    lowerBounds,
                                                    upperBounds,
                                                    maxNumberVertices,
                                                    limitVertices,
                                                    Volume);
    if (matrixArray == NULL)
            goto bail;
    
    Vertices = DSMatrixArrayMatrix(matrixArray, 0);

    /* The Vertices matrix has dimensions [Nr. Vertices x Independent Variables] */
    DSCaseFacetsForVertices(aCase, Vertices, maxNumberVertices, limitVertices, Volume);
    
    if (return_vertices_matrix == true){
            Volume->vertices = Vertices;
            dsCaseVolumeSetAverage(lowerBounds, Volume);
    }else{
            DSMatrixArrayFree(matrixArray);
    }
    
    
    
//    dsCallqHull();
    
    DSVariablePool *centroid = NULL;
    centroid = DSCaseCentroid_qhull(aCase,
                                    lowerBounds,
                                    upperBounds,
                                    maxNumberVertices,
                                    limitVertices);
    
bail:
    return Volume;

}

extern DSVariablePool * DSCaseCentroid_qhull(const DSCase *aCase,
                                             const DSVariablePool *lowerBounds,
                                             const DSVariablePool *upperBounds,
                                             const long int maxNumberVertices,
                                             const bool limitVertices){
    
    DSVariablePool *centroid = NULL;
    DSMatrixArray *matrixArray = NULL;
    DSMatrix * Vertices = NULL;
    DSCaseVolume *Volume = NULL;
    coordT *points = NULL;

     
     if (aCase == NULL) {
             DSError(M_DS_CASE_NULL, A_DS_ERROR);
             goto bail;
     }
     Volume = DSSecureMalloc(sizeof(DSCaseVolume));
     matrixArray = DSCaseNDVertexEnumerationVertices(aCase,
                                                     lowerBounds,
                                                     upperBounds,
                                                     maxNumberVertices,
                                                     limitVertices,
                                                     Volume);
     if (matrixArray == NULL)
             goto bail;
     
     Vertices = DSMatrixArrayMatrix(matrixArray, 0);
    
    /*
     The next steps for the calculation include:
     
     1. Transform Matrix [numpoints x dim] into vector[numpoints x (dim + 1)] for the delaunay trianqulation
     
     2. Call delaunay triangulation.
     
     3. Call function within the macros FORALLfacets & FOREACHvertex_ that calculates centroid and updates the total volume.
     
     */
    
    points = DSMatrixToArray(Vertices, true);
    centroid = DSVariablePoolCopy(DSCaseXi(aCase));
    DSVariablePoolSetReadWriteAdd((DSVariablePool *) centroid);
    dsCaseCalculateCentroid_qhull(points,
                                  DSMatrixColumns(Vertices)+1,
                                  DSMatrixRows(Vertices),
                                  centroid);
    
    DSMatrixArrayFree(matrixArray);
    DSSecureFree(Volume);
    DSSecureFree(points);

bail:
    
    return centroid;
}

void dsCaseCalculateCentroid_qhull(coordT *points,
                                   DSUInteger dim,
                                   DSUInteger numpoints,
                                   DSVariablePool *centroid){
    

        boolT ismalloc= False;    /* True if qhull should free points in qh_freeqhull() or reallocation */
        char flags[250];          /* option flags for qhull, see qh-quick.htm */
        FILE *errfile= stderr;    /* error messages from qhull code */
        int exitcode;             /* 0 if no error from qhull */
        facetT *facet;            /* set by FORALLfacets */
        int curlong, totlong;     /* memory remaining after qh_memfreeshort, used if !qh_NOmem  */
        int i=0;
        vertexT *vertex, **vertexp;
        double total_volume = 0.0, facet_volume = 0.0;
        coordT *case_centroid = NULL;
        coordT *facet_centroid = NULL;
        coordT *facet_points = NULL;
        int vertex_count = 0;
        const DSVariable ** variables;
    
        qhT qh_qh;    /* Create a new instance of Qhull (qhB) */
        qhT *qh = &qh_qh;

        QHULL_LIB_CHECK
        qh_zero(qh, errfile);
        sprintf(flags, "qhull s Tcv i Fa FA d QJ Pp"); // for delaunay computations
//        sprintf(flags, "qhull s i FA QJ"); // for volume computations

        qh_init_A(qh, stdin, stdout, stderr, 0, NULL);
        exitcode = setjmp(qh->errexit);
        fflush(NULL);
        qh->NOerrexit= False;
        if (!exitcode) {
                qh_initflags(qh, flags);
                qh_setdelaunay(qh, dim, numpoints, points);
                qh_init_B(qh, points, numpoints, dim, ismalloc);
                qh_qhull(qh);
                qh_check_output(qh);
                fflush(NULL);
//                qh_produce_output(qh);  /* delete this line to help avoid io_r.c */
                qh_prepare_output(qh);
                if (qh->VERIFYoutput && !qh->FORCEoutput && !qh->STOPadd && !qh->STOPcone && !qh->STOPpoint)
                  qh_check_points(qh);
                fflush(NULL);
                
                /* Define an array containing the case centroid */
                case_centroid = DSSecureCalloc(sizeof(coordT), dim-1);
                
                FORALLfacets {
                    if(facet->good){
                        /* 1. Initialice point array that will contain the vertices of the facet */
                        facet_points = DSSecureMalloc(sizeof(coordT)*(dim*(dim-1)));
                        /* 2. Initialice a variable pool that will contain the centroid of this facet */
                        facet_centroid = DSSecureCalloc(sizeof(coordT), dim-1);
                        vertex_count = 0;
                                FOREACHvertex_(facet->vertices){
                                    
                                    for (i=0; i<(dim-1); i++){
                                        /* Fill out the facet_points array with information of each vertex*/
                                        facet_points[i + vertex_count*(dim-1)] = vertex->point[i];
                                        
                                        /* Update facet_centroid array containing the centroid of the facet*/
                                        facet_centroid[i] += vertex->point[i];
                                    }
                                vertex_count++;
                                }
                        
                        /* Get the volume of the respective facet. Update the total volume */
                        facet_volume = DSVolumefromArray_qhull(facet_points, dim, dim-1);
                        total_volume += facet_volume;

                        for(i=0; i<(dim-1); i++){
                            /* scale the centroid of the facet by dividing by the number of points and multiplying by the volume of the facet */
                            facet_centroid[i] = (facet_centroid[i]*facet_volume)/(dim);
                            /* update the centroid case by adding the facet centroid to the case centroid*/
                            case_centroid[i] += facet_centroid[i];
                        }
                    
                        if (facet_centroid != NULL)
                            DSSecureFree(facet_centroid);
                        if (facet_points != NULL)
                            DSSecureFree(facet_points);
                    }
                }
            /* Scale the case centroid */
            
            variables = DSVariablePoolAllVariables(centroid);
            for(i=0; i<(dim-1); i++){
                
                case_centroid[i] = case_centroid[i]/total_volume;
                DSVariablePoolSetValueForVariableWithName(centroid,
                                                          DSVariableName(variables[i]),
                                                          case_centroid[i]);
            }
            
            if (case_centroid !=NULL)
                DSSecureFree(case_centroid);
        }
        qh->NOerrexit= True;
        #ifdef qh_NOmem
          qh_freeqhull(qh, qh_ALL);
        #else
          qh_freeqhull(qh, !qh_ALL);
          qh_memfreeshort(qh, &curlong, &totlong);
          if (curlong || totlong)
            fprintf(stderr, "qhull warning (user_eg2, run 2): did not free %d bytes of long memory (%d pieces)\n",
                 totlong, curlong);
        #endif
    
}


void dsCallqHull(void){
    
    int DIM = 3;
    int TOTpoints = 8;

    int dim= DIM;             /* dimension of points */
    int numpoints = TOTpoints;            /* number of points */
    coordT points[DIM*TOTpoints]; /* array of coordinates for each point */
    boolT ismalloc= False;    /* True if qhull should free points in qh_freeqhull() or reallocation */
    char flags[250];          /* option flags for qhull, see qh-quick.htm */
    FILE *outfile= stdout;    /* output from qh_produce_output()
                                 use NULL to skip qh_produce_output() */
    FILE *errfile= stderr;    /* error messages from qhull code */
    int exitcode;             /* 0 if no error from qhull */
    facetT *facet;            /* set by FORALLfacets */
    int curlong, totlong;     /* memory remaining after qh_memfreeshort, used if !qh_NOmem  */
    int i=0;
    double value;
    vertexT *vertex, **vertexp;

    qhT qh_qh;    /* Create a new instance of Qhull (qhB) */
    qhT *qh = &qh_qh;

    QHULL_LIB_CHECK
    qh_zero(qh, errfile);
//    sprintf(flags, "qhull s Tcv i Fa FA d QJ"); // for delaunay computations
    sprintf(flags, "qhull s i FA QJ"); // for volume computations

    
    qh_init_A(qh, stdin, stdout, stderr, 0, NULL);
    exitcode = setjmp(qh->errexit);
    fflush(NULL);
    qh->NOerrexit= False;
    if (!exitcode) {
        qh_initflags(qh, flags);
        
//        // Define the array. Vertices of case 15 [-3 to 3]
//        points[0] = -1.00031292439;
//        points[1] = -1.875;
//        points[2] = 0;
//
//        points[3] = -3.0;
//        points[4] = -1.875;
//        points[5] = 0;
//
//        points[6] = -3.0;
//        points[7] = -2.0;
//        points[8] = 0;
//
//        points[9] = -0.50031286478;
//        points[10] = -2.0;
//        points[11] = 0;
        
//        // Define the array. Vertices of case 11 [-4 to 4]
//        points[0] = -1.0;
//        points[1] = 4.0;
//        points[2] = 0;
//
//        points[3] = -4.0;
//        points[4] = 4.0;
//        points[5] = 0;
//
//        points[6] = -4.0;
//        points[7] = -1.87512516975;
//        points[8] = 0;
//
//        points[9] = -1.0;
//        points[10] = -1.87512516975;
//        points[11] = 0;
        
        
//        // Define the array. Vertices of case 1 [-4 to 4] @ MOTIF1-3D. Ro = 100
//        points[0] = 10000.0 ;
//        points[1] = 10000.0 ;
//        points[2] = 10000.0;
//        points[3] = 0;
//
////         point 2
//        points[4] = 10000.0;
//        points[5] = 1.0;
//        points[6] = 10000.0;
//        points[7] = 0;
//
////         point 3
//        points[8] = 1.0;
//        points[9] = 1.0;
//        points[10] = 10000.0;
//        points[11] = 0;
//
////        point 4
//        points[12] = 1.0;
//        points[13] = 10000.0;
//        points[14] = 10000.0;
//        points[15] = 0;
//
////        point 5
//        points[16] = 10000.0;
//        points[17] = 10000.0;
//        points[18] = 0.01;
//        points[19] = 0;
//
////        point 6
//        points[20] = 10000.0;
//        points[21] = 1.0;
//        points[22] = 0.01;
//        points[23] = 0;
//
////        point 7
//        points[24] = 1.0;
//        points[25] = 1.0;
//        points[26] = 0.01;
//        points[27] = 0;
//
////        point 8
//        points[28] = 1.0;
//        points[29] = 10000.0;
//        points[30] = 0.01;
//        points[31] = 0;
        
        
        // Define the array. Vertices of case 1 [-4 to 4] @ MOTIF1-3D. Ro = 100
        points[0] = 4.0 ;
        points[1] = 4.0 ;
        points[2] = 4.0;
        
//         point 2
        points[3] = 4.0;
        points[4] = 0.0;
        points[5] = 4.0;
        
//         point 3
        points[6] = 0;
        points[7] = 0;
        points[8] = 4.0;
        
//        point 4
        points[9] = 0.0;
        points[10] = 4.0;
        points[11] = 4.0;
        
//        point 5
        points[12] = 4.0;
        points[13] = 4.0;
        points[14] = -2.0;
        
//        point 6
        points[15] = 4.0;
        points[16] = 0.0;
        points[17] = -2.0;
        
//        point 7
        points[18] = 0.0;
        points[19] = 0.0;
        points[20] = -2.0;
        
//        point 8
        points[21] = 0.0;
        points[22] = 4.0;
        points[23] = -2.0;
    

//        qh_setdelaunay(qh, dim, numpoints, points);
        qh_init_B(qh, points, numpoints, DIM, ismalloc);
        qh_qhull(qh);
        qh_check_output(qh);
        fflush(NULL);

        qh_produce_output(qh);  /* delete this line to help avoid io_r.c */
        if (qh->VERIFYoutput && !qh->FORCEoutput && !qh->STOPadd && !qh->STOPcone && !qh->STOPpoint)
          qh_check_points(qh);
        fflush(NULL);
        
        
        // Playing a bit
        
        printf("The total area of the facet is: %f. The total volume is: %f \n ",
               qh->totarea,
               qh->totvol);
        
        printf("The number of Delaunay regions is: %u \n", qh->num_good);
        printf("About to print good facets \n");
        FORALLfacets {
            if(facet->good){
                    
                    printf("Facet %d. Area = %f \n", i, facet->f.area);
                    printf("The maximum number of elements of the vertices is %d \n", facet->vertices->maxsize);
                
                    FOREACHvertex_(facet->vertices){
                            printf("The id for the vertex is %d \n", vertex->id);
                            printf("The coordinates of this vertex are: %f %f %f \n",
                                  vertex->point[0],
                                  vertex->point[1],
                                  vertex->point[2]);
                    }

                    i++;
            }
        }
            
    }
    qh->NOerrexit= True;
    #ifdef qh_NOmem
      qh_freeqhull(qh, qh_ALL);
    #else
      qh_freeqhull(qh, !qh_ALL);
      qh_memfreeshort(qh, &curlong, &totlong);
      if (curlong || totlong)
        fprintf(stderr, "qhull warning (user_eg2, run 2): did not free %d bytes of long memory (%d pieces)\n",
             totlong, curlong);
    #endif
    
    
    printf("Testing routine DSVolumefromFloatArray_qhull \n");
    float volume;
    volume = DSVolumefromArray_qhull(points, numpoints, DIM);
    printf("The volume calculated from the routine is %f \n", volume);
    printf("End of dsCallqHull \n");
    
}



void dsCaseVolumeSetAverage(const DSVariablePool * variables,
                            DSCaseVolume * Volume ){
    
    
    DSMatrix * vertices;
    vertices = Volume->vertices;
    
    
    if (vertices == NULL || variables == NULL || Volume == NULL){
            DSError("NULL pointer ", A_DS_ERROR);
            goto bail;
    }
    
    DSMatrix * averages = NULL;
    DSInteger i;
    DSVariablePool * aux = NULL;
    const DSVariable ** var;
    double value = 0.0;
    
    averages = DSMatrixAverage(vertices, false);
    
    if (DSMatrixColumns(averages) != DSVariablePoolNumberOfVariables(variables)){
            DSError("The number of dimensions is not consistent with the number of variables ", A_DS_ERROR);
            goto bail;
    }
    
    var = DSVariablePoolAllVariables(variables);
    aux = DSVariablePoolCopy(variables);
    
    for (i=0; i<DSVariablePoolNumberOfVariables(variables); i++){
        value = DSMatrixDoubleValue(averages, 0, i);
        DSVariablePoolSetValueForVariableWithName(aux,
                                                  DSVariableName(var[i]),
                                                  value);
    }
    
    Volume->average = aux;
    
bail:
    return;
}

extern double DSCaseVolumeGetVolume(const DSCaseVolume * Volume_str){
    
    double volume = 0;
    
    if (Volume_str == NULL) {
            DSError("Pointer is null!", A_DS_ERROR);
            goto bail;
    }
    
    volume = Volume_str->volume;
    
bail:
    return volume;
    
    
}

extern double DSCaseVolumeGetVertices(const DSCaseVolume * Volume_str){
    
    double vertices = 0;
    
    if (Volume_str == NULL) {
            DSError("Pointer is null!", A_DS_ERROR);
            goto bail;
    }
    
    vertices = Volume_str->nr_vertices;
    
bail:
    return vertices;
}

extern DSMatrix * DSCaseVolumeGetVerticesMatrix(const DSCaseVolume * Volume_str){
    
    DSMatrix * vertices = NULL;
    
    if (Volume_str == NULL) {
            DSError("Pointer is null!", A_DS_ERROR);
            goto bail;
    }
    
    vertices = Volume_str->vertices;
    
bail:
    return vertices;
}

extern DSVariablePool * DSCaseVolumeGetOperatingPoint2D(const DSCaseVolume *Volume_str){
    
    DSVariablePool *operatingPoint = NULL;
    
    
    if (Volume_str == NULL) {
            DSError("Pointer is null!", A_DS_ERROR);
            goto bail;
    }
    
    operatingPoint = Volume_str->average;
bail:
    return operatingPoint;
}


extern double DSCaseDimension(const DSCase *aCase,
                              const DSVariablePool *lowerBounds,
                              const DSVariablePool *upperBounds){
    
    /* This function calculates the dimensionality of the vertices of the polytope  */

    DSMatrixArray *matrixArray = NULL;
    DSMatrix * Vertices = NULL;
    DSMatrix *aux = NULL;
    double dimension = 0.0;
    
    if (aCase == NULL) {
            DSError(M_DS_CASE_NULL, A_DS_ERROR);
            goto bail;
    }
        
    matrixArray = DSCaseNDVertexEnumeration(aCase, lowerBounds, upperBounds);
        
    if (matrixArray == NULL) {
            goto bail;
    }
    if (DSMatrixArrayNumberOfMatrices(matrixArray) < 2) {
            goto bail;
    }
    
    Vertices = DSMatrixArrayMatrix(matrixArray, 0);
    
//    printf("Reporting from function DSCaseDimension(). for case %s. The vertices of this case are: \n ", aCase->caseIdentifier);
//    DSMatrixPrint(Vertices);

    if (DSMatrixColumns(Vertices) > DSMatrixRows(Vertices))
        aux = DSMatrixTranspose(Vertices);
    else
        aux = Vertices;
    
    if (aux != NULL)
        dimension = DSMatrixRank(aux);
    else
        goto bail;
        
    if (matrixArray != NULL)
        DSMatrixArrayFree(matrixArray);
    
bail:
    return dimension;

}


#if defined (__APPLE__) && defined (__MACH__)
#pragma mark Intersection of cases
#endif


extern const bool DSCaseIntersectionListIsValid(const DSUInteger numberOfCases, const DSCase *firstCase, ...)
{
        bool isValid = false;
        DSUInteger i;
        const DSCase ** cases = NULL;
        va_list ap;
        if (firstCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (numberOfCases == 0) {
                DSError(M_DS_WRONG ": Number of cases must be at least one", A_DS_WARN);
                goto bail;
        }
        va_start(ap, firstCase);
        cases = DSSecureMalloc(sizeof(DSCase *)*numberOfCases);
        cases[0] = firstCase;
        for (i = 1; i < numberOfCases; i++) {
                cases[i] = va_arg(ap, DSCase *);
                if (cases[i] == NULL) {
                        DSError(M_DS_CASE_NULL, A_DS_ERROR);
                        break;
                }
        }
        va_end(ap);
        if (i == numberOfCases)
                isValid = DSCaseIntersectionIsValid(numberOfCases, cases);
        DSSecureFree(cases);

bail:
        return isValid;
}

#if defined (__APPLE__) && defined (__MACH__)
#pragma mark Pseudocase with intersection of cases
#endif

///**
// *
// */
//static DSPseudoCase * dsPseudoCaseFromConcurrentCasisInSlice(const DSUInteger numberOfCases, const DSCase ** cases, )
//{
//        DSUInteger i;
//        DSPseudoCase * caseIntersection = NULL;
//        DSMatrix *U = NULL, *Zeta = NULL, *temp;
//        if (numberOfCases == 0) {
//                DSError(M_DS_WRONG ": Number of cases must be at least one", A_DS_ERROR);
//                goto bail;
//        }
//        if (cases == NULL) {
//                DSError(M_DS_NULL ": Array of cases is NULL", A_DS_ERROR);
//                goto bail;
//        }
//        for (i = 0; i < numberOfCases; i++) {
//                if (DSCaseHasSolution(cases[i]) == false)
//                        goto bail;
//        }
//        U = DSMatrixCopy(DSCaseU(cases[0]));
//        Zeta = DSMatrixCopy(DSCaseZeta(cases[0]));
//        for (i = 1; i < numberOfCases; i++) {
//                temp = DSMatrixAppendMatrices(U, DSCaseU(cases[i]), false);
//                DSMatrixFree(U);
//                U = temp;
//                temp = DSMatrixAppendMatrices(Zeta, DSCaseZeta(cases[i]), false);
//                DSMatrixFree(Zeta);
//                Zeta = temp;
//                if (U == NULL || Zeta == NULL)
//                        goto bail;
//        }
//        caseIntersection = DSSecureCalloc(1, sizeof(DSCase));
//        DSCaseXd(caseIntersection) = DSCaseXd(cases[0]);
//        DSCaseXi(caseIntersection) = DSCaseXi(cases[0]);
//        DSCaseU(caseIntersection) = U;
//        DSCaseZeta(caseIntersection) = Zeta;
//        U = NULL;
//        Zeta = NULL;
//bail:
//        if (U != NULL)
//                DSMatrixFree(U);
//        if (Zeta != NULL)
//                DSMatrixFree(Zeta);
//        return caseIntersection;
//}

/**
 * 
 */
extern DSPseudoCase * DSPseudoCaseFromIntersectionOfCases(const DSUInteger numberOfCases, const DSCase ** cases)
{
        DSUInteger i;
        DSPseudoCase * caseIntersection = NULL;
        DSMatrix *Cd = NULL, *Ci = NULL, *delta = NULL, *U = NULL, *Zeta = NULL, *temp = NULL;
        char caseIdentifier[1000];
        if (numberOfCases == 0) {
                DSError(M_DS_WRONG ": Number of cases must be at least one", A_DS_ERROR);
                goto bail;
        }
        if (cases == NULL) {
                DSError(M_DS_NULL ": Array of cases is NULL", A_DS_ERROR);
                goto bail;
        }
        for (i = 0; i < numberOfCases; i++) {
                if (cases[i] == NULL) {
                        DSError(M_DS_CASE_NULL, A_DS_ERROR);
                        goto bail;
                }
                if (DSCaseHasSolution(cases[i]) == false)
                        goto bail;
        }
        U = DSMatrixCopy(DSCaseU(cases[0]));
        Zeta = DSMatrixCopy(DSCaseZeta(cases[0]));
        // These checks were added because when analyzing blowing behavior of a single s-system matrices Ci, Cd and Delta do not exist.
        if (DSCaseCi(cases[0]) != NULL)
            Ci = DSMatrixCopy(DSCaseCi(cases[0]));
        if (DSCaseCd(cases[0]) != NULL)
            Cd = DSMatrixCopy(DSCaseCd(cases[0]));
        if (DSCaseDelta(cases[0]) != NULL)
            delta = DSMatrixCopy(DSCaseDelta(cases[0]));
        sprintf(caseIdentifier, "%s", DSCaseIdentifier(cases[0]));
        for (i = 1; i < numberOfCases; i++) {
                sprintf(caseIdentifier, "%s, %s", caseIdentifier, DSCaseIdentifier(cases[i]));
            
                if (U != NULL && DSCaseU(cases[i]) != NULL ){
                    temp = DSMatrixAppendMatrices(U, DSCaseU(cases[i]), false);
                    DSMatrixFree(U);
                    U = temp;
                }
            
                if (Zeta != NULL && DSCaseZeta(cases[i]) != NULL ){
                    temp = DSMatrixAppendMatrices(Zeta, DSCaseZeta(cases[i]), false);
                    DSMatrixFree(Zeta);
                    Zeta = temp;
                }
            
                if (Cd != NULL && DSCaseCd(cases[i]) != NULL ){
                    temp = DSMatrixAppendMatrices(Cd, DSCaseCd(cases[i]), false);
                    DSMatrixFree(Cd);
                    Cd = temp;
                }
            
                if (Ci != NULL && DSCaseCi(cases[i]) != NULL){
                    temp = DSMatrixAppendMatrices(Ci, DSCaseCi(cases[i]), false);
                    DSMatrixFree(Ci);
                    Ci = temp;
                }
            
                if (delta != NULL && DSCaseDelta(cases[i])){
                    temp = DSMatrixAppendMatrices(delta, DSCaseDelta(cases[i]), false);
                    DSMatrixFree(delta);
                    delta = temp;
                }
            
//                if (U == NULL || Zeta == NULL || Cd == NULL || Ci == NULL || delta == NULL)
//                        goto bail;
            
//               This constraint was relaxed to allow for analysis of subcases of single S-systems (without Cd, Ci oder     delta matrices).
            
                  if (U == NULL && Zeta == NULL && Cd == NULL && Ci == NULL && delta == NULL)
                        goto bail;
            
        }
        caseIntersection = DSSecureCalloc(1, sizeof(DSCase));
        caseIntersection->freeVariables = true;
        caseIntersection->Xd_a = DSVariablePoolCopy(DSCaseXd_a(cases[0]));
        caseIntersection->Xd = DSVariablePoolCopy(DSCaseXd(cases[0]));
        caseIntersection->Xi = DSVariablePoolCopy(DSCaseXi(cases[0]));
        DSCaseU(caseIntersection) = U;
        DSCaseZeta(caseIntersection) = Zeta;
        DSCaseCi(caseIntersection) = Ci;
        DSCaseCd(caseIntersection) = Cd;
        DSCaseDelta(caseIntersection) = delta;
        DSCaseId(caseIntersection) = strdup(caseIdentifier);
        U = NULL;
        Zeta = NULL;
        Cd = NULL;
        Ci = NULL;
        delta = NULL;
        caseIntersection->ssys = NULL;
bail:
        if (U != NULL)
                DSMatrixFree(U);
        if (Zeta != NULL)
                DSMatrixFree(Zeta);
        if (Ci != NULL)
                DSMatrixFree(Ci);
        if (Cd != NULL)
                DSMatrixFree(Cd);
        if (delta != NULL)
                DSMatrixFree(delta);
        return caseIntersection;
}

/**
 *
 */
extern DSPseudoCase * DSPseudoCaseFromIntersectionOfCasesExcludingSlice(const DSUInteger numberOfCases, const DSCase ** cases, const DSUInteger numberOfExceptions, const char ** exceptionVarNames)
{
        DSUInteger i, j, k, l, currentRow, rows, columns_i, columns_d, *indices_i, *indices_d;
        DSUInteger newXi, newXd, numberExtra_i, numberExtra_d;
        DSPseudoCase * caseIntersection = NULL;
        DSMatrix * Cd = NULL, *Ci = NULL, *delta = NULL;
        DSMatrix *U = NULL, *Zeta = NULL, *tempCd, *tempU, *tempZeta;
        DSVariablePool *Xd, *Xi, *Xd_a;
        char * name = NULL;
        const char ** variableNames_i, ** variableNames_d;
        if (numberOfCases == 0) {
                DSError(M_DS_WRONG ": Number of cases must be at least one", A_DS_ERROR);
                goto bail;
        }
        if (cases == NULL) {
                DSError(M_DS_NULL ": Array of cases is NULL", A_DS_ERROR);
                goto bail;
        }
        for (i = 0; i < numberOfCases; i++) {
                if (DSCaseHasSolution(cases[i]) == false)
                        goto bail;
        }
        Xd_a = DSVariablePoolAlloc();
        Xd = DSVariablePoolAlloc();
        Xi = DSVariablePoolAlloc();
        variableNames_i = DSVariablePoolAllVariableNames(DSCaseXi(cases[0]));
        variableNames_d = DSVariablePoolAllVariableNames(DSCaseXd(cases[0]));
        indices_d = DSSecureCalloc(numberOfExceptions, sizeof(DSUInteger));
        indices_i = DSSecureCalloc(numberOfExceptions, sizeof(DSUInteger));
        newXi = 0;
        newXd = 0;
        k = 0;
        l = 0;
        for (j = 0; j < numberOfExceptions; j++) {
                if (DSVariablePoolHasVariableWithName(DSCaseXd(cases[0]), exceptionVarNames[j]))
                        indices_d[newXd++] = DSVariablePoolIndexOfVariableWithName(DSCaseXd(cases[0]), exceptionVarNames[j]);
                if (DSVariablePoolHasVariableWithName(DSCaseXi(cases[0]), exceptionVarNames[j]))
                        indices_i[newXi++] = DSVariablePoolIndexOfVariableWithName(DSCaseXi(cases[0]), exceptionVarNames[j]);
                for (i = 0; i < numberOfCases; i++) {
                        if (DSVariablePoolHasVariableWithName(DSCaseXd(cases[i]), exceptionVarNames[j])) {
                                continue;
                        } else if (DSVariablePoolHasVariableWithName(DSCaseXi(cases[i]), exceptionVarNames[j])) {
                                continue;
                        } else {
                                DSError(M_DS_WRONG ": Case does not have variable to except", A_DS_ERROR);
                                DSSecureFree(indices_i);
                                DSSecureFree(indices_d);
                                goto bail;
                        }
                }
        }
//        k = 0;
        name = DSSecureCalloc(sizeof(char), 200);
        for (i = 0; i < DSVariablePoolNumberOfVariables(DSCaseXi(cases[0])); i++) {
                for (j = 0; j < newXi; j++) {
                        if (i == indices_i[j])
                                break;
                }
                if (j == newXi) {
                        DSVariablePoolAddVariableWithName(Xi, variableNames_i[i]);
                } else {
                        sprintf(name, "$%s_0", variableNames_i[i]);
                        DSVariablePoolAddVariableWithName(Xi, name);
                }
        }
        for (i = 0; i < DSVariablePoolNumberOfVariables(DSCaseXd(cases[0])); i++) {
                for (j = 0; j < newXd; j++) {
                        if (i == indices_d[j])
                                break;
                }
                if (j == newXd) {
                        DSVariablePoolAddVariableWithName(Xd, variableNames_d[i]);
                        if (DSVariablePoolHasVariableWithName(DSCaseXd_a(cases[0]), variableNames_d[i]))
                                DSVariablePoolAddVariableWithName(Xd_a, variableNames_d[i]);
                } else {
                        sprintf(name, "$%s_0", variableNames_d[i]);
                        DSVariablePoolAddVariableWithName(Xd, name);
                        if (DSVariablePoolHasVariableWithName(DSCaseXd_a(cases[0]), variableNames_d[i]))
                                DSVariablePoolAddVariableWithName(Xd_a, name);

                }
        }
        k = 0;
        for (i = 0; i < newXi*(numberOfCases-1); i++) {
                j = i % newXi;
                sprintf(name, "$%s_%i", variableNames_i[indices_i[j]], i / newXi + 1);
                DSVariablePoolAddVariableWithName(Xi, name);
        }
        for (i = 0; i < newXd*(numberOfCases-1); i++) {
                j = i % newXd;
                sprintf(name, "$%s_%i", variableNames_d[indices_d[j]], i / newXd + 1);
                if (DSVariablePoolHasVariableWithName(DSCaseXd_a(cases[0]), variableNames_d[indices_d[j]]))
                        DSVariablePoolAddVariableWithName(Xd_a, name);
                DSVariablePoolAddVariableWithName(Xd, name);
        }
        DSSecureFree(variableNames_i);
        DSSecureFree(variableNames_d);
        DSSecureFree(name);
        numberExtra_i = newXi*(numberOfCases-1);
        numberExtra_d = newXd*(numberOfCases-1);
        rows = 0;
        columns_i = DSMatrixColumns(DSCaseU(cases[0]))+numberExtra_i;
        columns_d = DSMatrixColumns(DSCaseCd(cases[0]))+numberExtra_d;
        for (i = 0; i < numberOfCases; i++) {
                rows += DSMatrixRows(DSCaseZeta(cases[i]));
        }
        Cd = DSMatrixCalloc(rows, columns_d);
        Ci = DSMatrixCalloc(rows, columns_i);
        delta = DSMatrixCalloc(rows, 1);
        U = DSMatrixCalloc(rows, columns_i);
        Zeta = DSMatrixCalloc(rows, 1);
        currentRow = 0;
        for (i = 0; i < numberOfCases; i++) {
                tempCd = DSCaseCd(cases[i]);
                tempU = DSCaseU(cases[i]);
                tempZeta = DSCaseZeta(cases[i]);
                for (j = 0; j < DSMatrixRows(tempZeta); j++) {
                        DSMatrixSetDoubleValue(Zeta, currentRow, 0, DSMatrixDoubleValue(tempZeta, j, 0));
                        DSMatrixSetDoubleValue(delta, currentRow, 0, DSMatrixDoubleValue(DSCaseDelta(cases[i]), j, 0));
                        for (k = 0; k < DSMatrixColumns(tempU); k++) {
                                DSMatrixSetDoubleValue(Ci, currentRow, k, DSMatrixDoubleValue(DSCaseCi(cases[i]), j, k));
                                DSMatrixSetDoubleValue(U, currentRow, k, DSMatrixDoubleValue(tempU, j, k));
                        }
                        for (k = 0; k < DSMatrixColumns(tempCd); k++) {
                                DSMatrixSetDoubleValue(Cd, currentRow, k, DSMatrixDoubleValue(DSCaseCd(cases[i]), j, k));
                        }
                        if (i > 0) {
                                for (k = 0; k < newXi; k++) {
                                        DSMatrixSetDoubleValue(U,
                                                               currentRow,
                                                               DSMatrixColumns(tempU)+newXi*(i-1)+k,
                                                               DSMatrixDoubleValue(U, currentRow, indices_i[k]));
                                        DSMatrixSetDoubleValue(U, currentRow, indices_i[k], 0.0f);
                                        DSMatrixSetDoubleValue(Ci,
                                                               currentRow,
                                                               DSMatrixColumns(DSCaseCi(cases[i]))+newXi*(i-1)+k,
                                                               DSMatrixDoubleValue(Ci, currentRow, indices_i[k]));
                                        DSMatrixSetDoubleValue(Ci, currentRow, indices_i[k], 0.0f);
                                }
                                for (k = 0; k < newXd; k++) {
                                        DSMatrixSetDoubleValue(Cd,
                                                               currentRow,
                                                               DSMatrixColumns(tempCd)+newXd*(i-1)+k,
                                                               DSMatrixDoubleValue(Cd, currentRow, indices_d[k]));
                                        DSMatrixSetDoubleValue(Cd, currentRow, indices_d[k], 0.0f);
                                }
                        }
                        currentRow++;
                }
        }
        caseIntersection = DSSecureCalloc(1, sizeof(DSCase));
        caseIntersection->freeVariables = true;
        caseIntersection->Xd_a = Xd_a;
        caseIntersection->Xd = Xd;
        caseIntersection->Xi = Xi;
        DSCaseU(caseIntersection) = U;
        DSCaseZeta(caseIntersection) = Zeta;
        DSCaseCi(caseIntersection) = Ci;
        DSCaseCd(caseIntersection) = Cd;
        DSCaseDelta(caseIntersection) = delta;
        caseIntersection->ssys = NULL;
        U = NULL;
        Zeta = NULL;
        Ci = NULL;
        Cd = NULL;
        delta = NULL;
        DSSecureFree(indices_i);
        DSSecureFree(indices_d);
        DSCaseRemoveRedundantBoundaries(caseIntersection);
bail:
        if (Cd != NULL)
                DSMatrixFree(Cd);
        if (Ci != NULL)
                DSMatrixFree(Ci);
        if (delta != NULL)
                DSMatrixFree(delta);
        if (U != NULL)
                DSMatrixFree(U);
        if (Zeta != NULL)
                DSMatrixFree(Zeta);
        return caseIntersection;
}

extern const bool DSCaseIntersectionIsValid(const DSUInteger numberOfCases, const DSCase **cases)
{
        bool isValid = false;
        DSPseudoCase * caseIntersection = NULL;
        caseIntersection = DSPseudoCaseFromIntersectionOfCases(numberOfCases, cases);
        if (caseIntersection == NULL)
                goto bail;
        isValid = DSCaseIsValid(caseIntersection, true);
        DSSecureFree(caseIntersection);
bail:
        return isValid;
}

extern const bool DSCaseIntersectionIsValidAtSlice(const DSUInteger numberOfCases, const DSCase **cases,  const DSVariablePool * lowerBounds, const DSVariablePool *upperBounds)
{
        bool isValid = false;
        DSPseudoCase *caseIntersection = NULL;
        caseIntersection = DSPseudoCaseFromIntersectionOfCases(numberOfCases, cases);
        if (caseIntersection == NULL)
                goto bail;
        isValid = DSCaseIsValidAtSlice(caseIntersection, lowerBounds, upperBounds, true);
        DSSecureFree(caseIntersection);
bail:
        return isValid;
}

extern const bool DSCaseIntersectionExceptSliceIsValid(const DSUInteger numberOfCases, const DSCase **cases, const DSUInteger numberOfExceptions, const char ** exceptionVarNames)
{
        bool isValid = false;
        DSPseudoCase *caseIntersection = NULL;
        caseIntersection = DSPseudoCaseFromIntersectionOfCasesExcludingSlice(numberOfCases, cases, numberOfExceptions, exceptionVarNames);
        if (caseIntersection == NULL)
                goto bail;
        isValid = DSCaseIsValid(caseIntersection, true);
        DSSecureFree(caseIntersection);
bail:
        return isValid;
}

extern const bool DSCaseIntersectionExceptSliceIsValidAtSlice(const DSUInteger numberOfCases, const DSCase **cases, const DSUInteger numberOfExceptions, const char ** exceptionVarNames, const DSVariablePool * lowerBounds, const DSVariablePool * upperBounds)
{
        bool isValid = false;
        DSPseudoCase *caseIntersection = NULL;
        caseIntersection = DSPseudoCaseFromIntersectionOfCasesExcludingSlice(numberOfCases, cases, numberOfExceptions, exceptionVarNames);
        if (caseIntersection == NULL)
                goto bail;
        isValid = DSCaseIsValidAtSlice(caseIntersection, lowerBounds, upperBounds, true);
        DSSecureFree(caseIntersection);
bail:
        return isValid;
}

extern DSVariablePool * DSCaseIntersectionExceptSliceValidParameterSet(const DSUInteger numberOfCases, const DSCase **cases, const DSUInteger numberOfExceptions, const char ** exceptionVarNames)
{
        DSPseudoCase *caseIntersection = NULL;
        DSVariablePool * variablePool = NULL;
        caseIntersection = DSPseudoCaseFromIntersectionOfCasesExcludingSlice(numberOfCases, cases, numberOfExceptions, exceptionVarNames);
        if (caseIntersection == NULL)
                goto bail;
        variablePool = DSCaseValidParameterSet(caseIntersection);
        DSSecureFree(caseIntersection);
bail:
        return variablePool;
}

extern DSVariablePool * DSCaseIntersectionExceptSliceValidParameterSetByOptimizingFunction(const DSUInteger numberOfCases, const DSCase **cases, const DSUInteger numberOfExceptions, const char ** exceptionVarNames, const char * function, bool minimize)
{
        DSPseudoCase *caseIntersection = NULL;
        DSVariablePool * variablePool = NULL;
        caseIntersection = DSPseudoCaseFromIntersectionOfCasesExcludingSlice(numberOfCases, cases, numberOfExceptions, exceptionVarNames);
        if (caseIntersection == NULL)
                goto bail;
        variablePool = DSCaseValidParameterSetByOptimizingFunction(caseIntersection, function, minimize);
        DSSecureFree(caseIntersection);
bail:
        return variablePool;
}

extern DSVariablePool * DSCaseIntersectionExceptSliceValidParameterSetWithConstraints(const DSUInteger numberOfCases, const DSCase **cases, const DSUInteger numberOfExceptions, const char ** exceptionVarNames, const char ** constraints, DSUInteger numberOfConstraints)
{
        DSPseudoCase *caseIntersection = NULL;
        DSVariablePool * variablePool = NULL;
        caseIntersection = DSPseudoCaseFromIntersectionOfCasesExcludingSlice(numberOfCases, cases, numberOfExceptions, exceptionVarNames);
        DSCaseAddConstraints(caseIntersection, constraints, numberOfConstraints);
        if (caseIntersection == NULL)
                goto bail;
        variablePool = DSCaseValidParameterSet(caseIntersection);
        DSSecureFree(caseIntersection);
bail:
        return variablePool;
}

extern DSVariablePool * DSCaseIntersectionExceptSliceValidParameterSetWithConstraintsByOptimizingFunction(const DSUInteger numberOfCases, const DSCase **cases, const DSUInteger numberOfExceptions, const char ** exceptionVarNames, const char ** constraints, DSUInteger numberOfConstraints, const char * function, bool minimize)
{
        DSPseudoCase *caseIntersection = NULL;
        DSVariablePool * variablePool = NULL;
        caseIntersection = DSPseudoCaseFromIntersectionOfCasesExcludingSlice(numberOfCases, cases, numberOfExceptions, exceptionVarNames);
        DSCaseAddConstraints(caseIntersection, constraints, numberOfConstraints);
        if (caseIntersection == NULL)
                goto bail;
        variablePool = DSCaseValidParameterSetByOptimizingFunction(caseIntersection, function, minimize);
        DSSecureFree(caseIntersection);
bail:
        return variablePool;
}

extern DSVariablePool * DSCaseIntersectionExceptSliceValidParameterSetAtSlice(const DSUInteger numberOfCases, const DSCase **cases, const DSUInteger numberOfExceptions, const char ** exceptionVarNames, const DSVariablePool * lowerBounds, const DSVariablePool * upperBounds)
{
        DSPseudoCase *caseIntersection = NULL;
        DSVariablePool * variablePool = NULL;
        caseIntersection = DSPseudoCaseFromIntersectionOfCasesExcludingSlice(numberOfCases, cases, numberOfExceptions, exceptionVarNames);
        if (caseIntersection == NULL)
                goto bail;
        variablePool = DSCaseValidParameterSetAtSlice(caseIntersection, lowerBounds, upperBounds);
        DSSecureFree(caseIntersection);
bail:
        return variablePool;
}

extern DSVariablePool * DSCaseIntersectionExceptSliceValidParameterSetAtSliceByOptimizingFunction(const DSUInteger numberOfCases, const DSCase **cases, const DSUInteger numberOfExceptions, const char ** exceptionVarNames, const DSVariablePool * lowerBounds, const DSVariablePool * upperBounds, const char * function, bool minimize)
{
        DSPseudoCase *caseIntersection = NULL;
        DSVariablePool * variablePool = NULL;
        caseIntersection = DSPseudoCaseFromIntersectionOfCasesExcludingSlice(numberOfCases, cases, numberOfExceptions, exceptionVarNames);
        if (caseIntersection == NULL)
                goto bail;
        variablePool = DSCaseValidParameterSetAtSliceByOptimizingFunction(caseIntersection, lowerBounds, upperBounds, function, minimize);
        DSSecureFree(caseIntersection);
bail:
        return variablePool;
}

extern DSVertices * DSCaseIntersectionVerticesForSlice(const DSUInteger numberOfCases, const DSCase **cases, const DSVariablePool * lowerBounds, const DSVariablePool *upperBounds, const DSUInteger numberOfVariables, const char ** variables)
{
        DSVertices * vertices = NULL;
        DSPseudoCase *caseIntersection = NULL;
        if (numberOfCases == 0) {
                DSError(M_DS_WRONG ": Number of cases must be at least one", A_DS_ERROR);
                goto bail;
        }
        if (cases == NULL) {
                DSError(M_DS_NULL ": Array of cases is NULL", A_DS_ERROR);
                goto bail;
        }
        if (lowerBounds == NULL && upperBounds == NULL) {
                DSError(M_DS_VAR_NULL ": Variable pool with variables to fix is NULL", A_DS_ERROR);
                goto bail;
        }
        if (DSVariablePoolNumberOfVariables(lowerBounds) != DSVariablePoolNumberOfVariables(upperBounds)) {
                DSError(M_DS_WRONG ": Number of variables to bound must match", A_DS_ERROR);
                goto bail;
        }
        caseIntersection = DSPseudoCaseFromIntersectionOfCases(numberOfCases, cases);
        if (cases == NULL)
                goto bail;
        if (numberOfVariables == 1) {
                vertices = DSCaseVerticesFor1DSlice(caseIntersection, lowerBounds, upperBounds, variables[0]);
        } else {
                vertices = DSCaseVerticesForSlice(caseIntersection, lowerBounds, upperBounds, numberOfVariables, variables);
        }
        DSSecureFree(caseIntersection);
bail:
        return vertices;
}

extern DSMatrixArray * DSCaseIntersectionFacesFor3DSliceAndConnectivity(const DSUInteger numberOfCases, const DSCase **cases, const DSVariablePool * lowerBounds, const DSVariablePool *upperBounds, const char * xVariable, const char *yVariable, const char *zVariable)
{
        DSMatrixArray * faces = NULL;
        DSPseudoCase *caseIntersection = NULL;
        if (numberOfCases == 0) {
                DSError(M_DS_WRONG ": Number of cases must be at least one", A_DS_ERROR);
                goto bail;
        }
        if (cases == NULL) {
                DSError(M_DS_NULL ": Array of cases is NULL", A_DS_ERROR);
                goto bail;
        }
        if (lowerBounds == NULL && upperBounds == NULL) {
                DSError(M_DS_VAR_NULL ": Variable pool with variables to fix is NULL", A_DS_ERROR);
                goto bail;
        }
        if (DSVariablePoolNumberOfVariables(lowerBounds) != DSVariablePoolNumberOfVariables(upperBounds)) {
                DSError(M_DS_WRONG ": Number of variables to bound must match", A_DS_ERROR);
                goto bail;
        }
        caseIntersection = DSPseudoCaseFromIntersectionOfCases(numberOfCases, cases);
        if (cases == NULL)
                goto bail;
        faces = DSCaseFacesFor3DSliceAndConnectivity(caseIntersection, lowerBounds, upperBounds, xVariable, yVariable, zVariable);
        DSSecureFree(caseIntersection);
bail:
        return faces;
}

extern const bool DSUnstableCaseConditionsAreValid(DSUnstableCase *uCase, const DSUInteger *bSignature)
{
    
            bool isValid = false;
            glp_prob *linearProblem = NULL;
            DSMatrix * U1 = NULL, *U2 = NULL, *U = NULL, * Zeta1 = NULL, *Zeta = NULL;
            DSCase *aCase = uCase->originalCase;
    
            if (uCase == NULL || aCase == NULL ) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
            }
    
            if (bSignature == NULL)
                goto bail;
    
            // populate matrices Cd_unstable, Ci_unstable and delta_unstable
            dsUnstableCaseGetAdditionalConstraintMatrices(uCase, bSignature);
            
             if (uCase->Cd_unstable == NULL || uCase->Ci_unstable == NULL || uCase->delta_unstable == NULL ) {
                 goto bail;
             }
    
    
            if (DSCaseCd(aCase) != NULL && DSCaseCi(aCase) != NULL ){
                U1 = DSMatrixAppendMatrices(DSCaseCd(aCase), DSCaseCi(aCase), true);
                U2 = DSMatrixAppendMatrices(uCase->Cd_unstable, uCase->Ci_unstable, true);
                U = DSMatrixAppendMatrices(U1, U2, false);
                if (U1 != NULL)
                    DSMatrixFree(U1);
                if (U2 != NULL)
                    DSMatrixFree(U2);
                Zeta1 = DSCaseDelta(aCase);
                Zeta = DSMatrixAppendMatrices(Zeta1, uCase->delta_unstable, false);
            } else {
                 U = DSMatrixAppendMatrices(uCase->Cd_unstable, uCase->Ci_unstable, true);
                 Zeta = DSMatrixCopy(uCase->delta_unstable);
            }
            if (U != NULL && Zeta != NULL)
                    linearProblem = dsUnstableCaseLinearProblemForCaseValidity(U, Zeta, uCase, bSignature);
            if (U !=NULL)
                DSMatrixFree(U);
            if (Zeta != NULL)
                DSMatrixFree(Zeta);
    
            if (linearProblem != NULL) {
                glp_simplex(linearProblem, NULL);
                
                // -1e-14 original value // alt. value -1e-13
                if (glp_get_obj_val(linearProblem) <= -1e-13 && glp_get_prim_stat(linearProblem) == GLP_FEAS) {
                    isValid = true;
                    DSUnstableCaseDetermineBlowingBehavior(linearProblem, uCase, bSignature);
                }
                glp_delete_prob(linearProblem);
            }
        bail:
            return isValid;
    
}

extern glp_prob * dsUnstableCaseLinearProblemForMatrices(const DSMatrix *A_, const DSMatrix *B_,
                                                         DSUnstableCase *uCase, const DSUInteger * bSignature)
{
            glp_prob *linearProblem = NULL;
            int * ia = NULL, *ja = NULL;
            double *ar = NULL;
            DSUInteger i, numberOfXi, numberOfBoundaries;
            DSUInteger numberOfEqualities = uCase->Xd_e->numberOfVariables , numberOfKnifes = uCase->knifeEdge->numberOfVariables;
            DSUInteger numberOfBlowing = uCase->Xd_b->numberOfVariables;
            DSMatrixArray *lp_AB1 = NULL, *lp_AB2 = NULL;
            DSUInteger blowingIndex;
            char * name;
            const DSVariablePool *Xd = uCase->originalCase->Xd;
            DSMatrix *A, *B;
            DSUInteger numberOfOriginalBounds;
    
            if (uCase->originalCase->delta != NULL ){
                 numberOfOriginalBounds = DSMatrixRows(uCase->originalCase->delta) + numberOfBlowing;
            } else {
                 numberOfOriginalBounds = numberOfBlowing;
            }
    
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
            numberOfXi = DSMatrixColumns(A);
            numberOfBoundaries = DSMatrixRows(A);
    
        //    printf("-------reporting from dsUnstableCaseLinearProblemForMatrices \n");
        //    printf("The matrix A is: \n");
        //    DSMatrixPrint(A);
        //    printf("The matrix B is: \n");
        //    DSMatrixPrint(B);
    
            ia = DSMatrixRowsForGLPK(A);
            ja = DSMatrixColumnsForGLPK(A);
            ar = DSMatrixDataForGLPK(A);
    
            glp_add_rows(linearProblem, numberOfBoundaries);
            glp_add_cols(linearProblem, numberOfXi);
            glp_set_obj_dir(linearProblem, GLP_MIN);
            glp_load_matrix(linearProblem, numberOfBoundaries*numberOfXi, ia, ja, ar);
    
            // set bounds for original boundaries and additional constraints coming from blow up/down assumption.
            if (numberOfOriginalBounds != 0){
                for (i = 0; i < numberOfOriginalBounds; i++) {
                    glp_set_row_bnds(linearProblem, i+1, GLP_UP, 0.0,
                                     DSMatrixDoubleValue(B, i, 0));
                }
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
    
            // set bounds for all variables:
            for (i = 0; i < numberOfXi; i++)
                glp_set_col_bnds(linearProblem, i+1, GLP_DB, -3.0, 3.0);

    
            // set bounds for dependet variables free. This is useful for cases where the solution of a dependent pool depends on a Xd_b
            DSUInteger nr_dependent = uCase->originalCase->Xd->numberOfVariables;
            for (i = 0; i < nr_dependent; i++)
                glp_set_col_bnds(linearProblem, i+1, GLP_FR, 0.0, 0.0);
    
    
            //set value for blowing variables. Loop over number of Xd_b
            for(i=0; i<numberOfBlowing; i++){
                name = DSVariableName(DSVariablePoolVariableAtIndex(uCase->Xd_b, i));
                blowingIndex = DSVariablePoolIndexOfVariableWithName(Xd, name) ;
                if (bSignature[i] == 1){
                    glp_set_col_bnds(linearProblem, blowingIndex+1, GLP_FX, 12.0, 12.0);
                }else{
                    glp_set_col_bnds(linearProblem, blowingIndex+1, GLP_FX, -12.0, -12.0);
                }
            }
    
            if (ia != NULL)
                DSSecureFree(ia);
            if (ja != NULL)
                DSSecureFree(ja);
            if (ar != NULL)
                DSSecureFree(ar);
            if (uCase->Xd_e->numberOfVariables != 0)
                DSMatrixArrayFree(lp_AB1);
            DSMatrixArrayFree(lp_AB2);
        bail:
            return linearProblem;
}

extern glp_prob * dsUnstableCaseLinearProblemForCaseValidity(const DSMatrix * U,
                                                             const DSMatrix *zeta,
                                                             DSUnstableCase *uCase,
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
            glp_set_col_bnds(linearProblem, glp_get_num_cols(linearProblem), GLP_LO, -1.0, 0.0);
            glp_set_obj_coef(linearProblem, glp_get_num_cols(linearProblem), 1.0);
    
            DSMatrixFree(coefficients);
            if (slacks != NULL)
                DSMatrixFree(slacks);
        bail:
            return linearProblem;
}
