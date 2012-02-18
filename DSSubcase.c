//
//  DSSubcase.c
//  DesignSpaceToolboxV2
//
//  Created by Jason Lomnitz on 9/20/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#include <stdio.h>
#include "DSErrors.h"
#include "DSMemoryManager.h"
#include "DSExpression.h"
#include "DSSubcase.h"
#include "DSMatrix.h"
#include "DSMatrixArray.h"
#include "DSSSystem.h"
#include "DSCase.h"
#include "DSGMASystem.h"
#include "DSDesignSpace.h"
#include "DSDesignSpaceStack.h"

extern DSMatrix * DSSubcaseProblematicEquations(const DSCase * aCase)
{
        DSMatrix *problematic = NULL;
        bool isUnderdetermined = false;
        DSMatrix *nullspace = NULL, *A = NULL;
        DSUInteger i, j;
        double firstValue = NAN, current;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (DSCaseHasSolution(aCase) == true) {
                goto bail;
        }
        A = DSSSystemA(aCase->ssys);
        nullspace = DSMatrixLeftNullspace(A);
        DSMatrixFree(A);
        if (nullspace == NULL)
                goto bail;
        
        isUnderdetermined = true;
        problematic = DSMatrixCalloc(DSMatrixRows(nullspace), DSMatrixColumns(nullspace));
        for (i = 0; i < DSMatrixColumns(nullspace); i++) {
                firstValue = NAN;
                for (j = 0; j < DSMatrixRows(nullspace); j++) {
                        current = DSMatrixDoubleValue(nullspace, j, i); 
                        if (fabs(current) < 1E-14)
                                continue;
                        DSMatrixSetDoubleValue(problematic, j, i, 1.0);
                        if (isnan(firstValue) == true) {
                                firstValue = current;
                        } else if (fabs(current - firstValue) >= 1E-14) {
                                isUnderdetermined = false;
                                break;
                        }
                }
                if (j != DSMatrixRows(nullspace))
                        break;
        }
        if (isUnderdetermined == false) {
                DSMatrixFree(problematic);
                problematic = NULL;
        }
        DSMatrixFree(nullspace);
bail:
        return problematic;
}

extern DSMatrixArray * DSSubcaseProblematicTerms(const DSCase *aCase, const DSMatrix * dependentEquations)
{
        DSMatrixArray *dependentTerms = NULL;
        DSMatrix *G, *H, *g, *h, *termMatrix, *nullspace, *coefficients;
        DSUInteger i, j,k, *dependent, numDependent;
        double value;
        double sign;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (DSCaseHasSolution(aCase) == true) {
                goto bail;
        }
        if (dependentEquations == NULL) {
                DSError(M_DS_MAT_NULL ": Array of dependent equations must not be null", A_DS_ERROR);
                goto bail;
        }
        dependentTerms = DSMatrixArrayAlloc();
        dependent = DSSecureCalloc(sizeof(DSUInteger), DSSSystemNumberOfEquations(aCase->ssys));
        G = DSSSystemG(aCase->ssys);
        H = DSSSystemH(aCase->ssys);
        for (i = 0; i < DSMatrixColumns(dependentEquations); i++) {
                numDependent = 0;
                sign = 1.0;
                for (j = 0; j < DSMatrixRows(dependentEquations); j++) {
                        if (DSMatrixDoubleValue(dependentEquations, j, i) == 1.0) {
                                dependent[numDependent++] = j;
                        }
                }
                g = DSMatrixSubMatrixIncludingRows(G, numDependent, dependent);
                h = DSMatrixSubMatrixIncludingRows(H, numDependent, dependent);
                termMatrix = DSMatrixAppendMatrices(g, h, false);
                DSMatrixFree(g);
                DSMatrixFree(h);
                nullspace = DSMatrixLeftNullspace(termMatrix);
                coefficients = DSMatrixCalloc(numDependent, DSMatrixColumns(nullspace));
                for (j = 0; j < DSMatrixRows(nullspace); j++) {
                        for (k = 0; k < DSMatrixColumns(nullspace); k++) {
                                value = DSMatrixDoubleValue(nullspace, j, k);
                                if (fabs(value) <= 1E-14) {
                                        DSMatrixSetDoubleValue(nullspace, j, k, 0.0);
                                        continue;
                                }
                                DSMatrixSetDoubleValue(nullspace, j, k, copysign(1.0, value));
                                if (j / numDependent == 0)
                                        DSMatrixSetDoubleValue(coefficients, j % numDependent, k, DSMatrixDoubleValue(DSSSystemAlpha(aCase->ssys), dependent[j % numDependent], 0));
                                else
                                        DSMatrixSetDoubleValue(coefficients, j % numDependent, k, -DSMatrixDoubleValue(DSSSystemBeta(aCase->ssys), dependent[j % numDependent], 0));
                                sign *= -1.0;
                        }
                }
                DSMatrixFree(termMatrix);
                DSMatrixFree(nullspace);
                DSMatrixArrayAddMatrix(dependentTerms, coefficients);
        }
        DSMatrixFree(G);
        DSMatrixFree(H);
        DSSecureFree(dependent);
bail:
        return dependentTerms;
}

extern DSMatrixArray * DSSubcaseCoefficientsOfInterest(const DSCase * aCase, const DSMatrixArray * problematicTerms)
{
        DSMatrixArray * coefficientArray = NULL;
        DSMatrix *problematic = NULL;
        DSUInteger i, j;
        double min, value;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (DSCaseHasSolution(aCase) == true) {
                goto bail;
        }
        if (problematicTerms == NULL) {
                DSError(M_DS_MAT_NULL ": Array of problematic terms is null", A_DS_ERROR);
                goto bail;
        }
        coefficientArray = DSMatrixArrayAlloc();
        for (i = 0; i < DSMatrixArrayNumberOfMatrices(problematicTerms); i++) {
                problematic = DSMatrixLeftNullspace(DSMatrixArrayMatrix(problematicTerms, i));
                if (problematic == NULL)
                        continue;
                DSMatrixRoundToSignificantFigures(problematic, 14);
                min = INFINITY;
                for (j = 0; j < DSMatrixRows(problematic); j++) {
                        value = DSMatrixDoubleValue(problematic, j, 0);
                        if (fabs(value) == 0)
                                continue;
                        min = ((fabs(value) <= fabs(min)) ? value : min);
                }
                for (j = 0; j < DSMatrixRows(problematic); j++) {
                        value = DSMatrixDoubleValue(problematic, j, 0);
                        if (value == 0)
                                continue;
                        DSMatrixSetDoubleValue(problematic, j, 0, value/min);
                }
                DSMatrixArrayAddMatrix(coefficientArray, problematic);
        }
bail:
        return coefficientArray;
}
/*
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
*/
//static DSDesignSpace * dsSubcaseCreateUniqueSystemSubcase(const DSCase *aCase, const DSGMASystem * modifiedGMA, const DSMatrix  * problematicEquations, const DSExpression ** augmentedEquations)
//{
//        DSDesignSpace * ds = NULL;
//        DSUInteger i, j;
//        DSUInteger * equationIndex = NULL;
//        char **equations;
//        DSExpression **caseEquations;
//        DSExpression *newEquations;
//        if (aCase == NULL) {
//                DSError(M_DS_CASE_NULL, A_DS_ERROR);
//                goto bail;
//        }
//        if (modifiedGMA == NULL) {
//                DSError(M_DS_GMA_NULL, A_DS_ERROR);
//                goto bail;
//        }
//        if (problematicEquations == NULL) {
//                DSError(M_DS_MAT_NULL, A_DS_ERROR);
//                goto bail;
//        }
//        if (augmentedEquations == NULL) {
//                DSError(M_DS_NULL ": Augmented equations not found", A_DS_ERROR);
//                goto bail;
//        }
//        equationIndex = DSSecureMalloc(sizeof(DSUInteger)*DSMatrixColumns(problematicEquations));
//        for (i = 0; i < DSMatrixColumns(problematicEquations); i++) {
//                equationIndex[i] = DSMatrixRows(problematicEquations);
//                for (j = 0; j < DSMatrixRows(problematicEquations); j++) {
//                        if (DSMatrixDoubleValue(problematicEquations, j, i) == 0.0)
//                                continue;
//                        equationIndex[i] = j;
//                        break;
//                }
//        }
//        caseEquations = DSCaseEquations(aCase);
//        equations = DSSecureCalloc(sizeof(char *), DSCaseNumberOfEquations(aCase));
//        for (i = 0; i < DSCaseNumberOfEquations(aCase); i++) {
//                equations[i] = DSExpressionAsString(caseEquations[i]);
//                for (j = 0; j < DSMatrixColumns(problematicEquations); j++) {
//                        if (i != equationIndex[j])
//                                continue;
//                        newEquations = DSExpressionCopy(augmentedEquations[j]);
//                        DSSecureFree(equations[i]);
//                        equations[i] = DSExpressionAsString(newEquations);
//                        DSExpressionFree(newEquations);
//                }
//                //                printf("%i: %s\n", i, equations[i]);
//        }
//        ds = DSDesignSpaceByParsingStringsWithXi(DSGMASystemXd(modifiedGMA), DSGMASystemXi(modifiedGMA), equations, DSCaseNumberOfEquations(aCase));
//        for (i = 0; i < DSCaseNumberOfEquations(aCase); i++) {
//                DSSecureFree(equations[i]);
//                DSExpressionFree(caseEquations[i]);
//        }
//        DSSecureFree(equations);
//        DSSecureFree(caseEquations);
//bail:
//        if (equationIndex != NULL)
//                DSSecureFree(equationIndex);
//        return ds;
//}

static DSDesignSpace * dsSubcaseCreateUniqueSystemSubcase(const DSCase *aCase, const DSGMASystem * modifiedGMA, const DSMatrix  * problematicEquations, const DSExpression ** augmentedEquations)
{
        DSDesignSpace * ds = NULL;
        DSUInteger i, j;
        DSUInteger * equationIndex = NULL;
        char **equations;
        DSExpression **caseEquations;
        DSExpression *newEquations;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (modifiedGMA == NULL) {
                DSError(M_DS_GMA_NULL, A_DS_ERROR);
                goto bail;
        }
        if (problematicEquations == NULL) {
                DSError(M_DS_MAT_NULL, A_DS_ERROR);
                goto bail;
        }
        if (augmentedEquations == NULL) {
                DSError(M_DS_NULL ": Augmented equations not found", A_DS_ERROR);
                goto bail;
        }
        equationIndex = DSSecureMalloc(sizeof(DSUInteger)*DSMatrixColumns(problematicEquations));
        for (i = 0; i < DSMatrixColumns(problematicEquations); i++) {
                equationIndex[i] = DSMatrixRows(problematicEquations);
                for (j = 0; j < DSMatrixRows(problematicEquations); j++) {
                        if (DSMatrixDoubleValue(problematicEquations, j, i) == 0.0)
                                continue;
                        equationIndex[i] = j;
                        break;
                }
        }
        caseEquations = DSCaseEquations(aCase);
        equations = DSSecureCalloc(sizeof(char *), DSCaseNumberOfEquations(aCase));
        for (i = 0; i < DSCaseNumberOfEquations(aCase); i++) {
                equations[i] = DSExpressionAsString(caseEquations[i]);
                for (j = 0; j < DSMatrixColumns(problematicEquations); j++) {
                        if (i != equationIndex[j])
                                continue;
                        newEquations = DSExpressionCopy(augmentedEquations[j]);
                        DSSecureFree(equations[i]);
                        equations[i] = DSExpressionAsString(newEquations);
                        DSExpressionFree(newEquations);
                }
        }
        ds = DSDesignSpaceByParsingStringsWithXi(DSGMASystemXd(modifiedGMA), DSGMASystemXi(modifiedGMA), equations, DSCaseNumberOfEquations(aCase));
        for (i = 0; i < DSCaseNumberOfEquations(aCase); i++) {
                DSSecureFree(equations[i]);
                DSExpressionFree(caseEquations[i]);
        }
        DSSecureFree(equations);
        DSSecureFree(caseEquations);
bail:
        if (equationIndex != NULL)
                DSSecureFree(equationIndex);
        return ds;
}

extern void DSSubcaseDesignSpaceForUnderdeterminedCase(const DSCase * aCase, const DSDesignSpace * original)
{
        DSDesignSpace *subcases = NULL;
        DSGMASystem * temp = NULL;
        DSMatrix * problematicEquations = NULL;
        DSMatrixArray * problematicTerms = NULL;
        DSMatrixArray * coefficientArray = NULL;
        DSUInteger i, j, k, l;
        DSExpression **augmentedEquations;
        double value;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (original == NULL) {
                DSError(M_DS_DESIGN_SPACE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (DSCaseNumberOfEquations(aCase) != DSDesignSpaceNumberOfEquations(original)) {
                DSError(M_DS_WRONG ": Number of equation in design space must match number of equations in case", A_DS_ERROR);
                goto bail;
        }
        problematicEquations = DSSubcaseProblematicEquations(aCase);
        if (problematicEquations == NULL)
                goto bail;
        problematicTerms = DSSubcaseProblematicTerms(aCase, problematicEquations);
        if (problematicTerms == NULL)
                goto bail;
        coefficientArray = DSSubcaseCoefficientsOfInterest(aCase, problematicTerms);
        if (coefficientArray == NULL)
                goto bail;
        if (DSMatrixArrayNumberOfMatrices(problematicTerms) != DSMatrixArrayNumberOfMatrices(coefficientArray))
                goto bail;
        temp = DSGMASystemCopy(original->gma);
        augmentedEquations = DSSecureCalloc(sizeof(DSExpression *), DSMatrixColumns(problematicEquations));
        for (i = 0; i < DSMatrixColumns(problematicEquations); i++) {
                l = 0;
                for (j = 0; j < DSMatrixRows(problematicEquations); j++) {
                        if (DSMatrixDoubleValue(problematicEquations, j, i) == 0)
                                continue;
                        for (k = 0; k < DSMatrixColumns(DSGMASystemAlpha(temp)); k++) {
                                value = DSMatrixArrayDoubleWithIndices(coefficientArray, i, l, 0);
                                if (k+1 == aCase->signature[2*j])
                                        value = 0.0;
                                DSMatrixSetDoubleValue((DSMatrix *)DSGMASystemAlpha(temp), j, k, DSMatrixDoubleValue(DSGMASystemAlpha(temp), j, k)*value);
                        }
                        for (k = 0; k < DSMatrixColumns(DSGMASystemBeta(temp)); k++) {
                                value = DSMatrixArrayDoubleWithIndices(coefficientArray, i, l, 0);
                                if (k+1 == aCase->signature[2*j+1])
                                        value = 0.0;
                                DSMatrixSetDoubleValue((DSMatrix *)DSGMASystemBeta(temp), j, k, DSMatrixDoubleValue(DSGMASystemBeta(temp), j, k)*value);
                        }
                        l++;
                        augmentedEquations[i] = DSExpressionAddExpressions(augmentedEquations[i], DSGMASystemPositiveTermsForEquations(temp, j));
                        augmentedEquations[i] = DSExpressionAddExpressions(augmentedEquations[i], DSGMASystemNegativeTermsForEquations(temp, j));
                }
        }
        
        subcases = dsSubcaseCreateUniqueSystemSubcase(aCase, temp, problematicEquations, (const DSExpression **)augmentedEquations);
        if (subcases != NULL) {
                DSDesignSpaceAddConditions(subcases, aCase->Cd, aCase->Ci, aCase->delta);
                DSDesignSpaceStackPush(original->subcases, subcases);
        }
        for (i = 0; i < DSMatrixColumns(problematicEquations); i++)
                DSExpressionFree(augmentedEquations[i]);
        DSSecureFree(augmentedEquations);
        DSGMASystemFree(temp);
bail:
        if (problematicEquations != NULL)
                DSMatrixFree(problematicEquations);
        if (problematicTerms != NULL)
                DSMatrixArrayFree(problematicTerms);
        if (coefficientArray != NULL)
                DSMatrixArrayFree(coefficientArray);
        return;
}
