//
//  DSCyclicalCaseConstruction.c
//  DesignSpaceToolboxV2
//
//  Created by Jason Lomnitz on 8/7/13.
//
//

#include <stdio.h>
#include <string.h>
#include "DSCyclicalCase.h"
#include "DSMatrixArray.h"

extern void DSCyclicalCaseDesignSpaceCalculateSubCyclicalCases(DSDesignSpace *ds, DSDesignSpace * modifierDesignSpace, const DSUInteger * modifierSignature);
static DSUInteger dsCyclicalCaseNumberOfAugmentedSystems(const DSDesignSpace * original,
                                                         const DSMatrix * problematicEquations);
#if defined (__APPLE__) && defined (__MACH__)
#pragma mark Cyclical calculation functions -
#endif


extern DSMatrix * dsSubcaseProblematicEquations(const DSCase * aCase)
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

extern DSMatrixArray * dsSubcaseProblematicTerms(const DSCase *aCase, const DSMatrix * dependentEquations)
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
                if (numDependent == 0)
                        break;
                for (j = 0; j < DSMatrixRows(nullspace); j++) {
                        for (k = 0; k < DSMatrixColumns(nullspace); k++) {
                                value = DSMatrixDoubleValue(nullspace, j, k);
                                if (fabs(value) <= 1E-14) {
                                        DSMatrixSetDoubleValue(nullspace, j, k, 0.0);
                                        continue;
                                }
                                DSMatrixSetDoubleValue(nullspace, j, k, copysign(1.0, value));
                                if ((j / numDependent) == 0)
                                        DSMatrixSetDoubleValue(coefficients, j % numDependent, k, DSMatrixDoubleValue(DSSSystemAlpha(aCase->ssys), dependent[j % numDependent], 0));
                                else
                                        DSMatrixSetDoubleValue(coefficients, j % numDependent,
                                                               k,
                                                               DSMatrixDoubleValue(coefficients, j%numDependent, k)-DSMatrixDoubleValue(DSSSystemBeta(aCase->ssys), dependent[j % numDependent], 0));
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

extern DSMatrixArray * dsSubcaseCoefficientsOfInterest(const DSCase * aCase, const DSMatrixArray * problematicTerms)
{
        DSMatrixArray * coefficientArray = NULL;
        DSMatrix *problematic = NULL, *coefficients;
        DSUInteger i, j, k;
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
                for (k = 0; k < DSMatrixColumns(problematic); k++) {
                        min = INFINITY;
                        for (j = 0; j < DSMatrixRows(problematic); j++) {
                                value = DSMatrixDoubleValue(problematic, j, k);
                                if (fabs(value) == 0)
                                        continue;
                                min = ((fabs(value) <= fabs(min)) ? value : min);
                        }
                        for (j = 0; j < DSMatrixRows(problematic); j++) {
                                value = DSMatrixDoubleValue(problematic, j, k);
                                if (value == 0)
                                        continue;
                                DSMatrixSetDoubleValue(problematic, j, k, value/min);
                        }
                }
                coefficients = DSMatrixCalloc(DSMatrixRows(problematic), 1);
                for (j = 0; j < DSMatrixRows(problematic); j++) {
                        for (k = 0; k < DSMatrixColumns(problematic); k++) {
                                value = DSMatrixDoubleValue(coefficients, j, 0);
                                value += DSMatrixDoubleValue(problematic, j, k);
                                DSMatrixSetDoubleValue(coefficients, j, 0, value);
                        }
                }
                DSMatrixFree(problematic);
                DSMatrixArrayAddMatrix(coefficientArray, coefficients);
        }
bail:
        return coefficientArray;
}

static DSDesignSpace * dsSubcaseCreateUniqueSystemSubcase(const DSCase *aCase, const DSGMASystem * modifiedGMA, const DSMatrix  * problematicEquations, const DSExpression ** augmentedEquations)
{
        DSDesignSpace * ds = NULL;
        DSUInteger i, j;
        DSUInteger * equationIndex = NULL;
        char **equations, *temp_rhs, *temp_lhs;
        DSExpression **caseEquations;
        DSExpression  *eqLHS;
        DSVariablePool * Xda = NULL;
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
        caseEquations = DSCaseEquations(aCase);
        equations = DSSecureCalloc(sizeof(char *), DSCaseNumberOfEquations(aCase));
        Xda = DSVariablePoolAlloc();
        for (i = 0; i < DSVariablePoolNumberOfVariables(DSGMASystemXd_a(modifiedGMA)); i++)
                DSVariablePoolAddVariableWithName(Xda, DSVariableName(DSVariablePoolVariableAtIndex(DSGMASystemXd_a(modifiedGMA), i)));
        equationIndex = DSSecureMalloc(sizeof(DSUInteger)*DSMatrixColumns(problematicEquations));
        for (i = 0; i < DSCaseNumberOfEquations(aCase); i++)
                equations[i] = DSExpressionAsString(caseEquations[i]);
        for (i = 0; i < DSMatrixColumns(problematicEquations); i++) {
                for (j = 0; j < DSMatrixRows(problematicEquations); j++) {
                        if (DSMatrixDoubleValue(problematicEquations, j, i) == 0.0)
                                continue;
                        equationIndex[i] = j;
                        break;
                }
        }
        for (i = 0; i < DSCaseNumberOfEquations(aCase); i++) {
                for (j = 0; j < DSMatrixColumns(problematicEquations); j++) {
                        if (i != equationIndex[j])
                                continue;
                        DSSecureFree(equations[i]);
                        eqLHS = DSExpressionEquationLHSExpression(caseEquations[i]);
                        temp_lhs = DSExpressionAsString(eqLHS);
                        temp_rhs = DSExpressionAsString(augmentedEquations[j]);
                        DSExpressionFree(eqLHS);
                        equations[i] = DSSecureCalloc(sizeof(char),strlen(temp_rhs) + strlen(temp_lhs)+4);
                        equations[i] = strcpy(equations[i], temp_lhs);
                        equations[i] = strcat(equations[i], " = ");
                        equations[i] = strcat(equations[i], temp_rhs);
                        DSSecureFree(temp_rhs);
                        DSSecureFree(temp_lhs);
                }
        }
        ds = DSDesignSpaceByParsingStringsWithXi(equations, Xda, DSGMASystemXi(modifiedGMA), DSCaseNumberOfEquations(aCase));
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

#if defined (__APPLE__) && defined (__MACH__)
#pragma mark - Multiple augmented systems to obtain appropriate dynamics -
#endif

static void dsAddConstraintsForSubdominantDecays(DSDesignSpace * subcase, const DSCase * aCase, const DSDesignSpace * original, DSMatrix * problematicEquations, DSMatrixArray * coefficientArray, DSUInteger * subdominantDecays, DSUInteger * subdominantDecayTerms)
{
        DSGMASystem * gma = NULL;
        DSUInteger i, j, k, m, l, index = 0;
        DSUInteger numberOfConditions = 0, numberOfXd, numberOfXi;
        DSMatrix * Cd, *Ci, *delta;
        double subCoefficient = 0, coefficient, value;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (original == NULL) {
                DSError(M_DS_DESIGN_SPACE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (problematicEquations == NULL) {
                DSError(M_DS_MAT_NULL, A_DS_ERROR);
                goto bail;
        }
        if (subdominantDecays == NULL) {
                DSError(M_DS_NULL ": Array indicating the subdominant decay equation is NULL", A_DS_ERROR);
                goto bail;
        }
        if (subdominantDecayTerms == NULL) {
                DSError(M_DS_NULL ": Array indicating the subdominant decay terms is NULL", A_DS_ERROR);
                goto bail;
        }
        gma = DSGMASystemCopy(DSDesignSpaceGMASystem(original));
        numberOfXd = DSVariablePoolNumberOfVariables(DSGMASystemXd(gma));
        numberOfXi = DSVariablePoolNumberOfVariables(DSGMASystemXi(gma));
        for (i = 0; i < DSMatrixColumns(problematicEquations); i++) {
                for (j = 0; j < DSMatrixRows(problematicEquations); j++) {
                        if ((DSUInteger)DSMatrixDoubleValue(problematicEquations, j, i) == 0)
                                continue;
                        if (subdominantDecays[i] == j) {
                                numberOfConditions += (DSDesignSpaceSignature(original)[j*2+1]-2);
                                continue;
                        }
                        numberOfConditions += (DSDesignSpaceSignature(original)[j*2+1]-1);
                }
        }
        if (numberOfConditions == 0) {
                goto bail;
        }
        Cd = DSMatrixAlloc(numberOfConditions, DSVariablePoolNumberOfVariables(DSGMASystemXd(gma)));
        Ci = DSMatrixAlloc(numberOfConditions, DSVariablePoolNumberOfVariables(DSGMASystemXi(gma)));
        delta =DSMatrixAlloc(numberOfConditions, 1);
        index = 0;

        for (i = 0; i < DSMatrixColumns(problematicEquations); i++) {
                k = 0;
                for (j = 0; j < DSMatrixRows(problematicEquations); j++) {
                        if (DSMatrixDoubleValue(problematicEquations, j, i) == 0)
                                continue;
                        if (j == subdominantDecays[i]) {
                                subCoefficient = DSMatrixArrayDoubleWithIndices(coefficientArray, i, k, 0);
                        }
                        k++;
                }
                l = 0;
                for (j = 0; j < DSMatrixRows(problematicEquations); j++) {
                        if ((DSUInteger)DSMatrixDoubleValue(problematicEquations, j, i) == 0)
                                continue;
                        coefficient = DSMatrixArrayDoubleWithIndices(coefficientArray, i, l, 0);
                        for (k = 0; k < DSDesignSpaceSignature(original)[j*2+1];  k++) {
                                if (k+1 == DSCaseSignature(aCase)[j*2+1])
                                        continue;
                                if (k == subdominantDecayTerms[i] && j == subdominantDecays[i])
                                        continue;
                                value = log10(subCoefficient/coefficient*DSMatrixDoubleValue(DSGMASystemBeta(gma), subdominantDecays[i], subdominantDecayTerms[i])
                                              /DSMatrixDoubleValue(DSGMASystemBeta(gma), j, k));
                                DSMatrixSetDoubleValue(delta, index, 0, value);
                                for (m = 0; m < numberOfXd; m++) {
                                        value = DSMatrixArrayDoubleWithIndices(DSGMASystemHd(gma), subdominantDecays[i], subdominantDecayTerms[i], m);
                                        value -= DSMatrixArrayDoubleWithIndices(DSGMASystemHd(gma), j, k, m);
                                        DSMatrixSetDoubleValue(Cd, index, m, value);
                                }
                                for (m = 0; m < numberOfXi; m++) {
                                        value = DSMatrixArrayDoubleWithIndices(DSGMASystemHi(gma), subdominantDecays[i], subdominantDecayTerms[i], m);
                                        value -= DSMatrixArrayDoubleWithIndices(DSGMASystemHi(gma), j, k, m);
                                        DSMatrixSetDoubleValue(Ci, index, m, value);
                                }
                                index++;
                        }
                        l++;
                }
        }
        DSDesignSpaceAddConditions(subcase, Cd, Ci, delta);
        DSMatrixFree(Cd);
        DSMatrixFree(Ci);
        DSMatrixFree(delta);
bail:
        if (gma != NULL)
                DSGMASystemFree(gma);
        return;
}

static DSDesignSpace * dsCyclicalCaseCreateUniqueAugmentedSystem(const DSCase *aCase, const DSGMASystem * modifiedGMA, const DSMatrix  * problematicEquations, const DSExpression ** augmentedEquations, const DSUInteger * subdominantDecays)
{
        DSDesignSpace * ds = NULL;
        DSUInteger i, j;
        char **equations, *temp_rhs, *temp_lhs;
        DSExpression **caseEquations;
        DSExpression  *eqLHS;
        DSVariablePool * Xda = NULL;
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
        caseEquations = DSCaseEquations(aCase);
        equations = DSSecureCalloc(sizeof(char *), DSCaseNumberOfEquations(aCase));
        Xda = DSVariablePoolAlloc();
        for (i = 0; i < DSVariablePoolNumberOfVariables(DSGMASystemXd_a(modifiedGMA)); i++)
                DSVariablePoolAddVariableWithName(Xda, DSVariableName(DSVariablePoolVariableAtIndex(DSGMASystemXd_a(modifiedGMA), i)));
        for (i = 0; i < DSCaseNumberOfEquations(aCase); i++) {
                equations[i] = DSExpressionAsString(caseEquations[i]);
        }
        for (j = 0; j < DSMatrixColumns(problematicEquations); j++) {
                i = subdominantDecays[j];
                DSSecureFree(equations[i]);
                eqLHS = DSExpressionEquationLHSExpression(caseEquations[i]);
                temp_lhs = DSExpressionAsString(eqLHS);
                temp_rhs = DSExpressionAsString(augmentedEquations[j]);
                DSExpressionFree(eqLHS);
                equations[i] = NULL;
                asprintf(&equations[i], "%s = %s", temp_lhs, temp_rhs);
                DSSecureFree(temp_rhs);
                DSSecureFree(temp_lhs);
        }
        ds = DSDesignSpaceByParsingStringsWithXi(equations, Xda, DSGMASystemXi(modifiedGMA), DSCaseNumberOfEquations(aCase));
        for (i = 0; i < DSCaseNumberOfEquations(aCase); i++) {
                DSSecureFree(equations[i]);
                DSExpressionFree(caseEquations[i]);
        }
        DSSecureFree(equations);
        DSSecureFree(caseEquations);
bail:
        return ds;
}

static DSDesignSpace * dsCyclicalCaseAugmentedSystemForSubdominantDecays(const DSCase * aCase,
                                                                         const DSDesignSpace * original,
                                                                         DSMatrix * problematicEquations,
                                                                         const DSMatrixArray * problematicTerms,
                                                                         const DSMatrixArray * coefficientArray,
                                                                         const DSUInteger *subdominantDecaySpecies,
                                                                         const DSUInteger *subdominantDecayTerm)
{
        DSDesignSpace * augmentedSystem = NULL, *modifierDesignSpace;
        DSGMASystem * gma = NULL, *subcaseGMA;
        DSExpression ** augmentedEquations = NULL;
        DSUInteger i, j, k, l,*signature;
        double subCoefficient, value;
        char ** subcaseEquations;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (original == NULL) {
                DSError(M_DS_DESIGN_SPACE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (problematicEquations == NULL) {
                DSError(M_DS_MAT_NULL, A_DS_ERROR);
                goto bail;
        }
        if (problematicTerms == NULL) {
                DSError(M_DS_MAT_NULL ": Problematic term matrix array is NULL", A_DS_ERROR);
                goto bail;
        }
        if (coefficientArray == NULL) {
                DSError(M_DS_MAT_NULL ": Coefficient matrix array is NULL", A_DS_ERROR);
                goto bail;
        }
        if (subdominantDecaySpecies == NULL) {
                DSError(M_DS_NULL ": Array indicating the subdominant decay equation is NULL", A_DS_ERROR);
                goto bail;
        }
        gma = DSGMASystemCopy(DSDesignSpaceGMASystem(original));
        augmentedEquations = DSSecureCalloc(sizeof(DSExpression *), DSMatrixColumns(problematicEquations));
        subcaseEquations =  DSSecureCalloc(sizeof(char *), DSGMASystemNumberOfEquations(gma));
        signature = DSSecureMalloc(sizeof(DSUInteger)*DSGMASystemNumberOfEquations(gma)*2);
        for (i = 0; i < DSMatrixColumns(problematicEquations); i++) {
                l = 0;
                for (j = 0; j < DSMatrixRows(problematicEquations); j++) {
                        if (DSMatrixDoubleValue(problematicEquations, j, i) == 0)
                                continue;
                        if (j == subdominantDecaySpecies[i]) {
                                subCoefficient = DSMatrixArrayDoubleWithIndices(coefficientArray, i, l, 0);
                        }
                        l++;
                }
                l = 0;
//                positiveTerms = 0;
//                negativeTerms = 0;
                for (j = 0; j < DSMatrixRows(problematicEquations); j++) {
                        if (DSMatrixDoubleValue(problematicEquations, j, i) == 0)
                                continue;
                        value = DSMatrixArrayDoubleWithIndices(coefficientArray, i, l, 0);
//                        positiveTerms += DSDesignSpaceSignature(original)[2*j]-1;
                        for (k = 0; k < DSMatrixColumns(DSGMASystemAlpha(gma)); k++) {
                                if (k+1 == aCase->signature[2*j]) {
                                        DSMatrixSetDoubleValue((DSMatrix *)DSGMASystemAlpha(gma), j, k, 0.0f);
                                } else {
                                        DSMatrixSetDoubleValue((DSMatrix *)DSGMASystemAlpha(gma), j, k,
                                                               DSMatrixDoubleValue(DSGMASystemAlpha(gma), j, k)*value/subCoefficient);
                                }

                        }
                        l++;
                        augmentedEquations[i] = DSExpressionAddExpressions(augmentedEquations[i], DSGMASystemPositiveTermsForEquations(gma, j));
                }
                j = subdominantDecaySpecies[i];
//                negativeTerms = DSDesignSpaceSignature(original)[2*j+1]-1;
                for (k = 0; k < DSMatrixColumns(DSGMASystemBeta(gma)); k++) {
                        if (k+1 == aCase->signature[2*j+1]) {
                                DSMatrixSetDoubleValue((DSMatrix *)DSGMASystemBeta(gma), j, k, 0.0f);
                        } else {
                                DSMatrixSetDoubleValue((DSMatrix *)DSGMASystemBeta(gma), j, k,
                                                       DSMatrixDoubleValue(DSGMASystemBeta(gma), j, k));
                        }
                }
//                printf("%i %i ", positiveTerms, negativeTerms);
                augmentedEquations[i] = DSExpressionAddExpressions(augmentedEquations[i], DSGMASystemNegativeTermsForEquations(gma, j));
//                if (positiveTerms == 0 || negativeTerms == 0)
//                        DSExpressionPrint(augmentedEquations[i]);
        }
//        printf("\n");
        augmentedSystem = dsCyclicalCaseCreateUniqueAugmentedSystem(aCase,
                                                                    gma,
                                                                    problematicEquations,
                                                                    (const DSExpression **)augmentedEquations,
                                                                    subdominantDecaySpecies);
        DSDesignSpaceAddConditions(augmentedSystem, DSCaseCd(aCase), DSCaseCi(aCase), DSCaseDelta(aCase));
        dsAddConstraintsForSubdominantDecays(augmentedSystem, aCase, original, problematicEquations, coefficientArray, subdominantDecaySpecies, subdominantDecayTerm);
        for (i = 0; i < DSGMASystemNumberOfEquations(gma); i++) {
                subcaseEquations[i] = DSExpressionAsString(DSDesignSpaceEquations(original)[i]);
                signature[2*i] = DSCaseSignature(aCase)[2*i];
                signature[2*i+1] = DSCaseSignature(aCase)[2*i+1];
        }
        for (i = 0; i < DSMatrixColumns(problematicEquations); i++) {
                for (j = 0; j < DSMatrixRows(problematicEquations); j++) {
                        if (DSMatrixDoubleValue(problematicEquations, j, i) == 0)
                                continue;
                        signature[2*j] = 0;
                        signature[2*j+1] = 0;
                        DSSecureFree(subcaseEquations[j]);
                        subcaseEquations[j] = DSExpressionAsString(DSCaseEquations(aCase)[j]);
                }
                j =subdominantDecaySpecies[i];
                DSSecureFree(subcaseEquations[j]);
                subcaseEquations[j] = DSExpressionAsString(DSDesignSpaceEquations(augmentedSystem)[j]);
        }
        modifierDesignSpace = DSDesignSpaceByParsingStringsWithXi(subcaseEquations, DSGMASystemXd_a(gma), DSGMASystemXi(gma), DSDesignSpaceNumberOfEquations(original));
        modifierDesignSpace->Cd = DSMatrixCopy(augmentedSystem->Cd);
        modifierDesignSpace->Ci = DSMatrixCopy(augmentedSystem->Ci);
        modifierDesignSpace->delta = DSMatrixCopy(augmentedSystem->delta);
//        DSCyclicalCaseDesignSpaceCalculateSubCyclicalCases(augmentedSystem, modifierDesignSpace, signature);
//        DSDesignSpaceCalculateCyclicalCases(augmentedSystem);
//        printf("Original Case\n");
//        DSCasePrintEquations(aCase);
//        printf("New GMA\n");
//        DSGMASystemPrintEquations(DSDesignSpaceGMASystem(augmentedSystem));
//        printf("Subcase Equations\n");
//        for (i = 0; i < DSGMASystemNumberOfEquations(gma); i++) {
//                printf("%i %i %s\n", signature[2*i], signature[2*i+1], subcaseEquations[i]);
//        }
        for (i = 0; i < DSMatrixColumns(problematicEquations); i++) {
                DSExpressionFree(augmentedEquations[i]);
        }
        DSSecureFree(augmentedEquations);
        DSGMASystemFree(gma);
        DSDesignSpaceFree(modifierDesignSpace);
bail:
        return augmentedSystem;
}

extern void DSCyclicalCaseDesignSpaceCalculateSubCyclicalCase(DSDesignSpace *ds, DSCase * aCase, const DSDesignSpace * modifierDS)
{
        DSUInteger caseNumber;
        char * string = NULL;
        DSCyclicalCase * aSubcase = NULL;
        if (ds == NULL) {
                DSError(M_DS_DESIGN_SPACE_NULL, A_DS_ERROR);
                goto bail;
        }
        string = DSSecureCalloc(sizeof(char), 100);
        caseNumber = DSCaseNumber(aCase);
        sprintf(string, "%i", caseNumber);
        if (DSDictionaryValueForName(DSDesignSpaceCyclicalCaseDictionary(ds), string) == NULL) {
                aSubcase = DSCyclicalCaseForCaseInDesignSpace(modifierDS, aCase);
                if (aSubcase != NULL) {
                        DSDictionaryAddValueWithName((DSDictionary *)DSDesignSpaceCyclicalCaseDictionary(ds), string, aSubcase);
                }
        }
        if (string != NULL)
                DSSecureFree(string);
bail:
        return;
}

extern void DSCyclicalCaseDesignSpaceCalculateSubCyclicalCases(DSDesignSpace *ds, DSDesignSpace * modifierDesignSpace, const DSUInteger * modifierSignature)
{
        DSUInteger i, j, numberOfCases;
        DSCase * aCase = NULL;
        if (ds == NULL) {
                DSError(M_DS_DESIGN_SPACE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (modifierDesignSpace == NULL) {
                DSError(M_DS_DESIGN_SPACE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (modifierSignature == NULL) {
                DSError(M_DS_NULL, A_DS_ERROR);
                goto bail;
        }
        numberOfCases = DSDesignSpaceNumberOfCases(ds);
        if (numberOfCases == 0) {
                goto bail;
        }
        for (i = 0; i < numberOfCases; i++) {
                aCase = DSDesignSpaceCaseWithCaseNumber(ds, i+1);
                for (j = 0; j < DSCaseNumberOfEquations(aCase); j++) {
                        if (modifierSignature[2*j] != 0)
                                aCase->signature[2*j] = modifierSignature[2*j];
                        if (modifierSignature[2*j+1] != 0)
                                aCase->signature[2*j+1] = modifierSignature[2*j+1];
                }
                DSCyclicalCaseDesignSpaceCalculateSubCyclicalCase(ds, aCase, modifierDesignSpace);
                DSCaseFree(aCase);
        }
bail:
        return;
}

static DSUInteger dsCyclicalCaseNumberOfAugmentedSystems(const DSDesignSpace * original,
                                                         const DSMatrix * problematicEquations)
{
        DSUInteger i, j, max = 0, count;
        if (original == NULL) {
                DSError(M_DS_DESIGN_SPACE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (problematicEquations == NULL) {
                DSError(M_DS_MAT_NULL ": Matrix of problematic equations is NULL", A_DS_ERROR);
                goto bail;
        }
        max = 1;
        for (i = 0; i < DSMatrixColumns(problematicEquations); i++) {
                count = 0;
                for (j = 0; j < DSMatrixRows(problematicEquations); j++) {
                        if (DSMatrixDoubleValue(problematicEquations, j, i) == 0.0f)
                                continue;
                        count += (DSDesignSpaceSignature(original)[j*2+1] > 1);
                }
                max *= count;
        }
bail:
        return max;
}

static DSStack * dsCyclicalCaseCreateAugmentedSystems(const DSCase * aCase,
                                                      const DSDesignSpace * original,
                                                      DSMatrix * problematicEquations,
                                                      const DSMatrixArray * problematicTerms,
                                                      const DSMatrixArray * coefficientArray)
{
        DSStack * augmentedSystemsStack = NULL;
        DSUInteger i, j, k, current, index, max, *numberOfequations, numberOfTerms;
        DSUInteger * decayEquations = NULL;
        DSUInteger * decayTerms;
        DSDesignSpace * subcase;
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
        if (problematicEquations == NULL)
                goto bail;
        if (problematicTerms == NULL)
                goto bail;
        if (coefficientArray == NULL)
                goto bail;
        if (DSMatrixArrayNumberOfMatrices(problematicTerms) != DSMatrixArrayNumberOfMatrices(coefficientArray))
                goto bail;
        decayEquations = DSSecureCalloc(sizeof(DSUInteger), DSMatrixColumns(problematicEquations));
        decayTerms = DSSecureCalloc(sizeof(DSUInteger), DSMatrixColumns(problematicEquations));
        numberOfequations = DSSecureCalloc(sizeof(DSUInteger), DSMatrixColumns(problematicEquations));
        max = dsCyclicalCaseNumberOfAugmentedSystems(original, problematicEquations);
        augmentedSystemsStack = DSStackAlloc();
        for (i = 0; i < DSMatrixColumns(problematicEquations); i++) {
                for (j = 0; j < DSMatrixRows(problematicEquations); j++) {
                        if ((DSUInteger)DSMatrixDoubleValue(problematicEquations, j, i) == 0)
                                continue;
                        numberOfequations[i] += (DSDesignSpaceSignature(original)[j*2+1] > 1);
                }
        }
        for (i = 0; i < max; i++) {
                current = i;
                numberOfTerms = 0;
                for (j = 0; j < DSMatrixColumns(problematicEquations); j++) {
                        decayEquations[j] = current % numberOfequations[j];
                        index = 0;
                        for (k = 0; k < DSMatrixRows(problematicEquations); k++) {
                                if ((DSUInteger)DSMatrixDoubleValue(problematicEquations, k, j) == 0)
                                        continue;
                                if (DSDesignSpaceSignature(original)[k*2+1] > 1) {
                                        if (index == decayEquations[j])
                                                break;
                                        else
                                                index++;
                                }
                        }
                        decayEquations[j] = k;
                        numberOfTerms += DSDesignSpaceSignature(original)[k*2+1];
                        current = current / numberOfequations[j];
                }
                for (j = 0; j < numberOfTerms; j++) {
                        index = j;
                        for (k = 0; k < DSMatrixColumns(problematicEquations); k++) {
                                decayTerms[k] = index % DSDesignSpaceSignature(original)[decayEquations[k]*2+1];
                                if (DSCaseSignature(aCase)[decayEquations[k]*2+1] == decayTerms[k]+1)
                                        break;
                                index = index / DSDesignSpaceSignature(original)[decayEquations[k]*2+1];
                        }
                        if (k != DSMatrixColumns(problematicEquations))
                                continue;
                        subcase = dsCyclicalCaseAugmentedSystemForSubdominantDecays(aCase,
                                                                                    original,
                                                                                    problematicEquations,
                                                                                    problematicTerms,
                                                                                    coefficientArray,
                                                                                    decayEquations,
                                                                                    decayTerms);
                        DSStackPush(augmentedSystemsStack, subcase);
                }
        }
        DSSecureFree(decayEquations);
        DSSecureFree(decayTerms);
        DSSecureFree(numberOfequations);
bail:
        return augmentedSystemsStack;
}



#if defined (__APPLE__) && defined (__MACH__)
#pragma mark - Exposed function to generate the internal systems for cyclical cases -
#endif

extern DSStack * DSCyclicalCaseDesignSpacesForUnderdeterminedCase(const DSCase * aCase, const DSDesignSpace * original)
{
        DSStack *subcases = NULL;
        DSMatrix * problematicEquations = NULL;
        DSMatrixArray * problematicTerms = NULL;
        DSMatrixArray * coefficientArray = NULL;
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
        problematicEquations = dsSubcaseProblematicEquations(aCase);
        if (problematicEquations == NULL)
                goto bail;
        problematicTerms = dsSubcaseProblematicTerms(aCase, problematicEquations);
        if (problematicTerms == NULL)
                goto bail;
        coefficientArray = dsSubcaseCoefficientsOfInterest(aCase, problematicTerms);
        if (coefficientArray == NULL)
                goto bail;
        if (DSMatrixArrayNumberOfMatrices(problematicTerms) != DSMatrixArrayNumberOfMatrices(coefficientArray))
                goto bail;
        subcases = dsCyclicalCaseCreateAugmentedSystems(aCase,
                                                        original,
                                                        problematicEquations,
                                                        problematicTerms,
                                                        coefficientArray);
bail:
        if (problematicEquations != NULL)
                DSMatrixFree(problematicEquations);
        if (problematicTerms != NULL)
                DSMatrixArrayFree(problematicTerms);
        if (coefficientArray != NULL)
                DSMatrixArrayFree(coefficientArray);
        return subcases;
}

extern DSDesignSpace * DSCyclicalCaseInternalForUnderdeterminedCase(const DSCase * aCase, const DSDesignSpace * original)
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
        problematicEquations = dsSubcaseProblematicEquations(aCase);
        if (problematicEquations == NULL)
                goto bail;
        problematicTerms = dsSubcaseProblematicTerms(aCase, problematicEquations);
        if (problematicTerms == NULL)
                goto bail;
        coefficientArray = dsSubcaseCoefficientsOfInterest(aCase, problematicTerms);
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
                        value = DSMatrixArrayDoubleWithIndices(coefficientArray, i, l, 0);
                        for (k = 0; k < DSMatrixColumns(DSGMASystemAlpha(temp)); k++) {
                                if (k+1 == aCase->signature[2*j]) {
                                        DSMatrixSetDoubleValue((DSMatrix *)DSGMASystemAlpha(temp), j, k, 0.0f);
                                } else {
                                        DSMatrixSetDoubleValue((DSMatrix *)DSGMASystemAlpha(temp), j, k,
                                                               DSMatrixDoubleValue(DSGMASystemAlpha(temp), j, k)*value);
                                }
                                if (k+1 == aCase->signature[2*j+1]) {
                                        DSMatrixSetDoubleValue((DSMatrix *)DSGMASystemBeta(temp), j, k, 0.0f);
                                } else if (k < DSMatrixColumns(DSGMASystemBeta(temp))) {
                                        DSMatrixSetDoubleValue((DSMatrix *)DSGMASystemBeta(temp), j, k,
                                                               DSMatrixDoubleValue(DSGMASystemBeta(temp), j, k)*value);
                                }
                        }
                        l++;
                        augmentedEquations[i] = DSExpressionAddExpressions(augmentedEquations[i], DSGMASystemPositiveTermsForEquations(temp, j));
                        augmentedEquations[i] = DSExpressionAddExpressions(augmentedEquations[i], DSGMASystemNegativeTermsForEquations(temp, j));
                }
        }
        
        subcases = dsSubcaseCreateUniqueSystemSubcase(aCase, temp, problematicEquations, (const DSExpression **)augmentedEquations);
        if (subcases != NULL) {
                DSDesignSpaceAddConditions(subcases, aCase->Cd, aCase->Ci, aCase->delta);
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
        return subcases;
}
