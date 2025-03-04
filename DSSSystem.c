/**
 * \file DSSSystem.m
 * \brief Header file with functions for dealing with S-Systems.
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
#include <stdarg.h>
#include <math.h>
#include "DSTypes.h"
#include "DSErrors.h"
#include "DSMemoryManager.h"
#include "DSVariable.h"
#include "DSCase.h"
#include "DSSSystem.h"
#include "DSUnstableCase.h"
#include "DSExpression.h"
#include "DSExpressionTokenizer.h"
#include "DSSSystemGrammar.h"
#include "DSGMASystemParsingAux.h"
#include "DSMatrix.h"
#include "DSMatrixArray.h"
#include "DSGMASystem.h"
#include <gsl/gsl_matrix_double.h>

/**
 * \defgroup DSSSysACCESSORS
 *
 * \brief Internal S-System Accessor macros.
 * 
 * \details Used within DSSSystem.c to access the data within a S-System data type.
 * These macros are not to be used putside of this file, as they do not check the
 * data dor consistency and thus would not invoke the DSError function, making
 * it harder to trace errors.
 */
/*\{*/
#define DSSSysXi(x)                       ((x)->Xi)
#define DSSSysXd(x)                       ((x)->Xd)
#define DSSSysXd_a(x)                     ((x)->Xd_a)
#define DSSSysXd_t(x)                     ((x)->Xd_t)
#define DSSSysXd_b(x)                     ((x)->Xd_b)
#define DSSSysXd_a_c(x)                   ((x)->Xd_a_c)
#define DSSSysAlpha(x)                    ((x)->alpha)
#define DSSSysAlphaAdjusted(x)            ((x)->alpha_adjusted)
#define DSSSysBeta(x)                     ((x)->beta)
#define DSSSysBetaAdjusted(x)             ((x)->beta_adjusted)
#define DSSSysGd(x)                       ((x)->Gd)
#define DSSSysGi(x)                       ((x)->Gi)
#define DSSSysHd(x)                       ((x)->Hd)
#define DSSSysHi(x)                       ((x)->Hi)
#define DSSSysM(x)                        ((x)->M)

#define DSSSysIsSingular(x)               ((x)->isSingular)
#define DSSSysShouldFreeXd(x)             ((x)->shouldFreeXd)
#define DSSSysShouldFreeXi(x)             ((x)->shouldFreeXi)
/*\}*/

#if defined (__APPLE__) && defined (__MACH__)
#pragma  mark - Function Prototypes
#endif

static void dsSSystemSolveEquations(DSSSystem *ssys);


#if defined (__APPLE__) && defined (__MACH__)
#pragma  mark - Allocation, deallocation and initialization
#endif

/**
 * \brief Creates a empty S-System.
 * 
 * This function allocates the necessary memory space used by a S-System and
 * initializes it so that it is ready for processing.  The initialized GMA
 * has all its fields set to 0 or NULL.  This is interpreted as an empty GMA
 * and is necessary for parsing a set of equations.
 *
 * \return A DSSSystem pointer to the newly allocated GMASystem.
 */
static DSSSystem * DSSSystemAlloc(void)
{
        DSSSystem *sys  = NULL;
        sys  = DSSecureCalloc(sizeof(DSSSystem), 1);
        return sys ;
}

extern DSSSystem * DSSSystemCopy(const DSSSystem * original)
{
        DSSSystem * newSSys = NULL;
        if (original == NULL) {
                DSError(M_DS_NULL, A_DS_WARN);
                goto bail;
        }
        newSSys = DSSSystemAlloc();
        DSSSysXd(newSSys) = DSVariablePoolCopy(DSSSystemXd(original));
        DSSSysXd_t(newSSys) = DSVariablePoolCopy(DSSSystemXd_t(original));
        DSSSysXd_a(newSSys) = DSVariablePoolCopy(DSSSystemXd_a(original));
        if (DSSSysXd_b(original) != NULL)
            DSSSysXd_b(newSSys) = DSVariablePoolCopy(DSSSystemXd_b(original));
        if (DSSSysXd_a_c(original) != NULL)
            DSSSysXd_a_c(newSSys) = DSVariablePoolCopy(DSSSystemXd_a_c(original));
        DSSSysXi(newSSys) = DSVariablePoolCopy(DSSSystemXi(original));
        DSSSystemSetShouldFreeXd(newSSys, true);
        DSSSystemSetShouldFreeXi(newSSys, true);
        DSSSysGd(newSSys) = DSMatrixCopy(DSSSysGd(original));
        DSSSysHd(newSSys) = DSMatrixCopy(DSSSysHd(original));
        DSSSysGi(newSSys) = DSMatrixCopy(DSSSysGi(original));
        DSSSysHi(newSSys) = DSMatrixCopy(DSSSysHi(original));
        DSSSysAlpha(newSSys) = DSMatrixCopy(DSSSysAlpha(original));
        DSSSysBeta(newSSys) = DSMatrixCopy(DSSSysBeta(original));
        DSSSystemSetIsSingular(newSSys, DSSSystemIsSingular(original));
        DSSSystemSetIsUnstable(newSSys, DSSSystemIsUnstable(original));
        DSSSystemSetIsConserved(newSSys, DSSSystemIsConserved(original));
        DSSSystemSetAdjustCodominantStoichiometry(newSSys, DSSSystemAdjustCodominantStoichiometry(original));
        if (DSSSystemIsConserved(newSSys) == true){
            newSSys->numberOfConservations = original->numberOfConservations;
        }
        if (DSSSystemIsSingular(newSSys) == false) {
                DSSSysM(newSSys) = DSMatrixCopy(DSSSysM(original));
        }
        if (DSSSystemAdjustCodominantStoichiometry(newSSys) == true)
            DSSSysAlphaAdjusted(newSSys) = DSMatrixCopy(DSSSysAlphaAdjusted(original));
bail:
        return newSSys;
}

extern void DSSSystemFree(DSSSystem * sys)
{
        if (sys  == NULL) {
                DSError(M_DS_NULL ": S-System to free is NULL", A_DS_ERROR);
                goto bail;
        }
        if (DSSSystemShouldFreeXd(sys) == true) {
                DSVariablePoolSetReadWriteAdd(DSSSysXd(sys));
                DSVariablePoolFree(DSSSysXd(sys));
                if (DSSSysXd_t(sys) != NULL) {
                        DSVariablePoolSetReadWriteAdd(DSSSysXd_t(sys));
                        DSVariablePoolFree(DSSSysXd_t(sys));
                }
                if (DSSSysXd_a(sys) != NULL) {
                        DSVariablePoolSetReadWriteAdd(DSSSysXd_a(sys));
                        DSVariablePoolFree(DSSSysXd_a(sys));
                }
                if (DSSSysXd_b(sys) != NULL) {
                    DSVariablePoolSetReadWriteAdd(DSSSysXd_b(sys));
                    DSVariablePoolFree(DSSSysXd_b(sys));
                }
                if (DSSSysXd_a_c(sys) != NULL) {
                    DSVariablePoolSetReadWriteAdd(DSSSysXd_a_c(sys));
                    DSVariablePoolFree(DSSSysXd_a_c(sys));
                }
        }
        if (DSSSystemShouldFreeXi(sys) == true) {
                DSVariablePoolSetReadWriteAdd(DSSSysXi(sys));
                DSVariablePoolFree(DSSSysXi(sys));
        }
        if (DSSSysAlpha(sys) != NULL)
                DSMatrixFree(DSSSysAlpha(sys));
        if (DSSSysBeta(sys) != NULL)
                DSMatrixFree(DSSSysBeta(sys));
        if (DSSSysGd(sys) != NULL)
                DSMatrixFree(DSSSysGd(sys));
        if (DSSSysGi(sys) != NULL)
                DSMatrixFree(DSSSysGi(sys));
        if (DSSSysHd(sys) != NULL)
                DSMatrixFree(DSSSysHd(sys));
        if (DSSSysHi(sys) != NULL)
                DSMatrixFree(DSSSysHi(sys));
        if (DSSSysM(sys) != NULL)
                DSMatrixFree(DSSSysM(sys));
        if (DSSSysAlphaAdjusted(sys) != NULL)
                DSMatrixFree(DSSSysAlphaAdjusted(sys));
        DSSecureFree(sys);
bail:
        return;
}

#if defined (__APPLE__) && defined (__MACH__)
#pragma  mark - Factory methods
#endif



#if defined (__APPLE__) && defined (__MACH__)
#pragma  mark Internal Parsing Functions
#endif

static gma_parseraux_t * dsSSystemParseStringToTermList(const char * string)
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
        parser = DSSSystemParserAlloc(DSSecureMalloc);
        root = DSGMAParserAuxAlloc();
        parser_aux = root;
        current = tokens;
        while (current != NULL) {
                if (DSExpressionTokenType(current) == DS_EXPRESSION_TOKEN_START) {
                        current = DSExpressionTokenNext(current);
                        continue;
                }
                DSSSystemParser(parser, 
                                DSExpressionTokenType(current),
                                current,
                                ((void**)&parser_aux));
                current = DSExpressionTokenNext(current);
        }
        DSSSystemParser(parser, 
                          0, 
                          NULL,
                          ((void **)&parser_aux));
        DSSSystemParserFree(parser, DSSecureFree);
        DSExpressionTokenFree(tokens);
        if (DSGMAParserAuxParsingFailed(root) == true) {
                DSGMAParserAuxFree(root);
                root = NULL;
        }
bail:
        return root;
}

static gma_parseraux_t ** dsSSysTermListForAllStrings(char * const * const strings, const DSUInteger numberOfEquations)
{
        DSUInteger i;
        gma_parseraux_t **aux = NULL;
        DSExpression *expr;
        char *aString;
        bool failed = false;
        aux = DSSecureCalloc(sizeof(gma_parseraux_t *), numberOfEquations);
        for (i = 0; i < numberOfEquations; i++) {
                if (strings[i] == NULL) {
                        DSError(M_DS_WRONG ": String to parse is NULL", A_DS_ERROR);
                        failed = true;
                        break;
                }
                if (strlen(strings[i]) == 0) {
                        DSError(M_DS_WRONG ": String to parse is empty", A_DS_ERROR);
                        failed = true;
                        break;
                }
                expr = DSExpressionByParsingString(strings[i]);
                if (expr != NULL) {
                        aString = DSExpressionAsString(expr);
                        aux[i] = dsSSystemParseStringToTermList(aString);
                        DSSecureFree(aString);
                        DSExpressionFree(expr);
                }
                if (aux[i] == NULL) {
                        DSError(M_DS_PARSE ": Expression not in S-System format", A_DS_ERROR);
                        failed = true;
                        break;
                }
        }
        if (failed == true) {
                for (i = 0; i < numberOfEquations; i++)
                        if (aux[i] != NULL)
                                DSGMAParserAuxFree(aux[i]);
                DSSecureFree(aux);
                aux = NULL;
        }
bail:
        return aux;
}

static DSVariablePool * dsSSystemIdentifyIndependentVariables(const DSVariablePool * const Xd, gma_parseraux_t ** aux, const DSUInteger numberOfEquations)
{
        DSVariablePool * Xi = NULL;
        gma_parseraux_t *current = NULL;
        DSUInteger i, j;
        const char *name;
        if (aux == NULL) {
                DSError(M_DS_NULL ": Parser auxiliary data is NULL", A_DS_ERROR);
                goto bail;
        }
        if (Xd == NULL) {
                DSError(M_DS_NULL ": Independent variables are NULL", A_DS_ERROR);
                goto bail;
        }
        if (numberOfEquations == 0) {
                DSError(M_DS_WRONG ": No equations to parse", A_DS_WARN);
                goto bail;
        }
        if (DSVariablePoolNumberOfVariables(Xd) != numberOfEquations) {
                DSError(M_DS_WRONG ": Number of independent variables does not match number of equations", A_DS_ERROR);
                goto bail;
        }
        Xi = DSVariablePoolAlloc();
        for (i = 0; i < numberOfEquations; i++) {
                current = aux[i];
                while (current) {
                        for (j = 0; j < DSGMAParserAuxNumberOfBases(current); j++) {
                                if (DSGMAParserAuxBaseAtIndexIsVariable(current, j) == false)
                                        continue;
                                name = DSGMAParserAuxVariableAtIndex(current, j);
                                if (DSVariablePoolHasVariableWithName(Xd, name) == true)
                                        continue;
                                if (DSVariablePoolHasVariableWithName(Xi, name) == true)
                                        continue;
                                DSVariablePoolAddVariableWithName(Xi, name);
                        }
                        current = DSGMAParserAuxNextNode(current);
                }
        }
bail:
        return Xi;
}

static void dsSSystemInitializeMatrices(DSSSystem *sys)
{
        DSUInteger numberOfEquations, numberOfXd, numberOfXi;
        numberOfEquations = DSVariablePoolNumberOfVariables(DSSSysXd(sys));
        numberOfXd = numberOfEquations;
        numberOfXi = DSVariablePoolNumberOfVariables(DSSSysXi(sys));
        DSSSysAlpha(sys) = DSMatrixCalloc(numberOfEquations, 1);
        DSSSysBeta(sys) = DSMatrixCalloc(numberOfEquations, 1);
        DSSSysGd(sys) = DSMatrixCalloc(numberOfEquations, numberOfXd);
        if (numberOfXi != 0)
                DSSSysGi(sys) = DSMatrixCalloc(numberOfEquations, numberOfXi);
        DSSSysHd(sys) = DSMatrixCalloc(numberOfEquations, numberOfXd);
        if (numberOfXi != 0)
                DSSSysHi(sys) = DSMatrixCalloc(numberOfEquations, numberOfXi);
}

static void dsSSysProcessPositiveExponentBasePairs(DSSSystem *sys, gma_parseraux_t *current, DSUInteger equation)
{
        DSUInteger j;
        const char *varName;
        double currentValue;
        for (j = 0; j < DSGMAParserAuxNumberOfBases(current); j++) {
                if (DSGMAParserAuxBaseAtIndexIsVariable(current, j) == false) { 
                        DSMatrixSetDoubleValue(DSSSysAlpha(sys), 
                                               equation, 0, 
                                               DSGMAParseAuxsConstantBaseAtIndex(current, j));
                        continue;
                }
                varName = DSGMAParserAuxVariableAtIndex(current, j);
                if (DSVariablePoolHasVariableWithName(DSSSysXd(sys), varName) == true) {
                        currentValue = DSMatrixDoubleValue(DSSSysGd(sys), equation, DSVariablePoolIndexOfVariableWithName(DSSSysXd(sys),
                                                                                                                          varName));
                        DSMatrixSetDoubleValue(DSSSysGd(sys),
                                               equation,
                                               DSVariablePoolIndexOfVariableWithName(DSSSysXd(sys),
                                                                                     varName),
                                               currentValue+DSGMAParserAuxExponentAtIndex(current, j));
                } else if (DSVariablePoolHasVariableWithName(DSSSysXi(sys), varName) == true) {
                        currentValue = DSMatrixDoubleValue(DSSSysHi(sys), equation, DSVariablePoolIndexOfVariableWithName(DSSSysXi(sys),
                                                                                                                          varName));
                        DSMatrixSetDoubleValue(DSSSysGi(sys),
                                               equation,
                                               DSVariablePoolIndexOfVariableWithName(DSSSysXi(sys),
                                                                                     varName),
                                               currentValue+DSGMAParserAuxExponentAtIndex(current, j));
                }
        }
}

static void dsSSysProcessNegativeExponentBasePairs(DSSSystem *sys, gma_parseraux_t *current, DSUInteger equation)
{
        DSUInteger j;
        const char *varName;
        double currentValue;
        for (j = 0; j < DSGMAParserAuxNumberOfBases(current); j++) {
                if (DSGMAParserAuxBaseAtIndexIsVariable(current, j) == false) { 
                        DSMatrixSetDoubleValue(DSSSysBeta(sys), 
                                               equation, 0, 
                                               DSGMAParseAuxsConstantBaseAtIndex(current, j));
                        continue;
                }
                varName = DSGMAParserAuxVariableAtIndex(current, j);
                if (DSVariablePoolHasVariableWithName(DSSSysXd(sys), varName) == true) {
                        currentValue = DSMatrixDoubleValue(DSSSysHd(sys), equation, DSVariablePoolIndexOfVariableWithName(DSSSysXd(sys),
                                                                                                                          varName));
                        DSMatrixSetDoubleValue(DSSSysHd(sys),
                                               equation,
                                               DSVariablePoolIndexOfVariableWithName(DSSSysXd(sys),
                                                                                     varName),
                                               currentValue+DSGMAParserAuxExponentAtIndex(current, j));
                } else if (DSVariablePoolHasVariableWithName(DSSSysXi(sys), varName) == true) {
                        currentValue = DSMatrixDoubleValue(DSSSysHi(sys), equation, DSVariablePoolIndexOfVariableWithName(DSSSysXi(sys),
                                                                                                                          varName));
                        DSMatrixSetDoubleValue(DSSSysHi(sys),
                                               equation,
                                               DSVariablePoolIndexOfVariableWithName(DSSSysXi(sys),
                                                                                     varName),
                                               currentValue+DSGMAParserAuxExponentAtIndex(current, j));
                }
        }
}

#include <unistd.h> 

static void dsSSystemSolveEquations(DSSSystem *ssys)
{
        DSMatrix *M, *Ad;
        if (ssys == NULL) {
                DSError(M_DS_NULL ": S-System being modified is NULL", A_DS_ERROR);
                goto bail;
        }
        DSSSystemSetIsSingular(ssys, true);
        Ad = DSMatrixBySubstractingMatrix(DSSSysGd(ssys), DSSSysHd(ssys));
        M = DSMatrixInverse(Ad);
        if (M != NULL) {
                DSSSystemSetIsSingular(ssys, false);
                DSSSysM(ssys) = M;
        }
        DSMatrixFree(Ad);
bail:
        return;
}

extern void DSSSystemRecalculateSolution(DSSSystem * ssys)
{
        if (ssys == NULL) {
                DSError(M_DS_NULL ": S-System being modified is NULL", A_DS_ERROR);
                goto bail;
        }
        if (DSSSysM(ssys) != NULL) {
                DSMatrixFree(DSSSysM(ssys));
                DSSSysM(ssys) = NULL;
                dsSSystemSolveEquations(ssys);
        }
bail:
        return;
}

static void dsSSystemCreateSystemMatrices(DSSSystem *sys, gma_parseraux_t **aux)
{
        gma_parseraux_t *current;
        DSUInteger numberOfEquations;
        DSUInteger i;
        if (sys == NULL) {
                DSError(M_DS_NULL ": S-System being modified is NULL", A_DS_ERROR);
                goto bail;
        }
        if (aux == NULL) {
                DSError(M_DS_NULL ": Parser auxiliary data is NULL", A_DS_ERROR);
                goto bail;
        }
        if (DSSSysXd(sys) == NULL || DSSSysXi(sys) == NULL) {
                DSError(M_DS_WRONG ": S-System data is incomplete: Need Xi and Xd", A_DS_ERROR);
                goto bail;
        }
        numberOfEquations = DSVariablePoolNumberOfVariables(DSSSysXd(sys));
        dsSSystemInitializeMatrices(sys);
        for (i = 0; i < numberOfEquations; i++) {
                current = aux[i];
                while (current) {
                        switch (DSGMAParserAuxSign(current)) {
                                case AUX_SIGN_POSITIVE:
                                        DSMatrixSetDoubleValue(DSSSysAlpha(sys), 
                                                               i, 0, 
                                                               1.0);
                                        dsSSysProcessPositiveExponentBasePairs(sys, current, i);
                                        break;
                                case AUX_SIGN_NEGATIVE:
                                        DSMatrixSetDoubleValue(DSSSysBeta(sys), 
                                                               i, 0, 
                                                               1.0);
                                        dsSSysProcessNegativeExponentBasePairs(sys, current, i);
                                        break;
                                default:
                                        break;
                        }
                        current = DSGMAParserAuxNextNode(current);
                }
        }
        dsSSystemSolveEquations(sys);
bail:
        return;
}

DSSSystem * dsSSystemWithAlgebraicConstraints(const DSSSystem * originalSystem, DSVariablePool * newXd, DSUInteger numberDifferentialVariables, DSUInteger * differentialIndices, DSUInteger numberOfAlgebraicVariables, DSUInteger * algebraicIndices) {
        DSUInteger i, j, k;
        DSSSystem * collapsedSSystem = NULL;
        DSMatrix * oldAd, *oldAi, *oldB;
        DSMatrix *temp, *subM, *subAdAlgebraic, *subAd, *subB, *subAi;
        double value, factor;
        oldAd = DSSSystemAd(originalSystem);
        oldAi = DSSSystemAi(originalSystem);
        oldB = DSSSystemB(originalSystem);
        subAdAlgebraic = DSMatrixSubMatrixExcludingRowsAndColumns(oldAd, numberDifferentialVariables, numberDifferentialVariables, differentialIndices, differentialIndices);
        subAd = DSMatrixSubMatrixExcludingRowsAndColumns(oldAd, numberDifferentialVariables, numberOfAlgebraicVariables, differentialIndices, algebraicIndices);
        subAi = DSMatrixSubMatrixExcludingRows(oldAi, numberDifferentialVariables, differentialIndices);
        DSMatrixMultiplyByScalar(subAi, -1.0f);
        subB = DSMatrixSubMatrixExcludingRows(oldB, numberDifferentialVariables, differentialIndices);
        DSMatrixMultiplyByScalar(subAd, -1.0f);
        subM = DSMatrixInverse(subAdAlgebraic);
        DSMatrixFree(subAdAlgebraic);
        temp = DSMatrixByMultiplyingMatrix(subM, subAd);
        DSMatrixFree(subAd);
        subAd = temp;
        temp = DSMatrixByMultiplyingMatrix(subM, subAi);
        DSMatrixFree(subAi);
        subAi = temp;
        temp = DSMatrixByMultiplyingMatrix(subM, subB);
        DSMatrixFree(subB);
        subB = temp;
        collapsedSSystem = DSSSystemAlloc();
        DSSSysXd(collapsedSSystem) = newXd;
        DSSSysXd_t(collapsedSSystem) = DSVariablePoolCopy(newXd);
        DSSSysXd_a(collapsedSSystem) = DSVariablePoolAlloc();
        DSSSysXi(collapsedSSystem) = (DSVariablePool *)DSSSystemXi(originalSystem);
        DSSSystemSetShouldFreeXd(collapsedSSystem, true);
        DSSSystemSetShouldFreeXi(collapsedSSystem, false);
        DSSSysGd(collapsedSSystem) = DSMatrixSubMatrixExcludingRowsAndColumns(DSSSystemGd(originalSystem), numberOfAlgebraicVariables,numberOfAlgebraicVariables, algebraicIndices, algebraicIndices);
        DSSSysHd(collapsedSSystem) = DSMatrixSubMatrixExcludingRowsAndColumns(DSSSystemHd(originalSystem), numberOfAlgebraicVariables,numberOfAlgebraicVariables, algebraicIndices, algebraicIndices);
        DSSSysGi(collapsedSSystem) = DSMatrixSubMatrixExcludingRows(DSSSystemGi(originalSystem), numberOfAlgebraicVariables, algebraicIndices);
        DSSSysHi(collapsedSSystem) = DSMatrixSubMatrixExcludingRows(DSSSystemHi(originalSystem), numberOfAlgebraicVariables, algebraicIndices);
        DSSSysAlpha(collapsedSSystem) = DSMatrixSubMatrixExcludingRows(DSSSystemAlpha(originalSystem), numberOfAlgebraicVariables, algebraicIndices);
        DSSSysBeta(collapsedSSystem) = DSMatrixSubMatrixExcludingRows(DSSSystemBeta(originalSystem), numberOfAlgebraicVariables, algebraicIndices);
        for (i = 0; i < DSMatrixRows(DSSSysGd(collapsedSSystem)); i++) {
                for (j = 0; j < numberOfAlgebraicVariables; j++) {
                        factor = DSMatrixDoubleValue(DSSSystemGd(originalSystem), i, algebraicIndices[j]);
                        for (k = 0; k < DSMatrixColumns(DSSSysGd(collapsedSSystem)); k++) {
                                value = DSMatrixDoubleValue(DSSSysGd(collapsedSSystem),
                                                            i, k)+ factor*DSMatrixDoubleValue(subAd,
                                                                                              j, k);
                                DSMatrixSetDoubleValue(DSSSysGd(collapsedSSystem),
                                                       i, k, value);
                        }
                        for (k = 0; k < DSMatrixColumns(DSSSysGi(collapsedSSystem)); k++) {
                                value = DSMatrixDoubleValue(DSSSysGi(collapsedSSystem),
                                                            i, k)+ factor*DSMatrixDoubleValue(subAi,
                                                                                              j, k);
                                DSMatrixSetDoubleValue(DSSSysGi(collapsedSSystem),
                                                       i, k, value);
                        }
                        value = DSMatrixDoubleValue(DSSSysAlpha(collapsedSSystem), i, 0.0f)+ factor*DSMatrixDoubleValue(subB,
                                                                                                                        j, 0.0f);
                        DSMatrixSetDoubleValue(DSSSysAlpha(collapsedSSystem), i, 0.0f, value);
                        factor = DSMatrixDoubleValue(DSSSystemHd(originalSystem), i, algebraicIndices[j]);
                        for (k = 0; k < DSMatrixColumns(DSSSysHd(collapsedSSystem)); k++) {
                                value = DSMatrixDoubleValue(DSSSysHd(collapsedSSystem),
                                                            i, k)+ factor*DSMatrixDoubleValue(subAd,
                                                                                              j, k);
                                DSMatrixSetDoubleValue(DSSSysHd(collapsedSSystem),
                                                       i, k, value);
                        }
                        for (k = 0; k < DSMatrixColumns(DSSSysHi(collapsedSSystem)); k++) {
                                value = DSMatrixDoubleValue(DSSSysHi(collapsedSSystem),
                                                            i, k)+ factor*DSMatrixDoubleValue(subAi,
                                                                                              j, k);
                                DSMatrixSetDoubleValue(DSSSysHi(collapsedSSystem),
                                                       i, k, value);
                        }
                        value = DSMatrixDoubleValue(DSSSysBeta(collapsedSSystem), i, 0.0f)+ factor*DSMatrixDoubleValue(subB,
                                                                                                                    j, 0.0f);
                        DSMatrixSetDoubleValue(DSSSysBeta(collapsedSSystem), i, 0.0f, value);
                }
        }
        dsSSystemSolveEquations(collapsedSSystem);
        DSMatrixFree(subM);
        DSMatrixFree(subAd);
        DSMatrixFree(subAi);
        DSMatrixFree(subB);
        DSMatrixFree(oldAd);
        DSMatrixFree(oldAi);
        DSMatrixFree(oldB);
bail:
        return collapsedSSystem;
}

extern DSMatrixArray * DSSSystemSolvedAuxiliaryVariableMatrices(const DSSSystem * originalSSystem)
{
        DSMatrixArray * matrices = NULL;
        DSUInteger numberOfAlgebraicVaiables = 0, numberOfDifferentialVariables = 0;
        DSUInteger i, j, k;
        DSUInteger * algebraicIndices, *differentialIndices;
        const DSVariablePool * oldXd;
        DSVariablePool * newXd = NULL;
        char * name;
        DSMatrix * Ad_aa, * Ad_tt, *Ad_ta, *Ad_at;
        DSMatrix * Ai_a, *Ai_t, *B_a, *B_t, *M_a;
        DSMatrix * Ad_A, * Ai_A, * B_A, *temp;
        if (originalSSystem == NULL) {
                DSError(M_DS_SSYS_NULL, A_DS_ERROR);
                goto bail;
        }
        if (DSSSystemXd_a(originalSSystem) == NULL) {
                DSError(M_DS_VAR_NULL ": Xd_a variable pool is NULL", A_DS_ERROR);
                goto bail;
        }
        if (DSVariablePoolNumberOfVariables(DSSSystemXd_a(originalSSystem)) > DSVariablePoolNumberOfVariables(DSSSystemXd(originalSSystem))) {
                DSError(M_DS_WRONG ": Number of algebraic variables exceeds number of total variables", A_DS_ERROR);
                goto bail;
        }
        if (DSVariablePoolNumberOfVariables(DSSSystemXd_a(originalSSystem)) == 0) {
                goto bail;
        }
        oldXd = DSSSystemXd(originalSSystem);
        newXd = DSVariablePoolAlloc();
        for (i = 0; i < DSVariablePoolNumberOfVariables(oldXd); i++) {
                name = DSVariableName(DSVariablePoolVariableAtIndex(oldXd, i));
                if (DSVariablePoolHasVariableWithName(DSSSystemXd_a(originalSSystem), name) == true) {
                        continue;
                }
                DSVariablePoolAddVariableWithName(newXd, name);
        }
        if (DSVariablePoolNumberOfVariables(oldXd) - DSVariablePoolNumberOfVariables(newXd) != DSVariablePoolNumberOfVariables(DSSSystemXd_a(originalSSystem))) {
                DSError(M_DS_WRONG, A_DS_ERROR);
                DSVariablePoolFree(newXd);
                goto bail;
        }
        numberOfDifferentialVariables = DSVariablePoolNumberOfVariables(newXd);
        numberOfAlgebraicVaiables = DSVariablePoolNumberOfVariables(DSSSystemXd_a(originalSSystem));
        differentialIndices = DSSecureCalloc(sizeof(DSUInteger), numberOfDifferentialVariables);
        algebraicIndices = DSSecureCalloc(sizeof(DSUInteger), numberOfAlgebraicVaiables);
        for (i = 0, j = 0, k = 0; i < DSVariablePoolNumberOfVariables(oldXd); i++) {
                if (DSVariablePoolHasVariableWithName(newXd, DSVariableName(DSVariablePoolVariableAtIndex(oldXd, i))) == true) {
                        differentialIndices[j++] = i;
                } else {
                        algebraicIndices[k++] = i;
                }
        }
        Ad_aa = DSMatrixSubMatrixIncludingRowsAndColumns(DSSSystemAd(originalSSystem), numberOfAlgebraicVaiables, numberOfAlgebraicVaiables, algebraicIndices, algebraicIndices);
        Ad_ta = DSMatrixSubMatrixIncludingRowsAndColumns(DSSSystemAd(originalSSystem), numberOfDifferentialVariables, numberOfAlgebraicVaiables, differentialIndices, algebraicIndices);
        Ad_at = DSMatrixSubMatrixIncludingRowsAndColumns(DSSSystemAd(originalSSystem), numberOfAlgebraicVaiables, numberOfDifferentialVariables, algebraicIndices, differentialIndices);
        Ad_tt = DSMatrixSubMatrixIncludingRowsAndColumns(DSSSystemAd(originalSSystem), numberOfDifferentialVariables, numberOfDifferentialVariables, differentialIndices, differentialIndices);
        Ai_a = DSMatrixSubMatrixIncludingRows(DSSSystemAi(originalSSystem), numberOfAlgebraicVaiables, algebraicIndices);
        Ai_t = DSMatrixSubMatrixIncludingRows(DSSSystemAi(originalSSystem), numberOfDifferentialVariables, differentialIndices);
        B_a = DSMatrixSubMatrixIncludingRows(DSSSystemB(originalSSystem), numberOfAlgebraicVaiables, algebraicIndices);
        B_t = DSMatrixSubMatrixIncludingRows(DSSSystemB(originalSSystem), numberOfDifferentialVariables, differentialIndices);
        M_a = DSMatrixInverse(Ad_aa);
        
        temp = DSMatrixByMultiplyingMatrix(M_a, Ad_at);
        Ad_A = DSMatrixByMultiplyingMatrix(Ad_ta, temp);
        DSMatrixFree(temp);
        DSMatrixSubstractByMatrix(Ad_tt, Ad_A);
        DSMatrixFree(Ad_A);
        Ad_A = Ad_tt;
        
        temp = DSMatrixByMultiplyingMatrix(M_a, Ai_a);
        Ai_A = DSMatrixByMultiplyingMatrix(Ad_ta, temp);
        DSMatrixFree(temp);
        DSMatrixSubstractByMatrix(Ai_t, Ai_A);
        DSMatrixFree(Ai_A);
        Ai_A = Ai_t;

        temp = DSMatrixByMultiplyingMatrix(M_a, B_a);
        B_A = DSMatrixByMultiplyingMatrix(Ad_ta, temp);
        DSMatrixFree(temp);
        DSMatrixAddByMatrix(B_t, B_A);
        DSMatrixFree(B_A);
        B_A = B_t;
        
        DSMatrixFree(Ad_aa);
        DSMatrixFree(Ad_at);
        DSMatrixFree(Ad_ta);
        DSMatrixFree(Ad_tt);
        DSMatrixFree(Ai_a);
        DSMatrixFree(Ai_t);
        DSMatrixFree(B_a);
        DSMatrixFree(B_t);
        
        DSSecureFree(differentialIndices);
        DSSecureFree(algebraicIndices);
        matrices = DSMatrixArrayAlloc();
        DSMatrixArrayAddMatrix(matrices, Ad_A);
        DSMatrixArrayAddMatrix(matrices, Ai_A);
        DSMatrixArrayAddMatrix(matrices, B_A);
        
        DSMatrixFree(Ai_A);
        DSMatrixFree(Ad_A);
        DSMatrixFree(B_A);
bail:
        return matrices;
}

static void dsSSystemMatricesByPartitioningAuxiliaryVariables(const DSSSystem * ssystem,
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

static void dsSSystemSolutionForAlgebraicConstraints(const DSSSystem * ssystem,
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

static DSSSystem * dsSSystemDSSSystemByRemovingAlgebraicConstraintsInternal(const DSSSystem * originalSSystem,
                                                                            const DSMatrix * MAd_at,
                                                                            const DSMatrix * MAi_a,
                                                                            const DSMatrix * MB_a)
{
        DSSSystem * collapsedSSystem = NULL;
        DSMatrix *Kd, *Ki, *Kd_t, *Kd_a, *a, *LHS, *RHS, *temp;
        DSUInteger i, * indices = NULL, auxiliary_count, ii;
        const DSVariablePool * Xd_t, *Xd_a, *Xd;
        if (originalSSystem == NULL) {
                DSError(M_DS_SSYS_NULL, A_DS_ERROR);
                goto bail;
        }
        if (DSSSystemXd_a(originalSSystem) == NULL) {
                DSError(M_DS_VAR_NULL ": Xd_a variable pool is NULL", A_DS_ERROR);
                goto bail;
        }
        if (DSVariablePoolNumberOfVariables(DSSSystemXd_a(originalSSystem)) > DSVariablePoolNumberOfVariables(DSSSystemXd(originalSSystem))) {
                DSError(M_DS_WRONG ": Number of algebraic variables exceeds number of total variables", A_DS_ERROR);
                goto bail;
        }
        if (DSVariablePoolNumberOfVariables(DSSSystemXd_a(originalSSystem)) == 0) {
                collapsedSSystem = DSSSystemCopy(originalSSystem);
                goto bail;
        }
        if (DSVariablePoolNumberOfVariables(DSSSystemXd_t(originalSSystem)) == 0) {
                DSError(M_DS_WRONG ": System does not have dynamic variables.", A_DS_WARN);
                goto bail;
        }
        if (MAd_at == NULL || MAi_a == NULL || MB_a == NULL) {
                DSError(M_DS_MAT_NULL, A_DS_ERROR);
                goto bail;
        }
        Xd = DSSSystemXd(originalSSystem);
        Xd_t = DSSSystemXd_t(originalSSystem);
        Xd_a = DSSSystemXd_a(originalSSystem);
        indices = DSVariablePoolIndicesOfSubPool(Xd, Xd_a);
        auxiliary_count = DSVariablePoolNumberOfVariables(Xd_a);
        
        collapsedSSystem = DSSSystemAlloc();
        DSSSysXd(collapsedSSystem) = DSVariablePoolCopy(Xd_t);
        DSSSysXd_t(collapsedSSystem) = DSVariablePoolCopy(Xd_t);
        DSSSysXd_a(collapsedSSystem) = DSVariablePoolAlloc();
        if (DSSSysXd_b(originalSSystem) != NULL)
            DSSSysXd_b(collapsedSSystem) = DSVariablePoolCopy(DSSSystemXd_b(originalSSystem));
        DSSSysXi(collapsedSSystem) = DSVariablePoolCopy(DSSSystemXi(originalSSystem));
        DSSSystemSetShouldFreeXd(collapsedSSystem, true);
        DSSSystemSetShouldFreeXi(collapsedSSystem, true);
        DSSSystemSetIsUnstable(collapsedSSystem, DSSSystemIsUnstable(originalSSystem));
        DSSSystemSetIsConserved(collapsedSSystem, DSSSystemIsConserved(originalSSystem));
        for (i = 0; i < 2; i++) {
                if (i == 0) {
                        Kd = (DSMatrix *)DSMatrixSubMatrixExcludingRows(DSSSystemGd(originalSSystem),
                                                                        auxiliary_count,
                                                                        indices);
                        Ki = DSMatrixSubMatrixExcludingRows(DSSSystemGi(originalSSystem),
                                                            auxiliary_count,
                                                            indices);
                        a = DSMatrixSubMatrixExcludingRows(DSSSystemAlpha(originalSSystem), auxiliary_count, indices);
                } else {
                        Kd = (DSMatrix *)DSMatrixSubMatrixExcludingRows(DSSSystemHd(originalSSystem),
                                                                        auxiliary_count,
                                                                        indices);
                        Ki = DSMatrixSubMatrixExcludingRows(DSSSystemHi(originalSSystem),
                                                            auxiliary_count,
                                                            indices);
                        a = DSMatrixSubMatrixExcludingRows(DSSSystemBeta(originalSSystem), auxiliary_count, indices);
                }
                
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
                temp = DSMatrixByMultiplyingScalar(Kd_a, -1);
                DSMatrixFree(Kd_a);
                Kd_a = temp;
                // transform matrix LHS to logarithmic coordinates
                for (ii=0; ii<DSMatrixRows(LHS); ii++ )
                    DSMatrixSetDoubleValue(LHS, ii, 0, log10(DSMatrixDoubleValue(LHS, ii, 0)));
                RHS = DSMatrixByMultiplyingMatrix(Kd_a, MB_a);
                DSMatrixFree(Kd_a);
                a = DSMatrixByAddingMatrix(LHS, RHS);
                // transform a back into cartesian coordinates
                for (ii=0; ii<DSMatrixRows(a); ii++ )
                    DSMatrixSetDoubleValue(a, ii, 0, pow(10, DSMatrixDoubleValue(a, ii, 0)));
                DSMatrixFree(LHS);
                DSMatrixFree(RHS);
                if (i == 0) {
                        DSSSysGd(collapsedSSystem) = Kd_t;
                        DSSSysGi(collapsedSSystem) = Ki;
                        DSSSysAlpha(collapsedSSystem) = a;
                } else {
                        DSSSysHd(collapsedSSystem) = Kd_t;
                        DSSSysHi(collapsedSSystem) = Ki;
                        DSSSysBeta(collapsedSSystem) = a;
                }
        }
        DSSecureFree(indices);
        dsSSystemSolveEquations(collapsedSSystem);
bail:
        return collapsedSSystem;
}

extern DSSSystem * DSSSystemByRemovingAlgebraicConstraints(const DSSSystem * originalSSystem)
{
        DSSSystem * collapsedSystem = NULL;
        DSMatrix * M_a, * Ad_at = NULL, * Ad_aa = NULL, * Ai_a = NULL, * B_a = NULL;
        if (originalSSystem == NULL) {
                DSError(M_DS_SSYS_NULL, A_DS_ERROR);
                goto bail;
        }
        if (DSSSystemXd_a(originalSSystem) == NULL) {
                DSError(M_DS_VAR_NULL ": Xd_a variable pool is NULL", A_DS_ERROR);
                goto bail;
        }
        if (DSVariablePoolNumberOfVariables(DSSSystemXd_a(originalSSystem)) > DSVariablePoolNumberOfVariables(DSSSystemXd(originalSSystem))) {
                DSError(M_DS_WRONG ": Number of algebraic variables exceeds number of total variables", A_DS_ERROR);
                goto bail;
        }
        if (DSVariablePoolNumberOfVariables(DSSSystemXd_a(originalSSystem)) == 0) {
                collapsedSystem = DSSSystemCopy(originalSSystem);
                goto bail;
        }
        if (DSVariablePoolNumberOfVariables(DSSSystemXd_t(originalSSystem)) == 0) {
                DSError(M_DS_WRONG ": System does not have dynamic variables.", A_DS_WARN);
                goto bail;
        }
        dsSSystemMatricesByPartitioningAuxiliaryVariables(originalSSystem,
                                                          &Ad_at,&Ad_aa,&Ai_a,&B_a);
        if (Ad_aa == NULL) {
                goto bail;
        }
        M_a = DSMatrixInverse(Ad_aa);
        if (M_a == NULL) {
                goto bail;
        }
        dsSSystemSolutionForAlgebraicConstraints(originalSSystem,
                                                 M_a,
                                                 &Ad_at,
                                                 &Ai_a,
                                                 &B_a);
        DSMatrixFree(M_a);
        collapsedSystem = dsSSystemDSSSystemByRemovingAlgebraicConstraintsInternal(originalSSystem,
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

        return collapsedSystem;
}

extern DSSSystem * DSSSystemByRemovingAlgebraicConstraints_old(const DSSSystem * originalSSystem)
{
        DSUInteger numberOfAlgebraicVaiables = 0, numberOfDifferentialVariables = 0;
        DSUInteger i, j, k;
        DSUInteger * algebraicIndices, *differentialIndices;
        DSSSystem * collapsedSystem = NULL;
        const DSVariablePool * oldXd;
        DSVariablePool * newXd = NULL;
        char * name;
        if (originalSSystem == NULL) {
                DSError(M_DS_SSYS_NULL, A_DS_ERROR);
                goto bail;
        }
        if (DSSSystemXd_a(originalSSystem) == NULL) {
                DSError(M_DS_VAR_NULL ": Xd_a variable pool is NULL", A_DS_ERROR);
                goto bail;
        }
        if (DSVariablePoolNumberOfVariables(DSSSystemXd_a(originalSSystem)) > DSVariablePoolNumberOfVariables(DSSSystemXd(originalSSystem))) {
                DSError(M_DS_WRONG ": Number of algebraic variables exceeds number of total variables", A_DS_ERROR);
                goto bail;
        }
        if (DSVariablePoolNumberOfVariables(DSSSystemXd_a(originalSSystem)) == 0) {
                collapsedSystem = DSSSystemCopy(originalSSystem);
                goto bail;
        }
        oldXd = DSSSystemXd(originalSSystem);
        newXd = DSVariablePoolAlloc();
        for (i = 0; i < DSVariablePoolNumberOfVariables(oldXd); i++) {
                name = DSVariableName(DSVariablePoolVariableAtIndex(oldXd, i));
                if (DSVariablePoolHasVariableWithName(DSSSystemXd_a(originalSSystem), name) == true) {
                        continue;
                }
                DSVariablePoolAddVariableWithName(newXd, name);
        }
        if (DSVariablePoolNumberOfVariables(oldXd) - DSVariablePoolNumberOfVariables(newXd) != DSVariablePoolNumberOfVariables(DSSSystemXd_a(originalSSystem))) {
                DSError(M_DS_WRONG, A_DS_ERROR);
                DSVariablePoolFree(newXd);
                goto bail;
        }
        numberOfDifferentialVariables = DSVariablePoolNumberOfVariables(newXd);
        numberOfAlgebraicVaiables = DSVariablePoolNumberOfVariables(DSSSystemXd_a(originalSSystem));
        differentialIndices = DSSecureCalloc(sizeof(DSUInteger), numberOfDifferentialVariables);
        algebraicIndices = DSSecureCalloc(sizeof(DSUInteger), numberOfAlgebraicVaiables);
        for (i = 0, j = 0, k = 0; i < DSVariablePoolNumberOfVariables(oldXd); i++) {
                if (DSVariablePoolHasVariableWithName(newXd, DSVariableName(DSVariablePoolVariableAtIndex(oldXd, i))) == true) {
                        differentialIndices[j++] = i;
                } else {
                        algebraicIndices[k++] = i;
                }
        }
        collapsedSystem = dsSSystemWithAlgebraicConstraints(originalSSystem, newXd, numberOfDifferentialVariables, differentialIndices, numberOfAlgebraicVaiables, algebraicIndices);
        DSSecureFree(differentialIndices);
        DSSecureFree(algebraicIndices);
bail:
        return collapsedSystem;
}

#if defined (__APPLE__) && defined (__MACH__)
#pragma mark Public Functions
#endif


extern DSSSystem * DSSSystemByParsingStringList(char * const *  string, const DSVariablePool * const Xd_a, ...)
{
        DSSSystem *gma = NULL;
        DSUInteger numberOfStrings = 0;
        char const ** strings = NULL;
        const char * aString = NULL;
        if (string == NULL) {
                DSError(M_DS_NULL ": String to parse is NULL", A_DS_ERROR);
        }
        va_list ap;
	va_start(ap, Xd_a);
        strings = DSSecureCalloc(sizeof(char *), 1);
        strings[0] = (const char *)string;
        numberOfStrings++;
        aString = va_arg(ap, char *);
        while (aString != NULL) {
                strings = DSSecureRealloc(strings, sizeof(char *)*(numberOfStrings+1));
                strings[numberOfStrings++] = aString;
                aString = va_arg(ap, char *);
        }
        gma = DSSSystemByParsingStrings((char * const * )strings, Xd_a, numberOfStrings);
        DSSecureFree(strings);
bail:
        return gma;
}

extern DSSSystem * DSSSystemByParsingStrings(char * const * const strings, const DSVariablePool * const Xd_a, const DSUInteger numberOfEquations)
{
        DSSSystem * sys = NULL;
        gma_parseraux_t **aux = NULL;
        DSUInteger i, j;
        DSVariablePool *tempPool, * Xd, * Xda, * Xdt;
        char * variableName;
        DSExpression * expr, * lhs;
        if (strings == NULL) {
                DSError(M_DS_NULL ": Array of strings is NULL", A_DS_ERROR);
                goto bail;
        }
        if (numberOfEquations == 0) {
                DSError(M_DS_WRONG ": No equations to parse", A_DS_WARN);
                goto bail;
        }
        aux = dsSSysTermListForAllStrings(strings, numberOfEquations);
        if (aux == NULL)
                goto bail;
        Xd = DSVariablePoolAlloc();
        Xda = DSVariablePoolAlloc();
        Xdt = DSVariablePoolAlloc();
        for (i=0; i < numberOfEquations; i++) {
                expr = DSExpressionByParsingString(strings[i]);
                lhs = DSExpressionEquationLHSExpression(expr);
                if (DSExpressionType(lhs) == DS_EXPRESSION_TYPE_CONSTANT) {
                        // If different from 0, should substract rhs by lhs
                }
                tempPool = DSExpressionVariablesInExpression(lhs);
                if (DSVariablePoolNumberOfVariables(tempPool) == 1) {
                        variableName = DSVariableName(DSVariablePoolVariableAtIndex(tempPool, 0));
                        if (DSExpressionType(lhs) == DS_EXPRESSION_TYPE_VARIABLE) {
                                DSVariablePoolAddVariableWithName(Xda, variableName);
                        } else {
                                DSVariablePoolAddVariableWithName(Xdt, variableName);
                        }
                        if (DSVariablePoolHasVariableWithName(Xd, variableName) == false) {
                                DSVariablePoolAddVariableWithName(Xd, variableName);
                        }
                }
                DSExpressionFree(lhs);
                DSExpressionFree(expr);
                DSVariablePoolFree(tempPool);
        }
        if (Xd_a != NULL) {
                for (j = 0; j < DSVariablePoolNumberOfVariables(Xd_a); j++) {
                        variableName = DSVariableName(DSVariablePoolVariableAtIndex(Xd_a, j));
                        if (DSVariablePoolHasVariableWithName(Xd, variableName) == false) {
                                DSVariablePoolAddVariableWithName(Xd, variableName);
                                DSVariablePoolAddVariableWithName(Xda, variableName);
                        }
                }
        }
        if (DSVariablePoolNumberOfVariables(Xd) != numberOfEquations) {
                DSError(M_DS_WRONG ": Number of dependent variables does not match number of equations", A_DS_ERROR);
                goto bail;
        }
        
        sys = DSSSystemAlloc();
        DSSSysXd(sys) = Xd;
        DSVariablePoolSetReadWrite(DSSSysXd(sys));
        DSSSysXd_a(sys) = Xda;
        DSSSysXd_t(sys) = Xdt;
        DSVariablePoolSetReadWrite(DSSSysXd_a(sys));
        DSVariablePoolSetReadWrite(DSSSysXd_t(sys));
        DSSSysXi(sys) = dsSSystemIdentifyIndependentVariables(Xd, aux, numberOfEquations);
        DSSSystemSetShouldFreeXd(sys, true);
        DSSSystemSetShouldFreeXi(sys, true);
        DSVariablePoolSetReadWrite(DSSSysXi(sys));
        dsSSystemCreateSystemMatrices(sys, aux);
        for (i=0; i < numberOfEquations; i++)
                if (aux[i] != NULL)
                        DSGMAParserAuxFree(aux[i]);
        DSSecureFree(aux);
bail:
        return sys;
}

extern DSSSystem * DSSSystemFromGMAWithDominantTerms(const DSGMASystem * gma, const DSUInteger * termArray)
{
        return DSSSystemWithTermsFromGMA(gma, termArray);
}

extern DSSSystem * DSSSystemWithTermsFromGMA(const DSGMASystem * gma, const DSUInteger * termArray)
{
        DSSSystem *ssys = NULL;
        DSUInteger i, j, term1, term2, numberOfEquations, numberOfXi;
        if (gma == NULL) {
                DSError(M_DS_NULL ": Template GMA to make S-System is NULL", A_DS_ERROR);
                goto bail;
        }
        if (termArray == NULL) {
                DSError(M_DS_NULL ": Array of dominant terms is NULL", A_DS_ERROR);
                goto bail;
        }
        ssys = DSSSystemAlloc();
        DSSSysXd(ssys) = (DSVariablePool *)DSGMASystemXd(gma);
        DSSSysXi(ssys) = (DSVariablePool *)DSGMASystemXi(gma);
        DSSSysXd_a(ssys) = (DSVariablePool *)DSGMASystemXd_a(gma);
        DSSSysXd_t(ssys) = (DSVariablePool *)DSGMASystemXd_t(gma);
        DSSSystemSetShouldFreeXd(ssys, false);
        DSSSystemSetShouldFreeXi(ssys, false);
        dsSSystemInitializeMatrices(ssys);
        numberOfEquations = DSGMASystemNumberOfEquations(gma);
        numberOfXi = DSVariablePoolNumberOfVariables(DSSSysXi(ssys));
        for (i = 0; i < 2*numberOfEquations; i+=2) {
                term1 = termArray[i];
                term2 = termArray[i+1];
                if (term1 > DSGMASystemSignature(gma)[i] || term2 > DSGMASystemSignature(gma)[i+1])
                        break;
                DSMatrixSetDoubleValue(DSSSysAlpha(ssys), i/2, 0, 
                                       DSMatrixDoubleValue(DSGMASystemAlpha(gma), i/2, term1-1));
                DSMatrixSetDoubleValue(DSSSysBeta(ssys), i/2, 0, 
                                       DSMatrixDoubleValue(DSGMASystemBeta(gma), i/2, term2-1));
                for (j = 0; j < numberOfEquations; j++) {
                        DSMatrixSetDoubleValue(DSSSysGd(ssys), i/2, j, 
                                               DSMatrixArrayDoubleWithIndices(DSGMASystemGd(gma), i/2, term1-1, j));
                        DSMatrixSetDoubleValue(DSSSysHd(ssys), i/2, j, 
                                               DSMatrixArrayDoubleWithIndices(DSGMASystemHd(gma), i/2, term2-1, j));
                }
                for (j = 0; j < numberOfXi; j++) {
                        DSMatrixSetDoubleValue(DSSSysGi(ssys), i/2, j, 
                                               DSMatrixArrayDoubleWithIndices(DSGMASystemGi(gma), i/2, term1-1, j));
                        DSMatrixSetDoubleValue(DSSSysHi(ssys), i/2, j, 
                                               DSMatrixArrayDoubleWithIndices(DSGMASystemHi(gma), i/2, term2-1, j));
                }
        }
        if (i == 2*numberOfEquations) {
                dsSSystemSolveEquations(ssys);
        } else {
                DSSSystemFree(ssys);
                ssys = NULL;
        }
bail:
        return ssys;
}


static void dsSSystemAssignMatricesCyclicalCase(DSSSystem * ssys,
                                                const DSSSystem * originalSsys,
                                                DSGMASystem *gma, DSUInteger term1,
                                                DSUInteger term2, double factor, DSUInteger i,
                                                DSUInteger dependent_pool_index, bool set_special)

{
    
        DSUInteger j;
        DSUInteger numberOfEquations = DSGMASystemNumberOfEquations(gma);
        DSUInteger numberOfXi = DSVariablePoolNumberOfVariables(DSSSysXi(ssys));
    
    
    if (set_special == true){
    
        // If the outlet reaction does not correspond to the main variable and it is secondary (or not)
        // First we assign the original Alphas and Betas to the mainCycleVariable
    
        DSMatrixSetDoubleValue(DSSSysAlpha(ssys), i/2, 0,
                               DSMatrixDoubleValue(DSSSystemAlpha(originalSsys), i/2, 0));
        DSMatrixSetDoubleValue(DSSSysBeta(ssys), i/2, 0,
                               DSMatrixDoubleValue(DSSSystemBeta(originalSsys), i/2, 0));
    
        // Then we assign values for Alpha und Betas from the first equation to the equation defined by term2
        DSMatrixSetDoubleValue(DSSSysAlpha(ssys),dependent_pool_index, 0,
                               DSMatrixDoubleValue(DSGMASystemAlpha(gma), i/2, term1-1)*factor);
        DSMatrixSetDoubleValue(DSSSysBeta(ssys), dependent_pool_index, 0,
                               DSMatrixDoubleValue(DSGMASystemBeta(gma), i/2, term2-1)*factor);
        // Now assign Gd and Hd. Let's first assign values to the equation of the mainCycleVariable (from original S-System)
        for (j = 0; j < numberOfEquations; j++) {
            DSMatrixSetDoubleValue(DSSSysGd(ssys), i/2, j,
                                   DSMatrixDoubleValue(DSSSystemGd(originalSsys), i/2, j));
            DSMatrixSetDoubleValue(DSSSysHd(ssys), i/2, j,
                                   DSMatrixDoubleValue(DSSSystemHd(originalSsys),i/2, j));
            // And now let's assign the dominat term of the maincyclevariable equation to the corresponding outlet reaction
            DSMatrixSetDoubleValue(DSSSysGd(ssys),dependent_pool_index, j,
                                   DSMatrixArrayDoubleWithIndices(DSGMASystemGd(gma), i/2, term1-1, j));
            DSMatrixSetDoubleValue(DSSSysHd(ssys),dependent_pool_index, j,
                                   DSMatrixArrayDoubleWithIndices(DSGMASystemHd(gma), i/2, term2-1, j));
        }
        // Now assign Gi and Hi. Let's first assign values to the equation of the mainCycleVariable (from original S-System)
        for (j = 0; j < numberOfXi; j++) {
            DSMatrixSetDoubleValue(DSSSysGi(ssys), i/2, j,
                                   DSMatrixDoubleValue(DSSSystemGi(originalSsys), i/2, j));
            DSMatrixSetDoubleValue(DSSSysHi(ssys), i/2, j,
                                   DSMatrixDoubleValue(DSSSystemHi(originalSsys), i/2, j));
            
            // And now let's assign the dominat term of the maincyclevariable equation to the corresponding outlet reaction
            DSMatrixSetDoubleValue(DSSSysGi(ssys), dependent_pool_index, j,
                                   DSMatrixArrayDoubleWithIndices(DSGMASystemGi(gma), i/2, term1-1, j));
            DSMatrixSetDoubleValue(DSSSysHi(ssys),dependent_pool_index, j,
                                   DSMatrixArrayDoubleWithIndices(DSGMASystemHi(gma), i/2, term2-1, j));
        }
    
    } else {
        
        DSMatrixSetDoubleValue(DSSSysAlpha(ssys), i/2, 0,
                               DSMatrixDoubleValue(DSGMASystemAlpha(gma), i/2, term1-1));
        DSMatrixSetDoubleValue(DSSSysBeta(ssys), i/2, 0,
                               DSMatrixDoubleValue(DSGMASystemBeta(gma), i/2, term2-1));
        
        for (j = 0; j < numberOfEquations; j++) {
            DSMatrixSetDoubleValue(DSSSysGd(ssys), i/2, j,
                                   DSMatrixArrayDoubleWithIndices(DSGMASystemGd(gma), i/2, term1-1, j));
            DSMatrixSetDoubleValue(DSSSysHd(ssys), i/2, j,
                                   DSMatrixArrayDoubleWithIndices(DSGMASystemHd(gma), i/2, term2-1, j));
        }
        for (j = 0; j < numberOfXi; j++) {
            DSMatrixSetDoubleValue(DSSSysGi(ssys), i/2, j,
                                   DSMatrixArrayDoubleWithIndices(DSGMASystemGi(gma), i/2, term1-1, j));
            DSMatrixSetDoubleValue(DSSSysHi(ssys), i/2, j,
                                   DSMatrixArrayDoubleWithIndices(DSGMASystemHi(gma), i/2, term2-1, j));
        }
        
    }
}


DSSSystem * DSSSystemWithTermsFromGMACyclical(const DSDesignSpace * ds,
                                              const DSUInteger * termArray,
                                              DSuIntegerMatrix * three_digit,
                                              DSuIntegerMatrix *location)

{
    
    // unpack variables gma, mainCycleVariables, originalSsys, numberOfCycles
    DSUInteger *alreadyassigned = NULL, secondary_variable;
    DSUInteger dependent_pool_index;
    DSGMASystem * gma = ds->gma ;
    const DSUInteger * mainCycleVariables = ds->extensionData->mainCycleVariables, * secondary_variables;
    DSUInteger numberOfCycles = ds->extensionData->numberCycles;
    const DSSSystem *originalSsys = ds->extensionData->originalsSystem;
    DSSSystem *ssys = NULL;
    DSUInteger i, w, v, term1, term2, numberOfEquations, numberOfXi, ind_term, ind_eq, d1, d2;
    DSUInteger currentCycle = 0 ;
    const DSMatrix * beta_original = ds->extensionData->beta;
    double factor = 1.0, den = 1.0, num = 1.0;
    bool isMainCyclicalVariable = false, isSecondaryVariable = false, isMainVariable = false;
    bool alreadyassigned_bol = false;
    DSDesignSpace *previous, *original;
    
    if (gma == NULL) {
        DSError(M_DS_NULL ": Template GMA to make S-System is NULL", A_DS_ERROR);
        goto bail;
    }
    if (termArray == NULL) {
        DSError(M_DS_NULL ": Array of dominant terms is NULL", A_DS_ERROR);
        goto bail;
    }
    
//
//
//    // delete this after debugging
//    bool de_bug = false;
//    char aux[100];
//
//    sprintf(aux, "%s", "93_15");
//    if (ds->casePrefix != NULL){
//        if (strcmp(ds->casePrefix, aux) == 0){
//            if (DSCaseNumberForSignature(termArray, ds->gma) == 5){
//                sprintf(aux, "%s", "93_15_5");
//                de_bug = true;
//            }
//        }
//    }
//
//    if (de_bug == true){
//        previous = ds->extensionData->parent_ds;
//        original = previous->extensionData->parent_ds;
//
//        printf("****(DSSSystemWithTermsFromGMACyclical) The equations of the original ds are: \n");
//        for (i=0; i<5; i++){
//            printf("Equation %u: %s \n",i, DSExpressionAsString(DSGMASystemEquations(original->gma)[i]));
//        }
//        printf("The corresponding beta matrix is: \n");
//        DSMatrixPrint(original->gma->beta);
//
//        printf("The equations of the previous ds (93) ds are: ---------\n");
//        for (i=0; i<5; i++){
//            printf("Equation %u: %s \n",i, DSExpressionAsString(DSGMASystemEquations(previous->gma)[i]));
//        }
//        printf("The corresponding beta matrix is: \n");
//        DSMatrixPrint(previous->gma->beta);
//        printf("The correspinding H_L_term matrix is: \n");
//        DSuIntegerMatrixPrint(previous->extensionData->H_l_term);
//        printf("The corresponding H_L_eq matrix is: \n");
//        DSuIntegerMatrixPrint(previous->extensionData->H_l_eq);
//        printf("The main cyclical variables (%u) are: %u \n ", previous->extensionData->numberCycles,
//               previous->extensionData->mainCycleVariables[0]);
//        printf("The secondary cyclical variables (%u) are: %u \n", previous->extensionData->numberSecondaryVariables[0], previous->extensionData->allSecondaryVariables[0][0]);
//        printf("The original ssystem is: \n");
//        DSSSystemPrintEquations(previous->extensionData->originalsSystem);
//
//        printf("The equations of the current ds (93_15) are: ----------\n");
//        for (i=0; i<5; i++){
//            printf("Equation %u: %s \n",i, DSExpressionAsString(DSGMASystemEquations(ds->gma)[i]));
//        }
//        printf("The corresponding beta matrix is: \n");
//        DSMatrixPrint(ds->gma->beta);
//        printf("The correspinding H_L_term matrix is: \n");
//        DSuIntegerMatrixPrint(ds->extensionData->H_l_term);
//        printf("The corresponding H_L_eq matrix is: \n");
//        DSuIntegerMatrixPrint(ds->extensionData->H_l_eq);
//        printf("The main cyclical variables (%u) are: %u \n ", ds->extensionData->numberCycles,
//               ds->extensionData->mainCycleVariables[0]);
//        printf("The secondary cyclical variables (%u) are: %u \n", ds->extensionData->numberSecondaryVariables[0], ds->extensionData->allSecondaryVariables[0][0]);
//        printf("The original ssystem is: \n");
//        DSSSystemPrintEquations(ds->extensionData->originalsSystem);
//
//    }
//
    ssys = DSSSystemAlloc();
    DSSSysXd(ssys) = (DSVariablePool *)DSGMASystemXd(gma);
    DSSSysXi(ssys) = (DSVariablePool *)DSGMASystemXi(gma);
    DSSSysXd_a(ssys) = (DSVariablePool *)DSGMASystemXd_a(gma);
    DSSSysXd_t(ssys) = (DSVariablePool *)DSGMASystemXd_t(gma);
    DSSSystemSetShouldFreeXd(ssys, false);
    DSSSystemSetShouldFreeXi(ssys, false);
    dsSSystemInitializeMatrices(ssys);
    numberOfEquations = DSGMASystemNumberOfEquations(gma);
    numberOfXi = DSVariablePoolNumberOfVariables(DSSSysXi(ssys));
    
    if (numberOfCycles != 0)
        alreadyassigned = DSSecureMalloc(numberOfCycles*sizeof(DSUInteger));
    for(w = 0; w < numberOfCycles; w++)
        alreadyassigned[w] = numberOfEquations;
    
    
    for (i = 0; i < 2*numberOfEquations; i+=2) {
        term1 = termArray[i];
        term2 = termArray[i+1];
        if (term1 > DSGMASystemSignature(gma)[i] || term2 > DSGMASystemSignature(gma)[i+1])
            break;
        
        // define if equation with index i/2 is main or not.
        for(w = 0; w < numberOfCycles; w++)
            if( i/2 == mainCycleVariables[w])
                isMainCyclicalVariable = true;
        
        // If constructing equation for a main cyclical variable
        if (isMainCyclicalVariable == true){
            
//            if (de_bug == true)
//                printf("ismaincyclicalvariable is true for case %s \n", aux);
            
                    // define if outlet reaction (dependent_pool_index) is a secondary variable or main variable
                    isSecondaryVariable = false;
                    isMainVariable = false;
                    dependent_pool_index = DSuIntegerMatrixValue(ds->extensionData->H_l_eq, currentCycle, term2-1);
                    for(w = 0; w <numberOfCycles; w++){
                        if( dependent_pool_index == mainCycleVariables[w])
                            isMainVariable = true;
                        for (v = 0; v < ds->extensionData->numberSecondaryVariables[w]; v++){
                            secondary_variables = ds->extensionData->allSecondaryVariables[w];
                            secondary_variable = secondary_variables[v];
                            if( dependent_pool_index == secondary_variable)
                                isSecondaryVariable = true;
                        }
                    }
            
                    ind_eq = DSuIntegerMatrixValue(ds->extensionData->H_l_eq, currentCycle, term2-1);
                    ind_term = DSuIntegerMatrixValue(ds->extensionData->H_l_term, currentCycle, term2-1);
                    DSuIntegerMatrixSetValue(location, currentCycle, 0, ind_eq);
                    d1 = DSuIntegerMatrixValue(ds->extensionData->G_l_eq, currentCycle, term1-1);
                    DSuIntegerMatrixSetValue(three_digit, currentCycle, 0, d1 );
                    d2 = DSuIntegerMatrixValue(ds->extensionData->G_l_term, currentCycle, term1-1);
                    DSuIntegerMatrixSetValue(three_digit, currentCycle, 1, d2);
                    DSuIntegerMatrixSetValue(three_digit, currentCycle, 2, ind_term);
                    den = DSMatrixDoubleValue(beta_original, ind_eq, ind_term);
                    num = DSMatrixDoubleValue(DSGMASystemBeta(gma), i/2, term2-1);
                    factor = den/num;
//                    if (de_bug == true)
//                        printf("factor equals %f for case %s \n", factor, aux);
            
                    if ((isMainVariable == false && isSecondaryVariable == true ) ||
                        (isMainVariable == false && isSecondaryVariable == false)){
                        
                            // If the outlet reaction does not correspond to the main variable and it is secondary.
                            // First we assign the original Alphas and Betas to the mainCycleVariable
                        
//                            if (de_bug == true)
//                                printf("main variable is secondary for case %s \n", aux);
                        
                            dsSSystemAssignMatricesCyclicalCase(ssys,
                                                                originalSsys,
                                                                gma,term1,
                                                                term2, factor, i,
                                                                dependent_pool_index, true);
                        
                            alreadyassigned[currentCycle] = dependent_pool_index;
                        
                    }else{
                            //do the normal assignment for the main cyclical variable as output
//                            if (de_bug == true)
//                                printf("main variable is output for case %s \n", aux);
                        
                            dsSSystemAssignMatricesCyclicalCase(ssys,
                                                                originalSsys,
                                                                gma,term1,
                                                                term2, factor, i,
                                                                dependent_pool_index, false);
                    }
                    currentCycle++;
        }else{
            
//                    if (de_bug == true)
//                        printf("ismaincyclicalvariable is false for case %s \n", aux);
            
                    // The normal assignment, when we are NOT dealing with a main variable.
                    for(w = 0; w < numberOfCycles; w++)
                        if( i/2 == alreadyassigned[w])
                            alreadyassigned_bol = true;
            
                    if (alreadyassigned_bol == false)
                        dsSSystemAssignMatricesCyclicalCase(ssys,
                                                            originalSsys,
                                                            gma,term1,
                                                            term2, factor, i,
                                                            dependent_pool_index, false);
                    alreadyassigned_bol = false;
        }
        
        isMainCyclicalVariable = false;
        isSecondaryVariable = false;
    }
    
    if (i == 2*numberOfEquations) {
        dsSSystemSolveEquations(ssys);
    } else {
        DSSSystemFree(ssys);
        ssys = NULL;
    }
    
bail:
    if (alreadyassigned != NULL)
        DSSecureFree(alreadyassigned);
    
//    if (de_bug ==true)
//        printf("****(DSSSystemWithTermsFromGMACyclical) End of the function for case %s\n", aux);
    
    return ssys;
}

DSSSystem * DSSSystemWithTermsFromGMACyclical_beforefactoring(const DSDesignSpace * ds,
                                              const DSUInteger * termArray,
                                              DSuIntegerMatrix * three_digit,
                                              DSuIntegerMatrix *location)

{
    
    // unpack variables gma, mainCycleVariables, originalSsys, numberOfCycles
    DSUInteger *alreadyassigned = NULL, secondary_variable;
    DSUInteger dependent_pool_index;
    DSGMASystem * gma = ds->gma ;
    const DSUInteger * mainCycleVariables = ds->extensionData->mainCycleVariables, * secondary_variables;
    DSUInteger numberOfCycles = ds->extensionData->numberCycles;
    const DSSSystem *originalSsys = ds->extensionData->originalsSystem;
    DSSSystem *ssys = NULL;
    DSUInteger i, j, w, v, term1, term2, numberOfEquations, numberOfXi, ind_term, ind_eq, d1, d2;
    DSUInteger currentCycle = 0 ;
    const DSMatrix * beta_original = ds->extensionData->beta;
    double factor = 1.0, den = 1.0, num = 1.0;
    bool isMainCyclicalVariable = false, isSecondaryVariable = false, isMainVariable = false;
    bool alreadyassigned_bol = false;

    if (gma == NULL) {
        DSError(M_DS_NULL ": Template GMA to make S-System is NULL", A_DS_ERROR);
        goto bail;
    }
    if (termArray == NULL) {
        DSError(M_DS_NULL ": Array of dominant terms is NULL", A_DS_ERROR);
        goto bail;
    }
    
    ssys = DSSSystemAlloc();
    DSSSysXd(ssys) = (DSVariablePool *)DSGMASystemXd(gma);
    DSSSysXi(ssys) = (DSVariablePool *)DSGMASystemXi(gma);
    DSSSysXd_a(ssys) = (DSVariablePool *)DSGMASystemXd_a(gma);
    DSSSysXd_t(ssys) = (DSVariablePool *)DSGMASystemXd_t(gma);
    DSSSystemSetShouldFreeXd(ssys, false);
    DSSSystemSetShouldFreeXi(ssys, false);
    dsSSystemInitializeMatrices(ssys);
    numberOfEquations = DSGMASystemNumberOfEquations(gma);
    numberOfXi = DSVariablePoolNumberOfVariables(DSSSysXi(ssys));
    
    alreadyassigned = DSSecureMalloc(numberOfCycles*sizeof(DSUInteger));
    for(w = 0; w < numberOfCycles; w++)
            alreadyassigned[w] = numberOfEquations;
    

    for (i = 0; i < 2*numberOfEquations; i+=2) {
            term1 = termArray[i];
            term2 = termArray[i+1];
            if (term1 > DSGMASystemSignature(gma)[i] || term2 > DSGMASystemSignature(gma)[i+1])
                break;
        
            // define if equation with index i/2 is main or not.
            for(w = 0; w < numberOfCycles; w++)
                    if( i/2 == mainCycleVariables[w])
                        isMainCyclicalVariable = true;
                
            // If constructing equation for a main cyclical variable
            if (isMainCyclicalVariable == true){
                
                // define if outlet reaction (dependent_pool_index) is a secondary variable or main variable
                isSecondaryVariable = false;
                isMainVariable = false;
                dependent_pool_index = DSuIntegerMatrixValue(ds->extensionData->H_l_eq, currentCycle, term2-1);
                for(w = 0; w <numberOfCycles; w++){
                    if( dependent_pool_index == mainCycleVariables[w])
                        isMainVariable = true;
                    for (v = 0; v < ds->extensionData->numberSecondaryVariables[w]; v++){
                        secondary_variables = ds->extensionData->allSecondaryVariables[w];
                        secondary_variable = secondary_variables[v];
                        if( dependent_pool_index == secondary_variable)
                            isSecondaryVariable = true;
                    }
                }
                
                ind_eq = DSuIntegerMatrixValue(ds->extensionData->H_l_eq, currentCycle, term2-1);
                ind_term = DSuIntegerMatrixValue(ds->extensionData->H_l_term, currentCycle, term2-1);
                DSuIntegerMatrixSetValue(location, currentCycle, 0, ind_eq);
                d1 = DSuIntegerMatrixValue(ds->extensionData->G_l_eq, currentCycle, term1-1);
                DSuIntegerMatrixSetValue(three_digit, currentCycle, 0, d1 );
                d2 = DSuIntegerMatrixValue(ds->extensionData->G_l_term, currentCycle, term1-1);
                DSuIntegerMatrixSetValue(three_digit, currentCycle, 1, d2);
                DSuIntegerMatrixSetValue(three_digit, currentCycle, 2, ind_term);
                den = DSMatrixDoubleValue(beta_original, ind_eq, ind_term);
                num = DSMatrixDoubleValue(DSGMASystemBeta(gma), i/2, term2-1);
                factor = den/num;


                if ((isMainVariable == false && isSecondaryVariable == true) ||
                    (isMainVariable == false && isSecondaryVariable == false )){
                    
                    // If the outlet reaction does not correspond to the main variable and it is secondary (or not)
                    // First we assign the original Alphas and Betas to the mainCycleVariable
                    
                    DSMatrixSetDoubleValue(DSSSysAlpha(ssys), i/2, 0,
                                           DSMatrixDoubleValue(DSSSystemAlpha(originalSsys), i/2, 0));
                    DSMatrixSetDoubleValue(DSSSysBeta(ssys), i/2, 0,
                                           DSMatrixDoubleValue(DSSSystemBeta(originalSsys), i/2, 0));
                    
                    // Then we assign values for Alpha und Betas from the first equation to the equation defined by term2
                    DSMatrixSetDoubleValue(DSSSysAlpha(ssys),dependent_pool_index, 0,
                                           DSMatrixDoubleValue(DSGMASystemAlpha(gma), i/2, term1-1)*factor);
                    DSMatrixSetDoubleValue(DSSSysBeta(ssys), dependent_pool_index, 0,
                                           DSMatrixDoubleValue(DSGMASystemBeta(gma), i/2, term2-1)*factor);
                    // Now assign Gd and Hd. Let's first assign values to the equation of the mainCycleVariable (from original S-System)
                    for (j = 0; j < numberOfEquations; j++) {
                        DSMatrixSetDoubleValue(DSSSysGd(ssys), i/2, j,
                                               DSMatrixDoubleValue(DSSSystemGd(originalSsys), i/2, j));
                        DSMatrixSetDoubleValue(DSSSysHd(ssys), i/2, j,
                                               DSMatrixDoubleValue(DSSSystemHd(originalSsys),i/2, j));
                        // And now let's assign the dominat term of the maincyclevariable equation to the corresponding outlet reaction
                        DSMatrixSetDoubleValue(DSSSysGd(ssys),dependent_pool_index, j,
                                               DSMatrixArrayDoubleWithIndices(DSGMASystemGd(gma), i/2, term1-1, j));
                        DSMatrixSetDoubleValue(DSSSysHd(ssys),dependent_pool_index, j,
                                               DSMatrixArrayDoubleWithIndices(DSGMASystemHd(gma), i/2, term2-1, j));
                    }
                    // Now assign Gi and Hi. Let's first assign values to the equation of the mainCycleVariable (from original S-System)
                    for (j = 0; j < numberOfXi; j++) {
                        DSMatrixSetDoubleValue(DSSSysGi(ssys), i/2, j,
                                               DSMatrixDoubleValue(DSSSystemGi(originalSsys), i/2, j));
                        DSMatrixSetDoubleValue(DSSSysHi(ssys), i/2, j,
                                               DSMatrixDoubleValue(DSSSystemHi(originalSsys), i/2, j));
                        
                    // And now let's assign the dominat term of the maincyclevariable equation to the corresponding outlet reaction
                        DSMatrixSetDoubleValue(DSSSysGi(ssys), dependent_pool_index, j,
                                               DSMatrixArrayDoubleWithIndices(DSGMASystemGi(gma), i/2, term1-1, j));
                        DSMatrixSetDoubleValue(DSSSysHi(ssys),dependent_pool_index, j,
                                               DSMatrixArrayDoubleWithIndices(DSGMASystemHi(gma), i/2, term2-1, j));
                    }
                    
                    alreadyassigned[currentCycle] = dependent_pool_index;
                    
                }else{ //do the normal assignment for the main cyclical variable as output
                    
                    DSMatrixSetDoubleValue(DSSSysAlpha(ssys), i/2, 0,
                                           DSMatrixDoubleValue(DSGMASystemAlpha(gma), i/2, term1-1));
                    DSMatrixSetDoubleValue(DSSSysBeta(ssys), i/2, 0,
                                           DSMatrixDoubleValue(DSGMASystemBeta(gma), i/2, term2-1));
                    
                    for (j = 0; j < numberOfEquations; j++) {
                        DSMatrixSetDoubleValue(DSSSysGd(ssys), i/2, j,
                                               DSMatrixArrayDoubleWithIndices(DSGMASystemGd(gma), i/2, term1-1, j));
                        DSMatrixSetDoubleValue(DSSSysHd(ssys), i/2, j,
                                               DSMatrixArrayDoubleWithIndices(DSGMASystemHd(gma), i/2, term2-1, j));
                    }
                    for (j = 0; j < numberOfXi; j++) {
                        DSMatrixSetDoubleValue(DSSSysGi(ssys), i/2, j,
                                               DSMatrixArrayDoubleWithIndices(DSGMASystemGi(gma), i/2, term1-1, j));
                        DSMatrixSetDoubleValue(DSSSysHi(ssys), i/2, j,
                                               DSMatrixArrayDoubleWithIndices(DSGMASystemHi(gma), i/2, term2-1, j));
                    }
                }
            currentCycle++;
                
            }else{ // The normal assignment, when we are NOT dealing with a main variable.
                
                for(w = 0; w < numberOfCycles; w++)
                    if( i/2 == alreadyassigned[w])
                        alreadyassigned_bol = true;
                
                if (alreadyassigned_bol == false){
                    DSMatrixSetDoubleValue(DSSSysAlpha(ssys), i/2, 0,
                                           DSMatrixDoubleValue(DSGMASystemAlpha(gma), i/2, term1-1));
                    DSMatrixSetDoubleValue(DSSSysBeta(ssys), i/2, 0,
                                           DSMatrixDoubleValue(DSGMASystemBeta(gma), i/2, term2-1));
                    for (j = 0; j < numberOfEquations; j++) {
                        DSMatrixSetDoubleValue(DSSSysGd(ssys), i/2, j,
                                               DSMatrixArrayDoubleWithIndices(DSGMASystemGd(gma), i/2, term1-1, j));
                        DSMatrixSetDoubleValue(DSSSysHd(ssys), i/2, j,
                                               DSMatrixArrayDoubleWithIndices(DSGMASystemHd(gma), i/2, term2-1, j));
                    }
                    for (j = 0; j < numberOfXi; j++) {
                        DSMatrixSetDoubleValue(DSSSysGi(ssys), i/2, j,
                                               DSMatrixArrayDoubleWithIndices(DSGMASystemGi(gma), i/2, term1-1, j));
                        DSMatrixSetDoubleValue(DSSSysHi(ssys), i/2, j,
                                               DSMatrixArrayDoubleWithIndices(DSGMASystemHi(gma), i/2, term2-1, j));
                    }
                }
                alreadyassigned_bol = false;
            }
        
            isMainCyclicalVariable = false;
            isSecondaryVariable = false;
    }
    
    if (i == 2*numberOfEquations) {
            dsSSystemSolveEquations(ssys);
    } else {
            DSSSystemFree(ssys);
            ssys = NULL;
    }
    
bail:
    if (alreadyassigned != NULL)
        DSSecureFree(alreadyassigned);
    return ssys;
}

static DSSSystem * DSSSystemWithTermsFromGMACyclical_old(const DSDesignSpace * ds,
                                                         const DSUInteger * termArray,
                                                         DSuIntegerMatrix * three_digit,
                                                         DSuIntegerMatrix *location)

{
    // unpack variables gma, mainCycleVariables, originalSsys, numberOfCycles
    DSUInteger alreadyassigned, secondary_variable;
    DSUInteger dependent_pool_index;
    DSGMASystem * gma = ds->gma ;
    const DSUInteger * mainCycleVariables = ds->extensionData->mainCycleVariables, * secondary_variables;
    const DSUInteger numberOfCycles = ds->extensionData->numberCycles;
    const DSSSystem *originalSsys = ds->extensionData->originalsSystem;
    DSSSystem *ssys = NULL;
    DSUInteger i, j, w, v, term1, term2, numberOfEquations, numberOfXi, ind_term, ind_eq, d1, d2;
    DSUInteger currentCycle = 0 ;
    const DSMatrix * beta_original = ds->extensionData->beta;
    double factor = 1.0, den = 1.0, num = 1.0;
    bool isMainCyclicalVariable = false, isSecondaryVariable = false;
    
    if (gma == NULL) {
        DSError(M_DS_NULL ": Template GMA to make S-System is NULL", A_DS_ERROR);
        goto bail;
    }
    if (termArray == NULL) {
        DSError(M_DS_NULL ": Array of dominant terms is NULL", A_DS_ERROR);
        goto bail;
    }
    ssys = DSSSystemAlloc();
    DSSSysXd(ssys) = (DSVariablePool *)DSGMASystemXd(gma);
    DSSSysXi(ssys) = (DSVariablePool *)DSGMASystemXi(gma);
    DSSSysXd_a(ssys) = (DSVariablePool *)DSGMASystemXd_a(gma);
    DSSSysXd_t(ssys) = (DSVariablePool *)DSGMASystemXd_t(gma);
    DSSSystemSetShouldFreeXd(ssys, false);
    DSSSystemSetShouldFreeXi(ssys, false);
    dsSSystemInitializeMatrices(ssys);
    numberOfEquations = DSGMASystemNumberOfEquations(gma);
    numberOfXi = DSVariablePoolNumberOfVariables(DSSSysXi(ssys));
    
    for (i = 0; i < 2*numberOfEquations; i+=2) {
        term1 = termArray[i];
        term2 = termArray[i+1];
        if (term1 > DSGMASystemSignature(gma)[i] || term2 > DSGMASystemSignature(gma)[i+1])
            break;
        
        // define if equation with index i/2 is main or not.
        for(w = 0; w < numberOfCycles; w++)
            if( i/2 == mainCycleVariables[w])
                isMainCyclicalVariable = true;
        
        
        if (isMainCyclicalVariable == true){          // If constructing equation for a main cyclical variable
            
            dependent_pool_index = DSuIntegerMatrixValue(ds->extensionData->H_l_eq, currentCycle, term2-1);
            
            // define if outlet reaction (dependent_pool_index) is a secondary variable
            for(w = 0; w <numberOfCycles; w++){
                for (v = 0; v < ds->extensionData->numberSecondaryVariables[w]; v++){
                    secondary_variables = ds->extensionData->allSecondaryVariables[w];
                    secondary_variable = secondary_variables[v];
                    if( dependent_pool_index == secondary_variable){
                        isSecondaryVariable = true;
                    }
                }
            }
            
            if (DSMatrixArrayDoubleWithIndices(DSGMASystemHd(gma), i/2, term2-1, i/2) == 0.0 && isSecondaryVariable == true ){
                
                // If the outlet reaction does not correspond to the main variable and it is secondary.
                // First we assign the original Alphas and Betas to the mainCycleVariable
                
                ind_eq = DSuIntegerMatrixValue(ds->extensionData->H_l_eq, currentCycle, term2-1);
                ind_term = DSuIntegerMatrixValue(ds->extensionData->H_l_term, currentCycle, term2-1);
                DSuIntegerMatrixSetValue(location, currentCycle, 0, ind_eq);
                d1 = DSuIntegerMatrixValue(ds->extensionData->G_l_eq, currentCycle, term1-1);
                DSuIntegerMatrixSetValue(three_digit, currentCycle, 0, d1 );
                d2 = DSuIntegerMatrixValue(ds->extensionData->G_l_term, currentCycle, term1-1);
                DSuIntegerMatrixSetValue(three_digit, currentCycle, 1, d2);
                DSuIntegerMatrixSetValue(three_digit, currentCycle, 2, ind_term);
                den = DSMatrixDoubleValue(beta_original, ind_eq, ind_term);
                num = DSMatrixDoubleValue(DSGMASystemBeta(gma), i/2, term2-1);
                factor = den/num;
                
                DSMatrixSetDoubleValue(DSSSysAlpha(ssys), i/2, 0,
                                       DSMatrixDoubleValue(DSSSystemAlpha(originalSsys), i/2, 0));
                DSMatrixSetDoubleValue(DSSSysBeta(ssys), i/2, 0,
                                       DSMatrixDoubleValue(DSSSystemBeta(originalSsys), i/2, 0));
                
                // Then we assign values for Alpha und Betas from the first equation to the equation defined by term2
                DSMatrixSetDoubleValue(DSSSysAlpha(ssys),dependent_pool_index, 0,
                                       DSMatrixDoubleValue(DSGMASystemAlpha(gma), i/2, term1-1)*factor);
                DSMatrixSetDoubleValue(DSSSysBeta(ssys), dependent_pool_index, 0,
                                       DSMatrixDoubleValue(DSGMASystemBeta(gma), i/2, term2-1)*factor);
                // Now assign Gd and Hd. Let's first assign values to the equation of the mainCycleVariable (from original S-System)
                for (j = 0; j < numberOfEquations; j++) {
                    DSMatrixSetDoubleValue(DSSSysGd(ssys), i/2, j,
                                           DSMatrixDoubleValue(DSSSystemGd(originalSsys), i/2, j));
                    DSMatrixSetDoubleValue(DSSSysHd(ssys), i/2, j,
                                           DSMatrixDoubleValue(DSSSystemHd(originalSsys),i/2, j));
                    // And now let's assign the dominat term of the maincyclevariable equation to the corresponding outlet reaction
                    DSMatrixSetDoubleValue(DSSSysGd(ssys),dependent_pool_index, j,
                                           DSMatrixArrayDoubleWithIndices(DSGMASystemGd(gma), i/2, term1-1, j));
                    DSMatrixSetDoubleValue(DSSSysHd(ssys),dependent_pool_index, j,
                                           DSMatrixArrayDoubleWithIndices(DSGMASystemHd(gma), i/2, term2-1, j));
                }
                // Now assign Gi and Hi. Let's first assign values to the equation of the mainCycleVariable (from original S-System)
                for (j = 0; j < numberOfXi; j++) {
                    DSMatrixSetDoubleValue(DSSSysGi(ssys), i/2, j,
                                           DSMatrixDoubleValue(DSSSystemGi(originalSsys), i/2, j));
                    DSMatrixSetDoubleValue(DSSSysHi(ssys), i/2, j,
                                           DSMatrixDoubleValue(DSSSystemHi(originalSsys), i/2, j));
                    
                    // And now let's assign the dominat term of the maincyclevariable equation to the corresponding outlet reaction
                    DSMatrixSetDoubleValue(DSSSysGi(ssys), dependent_pool_index, j,
                                           DSMatrixArrayDoubleWithIndices(DSGMASystemGi(gma), i/2, term1-1, j));
                    DSMatrixSetDoubleValue(DSSSysHi(ssys),dependent_pool_index, j,
                                           DSMatrixArrayDoubleWithIndices(DSGMASystemHi(gma), i/2, term2-1, j));
                }
                alreadyassigned = dependent_pool_index;
            }else{ //do the normal assignment for the main cyclical variable as output
                DSMatrixSetDoubleValue(DSSSysAlpha(ssys), i/2, 0,
                                       DSMatrixDoubleValue(DSGMASystemAlpha(gma), i/2, term1-1));
                DSMatrixSetDoubleValue(DSSSysBeta(ssys), i/2, 0,
                                       DSMatrixDoubleValue(DSGMASystemBeta(gma), i/2, term2-1));
                ind_eq = DSuIntegerMatrixValue(ds->extensionData->H_l_eq, currentCycle, term2-1);
                ind_term = DSuIntegerMatrixValue(ds->extensionData->H_l_term, currentCycle, term2-1);
                DSuIntegerMatrixSetValue(location, currentCycle, 0, ind_eq);
                d1 = DSuIntegerMatrixValue(ds->extensionData->G_l_eq, currentCycle, term1-1);
                DSuIntegerMatrixSetValue(three_digit, currentCycle, 0, d1 );
                d2 = DSuIntegerMatrixValue(ds->extensionData->G_l_term, currentCycle, term1-1);
                DSuIntegerMatrixSetValue(three_digit, currentCycle, 1, d2);
                DSuIntegerMatrixSetValue(three_digit, currentCycle, 2, ind_term);
                
                for (j = 0; j < numberOfEquations; j++) {
                    DSMatrixSetDoubleValue(DSSSysGd(ssys), i/2, j,
                                           DSMatrixArrayDoubleWithIndices(DSGMASystemGd(gma), i/2, term1-1, j));
                    DSMatrixSetDoubleValue(DSSSysHd(ssys), i/2, j,
                                           DSMatrixArrayDoubleWithIndices(DSGMASystemHd(gma), i/2, term2-1, j));
                }
                for (j = 0; j < numberOfXi; j++) {
                    DSMatrixSetDoubleValue(DSSSysGi(ssys), i/2, j,
                                           DSMatrixArrayDoubleWithIndices(DSGMASystemGi(gma), i/2, term1-1, j));
                    DSMatrixSetDoubleValue(DSSSysHi(ssys), i/2, j,
                                           DSMatrixArrayDoubleWithIndices(DSGMASystemHi(gma), i/2, term2-1, j));
                }
            }
            currentCycle++;
        }else{ // The normal assignment, when we are not dealing with a main variable.
            if (i/2 != alreadyassigned){
                DSMatrixSetDoubleValue(DSSSysAlpha(ssys), i/2, 0,
                                       DSMatrixDoubleValue(DSGMASystemAlpha(gma), i/2, term1-1));
                DSMatrixSetDoubleValue(DSSSysBeta(ssys), i/2, 0,
                                       DSMatrixDoubleValue(DSGMASystemBeta(gma), i/2, term2-1));
                for (j = 0; j < numberOfEquations; j++) {
                    DSMatrixSetDoubleValue(DSSSysGd(ssys), i/2, j,
                                           DSMatrixArrayDoubleWithIndices(DSGMASystemGd(gma), i/2, term1-1, j));
                    DSMatrixSetDoubleValue(DSSSysHd(ssys), i/2, j,
                                           DSMatrixArrayDoubleWithIndices(DSGMASystemHd(gma), i/2, term2-1, j));
                }
                for (j = 0; j < numberOfXi; j++) {
                    DSMatrixSetDoubleValue(DSSSysGi(ssys), i/2, j,
                                           DSMatrixArrayDoubleWithIndices(DSGMASystemGi(gma), i/2, term1-1, j));
                    DSMatrixSetDoubleValue(DSSSysHi(ssys), i/2, j,
                                           DSMatrixArrayDoubleWithIndices(DSGMASystemHi(gma), i/2, term2-1, j));
                }
            }
        }
        isMainCyclicalVariable = false;
        isSecondaryVariable = false;
    }
    if (i == 2*numberOfEquations) {
        dsSSystemSolveEquations(ssys);
    } else {
        DSSSystemFree(ssys);
        ssys = NULL;
    }
bail:
    return ssys;
}


#if defined (__APPLE__) && defined (__MACH__)
#pragma mark - Getter Methods
#endif

extern const DSUInteger DSSSystemNumberOfEquations(const DSSSystem * ssys)
{
        DSUInteger numberOfEquations = 0;
        if (ssys == NULL) {
                DSError(M_DS_NULL ": S-System is NULL", A_DS_ERROR);
                goto bail;
        }
        numberOfEquations = DSVariablePoolNumberOfVariables(DSSSysXd(ssys));
bail:
        return numberOfEquations;
}

extern const DSUInteger DSSSystemNumberOfConservations(const DSSSystem * ssys)
{
    DSUInteger numberOfConservations = 0;
    if (ssys == NULL) {
        DSError(M_DS_NULL ": S-System is NULL", A_DS_ERROR);
        goto bail;
    }
    numberOfConservations = ssys->numberOfConservations;
bail:
    return numberOfConservations;
}


static void dsSSystemEquationAddPositiveTermToString(const DSSSystem *ssys, 
                                                       const DSUInteger equation,
                                                       char ** string, 
                                                       DSUInteger *length)
{
        DSUInteger i, numberOfXd, numberOfXi;
        DSMatrix *Gd, *Gi;
        DSMatrix *alpha;
        double value;
        char tempString[100];
        const char * name;
        if (ssys == NULL) {
                DSError(M_DS_SSYS_NULL, A_DS_ERROR);
                goto bail;
        }
        numberOfXd = DSVariablePoolNumberOfVariables(DSSSysXd(ssys));
        if (equation >= numberOfXd) {
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
        numberOfXi = DSVariablePoolNumberOfVariables(DSSSysXi(ssys));
        Gd = DSSSysGd(ssys);
        Gi = DSSSysGi(ssys);
        alpha = DSSSysAlpha(ssys);
        sprintf(tempString, "%lf", DSMatrixDoubleValue(alpha, equation, 0));
        if (*length-strlen(*string) < 100) {
                *length += 1000;
                *string = DSSecureRealloc(*string, sizeof(char)**length);
        }
        strncat(*string, tempString, *length-strlen(*string));
        for (i = 0; i < numberOfXd+numberOfXi; i++) {
                if (*length-strlen(*string) < 100) {
                        *length += 1000;
                        *string = DSSecureRealloc(*string, sizeof(char)**length);
                }
                if (i < numberOfXi) {
                        name = DSVariableName(DSVariablePoolAllVariables(DSSSysXi(ssys))[i]);
                        value = DSMatrixDoubleValue(Gi, equation, i);
                        
                } else {
                        name = DSVariableName(DSVariablePoolAllVariables(DSSSysXd(ssys))[i-numberOfXi]);
                        value = DSMatrixDoubleValue(Gd, equation, i-numberOfXi);
                }
                if (value == 0.0)
                        continue;
                if (value == 1.0)
                        sprintf(tempString, "*%s", name);
                else
                        sprintf(tempString, "*%s^%lf", name, value);
                strncat(*string, tempString, *length-strlen(*string));
        }
bail:
        return;
}

static void dsSSystemEquationAddNegativeTermToString(const DSSSystem *ssys, 
                                                       const DSUInteger equation, 
                                                       char ** string, 
                                                       DSUInteger *length)
{
        DSUInteger i, numberOfXd, numberOfXi;
        DSMatrix *Hd, *Hi;
        DSMatrix *beta;
        double value;
        char tempString[100];
        const char * name;
        if (ssys == NULL) {
                DSError(M_DS_SSYS_NULL, A_DS_ERROR);
                goto bail;
        }
        numberOfXd = DSVariablePoolNumberOfVariables(DSSSysXd(ssys));
        if (equation >= numberOfXd) {
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
        numberOfXi = DSVariablePoolNumberOfVariables(DSSSysXi(ssys));
        Hd = DSSSysHd(ssys);
        Hi = DSSSysHi(ssys);
        beta = DSSSysBeta(ssys);
        sprintf(tempString, "%lf", DSMatrixDoubleValue(beta, equation, 0));
        if (*length-strlen(*string) < 100) {
                *length += 1000;
                *string = DSSecureRealloc(*string, sizeof(char)**length);
        }
        strncat(*string, tempString, *length-strlen(*string));
        for (i = 0; i < numberOfXd+numberOfXi; i++) {
                if (*length-strlen(*string) < 100) {
                        *length += 1000;
                        *string = DSSecureRealloc(*string, sizeof(char)**length);
                }
                if (i < numberOfXi) {
                        name = DSVariableName(DSVariablePoolAllVariables(DSSSysXi(ssys))[i]);
                        value = DSMatrixDoubleValue(Hi, equation, i);
                        
                } else {
                        name = DSVariableName(DSVariablePoolAllVariables(DSSSysXd(ssys))[i-numberOfXi]);
                        value = DSMatrixDoubleValue(Hd, equation, i-numberOfXi);
                }
                if (value == 0.0)
                        continue;
                if (value == 1.0)
                        sprintf(tempString, "*%s", name);
                else
                        sprintf(tempString, "*%s^%lf", name, value);
                strncat(*string, tempString, *length-strlen(*string));
        }
bail:
        return;
}


extern DSExpression ** DSSSystemEquations(const DSSSystem *ssys)
{
        DSUInteger i, j, index, sum, numberOfEquations;
        DSExpression ** equations = NULL, *root, *lhs, *rhs;
        char *varName;
        if (ssys == NULL) {
                DSError(M_DS_SSYS_NULL, A_DS_ERROR);
                goto bail;
        }
        numberOfEquations = DSSSystemNumberOfEquations(ssys);
        if (numberOfEquations == 0) {
                DSError("S-System being accessed has no equations", A_DS_ERROR);
                goto bail;
        }
        equations = DSSecureCalloc(sizeof(DSExpression *), numberOfEquations);
        for (i = 0; i < numberOfEquations; i++) {
                varName = DSVariableName(DSVariablePoolVariableAtIndex(DSSSystemXd(ssys), i));
                root = dsExpressionAllocWithOperator('=');
                // Check if varName is algebraic
                if (DSVariablePoolHasVariableWithName(DSSSysXd_a(ssys), varName) == true) {
                        sum = 0;
                        for (j = 0; j < DSVariablePoolNumberOfVariables(DSSSystemXd(ssys)); j++) {
                                if (DSMatrixDoubleValue(DSSSystemHd(ssys), i, j) != 0) {
                                        sum++;
                                        index = j;
                                        if (DSMatrixDoubleValue(DSSSystemHd(ssys), i, j) != 1) {
                                                sum++;
                                        }
                                }
                        }
                        for (j = 0; j < DSVariablePoolNumberOfVariables(DSSSystemXi(ssys)); j++) {
                                if (DSMatrixDoubleValue(DSSSystemHi(ssys), i, j) != 0) {
                                        sum++;
                                }
                        }
                        // Check if matrices DSSSystemAlpha or DSSSystemBeta are equal to zero. If not, increase sum
                        if( DSMatrixDoubleValue(DSSSystemAlpha(ssys), i, 0) != 0 || DSMatrixDoubleValue(DSSSystemBeta(ssys), i, 0) != 0 )
                                sum++;
                        if (sum == 1 && index == DSVariablePoolIndexOfVariableWithName(DSSSystemXd(ssys), varName)) {
                                lhs = dsExpressionAllocWithVariableName(varName);
                                rhs = DSExpressionFromPowerlawInMatrixForm(i,
                                                                           DSSSystemGd(ssys),
                                                                           DSSSystemXd(ssys), DSSSystemGi(ssys), DSSSystemXi(ssys), DSSSystemAlpha(ssys));
                        } else {
                                lhs = dsExpressionAllocWithConstant(0.0);
                                rhs = DSExpressionSubstractExpressions(DSExpressionFromPowerlawInMatrixForm(i,
                                                                                                      DSSSystemGd(ssys),
                                                                                                      DSSSystemXd(ssys),
                                                                                                      DSSSystemGi(ssys),
                                                                                                      DSSSystemXi(ssys),
                                                                                                      DSSSystemAlpha(ssys)),
                                                                       DSExpressionFromPowerlawInMatrixForm(i,
                                                                                                      DSSSystemHd(ssys),
                                                                                                      DSSSystemXd(ssys),
                                                                                                      DSSSystemHi(ssys),
                                                                                                      DSSSystemXi(ssys),
                                                                                                      DSSSystemBeta(ssys)));
                        }
                } else {
                        lhs = dsExpressionAllocWithOperator('.');
                        DSExpressionAddBranch(lhs, dsExpressionAllocWithVariableName(varName));
                        rhs = DSExpressionSubstractExpressions(DSExpressionFromPowerlawInMatrixForm(i,
                                                                                              DSSSystemGd(ssys),
                                                                                              DSSSystemXd(ssys),
                                                                                              DSSSystemGi(ssys),
                                                                                              DSSSystemXi(ssys),
                                                                                              DSSSystemAlpha(ssys)),
                                                         DSExpressionFromPowerlawInMatrixForm(i,
                                                                                              DSSSystemHd(ssys),
                                                                                              DSSSystemXd(ssys),
                                                                                              DSSSystemHi(ssys),
                                                                                              DSSSystemXi(ssys),
                                                                                              DSSSystemBeta(ssys)));
                        
                }
                DSExpressionAddBranch(root, lhs);
                DSExpressionAddBranch(root, rhs);
                equations[i] = root;
        }
bail:
        return equations;
}

static void dsSSystemSolutionToString(const DSSSystem *ssys, 
                                      const DSUInteger equation, 
                                      char ** string, 
                                      DSUInteger *length, const bool inLog)
{
        DSUInteger i, numberOfXd, numberOfXi;
        DSMatrix *MAi, *MB, *B, *Ai;
        char tempString[100] = "\0";
        const char *name;
        double value;
        if (ssys == NULL) {
                DSError(M_DS_SSYS_NULL, A_DS_ERROR);
                goto bail;
        }
        numberOfXd = DSVariablePoolNumberOfVariables(DSSSysXd(ssys));
        if (equation >= numberOfXd) {
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
        if (DSSSystemHasSolution(ssys) == false) {
                goto bail;
        }
        numberOfXi = DSVariablePoolNumberOfVariables(DSSSysXi(ssys));
        B = DSSSystemB(ssys);
        Ai = DSSSystemAi(ssys);
        if (numberOfXi != 0) {
                MAi = DSMatrixByMultiplyingMatrix(DSSSysM(ssys), Ai);
                DSMatrixFree(Ai);
        }
        MB = DSMatrixByMultiplyingMatrix(DSSSysM(ssys), B);
        DSMatrixFree(B);
        if (inLog == true) 
                sprintf(tempString, "%lf", DSMatrixDoubleValue(MB, equation, 0));
        else
                sprintf(tempString, "10^%lf", DSMatrixDoubleValue(MB, equation, 0));
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
                name = DSVariableName(DSVariablePoolAllVariables(DSSSysXi(ssys))[i]);
                value = -DSMatrixDoubleValue(MAi, equation, i);
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
        if (numberOfXi != 0)
                DSMatrixFree(MAi);
        DSMatrixFree(MB);
bail:
        return;
}


static void dsUnstableCaseProcessPseudoInverse(DSMatrix *pInverse, const DSSSystem *ssys)
{
    // This step is necessary to obtain correct equations for anomalous cases of type I. The preprocessing consinst of two steps. In a first step, all rows of the pseudoinverse (except for blowing indices) are normalized with the values in the diagonal. Then, rows corresponding to the blowing variables are set to zero.
    DSUInteger *indices, i, j, n;
    DSMatrix *Ad = DSSSystemAd(ssys);
    DSUInteger numberBlowing = DSVariablePoolNumberOfVariables(ssys->Xd_b);
    double factor;
    bool normalize;
    
    // Set rows to zero that are blowing variables.
    indices = DSVariablePoolIndicesOfSubPool(ssys->Xd_t, ssys->Xd_b);
    
    // set rows contained in indices to zero.
    for (i=0; i<numberBlowing; i++){
        for(j=0; j<DSMatrixColumns(pInverse); j++){
            DSMatrixSetDoubleValue(pInverse, indices[i], j, 0.0);
        }
    }
    
//    // Normalize rows with the values in the diagonal if row of Ad contains a blowing variable. Skip if diagonal is zero.
//    for (i=0; i<DSMatrixRows(pInverse); i++){
//            normalize = false;
//            for(n=0; n<numberBlowing; n++)
//                if (DSMatrixDoubleValue(Ad, i, indices[n]) != 0.0)
//                normalize = true;
//            factor = fabs(DSMatrixDoubleValue(pInverse, i, i));
//            if (factor != 0.0 && normalize == true){
//                    for(j=0; j<DSMatrixColumns(pInverse); j++){
//                        DSMatrixSetDoubleValue(pInverse, i, j, DSMatrixDoubleValue(pInverse, i, j)/factor);
//                    }
//            }
//    }

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


static void dsuSSystemSolutionToString(const DSSSystem *ssys,
                                      const DSUInteger equation,
                                      char ** string,
                                      DSUInteger *length, const bool inLog)
{
    DSUInteger i, numberOfXd, numberOfXi, n, col, row;
    DSMatrix *MAi, *MB, *MB_blow, *MB_original, *B, *Ai, *Ad, *Ad2, *w, *pInvA, *Identity, *I_pInvA;
    const DSMatrix *pInverse;
    char tempString[100] = "\0";
    const char *name;
    double value;
    
    DSUInteger numberOfBlowing = DSVariablePoolNumberOfVariables(DSSSysXd_b(ssys));
    DSUInteger *indices = DSVariablePoolIndicesOfSubPool(DSSSysXd(ssys), DSSSysXd_b(ssys));

    
    if (ssys == NULL) {
        DSError(M_DS_SSYS_NULL, A_DS_ERROR);
        goto bail;
    }
    numberOfXd = DSVariablePoolNumberOfVariables(DSSSysXd(ssys));
    if (equation >= numberOfXd) {
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
    
    numberOfXi = DSVariablePoolNumberOfVariables(DSSSysXi(ssys));
    B = DSSSystemB(ssys);
    Ai = DSSSystemAi(ssys);
    Ad = DSSSystemAd(ssys);
    Ad2 = DSSSystemAd(ssys);
    
    // Set both rows and colums of blowing variables to zero.
    for(n=0; n<numberOfBlowing; n++){
        for(col=0; col<DSMatrixColumns(Ad2); col++)
            DSMatrixSetDoubleValue(Ad2, indices[n], col, 0.0);
        for(row=0; row<DSMatrixRows(Ad2); row++)
            DSMatrixSetDoubleValue(Ad2, row, indices[n], 0.0);
    }
    
    pInverse = DSMatrixPseudoInverse(Ad2);
    DSMatrixFree(Ad2);
    pInvA = DSMatrixByMultiplyingMatrix(pInverse, Ad);
    Identity = DSMatrixIdentity(DSMatrixRows(pInvA));
    I_pInvA = DSMatrixBySubstractingMatrix(Identity, pInvA);
    
    DSMatrixFree(Ad);
    DSMatrixFree(pInvA);
    DSMatrixFree(Identity);
    
    w = DSMatrixCalloc(DSMatrixRows(I_pInvA), 1);
    for (i=0; i<numberOfBlowing; i++){
        DSMatrixSetDoubleValue(w, indices[i], 0,
                               log10(DSVariablePoolValueForVariableWithName(ssys->Xd_b,   DSVariableName(DSVariablePoolAllVariables(DSSSysXd_b(ssys))[i]))));
    }

    MB_blow = DSMatrixByMultiplyingMatrix(I_pInvA, w);
    dsUnstableCaseProcessMB_blow(MB_blow);

    
    DSMatrixFree(w);
    DSMatrixFree(I_pInvA);
    
    if (numberOfXi != 0) {
        MAi = DSMatrixByMultiplyingMatrix(pInverse, Ai);
        DSMatrixFree(Ai);
    }
    MB_original = DSMatrixByMultiplyingMatrix(pInverse, B);
    MB = DSMatrixByAddingMatrix(MB_original, MB_blow);
    
    DSMatrixFree(B);
    DSMatrixFree(MB_original);
    DSMatrixFree(MB_blow);
    if (inLog == true)
        sprintf(tempString, "%lf", DSMatrixDoubleValue(MB, equation, 0));
    else
        sprintf(tempString, "10^%lf", DSMatrixDoubleValue(MB, equation, 0));
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
        name = DSVariableName(DSVariablePoolAllVariables(DSSSysXi(ssys))[i]);
        value = -DSMatrixDoubleValue(MAi, equation, i);
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
    if (numberOfXi != 0)
        DSMatrixFree(MAi);
    DSMatrixFree(MB);
    DSMatrixFree((DSMatrix *)pInverse);
bail:
    return;
}

extern DSExpression ** DSSSystemSolution(const DSSSystem *ssys)
{
        DSUInteger i, numberOfEquations, length;
        DSExpression ** solution = NULL;
        char *tempString, * equationString, *varName;
        if (ssys == NULL) {
                DSError(M_DS_SSYS_NULL, A_DS_ERROR);
                goto bail;
        }
        numberOfEquations = DSSSystemNumberOfEquations(ssys);
        if (numberOfEquations == 0) {
                DSError("S-System being accessed has no equations", A_DS_ERROR);
                goto bail;
        }
        if (DSSSystemHasSolution(ssys) == false) {
                goto bail;
        }
        solution = DSSecureCalloc(sizeof(DSExpression *), numberOfEquations);
        length = 1000;
        tempString = DSSecureCalloc(sizeof(char), length);
        for (i = 0; i < numberOfEquations; i++) {
                tempString[0] = '\0';
                dsSSystemSolutionToString(ssys, i, &tempString, &length, false);
                if (strlen(tempString) == 0)
                        break;
                varName = DSVariableName(DSVariablePoolVariableAtIndex(DSSSystemXd(ssys), i));
                equationString = DSSecureCalloc(
                                                sizeof(char),
                                                strlen(tempString)+strlen(varName)+4);
                equationString = strcpy(equationString, varName);
                equationString = strcat(equationString, " = ");
                equationString = strcat(equationString, tempString);
                solution[i] = DSExpressionByParsingString(equationString);
                DSSecureFree(equationString);
        }
        DSSecureFree(tempString);
bail:
        return solution;
}

extern DSExpression ** DSuSSystemSolution(const DSSSystem *ssys)
{
    DSUInteger i, numberOfEquations, length;
    DSExpression ** solution = NULL;
    char *tempString, * equationString, *varName;
    const DSSSystem * ssys_no_algebraic;
    
    if (ssys == NULL) {
        DSError(M_DS_SSYS_NULL, A_DS_ERROR);
        goto bail;
    }
    if (DSSSysXd_b(ssys) == NULL) {
        DSError("S-System being accessed has no blowing variables", A_DS_ERROR);
        goto bail;
    }
    if ( DSVariablePoolNumberOfVariables(DSSSystemXd_t( ssys)) != 0 ){
        ssys_no_algebraic = DSSSystemByRemovingAlgebraicConstraints(ssys);
    } else {
        ssys_no_algebraic = ssys;
    }
    numberOfEquations = DSSSystemNumberOfEquations(ssys_no_algebraic);
    if (numberOfEquations == 0) {
        DSError("S-System being accessed has no equations", A_DS_ERROR);
        goto bail;
    }
    solution = DSSecureCalloc(sizeof(DSExpression *), numberOfEquations);
    length = 1000;
    tempString = DSSecureCalloc(sizeof(char), length);
    for (i = 0; i < numberOfEquations; i++) {
        tempString[0] = '\0';
        dsuSSystemSolutionToString(ssys_no_algebraic, i, &tempString, &length, false);
        if (strlen(tempString) == 0)
            break;
        varName = DSVariableName(DSVariablePoolVariableAtIndex(DSSSystemXd(ssys_no_algebraic), i));
        equationString = DSSecureCalloc(
                                        sizeof(char),
                                        strlen(tempString)+strlen(varName)+4);
        equationString = strcpy(equationString, varName);
        equationString = strcat(equationString, " = ");
        equationString = strcat(equationString, tempString);
        solution[i] = DSExpressionByParsingString(equationString);
        DSSecureFree(equationString);
    }
    DSSecureFree(tempString);
    
    if (DSVariablePoolNumberOfVariables(DSSSystemXd_t( ssys))  != 0 ){
        DSSSystemFree((DSSSystem *)ssys_no_algebraic);
    }
bail:
    return solution;
}

extern DSExpression ** DSSSystemLogarithmicSolution(const DSSSystem *ssys)
{
        DSUInteger i, numberOfEquations, length;
        DSExpression ** solution = NULL;
        char *tempString, * equationString, *varName;
        if (ssys == NULL) {
                DSError(M_DS_SSYS_NULL, A_DS_ERROR);
                goto bail;
        }
        numberOfEquations = DSSSystemNumberOfEquations(ssys);
        if (numberOfEquations == 0) {
                DSError("S-System being accessed has no equations", A_DS_ERROR);
                goto bail;
        }
        if (DSSSystemHasSolution(ssys) == false) {
                goto bail;
        }
        solution = DSSecureCalloc(sizeof(DSExpression *), numberOfEquations);
        length = 1000;
        tempString = DSSecureCalloc(sizeof(char), length);
        for (i = 0; i < numberOfEquations; i++) {
                tempString[0] = '\0';
                dsSSystemSolutionToString(ssys, i, &tempString, &length, true);
                if (strlen(tempString) == 0)
                        break;
                varName = DSVariableName(DSVariablePoolVariableAtIndex(DSSSystemXd(ssys), i));
                equationString = DSSecureCalloc(sizeof(char),
                                                strlen(tempString)+strlen(varName)+10);
                equationString = strcpy(equationString, "log(");
                equationString = strcat(equationString, varName);
                equationString = strcat(equationString, ") = ");
                equationString = strcat(equationString, tempString);
                solution[i] = DSExpressionByParsingString(equationString);
                DSSecureFree(equationString);
        }
        DSSecureFree(tempString);
bail:
        return solution;
}

extern DSExpression ** DSuSSystemLogarithmicSolution(const DSSSystem *ssys)
{
    DSUInteger i, numberOfEquations, length;
    DSExpression ** solution = NULL;
    char *tempString, * equationString, *varName;
    const DSSSystem *ssys_no_algebraic;
    if (ssys == NULL) {
        DSError(M_DS_SSYS_NULL, A_DS_ERROR);
        goto bail;
    }
    if (DSSSysXd_b(ssys) == NULL) {
        DSError("S-System being accessed has no blowing variables", A_DS_ERROR);
        goto bail;
    }
    if (DSVariablePoolNumberOfVariables(DSSSystemXd_t( ssys))  != 0 ){
        ssys_no_algebraic = DSSSystemByRemovingAlgebraicConstraints(ssys);
    } else {
        ssys_no_algebraic = ssys;
    }
    numberOfEquations = DSSSystemNumberOfEquations(ssys_no_algebraic);
    if (numberOfEquations == 0) {
        DSError("S-System being accessed has no equations", A_DS_ERROR);
        goto bail;
    }
    solution = DSSecureCalloc(sizeof(DSExpression *), numberOfEquations);
    length = 1000;
    tempString = DSSecureCalloc(sizeof(char), length);
    for (i = 0; i < numberOfEquations; i++) {
        tempString[0] = '\0';
        dsuSSystemSolutionToString(ssys_no_algebraic, i, &tempString, &length, true);
        if (strlen(tempString) == 0)
            break;
        varName = DSVariableName(DSVariablePoolVariableAtIndex(DSSSystemXd(ssys_no_algebraic), i));
        equationString = DSSecureCalloc(sizeof(char),
                                        strlen(tempString)+strlen(varName)+10);
        equationString = strcpy(equationString, "log(");
        equationString = strcat(equationString, varName);
        equationString = strcat(equationString, ") = ");
        equationString = strcat(equationString, tempString);
        solution[i] = DSExpressionByParsingString(equationString);
        DSSecureFree(equationString);
    }
    DSSecureFree(tempString);
    if (DSVariablePoolNumberOfVariables(DSSSystemXd_t( ssys))  != 0 ){
        DSSSystemFree((DSSSystem *)ssys_no_algebraic);
    }
bail:
    return solution;
}

extern const DSMatrix * DSSSystemAlpha(const DSSSystem * ssys)
{
        DSMatrix *matrix = NULL;
        if (ssys == NULL) {
                DSError(M_DS_NULL ": S-System is NULL", A_DS_ERROR);
                goto bail;
        }
        matrix = DSSSysAlpha(ssys);
bail:
        return matrix;
}

extern const DSMatrix * DSSSystemAlphaAdjusted(const DSSSystem * ssys)
{
        DSMatrix *matrix = NULL;
        if (ssys == NULL) {
                DSError(M_DS_NULL ": S-System is NULL", A_DS_ERROR);
                goto bail;
        }
        matrix = DSSSysAlphaAdjusted(ssys);
bail:
        return matrix;
}

extern const DSMatrix * DSSSystemBeta(const DSSSystem * ssys)
{
        DSMatrix *matrix = NULL;
        if (ssys == NULL) {
                DSError(M_DS_NULL ": S-System is NULL", A_DS_ERROR);
                goto bail;
        }
        matrix = DSSSysBeta(ssys);
bail:
        return matrix;
}

extern const DSMatrix * DSSSystemGd(const DSSSystem * ssys)
{
        DSMatrix *matrix = NULL;
        if (ssys == NULL) {
                DSError(M_DS_NULL ": S-System is NULL", A_DS_ERROR);
                goto bail;
        }
        matrix = DSSSysGd(ssys);
bail:
        return matrix;
}

extern const DSMatrix * DSSSystemGi(const DSSSystem * ssys)
{
        DSMatrix *matrix = NULL;
        if (ssys == NULL) {
                DSError(M_DS_NULL ": S-System is NULL", A_DS_ERROR);
                goto bail;
        }
        matrix = DSSSysGi(ssys);
bail:
        return matrix;
}

extern const DSMatrix * DSSSystemHd(const DSSSystem * ssys)
{
        DSMatrix *matrix = NULL;
        if (ssys == NULL) {
                DSError(M_DS_NULL ": S-System is NULL", A_DS_ERROR);
                goto bail;
        }
        matrix = DSSSysHd(ssys);
bail:
        return matrix;
}

extern const DSMatrix * DSSSystemHi(const DSSSystem * ssys)
{
        DSMatrix *matrix = NULL;
        if (ssys == NULL) {
                DSError(M_DS_NULL ": S-System is NULL", A_DS_ERROR);
                goto bail;
        }
        matrix = DSSSysHi(ssys);
bail:
        return matrix;
}

extern const DSVariablePool * DSSSystemXd(const DSSSystem * ssys)
{
        DSVariablePool *pool = NULL;
        if (ssys == NULL) {
                DSError(M_DS_NULL ": S-System is NULL", A_DS_ERROR);
                goto bail;
        }
        pool = DSSSysXd(ssys);
bail:
        return pool;
}

extern const DSVariablePool * DSSSystemXd_a(const DSSSystem * ssys)
{
        DSVariablePool *pool = NULL;
        if (ssys == NULL) {
                DSError(M_DS_NULL ": S-System is NULL", A_DS_ERROR);
                goto bail;
        }
        pool = DSSSysXd_a(ssys);
bail:
        return pool;
}

extern const DSVariablePool * DSSSystemXd_t(const DSSSystem * const ssys)
{
        DSVariablePool *pool = NULL;
        if (ssys == NULL) {
                DSError(M_DS_NULL ": S-System is NULL", A_DS_ERROR);
                goto bail;
        }
        pool = DSSSysXd_t(ssys);
bail:
        return pool;
}

extern const DSVariablePool * DSSSystemXd_b(const DSSSystem * const ssys)
{
    DSVariablePool *pool = NULL;
    if (ssys == NULL) {
        DSError(M_DS_NULL ": S-System is NULL", A_DS_ERROR);
        goto bail;
    }
    pool = DSSSysXd_b(ssys);
bail:
    return pool;
}

extern const DSVariablePool * DSSSystemXd_a_c(const DSSSystem * const ssys)
{
    DSVariablePool *pool = NULL;
    if (ssys == NULL) {
        DSError(M_DS_NULL ": S-System is NULL", A_DS_ERROR);
        goto bail;
    }
    pool = DSSSysXd_a_c(ssys);
bail:
    return pool;
}

extern const DSVariablePool * DSSSystemXi(const DSSSystem * ssys)
{
        DSVariablePool *pool = NULL;
        if (ssys == NULL) {
                DSError(M_DS_NULL ": S-System is NULL", A_DS_ERROR);
                goto bail;
        }
        pool = DSSSysXi(ssys);
bail:
        return pool;
}

extern const DSMatrix * DSSSystemM(const DSSSystem * ssys)
{
        DSMatrix *M = NULL;
        if (ssys == NULL) {
                DSError(M_DS_NULL ": S-System is NULL", A_DS_ERROR);
                goto bail;
        }
        M = DSSSysM(ssys);
bail:
        return M;
}

extern DSMatrix * DSSSystemM_a(const DSSSystem * ssys)
{
        DSMatrix *Qd_a, *M_a = NULL;
        if (ssys == NULL) {
                DSError(M_DS_NULL ": S-System is NULL", A_DS_ERROR);
                goto bail;
        }
        Qd_a = DSSSystemQd_a(ssys);
        M_a = DSMatrixInverse(Qd_a);
bail:
        return M_a;
}

extern DSMatrix * DSSSystemAd(const DSSSystem * ssys)
{
        DSMatrix *Ad = NULL;
        if (ssys == NULL) {
                DSError(M_DS_NULL ": S-System is NULL", A_DS_ERROR);
                goto bail;
        }
        if (DSVariablePoolNumberOfVariables(DSSSysXd(ssys)) == 0)
                goto bail;
        if (DSSSysGd(ssys) == NULL || DSSSysHd(ssys) == NULL) {
                DSError(M_DS_MAT_NULL ": Gd/hd matrix is null", A_DS_ERROR);
                goto bail;
        }
        Ad = DSMatrixBySubstractingMatrix(DSSSysGd(ssys), DSSSysHd(ssys));
bail:
        return Ad;
}

extern DSMatrix * DSSSystemQd_a(const DSSSystem * ssys)
{
        DSMatrix *Ad_a = NULL;
        DSMatrix *Gd_a = NULL, *Hd_a = NULL;
        DSUInteger i, index, numberOfAuxiliary, * auxiliary_indices = NULL;
        if (ssys == NULL) {
                DSError(M_DS_NULL ": S-System is NULL", A_DS_ERROR);
                goto bail;
        }
        numberOfAuxiliary =DSVariablePoolNumberOfVariables(DSSSysXd_a(ssys));
        if (numberOfAuxiliary == 0)
                goto bail;
        auxiliary_indices = DSSecureCalloc(sizeof(DSUInteger *), numberOfAuxiliary);
        index = 0;
        for (i = 0; i < DSVariablePoolNumberOfVariables(DSSSysXd(ssys)); i++) {
                if (DSVariablePoolHasVariableWithName(DSSSysXd_a(ssys), DSVariableName(DSVariablePoolVariableAtIndex(DSSSysXd(ssys), i)))==true)
                        auxiliary_indices[index++] = i;
        }
        
        Gd_a = DSMatrixSubMatrixIncludingRowsAndColumns(DSSSysGd(ssys), numberOfAuxiliary, numberOfAuxiliary, auxiliary_indices, auxiliary_indices);
        Hd_a = DSMatrixSubMatrixIncludingRowsAndColumns(DSSSysHd(ssys), numberOfAuxiliary, numberOfAuxiliary, auxiliary_indices, auxiliary_indices);
        if (Gd_a == NULL || Hd_a == NULL) {
                DSError(M_DS_MAT_NULL ": Gd/hd matrix is null", A_DS_ERROR);
                goto bail;
        }
        Ad_a = DSMatrixBySubstractingMatrix(Gd_a, Hd_a);
bail:
        if (Gd_a != NULL)
                DSMatrixFree(Gd_a);
        if (Hd_a != NULL)
                DSMatrixFree(Hd_a);
        if (auxiliary_indices != NULL)
                DSSecureFree(auxiliary_indices);
        return Ad_a;
}

extern DSMatrix * DSSSystemQd_t(const DSSSystem * ssys)
{
        DSMatrix *Ad_a = NULL;
        DSMatrix *Gd_a = NULL, *Hd_a = NULL;
        DSUInteger i, index_a, index_t, numberOfAuxiliary, numberOfTime, * auxiliary_indices = NULL, *time_indices = NULL;
        if (ssys == NULL) {
                DSError(M_DS_NULL ": S-System is NULL", A_DS_ERROR);
                goto bail;
        }
        numberOfAuxiliary =DSVariablePoolNumberOfVariables(DSSSysXd_a(ssys));
        numberOfTime =DSVariablePoolNumberOfVariables(DSSSysXd_t(ssys));
        if (numberOfAuxiliary == 0)
                goto bail;
        if (numberOfTime == 0)
                goto bail;
        auxiliary_indices = DSSecureCalloc(sizeof(DSUInteger *), numberOfAuxiliary);
        time_indices = DSSecureCalloc(sizeof(DSUInteger *), numberOfTime);
        index_a = 0;
        index_t = 0;
        for (i = 0; i < DSVariablePoolNumberOfVariables(DSSSysXd(ssys)); i++) {
                if (DSVariablePoolHasVariableWithName(DSSSysXd_a(ssys), DSVariableName(DSVariablePoolVariableAtIndex(DSSSysXd(ssys), i)))==true)
                        auxiliary_indices[index_a++] = i;
                else
                        time_indices[index_t++] = i;
        }
        
        Gd_a = DSMatrixSubMatrixIncludingRowsAndColumns(DSSSysGd(ssys), numberOfAuxiliary, numberOfTime, auxiliary_indices, time_indices);
        Hd_a = DSMatrixSubMatrixIncludingRowsAndColumns(DSSSysHd(ssys), numberOfAuxiliary, numberOfTime, auxiliary_indices, time_indices);
        if (Gd_a == NULL || Hd_a == NULL) {
                DSError(M_DS_MAT_NULL ": Gd/hd matrix is null", A_DS_ERROR);
                goto bail;
        }
        Ad_a = DSMatrixBySubstractingMatrix(Gd_a, Hd_a);
bail:
        if (Gd_a != NULL)
                DSMatrixFree(Gd_a);
        if (Hd_a != NULL)
                DSMatrixFree(Hd_a);
        if (auxiliary_indices != NULL)
                DSSecureFree(auxiliary_indices);
        if (time_indices != NULL)
                DSSecureFree(time_indices);
        return Ad_a;
}

extern DSMatrix * DSSSystemQi_a(const DSSSystem * ssys)
{
        DSMatrix *Ai_a = NULL;
        DSMatrix *Gi_a = NULL, *Hi_a = NULL;
        DSUInteger i, index, numberOfAuxiliary, * auxiliary_indices = NULL;
        if (ssys == NULL) {
                DSError(M_DS_NULL ": S-System is NULL", A_DS_ERROR);
                goto bail;
        }
        numberOfAuxiliary =DSVariablePoolNumberOfVariables(DSSSysXd_a(ssys));
        if (numberOfAuxiliary == 0)
                goto bail;
        auxiliary_indices = DSSecureCalloc(sizeof(DSUInteger *), numberOfAuxiliary);
        index = 0;
        for (i = 0; i < DSVariablePoolNumberOfVariables(DSSSysXd(ssys)); i++) {
                if (DSVariablePoolHasVariableWithName(DSSSysXd_a(ssys), DSVariableName(DSVariablePoolVariableAtIndex(DSSSysXd(ssys), i)))==true)
                        auxiliary_indices[index++] = i;
        }
        Gi_a = DSMatrixSubMatrixIncludingRows(DSSSysGi(ssys), numberOfAuxiliary, auxiliary_indices);
        Hi_a = DSMatrixSubMatrixIncludingRows(DSSSysHi(ssys), numberOfAuxiliary, auxiliary_indices);
        if (Gi_a == NULL || Hi_a == NULL) {
                DSError(M_DS_MAT_NULL ": Gd/hd matrix is null", A_DS_ERROR);
                goto bail;
        }
        Ai_a = DSMatrixBySubstractingMatrix(Gi_a, Hi_a);
bail:
        if (Gi_a != NULL)
                DSMatrixFree(Gi_a);
        if (Hi_a != NULL)
                DSMatrixFree(Hi_a);
        if (auxiliary_indices != NULL)
                DSSecureFree(auxiliary_indices);
        return Ai_a;
}

extern DSMatrix * DSSSystemQB_a(const DSSSystem * ssys)
{
        DSMatrix *B_a = NULL, *B;
        DSUInteger i, index, numberOfAuxiliary, * auxiliary_indices = NULL;
        if (ssys == NULL) {
                DSError(M_DS_NULL ": S-System is NULL", A_DS_ERROR);
                goto bail;
        }
        numberOfAuxiliary =DSVariablePoolNumberOfVariables(DSSSysXd_a(ssys));
        if (numberOfAuxiliary == 0)
                goto bail;
        auxiliary_indices = DSSecureCalloc(sizeof(DSUInteger *), numberOfAuxiliary);
        index = 0;
        for (i = 0; i < DSVariablePoolNumberOfVariables(DSSSysXd(ssys)); i++) {
                if (DSVariablePoolHasVariableWithName(DSSSysXd_a(ssys), DSVariableName(DSVariablePoolVariableAtIndex(DSSSysXd(ssys), i)))==true)
                        auxiliary_indices[index++] = i;
        }
        B = DSSSystemB(ssys);
        B_a = DSMatrixSubMatrixIncludingRows(B, numberOfAuxiliary, auxiliary_indices);
        if (B_a == NULL) {
                DSError(M_DS_MAT_NULL ": B matrix is null", A_DS_ERROR);
                goto bail;
        }
        DSMatrixFree(B);
bail:
        if (auxiliary_indices != NULL)
                DSSecureFree(auxiliary_indices);
        return B_a;
}

extern DSMatrix * DSSSystemAd_a(const DSSSystem * ssys)
{
        DSMatrix *Ad_a = NULL;
        DSMatrix *Gd_a = NULL, *Hd_a = NULL;
        DSUInteger i, index, numberOfAuxiliary, * auxiliary_indices = NULL;
        if (ssys == NULL) {
                DSError(M_DS_NULL ": S-System is NULL", A_DS_ERROR);
                goto bail;
        }
        numberOfAuxiliary =DSVariablePoolNumberOfVariables(DSSSysXd_a(ssys));
        if (numberOfAuxiliary == 0)
                goto bail;
        auxiliary_indices = DSSecureCalloc(sizeof(DSUInteger *), numberOfAuxiliary);
        index = 0;
        for (i = 0; i < DSVariablePoolNumberOfVariables(DSSSysXd(ssys)); i++) {
                if (DSVariablePoolHasVariableWithName(DSSSysXd_a(ssys), DSVariableName(DSVariablePoolVariableAtIndex(DSSSysXd(ssys), i)))==true)
                        auxiliary_indices[index++] = i;
        }
        
        Gd_a = DSMatrixSubMatrixIncludingColumns(DSSSysGd(ssys), numberOfAuxiliary, auxiliary_indices);
        Hd_a = DSMatrixSubMatrixIncludingColumns(DSSSysHd(ssys), numberOfAuxiliary, auxiliary_indices);
        if (Gd_a == NULL || Hd_a == NULL) {
                DSError(M_DS_MAT_NULL ": Gd/hd matrix is null", A_DS_ERROR);
                goto bail;
        }
        Ad_a = DSMatrixBySubstractingMatrix(Gd_a, Hd_a);
bail:
        if (Gd_a != NULL)
                DSMatrixFree(Gd_a);
        if (Hd_a != NULL)
                DSMatrixFree(Hd_a);
        if (auxiliary_indices != NULL)
                DSSecureFree(auxiliary_indices);
        return Ad_a;
}

extern DSMatrix * DSSSystemAd_t(const DSSSystem * ssys)
{
        DSMatrix *Ad_t = NULL;
        DSMatrix *Gd_t = NULL, *Hd_t = NULL;
        DSUInteger i, index, numberOfAuxiliary, * auxiliary_indices = NULL;
        if (ssys == NULL) {
                DSError(M_DS_NULL ": S-System is NULL", A_DS_ERROR);
                goto bail;
        }
        numberOfAuxiliary =DSVariablePoolNumberOfVariables(DSSSysXd_a(ssys));
        if (numberOfAuxiliary == 0)
                goto bail;
        auxiliary_indices = DSSecureCalloc(sizeof(DSUInteger *), numberOfAuxiliary);
        index = 0;
        for (i = 0; i < DSVariablePoolNumberOfVariables(DSSSysXd(ssys)); i++) {
                if (DSVariablePoolHasVariableWithName(DSSSysXd_a(ssys), DSVariableName(DSVariablePoolVariableAtIndex(DSSSysXd(ssys), i)))==true)
                        auxiliary_indices[index++] = i;
        }
        Gd_t = DSMatrixSubMatrixExcludingColumns(DSSSysGd(ssys), numberOfAuxiliary, auxiliary_indices);
        Hd_t = DSMatrixSubMatrixExcludingColumns(DSSSysHd(ssys), numberOfAuxiliary, auxiliary_indices);
        if (Gd_t == NULL || Hd_t == NULL) {
                DSError(M_DS_MAT_NULL ": Gd/hd matrix is null", A_DS_ERROR);
                goto bail;
        }
        Ad_t = DSMatrixBySubstractingMatrix(Gd_t, Hd_t);
bail:
        if (Gd_t != NULL)
                DSMatrixFree(Gd_t);
        if (Hd_t != NULL)
                DSMatrixFree(Hd_t);
        if (auxiliary_indices != NULL)
                DSSecureFree(auxiliary_indices);
        return Ad_t;
}

extern DSMatrix * DSSSystemAi(const DSSSystem * ssys)
{
        DSMatrix *Ai = NULL;
        if (ssys == NULL) {
                DSError(M_DS_NULL ": S-System is NULL", A_DS_ERROR);
                goto bail;
        }
        if (DSVariablePoolNumberOfVariables(DSSSysXi(ssys)) == 0)
                goto bail;
        if (DSSSysGi(ssys) == NULL || DSSSysHi(ssys) == NULL) {
                DSError(M_DS_MAT_NULL ": Gi/hi matrix is null", A_DS_ERROR);
                goto bail;
        }
        Ai = DSMatrixBySubstractingMatrix(DSSSysGi(ssys), DSSSysHi(ssys));
bail:
        return Ai;
}

extern DSMatrix * DSSSystemB(const DSSSystem * ssys)
{
        DSMatrix *B = NULL;
        DSUInteger i;
        if (ssys == NULL) {
                DSError(M_DS_MAT_NULL ": B is NULL", A_DS_ERROR);
                goto bail;
        }
        if (DSSSysAlpha(ssys) == NULL || DSSSysBeta(ssys) == NULL) {
                DSError(M_DS_MAT_NULL ": Alpha/beta matrix is NULL", A_DS_ERROR);
                goto bail;
        }
        if (DSMatrixRows(DSSSysAlpha(ssys)) != DSMatrixRows(DSSSysBeta(ssys))) {
                DSError(M_DS_WRONG ": S-System alpha/beta matrix rows do not match", A_DS_ERROR);
                goto bail;
        }
        B = DSMatrixAlloc(DSMatrixRows(DSSSysBeta(ssys)), 1);
        for (i = 0; i < DSMatrixRows(B); i++) {
                DSMatrixSetDoubleValue(B, i, 0,
                                       log10(DSMatrixDoubleValue(DSSSysBeta(ssys), i, 0)/DSMatrixDoubleValue(DSSSysAlpha(ssys), i, 0)));
        }
bail:
        return B;
}

extern DSMatrix * DSSSystemA(const DSSSystem * ssys)
{
        DSMatrix *A = NULL, *G, *H;
        if (ssys == NULL) {
                DSError(M_DS_NULL ": S-System is NULL", A_DS_ERROR);
                goto bail;
        }
        G = DSSSystemG(ssys);
        H = DSSSystemH(ssys);
        A = DSMatrixBySubstractingMatrix(G, H);
        DSMatrixFree(G);
        DSMatrixFree(H);
bail:
        return A;       
}

extern DSMatrix * DSSSystemG(const DSSSystem *ssys)
{
        DSMatrix *G = NULL;
        if (ssys == NULL) {
                DSError(M_DS_NULL ": S-System is NULL", A_DS_ERROR);
                goto bail;
        }
//        if (DSMatrixRows(DSSSysGd(ssys))  != DSMatrixRows(DSSSysGi(ssys)))
//            goto bail;
        if (DSVariablePoolNumberOfVariables(DSSSystemXi(ssys)) != 0)
                G = DSMatrixAppendMatrices(DSSSysGd(ssys), DSSSysGi(ssys), true);
        else
                G = DSMatrixCopy(DSSSysGd(ssys));
bail:
        return G;
}

extern DSMatrix * DSSSystemH(const DSSSystem *ssys)
{
        DSMatrix *H = NULL;
        if (ssys == NULL) {
                DSError(M_DS_NULL ": S-System is NULL", A_DS_ERROR);
                goto bail;
        }
        if (DSVariablePoolNumberOfVariables(DSSSystemXi(ssys)) != 0)
                H = DSMatrixAppendMatrices(DSSSysHd(ssys), DSSSysHi(ssys), true);
        else
                H = DSMatrixCopy(DSSSysHd(ssys));
bail:
        return H;        
}


extern const bool DSSSystemHasSolution(const DSSSystem * ssys)
{
        bool hasSolution = false;
        if (ssys == NULL) {
                DSError(M_DS_NULL ": S-System is NULL", A_DS_ERROR);
                goto bail;
        }
        if (DSSSysM(ssys) != NULL)
                hasSolution = true;
bail:
        return hasSolution;
}

extern bool DSSSystemIsSingular(const DSSSystem *ssys)
{
        bool isSingular = false;
        if (ssys == NULL) {
                DSError(M_DS_SSYS_NULL, A_DS_ERROR);
                goto bail;
        }
        isSingular = (ssys->modifierFlags & DS_SSYSTEM_FLAG_SINGULAR) ? true : false;
bail:
        return isSingular;
}

extern bool DSSSystemIsConserved(const DSSSystem *ssys)
{
    bool isConserved = false;
    if (ssys == NULL) {
        DSError(M_DS_SSYS_NULL, A_DS_ERROR);
        goto bail;
    }
    isConserved = (ssys->modifierFlags & DS_SSYSTEM_FLAG_CONSERVED) ? true : false;
bail:
    return isConserved;
}

extern bool DSSSystemIsUnstable(const DSSSystem *ssys)
{
    bool isUnstable = false;
    if (ssys == NULL) {
        DSError(M_DS_SSYS_NULL, A_DS_ERROR);
        goto bail;
    }
    isUnstable = (ssys->modifierFlags & DS_SSYSTEM_FLAG_UNSTABLE) ? true : false;
bail:
    return isUnstable;
}

extern bool DSSSystemIsFalseBlowing(const DSSSystem *ssys)
{
    bool isFalseBlowing = false;
    if (ssys == NULL) {
        DSError(M_DS_SSYS_NULL, A_DS_ERROR);
        goto bail;
    }
    
       DSMatrix *Ai = NULL, *Ad = NULL, *b = NULL, *Ai_neg = NULL;
        DSMatrix *A = NULL, *Augmented = NULL;
        DSMatrix *transA = NULL, *transAugmented = NULL;
        DSSSystem *collapsedSystem = NULL;
        DSUInteger rank_A, rank_Augmented;
            
        if( DSSSystemXd_a(ssys) != 0){
            collapsedSystem = DSSSystemByRemovingAlgebraicConstraints(ssys);
            Ad = DSSSystemAd(collapsedSystem);
            Ai = DSSSystemAi(collapsedSystem);
            b = DSSSystemB(collapsedSystem);
            DSSSystemFree(collapsedSystem);
        }else {
            Ad = DSSSystemAd(ssys);
            Ai = DSSSystemAi(ssys);
            b = DSSSystemB(ssys);
        }
        
        //modify sign of matrix Ai
        Ai_neg = DSMatrixByMultiplyingScalar(Ai, -1.0);
        
        // construct matrix A and Augmented
        A = DSMatrixAppendMatrices(Ad, Ai_neg, true);
        Augmented = DSMatrixAppendMatrices(A, b, true);
        
        
        if (DSMatrixRows(Ad) < DSMatrixColumns(Ad)){
            transA = DSMatrixTranspose(Ad);
            rank_A = DSMatrixRank(transA);
            if (transA != NULL)
                DSMatrixFree(transA);
        }else{
            rank_A = DSMatrixRank(Ad);
        }
        
        if (DSMatrixRows(Augmented) < DSMatrixColumns(Augmented)){
            transAugmented = DSMatrixTranspose(Augmented);
            rank_Augmented = DSMatrixRank(transAugmented);
            if (transAugmented != NULL)
                DSMatrixFree(transAugmented);
        }else{
            rank_Augmented = DSMatrixRank(Augmented);
        }
    
        // find out if this case was a false blow up
        if (rank_Augmented != DSMatrixRows(Ad))
                isFalseBlowing = true;
    
        
        // delete variables
        if (Ai != NULL)
            DSMatrixFree(Ai);
        if (Ai_neg != NULL)
            DSMatrixFree(Ai_neg);
        if (Ad != NULL)
            DSMatrixFree(Ad);
        if (b != NULL)
            DSMatrixFree(b);
        if (A != NULL)
            DSMatrixFree(A);
        if (Augmented != NULL)
            DSMatrixFree(Augmented);
    
bail:
    return isFalseBlowing;
}

extern bool DSSSystemAdjustCodominantStoichiometry(const DSSSystem *ssys){
    
    bool AdjustCodominantStoichiometry = false;
    if (ssys == NULL) {
        DSError(M_DS_SSYS_NULL, A_DS_ERROR);
        goto bail;
    }
    AdjustCodominantStoichiometry = (ssys->modifierFlags & DS_SSYSTEM_FLAG_ADJUST_CODOMINANT_STOICHIOMETRY) ? true : false;
bail:
    return AdjustCodominantStoichiometry;
    
}

extern bool DSSSystemShouldFreeXd(const DSSSystem *ssys)
{
        bool shouldFree = false;
        if (ssys == NULL) {
                DSError(M_DS_SSYS_NULL, A_DS_ERROR);
                goto bail;
        }
        shouldFree = (ssys->modifierFlags & DS_SSYSTEM_FLAG_FREE_XD) ? true : false;
bail:
        return shouldFree;
}

extern bool DSSSystemShouldFreeXi(const DSSSystem *ssys)
{
        bool shouldFree = false;
        if (ssys == NULL) {
                DSError(M_DS_SSYS_NULL, A_DS_ERROR);
                goto bail;
        }
        shouldFree = (ssys->modifierFlags & DS_SSYSTEM_FLAG_FREE_XI) ? true : false;
bail:
        return shouldFree;
}


extern void DSSSystemSetIsSingular(DSSSystem *ssys, bool isSingular)
{
        unsigned char newFlag;
        if (ssys == NULL) {
                DSError(M_DS_SSYS_NULL, A_DS_ERROR);
                goto bail;
        }
        newFlag = ssys->modifierFlags & ~DS_SSYSTEM_FLAG_SINGULAR;
        ssys->modifierFlags = (isSingular ? DS_SSYSTEM_FLAG_SINGULAR : 0) | newFlag;
bail:
        return;
}

extern void DSSSystemSetIsConserved(DSSSystem *ssys, bool isConserved)
{
    unsigned char newFlag;
    if (ssys == NULL) {
        DSError(M_DS_SSYS_NULL, A_DS_ERROR);
        goto bail;
    }
    newFlag = ssys->modifierFlags & ~DS_SSYSTEM_FLAG_CONSERVED;
    ssys->modifierFlags = (isConserved ? DS_SSYSTEM_FLAG_CONSERVED : 0) | newFlag;
bail:
    return;
}

extern void DSSSystemSetIsUnstable(DSSSystem *ssys, bool isUnstable)
{
    unsigned char newFlag;
    if (ssys == NULL) {
        DSError(M_DS_SSYS_NULL, A_DS_ERROR);
        goto bail;
    }
    newFlag = ssys->modifierFlags & ~DS_SSYSTEM_FLAG_UNSTABLE;
    ssys->modifierFlags = (isUnstable ? DS_SSYSTEM_FLAG_UNSTABLE : 0) | newFlag;
bail:
    return;
}

extern void DSSSystemSetAdjustCodominantStoichiometry(DSSSystem *ssys, bool AdjustStoichiometry){
    
    unsigned char newFlag;
    if (ssys == NULL) {
        DSError(M_DS_SSYS_NULL, A_DS_ERROR);
        goto bail;
    }
    newFlag = ssys->modifierFlags & ~DS_SSYSTEM_FLAG_ADJUST_CODOMINANT_STOICHIOMETRY;
    ssys->modifierFlags = (AdjustStoichiometry ? DS_SSYSTEM_FLAG_ADJUST_CODOMINANT_STOICHIOMETRY : 0) | newFlag;
bail:
    return;

}


extern void DSSSystemSetShouldFreeXd(DSSSystem *ssys, bool shouldFreeXd)
{
        unsigned char newFlag;
        if (ssys == NULL) {
                DSError(M_DS_SSYS_NULL, A_DS_ERROR);
                goto bail;
        }
        newFlag = ssys->modifierFlags & ~DS_SSYSTEM_FLAG_FREE_XD;
        ssys->modifierFlags = (shouldFreeXd ? DS_SSYSTEM_FLAG_FREE_XD : 0) | newFlag;
bail:
        return;
}

extern void DSSSystemSetShouldFreeXi(DSSSystem *ssys, bool shouldFreeXi)
{
        unsigned char newFlag;
        if (ssys == NULL) {
                DSError(M_DS_SSYS_NULL, A_DS_ERROR);
                goto bail;
        }
        newFlag = ssys->modifierFlags & ~DS_SSYSTEM_FLAG_FREE_XI;
        ssys->modifierFlags = (shouldFreeXi ? DS_SSYSTEM_FLAG_FREE_XI : 0) | newFlag;
bail:
        return;
}


#if defined (__APPLE__) && defined (__MACH__)
#pragma mark - S-System functions
#endif

extern DSMatrix * DSSSystemSteadyStateValues(const DSSSystem *ssys, const DSVariablePool *Xi0)
{
        DSMatrix * steadyState = NULL;
        DSMatrix *Xi = NULL;
        DSMatrix *MAi, *MAiXi, *B = NULL, *Ai = NULL;
        DSVariablePool *pool = NULL;
        DSUInteger i;
        const char *name;
        if (ssys == NULL) {
                DSError(M_DS_SSYS_NULL, A_DS_ERROR);
                goto bail;
        }
        if (Xi0 == NULL && DSVariablePoolNumberOfVariables(DSSSysXi(ssys)) != 0) {
                DSError(M_DS_VAR_NULL ": Xi0 variable pool is NULL", A_DS_ERROR);
                goto bail;
        }
        if (DSSSystemHasSolution(ssys) == false)
                goto bail;
        B = DSSSystemB(ssys);
        steadyState = DSMatrixByMultiplyingMatrix(DSSSysM(ssys), B);
        if (DSVariablePoolNumberOfVariables(DSSSysXi(ssys)) != 0) {
                pool = DSVariablePoolAlloc();
                for (i=0; i < DSVariablePoolNumberOfVariables(DSSSysXi(ssys)); i++) {
                        name =  DSVariableName(DSVariablePoolAllVariables(DSSSysXi(ssys))[i]);
                        DSVariablePoolAddVariableWithName(pool, name);
                        if (DSVariablePoolHasVariableWithName(Xi0, name) == false) {
                                DSMatrixFree(steadyState);
                                DSVariablePoolFree(pool);
                                steadyState = NULL;
                                goto bail;
                        }
                        DSVariablePoolSetValueForVariableWithName(pool,  name, DSVariableValue(DSVariablePoolVariableWithName(Xi0, name)));
                }
                Ai = DSSSystemAi(ssys);
                Xi = DSVariablePoolValuesAsVector(pool, false);
                DSMatrixApplyFunction(Xi, log10);
                MAi = DSMatrixByMultiplyingMatrix(DSSSysM(ssys), Ai);
                MAiXi = DSMatrixByMultiplyingMatrix(MAi, Xi);
                DSMatrixSubstractByMatrix(steadyState, MAiXi);
                DSMatrixFree(Ai);
                DSMatrixFree(Xi);
                DSMatrixFree(MAi);
                DSMatrixFree(MAiXi);
                DSVariablePoolFree(pool);
        }
bail:
        if (B != NULL)
                DSMatrixFree(B);
        return steadyState;
}

extern DSMatrix * DSuSSystemSteadyStateValues(const DSSSystem *ssys, const DSVariablePool *Xi0)
{
    DSMatrix * steadyState = NULL;
    DSMatrix *Xi = NULL;
    DSMatrix *MAi, *MAiXi, *B = NULL, *Ai = NULL, *Ad, *Ad2,  *pInvA, *Identity, *I_pInvA, *w, *MB_blow;
    DSMatrix *MB_original;
    const DSMatrix *pInverse;
    DSVariablePool *pool = NULL;
    DSUInteger i;
    const char *name;
    const DSSSystem *ssys_no_algebraic;
    DSUInteger n, row, col;
    if (ssys == NULL) {
        DSError(M_DS_SSYS_NULL, A_DS_ERROR);
        goto bail;
    }
    if (Xi0 == NULL && DSVariablePoolNumberOfVariables(DSSSysXi(ssys)) != 0) {
        DSError(M_DS_VAR_NULL ": Xi0 variable pool is NULL", A_DS_ERROR);
        goto bail;
    }
    if (DSSSystemIsUnstable(ssys) == false){
        goto bail;
    }
    if ( DSVariablePoolNumberOfVariables(DSSSystemXd_t( ssys)) != 0 ){
        ssys_no_algebraic = DSSSystemByRemovingAlgebraicConstraints(ssys);
    } else {
        ssys_no_algebraic = ssys;
    }
    
    DSUInteger numberOfBlowing = DSVariablePoolNumberOfVariables(DSSSysXd_b(ssys_no_algebraic));
    DSUInteger *indices = DSVariablePoolIndicesOfSubPool(DSSSysXd(ssys_no_algebraic), DSSSysXd_b(ssys_no_algebraic));
    
    B = DSSSystemB(ssys_no_algebraic);
    Ad = DSSSystemAd(ssys_no_algebraic);
    Ad2 = DSSSystemAd(ssys_no_algebraic);
    
    // Set both rows and colums of blowing variables to zero.
    for(n=0; n<numberOfBlowing; n++){
        for(col=0; col<DSMatrixColumns(Ad2); col++)
        DSMatrixSetDoubleValue(Ad2, indices[n], col, 0.0);
        for(row=0; row<DSMatrixRows(Ad2); row++)
        DSMatrixSetDoubleValue(Ad2, row, indices[n], 0.0);
    }
    
    pInverse = DSMatrixPseudoInverse(Ad2);
    DSMatrixFree(Ad2);

//    pInverse = DSUnstableCaseGetSubSetPseudoInverse(ssys_no_algebraic->Xd_t, NULL, ssys_no_algebraic->Xd_b, Ad);
//    dsUnstableCaseProcessPseudoInverse(pInverse, ssys);
    
    pInvA = DSMatrixByMultiplyingMatrix(pInverse, Ad);
    Identity = DSMatrixIdentity(DSMatrixRows(pInvA));
    I_pInvA = DSMatrixBySubstractingMatrix(Identity, pInvA);
    
    DSMatrixFree(Ad);
    DSMatrixFree(pInvA);
    DSMatrixFree(Identity);
    w = DSMatrixCalloc(DSMatrixRows(I_pInvA), 1);
    
    for (i=0; i<numberOfBlowing; i++){
        DSMatrixSetDoubleValue(w, indices[i], 0,
                               log10(DSVariablePoolValueForVariableWithName(ssys_no_algebraic->Xd_b, DSVariableName(DSVariablePoolAllVariables(DSSSysXd_b(ssys_no_algebraic))[i]))));
    }
    MB_blow = DSMatrixByMultiplyingMatrix(I_pInvA, w);
    dsUnstableCaseProcessMB_blow(MB_blow);
    DSMatrixFree(w);
    MB_original = DSMatrixByMultiplyingMatrix(pInverse, B);
    steadyState = DSMatrixByAddingMatrix(MB_original, MB_blow);
    DSMatrixFree(MB_original);
    DSMatrixFree(MB_blow);
    DSMatrixFree(I_pInvA);
    
    
    if (DSVariablePoolNumberOfVariables(DSSSysXi(ssys_no_algebraic)) != 0) {
        pool = DSVariablePoolAlloc();
        for (i=0; i < DSVariablePoolNumberOfVariables(DSSSysXi(ssys_no_algebraic)); i++) {
            name =  DSVariableName(DSVariablePoolAllVariables(DSSSysXi(ssys_no_algebraic))[i]);
            DSVariablePoolAddVariableWithName(pool, name);
            if (DSVariablePoolHasVariableWithName(Xi0, name) == false) {
                DSMatrixFree(steadyState);
                DSVariablePoolFree(pool);
                steadyState = NULL;
                printf("bail variablepool has no variable \n");
                goto bail;
            }
            DSVariablePoolSetValueForVariableWithName(pool,  name, DSVariableValue(DSVariablePoolVariableWithName(Xi0, name)));
        }
        Ai = DSSSystemAi(ssys_no_algebraic);
        Xi = DSVariablePoolValuesAsVector(pool, false);
        DSMatrixApplyFunction(Xi, log10);
        MAi = DSMatrixByMultiplyingMatrix(pInverse, Ai);
        MAiXi = DSMatrixByMultiplyingMatrix(MAi, Xi);
        DSMatrixSubstractByMatrix(steadyState, MAiXi);
        DSMatrixFree(Ai);
        DSMatrixFree(Xi);
        DSMatrixFree(MAi);
        DSMatrixFree(MAiXi);
        DSVariablePoolFree(pool);
    }
    DSMatrixFree((DSMatrix *) pInverse);
    if ( DSVariablePoolNumberOfVariables(DSSSystemXd_t( ssys)) != 0 )
        DSSSystemFree((DSSSystem *)ssys_no_algebraic);
bail:
    if (B != NULL)
        DSMatrixFree(B);
    return steadyState;
}

extern DSMatrix * DSuSSystemSteadyStateValuesForConservedVariables(const DSSSystem *ssys, const DSVariablePool *Xi0){
    
            // This function should calculate steady states for variable pool Xd_a_c. Function DSSSystemSteadyStateValues()
            // was taken as template.
    
            DSMatrix * steadyState = NULL;
            DSMatrix *Xi = NULL;
            DSMatrix *MAi, *MAiXi, *B = NULL, *Ai = NULL, *M = NULL, *Ad = NULL;
            DSVariablePool *pool = NULL;
            DSUInteger i, n, *indices = NULL;
            const char *name;
            if (ssys == NULL) {
                DSError(M_DS_SSYS_NULL, A_DS_ERROR);
                goto bail;
            }
            if (Xi0 == NULL && DSVariablePoolNumberOfVariables(DSSSysXi(ssys)) != 0) {
                DSError(M_DS_VAR_NULL ": Xi0 variable pool is NULL", A_DS_ERROR);
                goto bail;
            }
            if (DSSSystemIsUnstable(ssys) == false){
                goto bail;
            }
    
            n = DSVariablePoolNumberOfVariables(DSSSysXd_a_c(ssys));
            indices = DSVariablePoolIndicesOfSubPool(DSSSysXd(ssys), DSSSysXd_a_c(ssys));
    
            if (indices == NULL){
                goto bail;
            }
    
            // let's get vector B
            B = DSMatrixSubMatrixIncludingRows(DSSSystemB(ssys), n, indices);
    
            // Now let's get matrix Ai
            Ai = DSMatrixSubMatrixIncludingRows(DSSSystemAi(ssys), n, indices) ;
    
            // Now let's get Matrix M
            Ad = DSMatrixSubMatrixIncludingRowsAndColumns(DSSSystemAd(ssys), n, n, indices, indices);
            M = DSMatrixInverse(Ad);
    
            steadyState = DSMatrixByMultiplyingMatrix(M, B);
            if (DSVariablePoolNumberOfVariables(DSSSysXi(ssys)) != 0) {
                pool = DSVariablePoolAlloc();
                for (i=0; i < DSVariablePoolNumberOfVariables(DSSSysXi(ssys)); i++) {
                    name =  DSVariableName(DSVariablePoolAllVariables(DSSSysXi(ssys))[i]);
                    DSVariablePoolAddVariableWithName(pool, name);
                    if (DSVariablePoolHasVariableWithName(Xi0, name) == false) {
                        DSMatrixFree(steadyState);
                        DSVariablePoolFree(pool);
                        steadyState = NULL;
                        goto bail;
                    }
                    DSVariablePoolSetValueForVariableWithName(pool,  name, DSVariableValue(DSVariablePoolVariableWithName(Xi0, name)));
                }
                
                Xi = DSVariablePoolValuesAsVector(pool, false);
                DSMatrixApplyFunction(Xi, log10);
                MAi = DSMatrixByMultiplyingMatrix(M, Ai);
                MAiXi = DSMatrixByMultiplyingMatrix(MAi, Xi);
                DSMatrixSubstractByMatrix(steadyState, MAiXi);
                DSMatrixFree(Ai);
                DSMatrixFree(Xi);
                DSMatrixFree(MAi);
                DSMatrixFree(MAiXi);
                DSVariablePoolFree(pool);
            }

        bail:
            if (B != NULL)
                DSMatrixFree(B);
            if (indices != NULL)
                DSSecureFree(indices);
            if (Ad != NULL)
                DSMatrixFree(Ad);
            if (M != NULL)
                DSMatrixFree(M);
            return steadyState;

}

extern DSMatrix * DSSSystemAuxiliaryVariablesForSteadyState(const DSSSystem *ssys, const DSVariablePool *Xdt0, const DSVariablePool *Xi0)
{
        DSMatrix * auxSolution = NULL;
        DSMatrix *Yi = NULL, *Yd_t = NULL;
        DSMatrix *M = NULL, *MA, *MAY, *B = NULL, *Ai = NULL, *Ad_t = NULL;
        DSVariablePool *pool = NULL;
        DSUInteger i;
        const char *name;
        if (ssys == NULL) {
                DSError(M_DS_SSYS_NULL, A_DS_ERROR);
                goto bail;
        }
        if (Xi0 == NULL && DSVariablePoolNumberOfVariables(DSSSysXi(ssys)) != 0) {
                DSError(M_DS_VAR_NULL ": Xi0 variable pool is NULL", A_DS_ERROR);
                goto bail;
        }
        if (Xdt0 == NULL && DSVariablePoolNumberOfVariables(DSSSysXd_t(ssys)) != 0) {
                DSError(M_DS_VAR_NULL ": Xd_t0 variable pool is NULL", A_DS_ERROR);
                goto bail;
        }
        B = DSSSystemQB_a(ssys);
        M = DSSSystemM_a(ssys);
        auxSolution = DSMatrixByMultiplyingMatrix(M, B);
        DSMatrixFree(B);
        if (DSVariablePoolNumberOfVariables(DSSSysXi(ssys)) != 0) {
                pool = DSVariablePoolAlloc();
                for (i=0; i < DSVariablePoolNumberOfVariables(DSSSysXi(ssys)); i++) {
                        name =  DSVariableName(DSVariablePoolAllVariables(DSSSysXi(ssys))[i]);
                        DSVariablePoolAddVariableWithName(pool, name);
                        if (DSVariablePoolHasVariableWithName(Xi0, name) == false) {
                                DSMatrixFree(auxSolution);
                                auxSolution = NULL;
                                DSVariablePoolFree(pool);
                                goto bail;
                        }
                        DSVariablePoolSetValueForVariableWithName(pool,  name, DSVariableValue(DSVariablePoolVariableWithName(Xi0, name)));
                }
                Yi = DSVariablePoolValuesAsVector(pool, false);
                DSMatrixApplyFunction(Yi, log10);
                Ai = DSSSystemQi_a(ssys);
                MA = DSMatrixByMultiplyingMatrix(M, Ai);
                MAY = DSMatrixByMultiplyingMatrix(MA, Yi);
                DSMatrixSubstractByMatrix(auxSolution, MAY);
                DSMatrixFree(Yi);
                DSMatrixFree(Ai);
                DSMatrixFree(MA);
                DSMatrixFree(MAY);
                DSVariablePoolFree(pool);
        }
        if (DSVariablePoolNumberOfVariables(DSSSysXd_t(ssys)) != 0) {
                pool = DSVariablePoolAlloc();
                for (i=0; i < DSVariablePoolNumberOfVariables(DSSSysXd_t(ssys)); i++) {
                        name =  DSVariableName(DSVariablePoolAllVariables(DSSSysXd_t(ssys))[i]);
                        DSVariablePoolAddVariableWithName(pool, name);
                        if (DSVariablePoolHasVariableWithName(Xdt0, name) == false) {
                                DSMatrixFree(auxSolution);
                                auxSolution = NULL;
                                DSVariablePoolFree(pool);
                                goto bail;
                        }
                        DSVariablePoolSetValueForVariableWithName(pool,  name, DSVariableValue(DSVariablePoolVariableWithName(Xdt0, name)));
                }
                Yd_t = DSVariablePoolValuesAsVector(pool, false);
                DSMatrixApplyFunction(Yd_t, log10);
                Ad_t = DSSSystemQd_t(ssys);
                MA = DSMatrixByMultiplyingMatrix(M, Ad_t);
                MAY = DSMatrixByMultiplyingMatrix(MA, Yd_t);
                DSMatrixSubstractByMatrix(auxSolution, MAY);
                DSMatrixFree(Yd_t);
                DSMatrixFree(Ad_t);
                DSMatrixFree(MA);
                DSMatrixFree(MAY);
                DSVariablePoolFree(pool);
        }
bail:
        if (M != NULL)
                DSMatrixFree(M);
        return auxSolution;
}

extern double DSSSystemSteadyStateFunction(const DSSSystem *ssys, const DSVariablePool *Xi0, const char * function)
{
        DSMatrix *ss=NULL;
        DSUInteger i;
        DSVariablePool *pool = NULL;
        const char *name;
        double value = NAN;
        DSExpression *expr = NULL;
        if (ssys == NULL) {
                DSError(M_DS_SSYS_NULL, A_DS_ERROR);
                goto bail;
        }
        if (Xi0 == NULL && DSVariablePoolNumberOfVariables(DSSSysXi(ssys)) != 0) {
                DSError(M_DS_VAR_NULL ": Xi0 variable pool is NULL", A_DS_ERROR);
                goto bail;
        }
        if (DSSSystemHasSolution(ssys) == false)
                goto bail;
        ss = DSSSystemSteadyStateValues(ssys, Xi0);
        if (ss == NULL)
                goto bail;
        pool = DSVariablePoolAlloc();
        for (i = 0; i < DSVariablePoolNumberOfVariables(DSSSysXd(ssys)); i++) {
                name = DSVariableName(DSVariablePoolAllVariables(DSSSysXd(ssys))[i]);
                DSVariablePoolAddVariableWithName(pool, name);
                DSVariablePoolSetValueForVariableWithName(pool, name, pow(10, DSMatrixDoubleValue(ss, i, 0)));
        }
        for (i = 0; i < DSVariablePoolNumberOfVariables(Xi0); i++) {
                name = DSVariableName(DSVariablePoolAllVariables(Xi0)[i]);
                DSVariablePoolAddVariableWithName(pool, name);
                DSVariablePoolSetValueForVariableWithName(pool, name, DSVariableValue(DSVariablePoolAllVariables(Xi0)[i]));
        }
        expr = DSExpressionByParsingString(function);
        if (expr != NULL)
                value = DSExpressionEvaluateWithVariablePool(expr, pool);
bail:
        if (ss != NULL)
                DSMatrixFree(ss);
        if (pool != NULL)
                DSVariablePoolFree(pool);
        if (expr != NULL)
                DSExpressionFree(expr);
        return value;   
}

extern DSMatrix * DSSSystemSteadyStateFluxForDependentVariables(const DSSSystem * ssys,
                                                                const DSVariablePool * Xd0,
                                                                const DSVariablePool * Xi0)
{
        DSMatrix * flux = NULL, *Xi = NULL, *ss=NULL;
        DSMatrix * gi0= NULL, *alpha = NULL;
        DSUInteger i;
        DSVariablePool *pool = NULL;
        const char *name;
        double value;
        if (ssys == NULL) {
                DSError(M_DS_SSYS_NULL, A_DS_ERROR);
                goto bail;
        }
        if (Xd0 == NULL && DSVariablePoolNumberOfVariables(DSSSysXd(ssys)) != 0) {
                DSError(M_DS_VAR_NULL ": Xd0 variable pool is NULL", A_DS_ERROR);
                goto bail;
        }
        if (Xi0 == NULL && DSVariablePoolNumberOfVariables(DSSSysXi(ssys)) != 0) {
                DSError(M_DS_VAR_NULL ": Xi0 variable pool is NULL", A_DS_ERROR);
                goto bail;
        }
        alpha = DSMatrixCopy(DSSSystemAlpha(ssys));
        DSMatrixApplyFunction(alpha, log10);
        ss = DSVariablePoolValuesAsVector(Xd0, false);
        DSMatrixApplyFunction(ss, log10);
        flux = DSMatrixByMultiplyingMatrix(DSSSystemGd(ssys), ss);
        if (DSVariablePoolNumberOfVariables(DSSSysXi(ssys)) != 0) {
                pool = DSVariablePoolAlloc();
                for (i=0; i < DSVariablePoolNumberOfVariables(DSSSysXi(ssys)); i++) {
                        name =  DSVariableName(DSVariablePoolAllVariables(DSSSysXi(ssys))[i]);
                        DSVariablePoolAddVariableWithName(pool, name);;
                        if (DSVariablePoolHasVariableWithName(Xi0, name) == false) {
                                DSError(M_DS_WRONG ": Variable Pool does not have independent variable", A_DS_ERROR);
                                DSMatrixFree(flux);
                                flux = NULL;
                                goto bail;
                        }
                        value = DSVariableValue(DSVariablePoolVariableWithName(Xi0, name));
                        DSVariablePoolSetValueForVariableWithName(pool, name, value);
                }
                Xi = DSVariablePoolValuesAsVector(pool, false);
                DSMatrixApplyFunction(Xi, log10);
                gi0 = DSMatrixByMultiplyingMatrix(DSSSystemGi(ssys), Xi);
                DSMatrixAddByMatrix(flux, gi0);
                DSMatrixFree(Xi);
                DSMatrixFree(gi0);
        }
        DSMatrixAddByMatrix(flux, alpha);
bail:
        if (alpha != NULL)
                DSMatrixFree(alpha);
        if (ss != NULL)
                DSMatrixFree(ss);
        if (pool != NULL)
                DSVariablePoolFree(pool);
        return flux;
}

extern DSMatrix * DSSSystemSteadyStateFlux(const DSSSystem *ssys, const DSVariablePool *Xi0)
{
        DSMatrix * flux = NULL, *Xi = NULL, *ss=NULL;
        DSMatrix * gi0= NULL, *alpha = NULL;
        DSUInteger i, index1;
        DSVariablePool *pool = NULL;
        const char *name;
        double value;
        if (ssys == NULL) {
                DSError(M_DS_SSYS_NULL, A_DS_ERROR);
                goto bail;
        }
        if (Xi0 == NULL && DSVariablePoolNumberOfVariables(DSSSysXi(ssys)) != 0) {
                DSError(M_DS_VAR_NULL ": Xi0 variable pool is NULL", A_DS_ERROR);
                goto bail;
        }
        if (DSSSystemHasSolution(ssys) == false)
                goto bail;
        alpha = DSMatrixCopy(DSSSystemAlpha(ssys));
        DSMatrixApplyFunction(alpha, log10);
        ss = DSSSystemSteadyStateValues(ssys, Xi0);
        flux = DSMatrixByMultiplyingMatrix(DSSSystemGd(ssys), ss);
        if (DSVariablePoolNumberOfVariables(DSSSysXi(ssys)) != 0) {
                pool = DSVariablePoolAlloc();
                for (i=0; i < DSVariablePoolNumberOfVariables(DSSSysXi(ssys)); i++) {
                        name =  DSVariableName(DSVariablePoolAllVariables(DSSSysXi(ssys))[i]);
                        DSVariablePoolAddVariableWithName(pool, name);;
                        if (DSVariablePoolHasVariableWithName(Xi0, name) == false) {
                                DSError(M_DS_WRONG ": Variable Pool does not have independent variable", A_DS_ERROR);
                                DSMatrixFree(flux);
                                flux = NULL;
                                DSVariablePoolFree(pool);
                                goto bail;
                        }
                        value = DSVariableValue(DSVariablePoolVariableWithName(Xi0, name));
                        DSVariablePoolSetValueForVariableWithName(pool, name, value);
                }
                Xi = DSVariablePoolValuesAsVector(pool, false);
                DSMatrixApplyFunction(Xi, log10);
                gi0 = DSMatrixByMultiplyingMatrix(DSSSystemGi(ssys), Xi);
                DSMatrixAddByMatrix(flux, gi0);
                DSMatrixFree(Xi);
                DSMatrixFree(gi0);
        }
        DSMatrixAddByMatrix(flux, alpha);
//        if (ssys->fluxDictionary != NULL) {
//                for (i = 0; i < DSVariablePoolNumberOfVariables(DSSSystemXd(ssys)); i++) {
//                        name = DSVariableName(DSVariablePoolVariableAtIndex(DSSSystemXd(ssys), i));
//                        altName = DSDictionaryValueForName(ssys->fluxDictionary, name);
//                        if (altName == NULL)
//                                continue;
//                        index1 = DSVariablePoolIndexOfVariableWithName(DSSSystemXd(ssys), name);
//                        index2 = DSVariablePoolIndexOfVariableWithName(DSSSystemXd(ssys), altName);
//                        DSMatrixSetDoubleValue(flux, index1, 0,
//                                               DSMatrixDoubleValue(flux, index2, 0));
//                }
//        }
        if (ssys->fluxDictionary != NULL) {
//                DSMatrixFree(flux);
                for (i = 0; i < DSVariablePoolNumberOfVariables(DSSSystemXd(ssys)); i++) {
                        name = DSVariableName(DSVariablePoolVariableAtIndex(DSSSystemXd(ssys), i));
                        DSVariablePoolAddVariableWithName(pool, name);
                        DSVariablePoolSetValueForVariableWithName(pool,
                                                                  name,
                                                                  pow(10, DSMatrixDoubleValue(ss, i, 0)));
                }
                for (i = 0; i < DSVariablePoolNumberOfVariables(DSSSystemXd(ssys)); i++) {
                        name = DSVariableName(DSVariablePoolVariableAtIndex(DSSSystemXd(ssys), i));
                        index1 = DSVariablePoolIndexOfVariableWithName(DSSSystemXd(ssys), name);
                        DSMatrixSetDoubleValue(flux, index1, 0,
                                               log10(DSExpressionEvaluateWithVariablePool(DSDictionaryValueForName(ssys->fluxDictionary,name), pool)));
                }


//                flux = DSSSystemSteadyStateValues(ssys, Xi0);
        }
bail:
        if (alpha != NULL)
                DSMatrixFree(alpha);
        if (ss != NULL)
                DSMatrixFree(ss);
        if (pool != NULL)
                DSVariablePoolFree(pool);
//    printf("Reporting from function DSSSystemSteadyStateFlux. The steady state flux matrix is: \n");
//    DSMatrixPrint(flux);
        return flux;
}

extern DSMatrix * DSSSystemSteadyStateFluxForConservedVariables(const DSSSystem *ssys, const DSVariablePool *Xi0)

    // This function is used to add fluxes from Xd_a_c (a subset of Xd_t in gma system). Function DSSSystemSteadyStateFlux is used as template.

{
    DSMatrix * flux = NULL, *Xi = NULL, *ss=NULL;
    DSMatrix * gi0= NULL, *alpha = NULL;
    DSUInteger i, index1;
    DSVariablePool *pool = NULL;
    const char *name;
    double value;
    if (ssys == NULL) {
        DSError(M_DS_SSYS_NULL, A_DS_ERROR);
        goto bail;
    }
    if (Xi0 == NULL && DSVariablePoolNumberOfVariables(DSSSysXi(ssys)) != 0) {
        DSError(M_DS_VAR_NULL ": Xi0 variable pool is NULL", A_DS_ERROR);
        goto bail;
    }
    if (DSSSystemHasSolution(ssys) == false || DSSSystemIsUnstable(ssys) == false)
        goto bail;
    alpha = DSMatrixCopy(DSSSystemAlpha(ssys));
    DSMatrixApplyFunction(alpha, log10);
    ss = DSSSystemSteadyStateValues(ssys, Xi0);
    flux = DSMatrixByMultiplyingMatrix(DSSSystemGd(ssys), ss);
    if (DSVariablePoolNumberOfVariables(DSSSysXi(ssys)) != 0) {
        pool = DSVariablePoolAlloc();
        for (i=0; i < DSVariablePoolNumberOfVariables(DSSSysXi(ssys)); i++) {
            name =  DSVariableName(DSVariablePoolAllVariables(DSSSysXi(ssys))[i]);
            DSVariablePoolAddVariableWithName(pool, name);;
            if (DSVariablePoolHasVariableWithName(Xi0, name) == false) {
                DSError(M_DS_WRONG ": Variable Pool does not have independent variable", A_DS_ERROR);
                DSMatrixFree(flux);
                flux = NULL;
                DSVariablePoolFree(pool);
                goto bail;
            }
            value = DSVariableValue(DSVariablePoolVariableWithName(Xi0, name));
            DSVariablePoolSetValueForVariableWithName(pool, name, value);
        }
        Xi = DSVariablePoolValuesAsVector(pool, false);
        DSMatrixApplyFunction(Xi, log10);
        gi0 = DSMatrixByMultiplyingMatrix(DSSSystemGi(ssys), Xi);
        DSMatrixAddByMatrix(flux, gi0);
        DSMatrixFree(Xi);
        DSMatrixFree(gi0);
    }
    DSMatrixAddByMatrix(flux, alpha);
    if (ssys->fluxDictionary != NULL) {
        for (i = 0; i < DSVariablePoolNumberOfVariables(DSSSystemXd(ssys)); i++) {
            name = DSVariableName(DSVariablePoolVariableAtIndex(DSSSystemXd(ssys), i));
            DSVariablePoolAddVariableWithName(pool, name);
            DSVariablePoolSetValueForVariableWithName(pool,
                                                      name,
                                                      pow(10, DSMatrixDoubleValue(ss, i, 0)));
        }
        for (i = 0; i < DSVariablePoolNumberOfVariables(DSSSystemXd(ssys)); i++) {
            name = DSVariableName(DSVariablePoolVariableAtIndex(DSSSystemXd(ssys), i));
            index1 = DSVariablePoolIndexOfVariableWithName(DSSSystemXd(ssys), name);
            DSMatrixSetDoubleValue(flux, index1, 0,
                                   log10(DSExpressionEvaluateWithVariablePool(DSDictionaryValueForName(ssys->fluxDictionary,name), pool)));
        }
    }
bail:
    if (alpha != NULL)
        DSMatrixFree(alpha);
    if (ss != NULL)
        DSMatrixFree(ss);
    if (pool != NULL)
        DSVariablePoolFree(pool);
    return flux;
}


extern DSMatrix * DSuSSystemSteadyStateFlux(const DSSSystem *ssys, const DSVariablePool *Xi0)
{
    DSMatrix * flux = NULL, *Xi = NULL, *ss=NULL;
    DSMatrix * gi0= NULL, *alpha = NULL;
    DSUInteger i, index1;
    DSVariablePool *pool = NULL;
    const DSSSystem * ssys_no_algebraic;
    
    const char *name;
    double value;
    if (ssys == NULL) {
        DSError(M_DS_SSYS_NULL, A_DS_ERROR);
        goto bail;
    }
    if (Xi0 == NULL && DSVariablePoolNumberOfVariables(DSSSysXi(ssys)) != 0) {
        DSError(M_DS_VAR_NULL ": Xi0 variable pool is NULL", A_DS_ERROR);
        goto bail;
    }
    if (DSSSystemIsUnstable(ssys) == false)
        goto bail;
    
    if ( DSVariablePoolNumberOfVariables(DSSSystemXd_t( ssys)) != 0 ){
        ssys_no_algebraic = DSSSystemByRemovingAlgebraicConstraints(ssys);
    } else {
        ssys_no_algebraic = ssys;
    }
    alpha = DSMatrixCopy(DSSSystemAlpha(ssys_no_algebraic));
    DSMatrixApplyFunction(alpha, log10);
    ss = DSuSSystemSteadyStateValues(ssys_no_algebraic, Xi0);
    flux = DSMatrixByMultiplyingMatrix(DSSSystemGd(ssys_no_algebraic), ss);
    if (DSVariablePoolNumberOfVariables(DSSSysXi(ssys_no_algebraic)) != 0) {
        pool = DSVariablePoolAlloc();
        for (i=0; i < DSVariablePoolNumberOfVariables(DSSSysXi(ssys_no_algebraic)); i++) {
            name =  DSVariableName(DSVariablePoolAllVariables(DSSSysXi(ssys_no_algebraic))[i]);
            DSVariablePoolAddVariableWithName(pool, name);;
            if (DSVariablePoolHasVariableWithName(Xi0, name) == false) {
                DSError(M_DS_WRONG ": Variable Pool does not have independent variable", A_DS_ERROR);
                DSMatrixFree(flux);
                flux = NULL;
                DSVariablePoolFree(pool);
                goto bail;
            }
            value = DSVariableValue(DSVariablePoolVariableWithName(Xi0, name));
            DSVariablePoolSetValueForVariableWithName(pool, name, value);
        }
        Xi = DSVariablePoolValuesAsVector(pool, false);
        DSMatrixApplyFunction(Xi, log10);
        gi0 = DSMatrixByMultiplyingMatrix(DSSSystemGi(ssys_no_algebraic), Xi);
        DSMatrixAddByMatrix(flux, gi0);
        DSMatrixFree(Xi);
        DSMatrixFree(gi0);
    }
    DSMatrixAddByMatrix(flux, alpha);
    if (ssys_no_algebraic->fluxDictionary != NULL) {
        for (i = 0; i < DSVariablePoolNumberOfVariables(DSSSystemXd(ssys_no_algebraic)); i++) {
            name = DSVariableName(DSVariablePoolVariableAtIndex(DSSSystemXd(ssys_no_algebraic), i));
            DSVariablePoolAddVariableWithName(pool, name);
            DSVariablePoolSetValueForVariableWithName(pool,
                                                      name,
                                                      pow(10, DSMatrixDoubleValue(ss, i, 0)));
        }
        for (i = 0; i < DSVariablePoolNumberOfVariables(DSSSystemXd(ssys_no_algebraic)); i++) {
            name = DSVariableName(DSVariablePoolVariableAtIndex(DSSSystemXd(ssys_no_algebraic), i));
            index1 = DSVariablePoolIndexOfVariableWithName(DSSSystemXd(ssys_no_algebraic), name);
            DSMatrixSetDoubleValue(flux, index1, 0,
                                   log10(DSExpressionEvaluateWithVariablePool(DSDictionaryValueForName(ssys_no_algebraic->fluxDictionary,name), pool)));
        }
    }
bail:
    if (alpha != NULL)
        DSMatrixFree(alpha);
    if (ss != NULL)
        DSMatrixFree(ss);
    if (pool != NULL)
        DSVariablePoolFree(pool);
    if ( DSVariablePoolNumberOfVariables(DSSSystemXd_t( ssys)) != 0 )
        DSSSystemFree((DSSSystem *)ssys_no_algebraic);
    
    return flux;
}

static void dsSSystemRouthArrayProcessZeroRoots(DSMatrix * routhMatrix, const DSUInteger row, const double threshold)
{
        DSUInteger i, order, newRowSize;
        bool rowEmpty;
        double value;
        if (routhMatrix == NULL) {
                DSError(M_DS_MAT_NULL, A_DS_ERROR);
                goto bail;
        }
        if (row >= DSMatrixRows(routhMatrix)) {
                DSError(M_DS_MAT_OUTOFBOUNDS, A_DS_ERROR);
                goto bail;
        }
        rowEmpty = true;
        for (i = 0; i < DSMatrixColumns(routhMatrix); i++) {
                if (fabs(DSMatrixDoubleValue(routhMatrix, row, i)) >= threshold) {
                        rowEmpty = false;
                }
        }
        if (rowEmpty == false) {
                DSMatrixSetDoubleValue(routhMatrix, row, 0, threshold);
                goto bail;
        }
        order = DSMatrixRows(routhMatrix)-row;
        newRowSize = order/2 + order % 2;
        for (i = 0; i < newRowSize; i++) {
                value = order*DSMatrixDoubleValue(routhMatrix, row-1, i);
                DSMatrixSetDoubleValue(routhMatrix, row, i, value);
                order -= 2;
        }
bail:
        return;
}

extern DSMatrix * DSSSystemRouthArrayForPoolTurnover(const DSSSystem *ssys, const DSMatrix * F, bool * hasImaginaryRoots)
{
        DSSSystem * reduced = NULL;
        DSMatrix * FA = NULL;
        DSMatrix * routhArray = NULL;
        DSMatrix * phi = NULL;
        DSMatrix * routhMatrix = NULL;
        DSMatrix *Ad;
        DSUInteger i, j;
        double value;
        double threshold = 1e-8;
        if (ssys == NULL) {
                DSError(M_DS_SSYS_NULL, A_DS_ERROR);
                goto bail;
        }
        if (F == NULL) {
                DSError(M_DS_MAT_NULL ": F matrix is NULL", A_DS_ERROR);
                goto bail;
        }
        if (DSVariablePoolNumberOfVariables(DSSSystemXd(ssys)) > 0) {
                reduced = DSSSystemByRemovingAlgebraicConstraints(ssys);
                ssys = reduced;
        }
        if (DSMatrixRows(F) != DSSSystemNumberOfEquations(ssys)) {
                DSError(M_DS_MAT_OUTOFBOUNDS, A_DS_ERROR);
                goto bail;
        }
        Ad =DSSSystemAd(ssys);
        FA = DSMatrixByMultiplyingMatrix(F, Ad);
        DSMatrixFree(Ad);
        phi = DSMatrixCharacteristicPolynomialCoefficients(FA);
        DSMatrixFree(FA);
        routhMatrix = DSMatrixCalloc(DSMatrixColumns(phi), DSMatrixColumns(phi));
//        routhArray = DSMatrixCalloc(DSMatrixRows(routhMatrix), 1);
        /* Make first row of routh matrix */
        for (i = 0; i < DSMatrixColumns(routhMatrix); i++) {
                value = 0.0f;
                if (2*i < DSMatrixColumns(phi))
                        value = DSMatrixDoubleValue(phi, 0, 2*i);
                DSMatrixSetDoubleValue(routhMatrix, 0, i, value);
        }
        for (i = 0; i < DSMatrixColumns(routhMatrix); i++) {
                value = 0.0f;
                if ((2*i)+1 < DSMatrixColumns(phi))
                        value = DSMatrixDoubleValue(phi, 0, (2*i)+1);
                DSMatrixSetDoubleValue(routhMatrix, 1, i, value);
        }
        if (hasImaginaryRoots != NULL) {
                *hasImaginaryRoots = false;
        }
        for (i = 2; i < DSMatrixRows(routhMatrix); i++) {
                for (j = 0; j < DSMatrixColumns(routhMatrix); j++) {
                        if (j == DSMatrixColumns(routhMatrix)-1) {
                                DSMatrixSetDoubleValue(routhMatrix, i, j, 0.0f);
                                continue;
                        }
                        value = DSMatrixDoubleValue(routhMatrix, i-1, 0);
                        value = (value*DSMatrixDoubleValue(routhMatrix, i-2, j+1)-DSMatrixDoubleValue(routhMatrix, i-2, 0)*DSMatrixDoubleValue(routhMatrix, i-1, j+1))/value;
//                        if (fabs(value) < threshold)
//                                value = 0.f;
                        DSMatrixSetDoubleValue(routhMatrix, i, j, value);
                        if (value == 0.f && j == 0) {
                                dsSSystemRouthArrayProcessZeroRoots(routhMatrix, i, threshold);
                                if (hasImaginaryRoots != NULL) {
                                        *hasImaginaryRoots = true;
                                }
                        }
//                        if (j == 0) {
//                                DSMatrixSetDoubleValue(routhArray, i, 0, DSMatrixDoubleValue(routhMatrix, i, j));
//                        }
                }
        }
        routhArray = DSMatrixSubMatrixIncludingColumnList(routhMatrix, 1, 0);
        DSMatrixFree(routhMatrix);
        DSMatrixFree(phi);
bail:
        if (reduced != NULL) {
                DSSSystemFree(reduced);
        }
        return routhArray;
}

extern DSMatrix * DSSSystemRouthArrayForSteadyState(const DSSSystem *ssys,
                                                    const DSVariablePool *Xd0,
                                                    const DSVariablePool *Xi0)
{
        DSMatrix * routhArray = NULL;
        DSMatrix * steadyState = NULL;
        DSMatrix * flux = NULL;
        DSMatrix * F = NULL;
        DSUInteger i;
        
        if (ssys == NULL) {
                DSError(M_DS_SSYS_NULL, A_DS_ERROR);
                goto bail;
        }
        if (Xd0 == NULL && DSVariablePoolNumberOfVariables(DSSSysXd(ssys)) != 0) {
                DSError(M_DS_VAR_NULL ": Xd0 variable pool is NULL", A_DS_ERROR);
                goto bail;
        }
        if (Xi0 == NULL && DSVariablePoolNumberOfVariables(DSSSysXi(ssys)) != 0) {
                DSError(M_DS_VAR_NULL ": Xi0 variable pool is NULL", A_DS_ERROR);
                goto bail;
        }
        steadyState = DSVariablePoolValuesAsVector(Xd0, false);
        flux = DSSSystemSteadyStateFluxForDependentVariables(ssys, Xd0, Xi0);
        F = DSMatrixIdentity(DSMatrixRows(flux));
        for (i = 0; i < DSMatrixColumns(F); i++) {
                DSMatrixSetDoubleValue(F,
                                       i,
                                       i,
                                       pow(10, DSMatrixDoubleValue(flux, i, 0))/pow(10,DSMatrixDoubleValue(steadyState, i, 0)));
        }
        routhArray = DSSSystemRouthArrayForPoolTurnover(ssys, F, NULL);
        DSMatrixFree(F);
bail:
        return routhArray;
}

extern DSMatrix * DSSSystemRouthArray(const DSSSystem *ssys, const DSVariablePool *Xi0, bool * hasImaginaryRoots)
{
        DSMatrix * routhArray = NULL;
        DSMatrix * steadyState = NULL;
        DSMatrix * flux = NULL;
        DSMatrix * F = NULL;
        DSUInteger i;
        
        if (ssys == NULL) {
                DSError(M_DS_SSYS_NULL, A_DS_ERROR);
                goto bail;
        }
        if (Xi0 == NULL && DSVariablePoolNumberOfVariables(DSSSysXi(ssys)) != 0) {
                DSError(M_DS_VAR_NULL ": Xi0 variable pool is NULL", A_DS_ERROR);
                goto bail;
        }
        if (DSSSystemHasSolution(ssys) == false)
                goto bail;
        steadyState = DSSSystemSteadyStateValues(ssys, Xi0);
        flux = DSSSystemSteadyStateFlux(ssys, Xi0);
        F = DSMatrixIdentity(DSMatrixRows(flux));
        for (i = 0; i < DSMatrixColumns(F); i++) {
                DSMatrixSetDoubleValue(F,
                                       i,
                                       i,
                                       pow(10, DSMatrixDoubleValue(flux, i, 0))/pow(10,DSMatrixDoubleValue(steadyState, i, 0)));
        }
        routhArray = DSSSystemRouthArrayForPoolTurnover(ssys, F, hasImaginaryRoots);
        DSMatrixFree(F);
        DSMatrixFree(steadyState);
        DSMatrixFree(flux);
bail:
        return routhArray;
}

extern DSUInteger DSSSystemNumberOfPositiveRootsForRouthArray(const DSMatrix *routhArray)
{
        DSUInteger positiveRoots = 0;
        DSUInteger i, length;
        double sign, baseSign = 1;

        if (routhArray == NULL) {
                DSError(M_DS_MAT_NULL ": Routh Array Matrix is Null", A_DS_ERROR);
                goto bail;
        }
        
        length = DSMatrixRows(routhArray);
        baseSign = (DSMatrixDoubleValue(routhArray, 0, 0) > 0) ? 1. : -1.;
        for (i = 1; i < length; i++) {
                sign = (DSMatrixDoubleValue(routhArray, i, 0) > 0) ? 1. : -1.;
                if (sign*baseSign < 0)
                        positiveRoots++;
                baseSign = sign;
        }
bail:
        return positiveRoots;
}

extern DSUInteger DSSSystemPositiveRootsForSteadyStateAndFlux(const DSSSystem *ssys,
                                                              const DSVariablePool *Xd0,
                                                              const DSVariablePool *Xi0,
                                                              const DSVariablePool *flux0)
{
        DSMatrix *steadyState, *flux, *F, * routhArray = NULL;
        char flux_name[100];
        DSUInteger i, positiveRoots = 0;
        if (ssys == NULL) {
                DSError(M_DS_SSYS_NULL, A_DS_ERROR);
                goto bail;
        }
        if (Xd0 == NULL && DSVariablePoolNumberOfVariables(DSSSysXd(ssys)) != 0) {
                DSError(M_DS_VAR_NULL ": Xd0 variable pool is NULL", A_DS_ERROR);
                goto bail;
        }
        if (Xi0 == NULL && DSVariablePoolNumberOfVariables(DSSSysXi(ssys)) != 0) {
                DSError(M_DS_VAR_NULL ": Xi0 variable pool is NULL", A_DS_ERROR);
                goto bail;
        }
        if (flux0 == NULL && DSVariablePoolNumberOfVariables(DSSSysXd(ssys)) != 0) {
                DSError(M_DS_VAR_NULL ": Flux variable pool is NULL", A_DS_ERROR);
                goto bail;
        }
        steadyState = DSVariablePoolValuesAsVector(Xd0, false);
        flux = DSMatrixAlloc(DSVariablePoolNumberOfVariables(DSSSystemXd(ssys)), 1);
        for (i = 0; i < DSVariablePoolNumberOfVariables(DSSSystemXd(ssys)); i++) {
                sprintf(flux_name, "V_%s", DSVariableName(DSVariablePoolVariableAtIndex(DSSSysXd(ssys), i)));
                DSMatrixSetDoubleValue(flux, i, 0, DSVariablePoolValueForVariableWithName(flux0, flux_name));
        }
        F = DSMatrixIdentity(DSMatrixRows(flux));
        for (i = 0; i < DSMatrixColumns(F); i++) {
                DSMatrixSetDoubleValue(F,
                                       i,
                                       i,
                                       pow(10, DSMatrixDoubleValue(flux, i, 0))/pow(10,DSMatrixDoubleValue(steadyState, i, 0)));
        }
        routhArray = DSSSystemRouthArrayForPoolTurnover(ssys, F, NULL);
        DSMatrixFree(F);
        DSMatrixFree(flux);
        DSMatrixFree(steadyState);
//        routhArray = DSSSystemRouthArrayForSteadyState(ssys, Xd0, Xi0);
        if (routhArray == NULL) {
                goto bail;
        }
        
        positiveRoots = DSSSystemNumberOfPositiveRootsForRouthArray(routhArray);
        DSMatrixFree(routhArray);
bail:
        return positiveRoots;
}

extern DSUInteger DSSSystemPositiveRootsForSteadyState(const DSSSystem *ssys,
                                                       const DSVariablePool *Xd0,
                                                       const DSVariablePool *Xi0)
{
        DSMatrix * routhArray = NULL;
        DSUInteger positiveRoots = 0;
        if (ssys == NULL) {
                DSError(M_DS_SSYS_NULL, A_DS_ERROR);
                goto bail;
        }
        if (Xd0 == NULL && DSVariablePoolNumberOfVariables(DSSSysXd(ssys)) != 0) {
                DSError(M_DS_VAR_NULL ": Xd0 variable pool is NULL", A_DS_ERROR);
                goto bail;
        }
        if (Xi0 == NULL && DSVariablePoolNumberOfVariables(DSSSysXi(ssys)) != 0) {
                DSError(M_DS_VAR_NULL ": Xi0 variable pool is NULL", A_DS_ERROR);
                goto bail;
        }
        routhArray = DSSSystemRouthArrayForSteadyState(ssys, Xd0, Xi0);
        if (routhArray == NULL) {
                goto bail;
        }
        positiveRoots = DSSSystemNumberOfPositiveRootsForRouthArray(routhArray);
        DSMatrixFree(routhArray);
bail:
        return positiveRoots;
}

extern DSUInteger DSSSystemPositiveRoots(const DSSSystem *ssys, const DSVariablePool *Xi0, bool * hasImaginaryRoots)
{
        DSMatrix * routhArray = NULL;
        DSUInteger positiveRoots = 0;
        if (ssys == NULL) {
                DSError(M_DS_SSYS_NULL, A_DS_ERROR);
                goto bail;
        }
        if (Xi0 == NULL && DSVariablePoolNumberOfVariables(DSSSysXi(ssys)) != 0) {
                DSError(M_DS_VAR_NULL ": Xi0 variable pool is NULL", A_DS_ERROR);
                goto bail;
        }
        if (DSSSystemHasSolution(ssys) == false)
                goto bail;
        routhArray = DSSSystemRouthArray(ssys, Xi0, hasImaginaryRoots);
        if (routhArray == NULL) {
                goto bail;
        }
        
        positiveRoots = DSSSystemNumberOfPositiveRootsForRouthArray(routhArray);
        if (positiveRoots == 3) {
                printf("%i\n", DSSSystemCharacteristicEquationCoefficientsNumberSignChanges(ssys, Xi0));
        }
        DSMatrixFree(routhArray);
bail:
        return positiveRoots;
}

extern DSUInteger DSSSystemRouthIndex(const DSSSystem *ssys, const DSVariablePool *Xi0)
{
        DSMatrix * routhArray = NULL;
        DSUInteger routhIndex = 0;
        DSUInteger i, length;
        double value, baseSign = 1;
        if (ssys == NULL) {
                DSError(M_DS_SSYS_NULL, A_DS_ERROR);
                goto bail;
        }
        if (Xi0 == NULL && DSVariablePoolNumberOfVariables(DSSSysXi(ssys)) != 0) {
                DSError(M_DS_VAR_NULL ": Xi0 variable pool is NULL", A_DS_ERROR);
                goto bail;
        }
        if (DSSSystemHasSolution(ssys) == false)
                goto bail;
        routhArray = DSSSystemRouthArray(ssys, Xi0, NULL);
        if (routhArray == NULL) {
                goto bail;
        }
                
        length = DSMatrixRows(routhArray);
        baseSign = (DSMatrixDoubleValue(routhArray, 0, 0) > 0) ? 1. : -1.;
        for (i = 0; i < length; i++) {
                value = DSMatrixDoubleValue(routhArray, i, 0);
                value *= baseSign;
                if (value < 0)
                        routhIndex += pow(2, i);
        }
        DSMatrixFree(routhArray);
bail:
        return routhIndex;
}

extern DSUInteger DSSSystemCharacteristicEquationCoefficientIndex(const DSSSystem *ssys, const DSVariablePool *Xi0)
{
        DSMatrix * coefficientArray = NULL;
        DSMatrix * F;
        DSMatrix * FA;
        DSMatrix * steadyState, * flux;
        DSUInteger Index = 0;
        DSUInteger i, length;
        double value, baseSign = 1;
        if (ssys == NULL) {
                DSError(M_DS_SSYS_NULL, A_DS_ERROR);
                goto bail;
        }
        if (Xi0 == NULL && DSVariablePoolNumberOfVariables(DSSSysXi(ssys)) != 0) {
                DSError(M_DS_VAR_NULL ": Xi0 variable pool is NULL", A_DS_ERROR);
                goto bail;
        }
        if (DSSSystemHasSolution(ssys) == false)
                goto bail;
        steadyState = DSSSystemSteadyStateValues(ssys, Xi0);
        flux = DSSSystemSteadyStateFlux(ssys, Xi0);
        F = DSMatrixIdentity(DSMatrixRows(flux));
        for (i = 0; i < DSMatrixColumns(F); i++) {
                DSMatrixSetDoubleValue(F,
                                       i,
                                       i,
                                       pow(10, DSMatrixDoubleValue(flux, i, 0))/pow(10,DSMatrixDoubleValue(steadyState, i, 0)));
        }
        FA = DSMatrixByMultiplyingMatrix(F, DSSSystemAd(ssys));
        coefficientArray = DSMatrixCharacteristicPolynomialCoefficients(FA);
        DSMatrixFree(steadyState);
        DSMatrixFree(flux);
        DSMatrixFree(F);
        DSMatrixFree(FA);
        if (coefficientArray == NULL) {
                goto bail;
        }
        
        length = DSMatrixRows(coefficientArray);
        baseSign = (DSMatrixDoubleValue(coefficientArray, 0, 0) > 0) ? 1. : -1.;
        for (i = 0; i < length; i++) {
                value = DSMatrixDoubleValue(coefficientArray, i, 0);
                value *= baseSign;
                if (value < 0)
                        Index += pow(2, i);
        }
        DSMatrixFree(coefficientArray);
bail:
        return Index;
}

extern DSUInteger DSSSystemCharacteristicEquationCoefficientsNumberSignChanges(const DSSSystem *ssys, const DSVariablePool *Xi0)
{
        DSMatrix * coefficientArray = NULL;
        DSMatrix * F;
        DSMatrix * FA;
        DSMatrix * steadyState, * flux;
        DSUInteger Index = 0;
        DSUInteger i, length;
        double value, baseSign = 1;
        if (ssys == NULL) {
                DSError(M_DS_SSYS_NULL, A_DS_ERROR);
                goto bail;
        }
        if (Xi0 == NULL && DSVariablePoolNumberOfVariables(DSSSysXi(ssys)) != 0) {
                DSError(M_DS_VAR_NULL ": Xi0 variable pool is NULL", A_DS_ERROR);
                goto bail;
        }
        if (DSSSystemHasSolution(ssys) == false)
                goto bail;
        steadyState = DSSSystemSteadyStateValues(ssys, Xi0);
        flux = DSSSystemSteadyStateFlux(ssys, Xi0);
        F = DSMatrixIdentity(DSMatrixRows(flux));
        for (i = 0; i < DSMatrixColumns(F); i++) {
                DSMatrixSetDoubleValue(F,
                                       i,
                                       i,
                                       pow(10, DSMatrixDoubleValue(flux, i, 0))/pow(10,DSMatrixDoubleValue(steadyState, i, 0)));
        }
        FA = DSMatrixByMultiplyingMatrix(F, DSSSystemAd(ssys));
        coefficientArray = DSMatrixCharacteristicPolynomialCoefficients(FA);
        DSMatrixFree(steadyState);
        DSMatrixFree(flux);
        DSMatrixFree(F);
        DSMatrixFree(FA);
        if (coefficientArray == NULL) {
                goto bail;
        }
        length = DSMatrixColumns(coefficientArray);
        baseSign = 1.;
        for (i = 1; i < length; i++) {
                value = DSMatrixDoubleValue(coefficientArray, 0, i);
                value *= baseSign;
                if (value < 0)
                        Index++;
                baseSign = (DSMatrixDoubleValue(coefficientArray, 0, i) > 0) ? 1. : -1.;
        }
        DSMatrixFree(coefficientArray);
bail:
        return Index;
}

extern double DSSSystemLogarithmicGain(const DSSSystem *ssys, const char *XdName, const char *XiName)
{
        double logGain = INFINITY;
        DSUInteger XdIndex = 0;
        DSUInteger XiIndex = 0;
        DSMatrix * L = NULL, * Ai = NULL;
        if (ssys == NULL) {
                DSError(M_DS_SSYS_NULL, A_DS_ERROR);
                goto bail;
        }
        if (XdName == NULL) {
                DSError(M_DS_NULL, A_DS_ERROR);
                goto bail;
        }
        if (XiName == NULL) {
                DSError(M_DS_NULL, A_DS_ERROR);
                goto bail;
        }
        if (DSVariablePoolHasVariableWithName(DSSSysXd(ssys), XdName) == false) {
                DSError(M_DS_WRONG, A_DS_ERROR);
                goto bail;
        } else {
                XdIndex=DSVariablePoolIndexOfVariableWithName(DSSSysXd(ssys), XdName);
        }
        if (DSVariablePoolHasVariableWithName(DSSSysXi(ssys), XiName) == false) {
                DSError(M_DS_WRONG, A_DS_ERROR);
                goto bail;
        } else {
                XiIndex = DSVariablePoolIndexOfVariableWithName(DSSSysXi(ssys), XiName);                
        }
        Ai = DSSSystemAi(ssys);
        if (Ai == NULL) {
                goto bail;
        }
        L = DSMatrixByMultiplyingMatrix(DSSSystemM(ssys), Ai);
        if (L == NULL) {
                DSError(M_DS_MAT_NULL, A_DS_ERROR);
                goto bail;
        }
        logGain = -DSMatrixDoubleValue(L, XdIndex, XiIndex);
        DSMatrixFree(L);
        DSMatrixFree(Ai);
bail:
        return logGain;
}

extern double DSuSSystemLogarithmicGain(const DSSSystem *ssys, const char *XdName, const char *XiName)
{
    double logGain = INFINITY;
    DSUInteger XdIndex = 0;
    DSUInteger XiIndex = 0;
    DSMatrix * L = NULL, * Ai = NULL;
    
    if (ssys == NULL) {
        DSError(M_DS_SSYS_NULL, A_DS_ERROR);
        goto bail;
    }
//    if ( DSVariablePoolNumberOfVariables(DSSSystemXd_t(ssys)) != 0 ){
//        ssys_no_algebraic = DSSSystemByRemovingAlgebraicConstraints(ssys);
//    } else {
//        ssys_no_algebraic = ssys;
//    }
    if (XdName == NULL) {
        DSError(M_DS_NULL, A_DS_ERROR);
        goto bail;
    }
    if (XiName == NULL) {
        DSError(M_DS_NULL, A_DS_ERROR);
        goto bail;
    }
    if (DSVariablePoolHasVariableWithName(DSSSysXd(ssys), XdName) == false) {
        DSError(M_DS_WRONG, A_DS_ERROR);
        goto bail;
    } else {
        XdIndex=DSVariablePoolIndexOfVariableWithName(DSSSysXd(ssys), XdName);
    }
    if (DSVariablePoolHasVariableWithName(DSSSysXi(ssys), XiName) == false) {
        DSError(M_DS_WRONG, A_DS_ERROR);
        goto bail;
    } else {
        XiIndex = DSVariablePoolIndexOfVariableWithName(DSSSysXi(ssys), XiName);
    }
    Ai = DSSSystemAi(ssys);
    if (Ai == NULL) {
        goto bail;
    }
    if (DSSSystemIsUnstable(ssys) == false ){
        goto bail;
    }
    
    DSMatrix *Ad = DSSSystemAd(ssys);
    const DSMatrix *pInverse = DSUnstableCaseGetSubSetPseudoInverse(ssys->Xd_t, NULL, ssys->Xd_b, Ad);
    L = DSMatrixByMultiplyingMatrix(pInverse, Ai);
    DSMatrixFree(Ad);
    DSMatrixFree((DSMatrix *)pInverse);
    
    if (L == NULL) {
        DSError(M_DS_MAT_NULL, A_DS_ERROR);
        goto bail;
    }
    
    logGain = -DSMatrixDoubleValue(L, XdIndex, XiIndex);
    DSMatrixFree(L);
    DSMatrixFree(Ai);

bail:
    return logGain;
}


#if defined (__APPLE__) && defined (__MACH__)
#pragma mark - Utility functions
#endif

static void dsSSystemPartitionSolutionMatrices(const DSSSystem * ssystem,
                                               const DSUInteger numberOfVariables,
                                               const DSUInteger * variablesToPartition,
                                               DSMatrix ** ADn,
                                               DSMatrix ** ADp,
                                               DSMatrix ** AIn,
                                               DSMatrix ** Bn,
                                               DSVariablePool ** yn,
                                               DSVariablePool ** yp)
{
        DSMatrix * tempMatrix, *Ad, *Ai, *B;
        DSUInteger i;
        char * name;
        if (ADn == NULL || ADp == NULL || AIn == NULL || Bn == NULL) {
                DSError(M_DS_NULL ": Matrix pointers to hold partitioned matrices cannot be null", A_DS_ERROR);
                goto bail;
        }
        *ADn = NULL;
        *ADp = NULL;
        *AIn = NULL;
        *Bn = NULL;
        if (ssystem == NULL) {
                DSError(M_DS_SSYS_NULL, A_DS_ERROR);
                goto bail;
        }
        if (variablesToPartition == NULL) {
                goto bail;
        }
        Ad = DSSSystemAd(ssystem);
        Ai = DSSSystemAi(ssystem);
        B = DSSSystemB(ssystem);
        tempMatrix = DSMatrixSubMatrixIncludingRows(Ad, numberOfVariables, variablesToPartition);
        *ADp = DSMatrixSubMatrixExcludingColumns(tempMatrix, numberOfVariables, variablesToPartition);
        *ADn = DSMatrixSubMatrixIncludingColumns(tempMatrix, numberOfVariables, variablesToPartition);
        DSMatrixFree(tempMatrix);
        *AIn = DSMatrixSubMatrixIncludingRows(Ai, numberOfVariables, variablesToPartition);
        *Bn = DSMatrixSubMatrixIncludingRows(B, numberOfVariables, variablesToPartition);
        *yn = DSVariablePoolAlloc();
        *yp = DSVariablePoolAlloc();
        DSMatrixFree(Ad);
        DSMatrixFree(Ai);
        DSMatrixFree(B);
        for (i = 0; i < numberOfVariables; i++) {
                name = DSVariableName(DSVariablePoolVariableAtIndex(DSSSystemXd(ssystem), variablesToPartition[i]));
                DSVariablePoolAddVariableWithName(*yn, name);
        }
        for (i = 0; i < DSVariablePoolNumberOfVariables(DSSSystemXd(ssystem)); i++) {
                name = DSVariableName(DSVariablePoolVariableAtIndex(DSSSystemXd(ssystem), i));
                if (DSVariablePoolHasVariableWithName(*yn, name) == false) {
                        DSVariablePoolAddVariableWithName(*yp, name);
                }
        }
bail:
        return;
}


static void dsSSystemSolutionOfPartitionedMatrices(const DSSSystem * ssystem,
                                                   const DSUInteger numberOfVariables,
                                                   const DSUInteger * partitionVariables,
                                                   DSMatrix ** LI,
                                                   DSMatrix **Lp,
                                                   DSMatrix **MBn,
                                                   DSVariablePool ** yn,
                                                   DSVariablePool ** yp)
{
        DSMatrix *ADn = NULL, * AIn = NULL, * ADp = NULL, * Bn = NULL, * Mn = NULL;
        if (LI == NULL || Lp == NULL || MBn == NULL) {
                DSError(M_DS_NULL ": Matrix pointers to hold partitioned matrices cannot be null", A_DS_ERROR);
                goto bail;
        }
        *LI = NULL;
        *Lp = NULL;
        *MBn = NULL;
        if (ssystem == NULL) {
                DSError(M_DS_SSYS_NULL, A_DS_ERROR);
                goto bail;
        }
        if (partitionVariables == NULL) {
                goto bail;
        }
        dsSSystemPartitionSolutionMatrices(ssystem, numberOfVariables, partitionVariables, &ADn, &ADp, &AIn, &Bn, yn, yp);
        if (ADn == NULL || ADp == NULL || AIn == NULL || Bn == NULL) {
                goto bail;
        }
        Mn = DSMatrixInverse(ADn);
        if (Mn == NULL) {
                goto bail;
        }
        *LI = DSMatrixByMultiplyingMatrix(Mn, AIn);
        *Lp = DSMatrixByMultiplyingMatrix(Mn, ADp);
        *MBn = DSMatrixByMultiplyingMatrix(Mn, Bn);
bail:
        if (ADn != NULL)
                DSMatrixFree(ADn);
        if (ADp != NULL)
                DSMatrixFree(ADp);
        if (AIn != NULL)
                DSMatrixFree(AIn);
        if (Bn != NULL)
                DSMatrixFree(Bn);
        if (Mn != NULL)
                DSMatrixFree(Mn);
        return;
}

static DSSSystem * dsSSystemForQuasiSteadyState(const DSSSystem * ssystem, DSUInteger numberOfVariables, const char ** variableNames)
{
        DSSSystem * newSSystem = NULL;
        DSMatrix *AD, * Gdp, *Gdn, * Hdp, *Hdn, *Gip, *Hip, * LI, *Lp, *MBn, *alpha, *beta, * tempMatrix;
        DSUInteger *variablesToPartition = NULL;
        DSVariablePool *yn, * Xd, *Xd_t, *Xd_a;
        DSUInteger i;
        char * name;
        if (ssystem == NULL) {
                DSError(M_DS_SSYS_NULL, A_DS_ERROR);
                goto bail;
        }
        if (numberOfVariables == 0) {
                DSError(M_DS_WRONG, A_DS_ERROR);
                goto bail;
        }
        if (variableNames == NULL) {
                DSError(M_DS_NULL, A_DS_ERROR);
                goto bail;
        }
        variablesToPartition = DSSecureMalloc(sizeof(DSUInteger)*numberOfVariables);
        Xd_t = DSVariablePoolAlloc();
        Xd_a = DSVariablePoolAlloc();
        for (i = 0; i < numberOfVariables; i++) {
                variablesToPartition[i] = DSVariablePoolIndexOfVariableWithName(DSSSystemXd(ssystem), variableNames[i]);
        }
        AD = DSSSystemAd(ssystem);
        Gdp = DSMatrixSubMatrixExcludingRowsAndColumns(DSSSystemGd(ssystem), numberOfVariables, numberOfVariables, variablesToPartition, variablesToPartition);
        Hdp = DSMatrixSubMatrixExcludingRowsAndColumns(DSSSystemHd(ssystem), numberOfVariables, numberOfVariables, variablesToPartition, variablesToPartition);
        tempMatrix = DSMatrixSubMatrixExcludingRows(DSSSystemGd(ssystem), numberOfVariables, variablesToPartition);
        Gdn = DSMatrixSubMatrixIncludingColumns(tempMatrix, numberOfVariables, variablesToPartition);
        DSMatrixFree(tempMatrix);
        tempMatrix = DSMatrixSubMatrixExcludingRows(DSSSystemHd(ssystem), numberOfVariables, variablesToPartition);
        Hdn = DSMatrixSubMatrixIncludingColumns(tempMatrix, numberOfVariables, variablesToPartition);
        DSMatrixFree(tempMatrix);
        Gip = DSMatrixSubMatrixExcludingRows(DSSSystemGi(ssystem), numberOfVariables, variablesToPartition);
        Hip = DSMatrixSubMatrixExcludingRows(DSSSystemHi(ssystem), numberOfVariables, variablesToPartition);
        alpha = DSMatrixSubMatrixExcludingRows(DSSSystemAlpha(ssystem), numberOfVariables, variablesToPartition);
        beta = DSMatrixSubMatrixExcludingRows(DSSSystemBeta(ssystem), numberOfVariables, variablesToPartition);
        dsSSystemSolutionOfPartitionedMatrices(ssystem, numberOfVariables, variablesToPartition, &LI, &Lp, &MBn, &yn, &Xd);
        tempMatrix = DSMatrixByMultiplyingMatrix(Gdn, Lp);
        DSMatrixSubstractByMatrix(Gdp, tempMatrix);
        DSMatrixFree(tempMatrix);
        tempMatrix = DSMatrixByMultiplyingMatrix(Hdn, Lp);
        DSMatrixSubstractByMatrix(Hdp, tempMatrix);
        DSMatrixFree(tempMatrix);
        tempMatrix = DSMatrixByMultiplyingMatrix(Gdn, LI);
        DSMatrixSubstractByMatrix(Gip, tempMatrix);
        DSMatrixFree(tempMatrix);
        tempMatrix = DSMatrixByMultiplyingMatrix(Hdn, LI);
        DSMatrixSubstractByMatrix(Hip, tempMatrix);
        DSMatrixFree(tempMatrix);
        tempMatrix = DSMatrixByMultiplyingMatrix(Gdn, MBn);
        for (i = 0; i < DSMatrixRows(tempMatrix); i++) {
                DSMatrixSetDoubleValue(tempMatrix, i, 0, pow(10., DSMatrixDoubleValue(tempMatrix, i, 0)));
        }
        DSMatrixAddByMatrix(alpha, tempMatrix);
        DSMatrixFree(tempMatrix);
        tempMatrix = DSMatrixByMultiplyingMatrix(Hdn, MBn);
        for (i = 0; i < DSMatrixRows(tempMatrix); i++) {
                DSMatrixSetDoubleValue(tempMatrix, i, 0, pow(10., DSMatrixDoubleValue(tempMatrix, i, 0)));
        }
        DSMatrixAddByMatrix(beta, tempMatrix);
        DSMatrixFree(tempMatrix);
        newSSystem = DSSSystemAlloc();
        DSSSysAlpha(newSSystem) = alpha;
        DSSSysGd(newSSystem) = Gdp;
        DSSSysGi(newSSystem) = Gip;
        DSSSysBeta(newSSystem) = beta;
        DSSSysHd(newSSystem) = Hdp;
        DSSSysHi(newSSystem) = Hip;
        for (i = 0; i < DSVariablePoolNumberOfVariables(DSSSysXd_t(ssystem)); i++) {
                name = DSVariableName(DSVariablePoolVariableAtIndex(DSSSysXd_t(ssystem), i));
                if (DSVariablePoolHasVariableWithName(Xd, name) == true)
                        DSVariablePoolAddVariableWithName(Xd_t, name);
        }
        for (i = 0; i < DSVariablePoolNumberOfVariables(DSSSysXd_a(ssystem)); i++) {
                name = DSVariableName(DSVariablePoolVariableAtIndex(DSSSysXd_a(ssystem), i));
                if (DSVariablePoolHasVariableWithName(Xd, name) == true)
                        DSVariablePoolAddVariableWithName(Xd_a, name);
        }
        DSSSysXd_a(newSSystem) = Xd_a;
        DSSSysXd_t(newSSystem) = Xd_t;
        DSSSysXd(newSSystem) = Xd;
        DSSSysXi(newSSystem) = DSVariablePoolCopy(DSSSysXi(ssystem));
        dsSSystemSolveEquations(newSSystem);
        DSMatrixFree(Gdn);
        DSMatrixFree(Hdn);
        DSMatrixFree(LI);
        DSMatrixFree(Lp);
        DSMatrixFree(MBn);
        DSVariablePoolFree(yn);
        DSSecureFree(variablesToPartition);
        DSSSystemPrint(newSSystem);
bail:
        return newSSystem;
}

extern DSSSystem * DSSSystemWithQuasiSteadyStates(const DSSSystem * ssystem, DSUInteger numberOfVariables, const char ** variableNames)
{
        return dsSSystemForQuasiSteadyState(ssystem, numberOfVariables, variableNames);
}

extern void DSSSystemPrint(const DSSSystem * ssys)
{
        int (*print)(const char *, ...);
        //        DSUInteger i;
        if (ssys == NULL) {
                DSError(M_DS_NULL ": S-System to print is NULL", A_DS_ERROR);
                goto bail;
        }
        if (DSPrintf == NULL)
                print = printf;
        else
                print = DSPrintf;
        print("\t  # Xd: %i\n\t  # Xi: %i\n\t     G: %ix%i\n\t     H: %ix%i\n\t Alpha: %ix1\n\t  Beta: %ix1\n\t   Sol: %s",
              DSVariablePoolNumberOfVariables(DSSSysXd(ssys)),
              DSVariablePoolNumberOfVariables(DSSSysXi(ssys)),
              DSMatrixRows(DSSSysGd(ssys)), 
              DSMatrixColumns(DSSSysGd(ssys))+((DSSSysGi(ssys) != NULL) ? DSMatrixColumns(DSSSysGi(ssys)) : 0),
              DSMatrixRows(DSSSysHd(ssys)), 
              DSMatrixColumns(DSSSysHd(ssys))+((DSSSysHi(ssys) != NULL) ? DSMatrixColumns(DSSSysHi(ssys)) : 0),
              DSMatrixRows(DSSSysAlpha(ssys)),
              DSMatrixRows(DSSSysBeta(ssys)),
              (DSSSystemHasSolution(ssys)  ? "YES" : "NO"));
        print("\n");
bail:
        return;
}


extern void DSSSystemPrintEquations(const DSSSystem *ssys)
{
        DSUInteger i;
        DSExpression ** equations = NULL;
        if (ssys == NULL) {
                DSError(M_DS_SSYS_NULL, A_DS_ERROR);
                goto bail;
        }
        equations = DSSSystemEquations(ssys);
        if (equations != NULL) {
                for (i= 0; i < DSSSystemNumberOfEquations(ssys); i++) {
                        DSExpressionPrint(equations[i]);
                        DSExpressionFree(equations[i]);
                }
                DSSecureFree(equations);
        }
bail:
        return;
}

extern void DSSSystemPrintSolution(const DSSSystem *ssys)
{
        DSUInteger i;
        DSExpression ** solution = NULL;
        if (ssys == NULL) {
                DSError(M_DS_SSYS_NULL, A_DS_ERROR);
                goto bail;
        }
        solution = DSSSystemSolution(ssys);
        if (solution != NULL) {
                for (i= 0; i < DSSSystemNumberOfEquations(ssys); i++) {
                        DSExpressionPrint(solution[i]);
                        DSExpressionFree(solution[i]);
                }
                DSSecureFree(solution);
        }
bail:
        return;
}

extern void DSSSystemPrintLogarithmicSolution(const DSSSystem *ssys)
{
        DSUInteger i;
        DSExpression ** solution = NULL;
        if (ssys == NULL) {
                DSError(M_DS_SSYS_NULL, A_DS_ERROR);
                goto bail;
        }
        solution = DSSSystemLogarithmicSolution(ssys);
        if (solution != NULL) {
                for (i= 0; i < DSSSystemNumberOfEquations(ssys); i++) {
                        DSExpressionPrint(solution[i]);
                        DSExpressionFree(solution[i]);
                }
                DSSecureFree(solution);
        }
bail:
        return;
}

extern void DSSSystemAdjustStoichiometryOfCodominantCase(DSSSystem *ssys){
    
    if (ssys == NULL) {
            DSError(M_DS_SSYS_NULL, A_DS_ERROR);
            goto bail;
    }
    
    if (DSSSystemAdjustCodominantStoichiometry(ssys) == false)
        goto bail;
    
    DSMatrix *alpha_adjusted = NULL;
    DSMatrixFree(ssys->alpha);
    alpha_adjusted = DSMatrixCopy(ssys->alpha_adjusted);
    ssys->alpha = alpha_adjusted;
    
bail:
    return;

}

#if defined(__APPLE__) && defined (__MACH__)
#pragma mark - Data Serialization
#endif

extern DSSSystemMessage * DSSSystemEncode(const DSSSystem * ssys)
{
        DSSSystemMessage * message = NULL;
        DSUInteger i;
        const DSVariablePool * X;
        if (ssys == NULL) {
                DSError(M_DS_SSYS_NULL, A_DS_ERROR);
                goto bail;
        }
        message = DSSecureMalloc(sizeof(DSSSystemMessage));
        dsssystem_message__init(message);
        message->alpha = DSMatrixEncode(DSSSystemAlpha(ssys));
        message->beta = DSMatrixEncode(DSSSystemBeta(ssys));
        message->gd = DSMatrixEncode(DSSSystemGd(ssys));
        message->hd = DSMatrixEncode(DSSSystemHd(ssys));
        message->gi = DSMatrixEncode(DSSSystemGi(ssys));
        message->hi = DSMatrixEncode(DSSSystemHi(ssys));
        message->modifierflag = ssys->modifierFlags;
        if (DSSSystemM(ssys) != NULL) {
                message->m = DSMatrixEncode(DSSSystemM(ssys));
        } else {
                message->m = NULL;
        }
        X = DSSSystemXd(ssys);
        message->n_xd = DSVariablePoolNumberOfVariables(X);
        message->xd = DSSecureMalloc(sizeof(char*)*message->n_xd);
        for (i = 0; i < DSVariablePoolNumberOfVariables(X); i++) {
                message->xd[i] = strdup(DSVariableName(DSVariablePoolVariableAtIndex(X, i)));
        }
        X = DSSSystemXi(ssys);
        message->n_xi = DSVariablePoolNumberOfVariables(X);
        message->xi = DSSecureMalloc(sizeof(char*)*message->n_xi);
        for (i = 0; i < DSVariablePoolNumberOfVariables(X); i++) {
                message->xi[i] = strdup(DSVariableName(DSVariablePoolVariableAtIndex(X, i)));
        }
        X = DSSSystemXd_a(ssys);
        message->n_xd_a = DSVariablePoolNumberOfVariables(X);
        message->xd_a = DSSecureMalloc(sizeof(char*)*message->n_xd_a);
        for (i = 0; i < DSVariablePoolNumberOfVariables(X); i++) {
                message->xd_a[i] = strdup(DSVariableName(DSVariablePoolVariableAtIndex(X, i)));
        }
        X = DSSSystemXd_t(ssys);
        message->n_xd_t = DSVariablePoolNumberOfVariables(X);
        message->xd_t = DSSecureMalloc(sizeof(char*)*message->n_xd_t);
        for (i = 0; i < DSVariablePoolNumberOfVariables(X); i++) {
                message->xd_t[i] = strdup(DSVariableName(DSVariablePoolVariableAtIndex(X, i)));
        }
        message->n_xd_t = DSVariablePoolNumberOfVariables(X);
        if (DSSSystemIsConserved(ssys) == true){
            message->has_numberofconservations = true;
            message->numberofconservations = ssys->numberOfConservations;
        } else{
            message->has_numberofconservations = false;
        }

bail:
        return message;
}

extern DSSSystem * DSSSystemFromSSystemMessage(const DSSSystemMessage * message)
{
        DSSSystem * ssystem = NULL;
        DSUInteger i;
        if (message == NULL) {
                printf("message is NULL\n");
                goto bail;
        }
        ssystem = DSSSystemAlloc();
        ssystem->alpha = DSMatrixFromMatrixMessage(message->alpha);
        ssystem->beta = DSMatrixFromMatrixMessage(message->beta);
        ssystem->Gd = DSMatrixFromMatrixMessage(message->gd);
        ssystem->Gi = DSMatrixFromMatrixMessage(message->gi);
        ssystem->Hd = DSMatrixFromMatrixMessage(message->hd);
        ssystem->Hi = DSMatrixFromMatrixMessage(message->hi);
        ssystem->modifierFlags = message->modifierflag;
        if (message->m != NULL) {
                ssystem->M = DSMatrixFromMatrixMessage(message->m);
        } else {
                ssystem->M = NULL;
        }
        ssystem->Xd = DSVariablePoolAlloc();
        ssystem->Xi = DSVariablePoolAlloc();
        ssystem->Xd_a = DSVariablePoolAlloc();
        ssystem->Xd_t = DSVariablePoolAlloc();
        for (i = 0; i < message->n_xd; i++) {
                DSVariablePoolAddVariableWithName(ssystem->Xd, message->xd[i]);
        }
        for (i = 0; i < message->n_xd_a; i++) {
                DSVariablePoolAddVariableWithName(ssystem->Xd_a, message->xd_a[i]);
        }
        for (i = 0; i < message->n_xd_t; i++) {
                DSVariablePoolAddVariableWithName(ssystem->Xd_t, message->xd_t[i]);
        }
        for (i = 0; i < message->n_xi; i++) {
                DSVariablePoolAddVariableWithName(ssystem->Xi, message->xi[i]);
        }
        DSSSystemSetShouldFreeXd(ssystem, true);
        DSSSystemSetShouldFreeXi(ssystem, true);
        if (message->has_numberofconservations == true){
            ssystem->numberOfConservations = message->numberofconservations;
            DSSSystemSetIsConserved(ssystem, true);
        }
    
bail:
        return ssystem;
}

extern DSSSystem * DSSSystemDecode(size_t length, const void * buffer)
{
        DSSSystem * ssystem = NULL;
        DSSSystemMessage * message;
        message = dsssystem_message__unpack(NULL, length, buffer);
        ssystem = DSSSystemFromSSystemMessage(message);
        dsssystem_message__free_unpacked(message, NULL);
bail:
        return ssystem;
}


