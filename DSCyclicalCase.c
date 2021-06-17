/**
 * \file DSCyclicalCase.c
 * \brief Implementation file with functions for dealing with subcases.
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
 * \todo Cleanup this file!!
 */
#include <string.h>
#include <stdio.h>
#include "DSCyclicalCase.h"

extern DSDesignSpace * DSCyclicalCaseInternalForUnderdeterminedCase(const DSCase * aCase, const DSDesignSpace * original);
extern DSDesignSpace * DSCyclicalCaseDesignSpacesForUnderdeterminedCase(const DSCase * aCase, const DSDesignSpace * original);

extern void DSCaseRemoveZeroBoundaries(DSCase *aCase);

#if defined (__APPLE__) && defined (__MACH__)
#pragma mark - Allocation, deallocation and initialization
#endif


extern DSCyclicalCase * DSCyclicalCaseForCaseInDesignSpace(const DSDesignSpace * ds, const DSCase * aCase)
{
        DSCyclicalCase * cyclicalCase = NULL;
        DSDesignSpace * subcase;
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
        cyclicalCase = DSSecureCalloc(sizeof(DSCyclicalCase), 1);
        subcase = DSCyclicalCaseDesignSpacesForUnderdeterminedCase(aCase, ds);
        if (subcase == NULL) {
                DSSecureFree(cyclicalCase);
                cyclicalCase = NULL;
                goto bail;
        }
        cyclicalCase->internalDesignspace = subcase;
        cyclicalCase->originalCase = DSCaseCopy(aCase);
bail:
        return cyclicalCase;
}

extern void DSCyclicalCaseFree(DSCyclicalCase * aSubcase)
{
        if (aSubcase == NULL) {
                DSError(M_DS_SUBCASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (aSubcase->internalDesignspace)
                DSDesignSpaceFree(aSubcase->internalDesignspace);
        if (aSubcase->originalCase != NULL)
                DSCaseFree(aSubcase->originalCase);
        DSSecureFree(aSubcase);
bail:
        return;
}

#if defined (__APPLE__) && defined (__MACH__)
#pragma mark - Getter
#endif


extern const DSVariablePool * DSCyclicalCaseXd(const DSCyclicalCase * cyclicalCase)
{
        const DSVariablePool * Xd = NULL;
        if (cyclicalCase == NULL) {
                DSError(M_DS_SUBCASE_NULL, A_DS_ERROR);
                goto bail;
        }
        Xd = DSSSystemXd(DSCaseSSystem(DSCyclicalCaseOriginalCase(cyclicalCase)));
bail:
        return Xd;
}

extern const DSVariablePool *  DSCyclicalCaseXi(const DSCyclicalCase * cyclicalCase)
{
        const DSVariablePool * Xi = NULL;
        if (cyclicalCase == NULL) {
                DSError(M_DS_SUBCASE_NULL, A_DS_ERROR);
                goto bail;
        }
        Xi = DSSSystemXi(DSCaseSSystem(DSCyclicalCaseOriginalCase(cyclicalCase)));
bail:
        return Xi;
}

extern const DSVariablePool * DSCyclicalCaseMainCycleVariables(const DSCyclicalCase *cyclicalCase,
                                                               DSUInteger cycle)
{
    
    DSVariablePool * X_main = NULL;
    const DSVariablePool * Xd = NULL;
    DSUInteger i;
    
    
    if (cyclicalCase == NULL) {
        DSError(M_DS_SUBCASE_NULL, A_DS_ERROR);
        goto bail;
    }
    if (cycle > DSCyclicalCaseNumberOfCycles(cyclicalCase)) {
        DSError("cycle index exceeds number of cycles", A_DS_ERROR);
        goto bail;
    }
    
    Xd = DSDesignSpaceXd(cyclicalCase->internalDesignspace);
    X_main = DSVariablePoolAlloc();
    i = cyclicalCase->internalDesignspace->extensionData->mainCycleVariables[cycle - 1];
    DSVariablePoolAddVariableWithName(X_main, DSVariablePoolVariableAtIndex(Xd, i)->name);
    
bail:
    return X_main;

}

extern const DSVariablePool * DSCyclicalCaseSecondaryCycleVariables(const DSCyclicalCase *cyclicalCase,
                                                                    DSUInteger cycle)
{
    
    DSVariablePool * X_sec = NULL;
    const DSVariablePool * Xd = NULL;
    DSUInteger i, ii, index;
    
    
    if (cyclicalCase == NULL) {
        DSError(M_DS_SUBCASE_NULL, A_DS_ERROR);
        goto bail;
    }
    if (cycle > DSCyclicalCaseNumberOfCycles(cyclicalCase)) {
        DSError("cycle index exceeds number of cycles", A_DS_ERROR);
        goto bail;
    }
    
    Xd = DSDesignSpaceXd(cyclicalCase->internalDesignspace);
    X_sec = DSVariablePoolAlloc();
    i = cyclicalCase->internalDesignspace->extensionData->numberSecondaryVariables[cycle - 1];
    
    for (ii=0; ii<i; ii++){
        index = cyclicalCase->internalDesignspace->extensionData->allSecondaryVariables[cycle-1][ii];
        DSVariablePoolAddVariableWithName(X_sec, DSVariablePoolVariableAtIndex(Xd, index)->name);
    }
    
bail:
    return X_sec;
    
}


extern const DSDesignSpace * DSCyclicalCaseInternalDesignSpace(const DSCyclicalCase * subcase)
{
        DSDesignSpace * ds = NULL;
        if (subcase == NULL) {
                DSError(M_DS_SUBCASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (subcase->internalDesignspace == NULL)
                goto bail;
        ds = subcase->internalDesignspace;
bail:
        return ds;
}

extern const DSCase * DSCyclicalCaseOriginalCase(const DSCyclicalCase * cyclicalCase)
{
        DSCase * aCase = NULL;
        if (cyclicalCase == NULL) {
                DSError(M_DS_SUBCASE_NULL, A_DS_ERROR);
                goto bail;
        }
        aCase = cyclicalCase->originalCase;
bail:
        return aCase;
}

extern const DSUInteger DSCyclicalCaseNumberOfValidSubcases(const DSCyclicalCase *cyclicalCase)
{
        DSDesignSpace * ds;
        DSUInteger numberOfValidSubcases = 0;
        if (cyclicalCase == NULL) {
                DSError(M_DS_CASE_NULL ": Cyclical case is null", A_DS_ERROR);
                goto bail;
        }
        ds = cyclicalCase->internalDesignspace;
                if (ds == NULL) {
                        DSError(M_DS_DESIGN_SPACE_NULL, A_DS_ERROR);
                        goto bail;
                }
                numberOfValidSubcases = DSDesignSpaceNumberOfValidCases(ds);
bail:
        return numberOfValidSubcases;
}

extern const DSUInteger DSCyclicalCaseNumberOfCycles(const DSCyclicalCase *cyclicalCase)
{
    DSDesignSpace * ds;
    DSUInteger numberOfCycles = 0;
    if (cyclicalCase == NULL) {
        DSError(M_DS_CASE_NULL ": Cyclical case is null", A_DS_ERROR);
        goto bail;
    }
    ds = cyclicalCase->internalDesignspace;
    if (ds == NULL) {
        DSError(M_DS_DESIGN_SPACE_NULL, A_DS_ERROR);
        goto bail;
    }
    numberOfCycles = ds->extensionData->numberCycles;
bail:
    return numberOfCycles;
}

extern const DSUInteger DSCyclicalCaseNumberOfValidBlowingSubcases(const DSCyclicalCase *cyclicalCase)
{
    DSDesignSpace * ds;
    DSUInteger numberOfValidBlowingSubcases = 0;
    if (cyclicalCase == NULL) {
        DSError(M_DS_CASE_NULL ": Cyclical case is null", A_DS_ERROR);
        goto bail;
    }
    ds = cyclicalCase->internalDesignspace;
    if (ds == NULL) {
        DSError(M_DS_DESIGN_SPACE_NULL, A_DS_ERROR);
        goto bail;
    }
//    numberOfValidBlowingSubcases = DSDictionaryCount(ds->unstableCases);
    numberOfValidBlowingSubcases = DSDesignSpaceNumberOfValidBlowingCases(ds, true);
    
bail:
    return numberOfValidBlowingSubcases;
}

extern const DSUInteger DSCyclicalCaseNumberOfSubcases(const DSCyclicalCase * cyclicalCase)
{
        DSUInteger  numberOfCases = 0;
        if (cyclicalCase == NULL) {
                DSError(M_DS_CASE_NULL ": Cyclical case is null", A_DS_ERROR);
                goto bail;
        }
        numberOfCases = DSDesignSpaceNumberOfCases(cyclicalCase->internalDesignspace);
//        for (i = 0; i < cyclicalCase->numberOfInternal; i++) {
//                numberOfCases += DSDesignSpaceNumberOfCases(cyclicalCase->internalDesignspaces[i]);
//        }
bail:
        return numberOfCases;
}

extern const DSUInteger DSCyclicalCaseDSDictionaryNumberOfValidSubCasesAndValidBlowingSubcases(const DSDictionary *cyclicalCases){
    
    DSUInteger i, valid_cyclical_cases = 0;
    DSCyclicalCase *cyclicalCase;
    
    if(cyclicalCases == NULL)
        goto bail;
    
    for(i = 0; i<DSDictionaryCount(cyclicalCases); i++ ){
        cyclicalCase = DSDictionaryValueForName(cyclicalCases, DSDictionaryNames(cyclicalCases)[i]);
        if (DSCyclicalCaseNumberOfValidSubcases(cyclicalCase) !=0 ||
            DSCyclicalCaseNumberOfValidBlowingSubcases(cyclicalCase) !=0)
//        if (DSCyclicalCaseNumberOfValidSubcases(cyclicalCase) != 0)
            valid_cyclical_cases++;
    }
    
bail:
    return valid_cyclical_cases;
}

extern DSCase * DSCyclicalCaseSubcaseWithCaseNumber(const DSCyclicalCase * cyclicalCase, const DSUInteger subcaseNumber)
{
        DSCase * aSubcase = NULL;
        if (cyclicalCase == NULL) {
                DSError(M_DS_CASE_NULL ": Cyclical case is null", A_DS_ERROR);
                goto bail;
        }
        aSubcase = DSDesignSpaceCaseWithCaseNumber(cyclicalCase->internalDesignspace, subcaseNumber);
//        DSCaseRemoveZeroBoundaries(aSubcase);
bail:
        return aSubcase;
}

extern const DSCyclicalCase * DSCyclicalCaseCyclicalSubcaseWithCaseNumber(const DSCyclicalCase * cyclicalCase, const DSUInteger subcaseNumber)
{
        const DSCyclicalCase * aSubcase = NULL;
        if (cyclicalCase == NULL) {
                DSError(M_DS_CASE_NULL ": Cyclical case is null", A_DS_ERROR);
                goto bail;
        }
        aSubcase = DSDesignSpaceCyclicalCaseWithCaseNumber(cyclicalCase->internalDesignspace, subcaseNumber);
bail:
        return aSubcase;
}

extern const DSUInteger DSCyclicalCaseNumberOfEquations(const DSCyclicalCase *cyclicalCase)
{
        return DSCaseNumberOfEquations(DSCyclicalCaseOriginalCase(cyclicalCase));
}

extern DSExpression ** DSCyclicalCaseEquations(const DSCyclicalCase *cyclicalCase)
{
        return DSCaseEquations(DSCyclicalCaseOriginalCase(cyclicalCase));
}

extern const DSUInteger DSCyclicalCaseNumberOfConditions(const DSCyclicalCase *cyclicalCase)
{
        return DSCaseNumberOfConditions(DSCyclicalCaseOriginalCase(cyclicalCase));
}

extern DSExpression ** DSCyclicalCaseConditions(const DSCyclicalCase *cyclicalCase)
{
        return DSCaseConditions(DSCyclicalCaseOriginalCase(cyclicalCase));
}

extern DSExpression ** DSCyclicalCaseLogarithmicConditions(const DSCyclicalCase *cyclicalCase)
{
        return DSCaseLogarithmicConditions(DSCyclicalCaseOriginalCase(cyclicalCase));
}

extern const DSUInteger DSCyclicalCaseNumberOfBoundaries(const DSCyclicalCase *cyclicalCase)
{
        return DSCaseNumberOfBoundaries(DSCyclicalCaseOriginalCase(cyclicalCase));
}

extern DSExpression ** DSCyclicalCaseBoundaries(const DSCyclicalCase *cyclicalCase)
{
        return DSCaseBoundaries(DSCyclicalCaseOriginalCase(cyclicalCase));
}

extern DSExpression ** DSCyclicalCaseLogarithmicBoundaries(const DSCyclicalCase *cyclicalCase)
{
        return DSCaseLogarithmicBoundaries(DSCyclicalCaseOriginalCase(cyclicalCase));
}

extern DSUInteger DSCyclicalCaseNumber(const DSCyclicalCase *cyclicalCase)
{
        return DSCaseNumber(DSCyclicalCaseOriginalCase(cyclicalCase));
}

extern const char * DSCyclicalCaseIdentifier(const DSCyclicalCase *cyclicalCase)
{
        return DSCaseIdentifier(DSCyclicalCaseOriginalCase(cyclicalCase));
}

extern const DSUInteger * DSCyclicalCaseSignature(const DSCyclicalCase *cyclicalCase)
{
        return DSCaseSignature(DSCyclicalCaseOriginalCase(cyclicalCase));
}

extern char * DSCyclicalCaseSignatureToString(const DSCyclicalCase *cyclicalCase)
{
    if (cyclicalCase->internalDesignspace->parent3DigitsSignature == NULL){
        return DSCaseSignatureToString(DSCyclicalCaseOriginalCase(cyclicalCase));
    }else{
                char temp[100];
                char * string = NULL;
                char space[50];
                DSUInteger i;
                DSCase *aCase = cyclicalCase->originalCase;
                DSUInteger *parentSig = cyclicalCase->internalDesignspace->parent3DigitsSignature;
                DSUInteger numberOfEquations;
        
                numberOfEquations =  DSCaseNumberOfEquations(aCase) + DSCaseNumberOfConservations(aCase) + DSCaseNumberOfInheritedConservations(aCase);
        
                string = DSSecureCalloc(sizeof(char), 5*numberOfEquations);
                strcpy(space, "  ");
                    for (i = 0; i < 3*numberOfEquations; i++) {
                        if (parentSig[i] >= 10){
                            if (parentSig[i] != 0.0){
                                sprintf(temp, "(%i)", parentSig[i]);
                            } else
                                continue;
                        } else {
                            if (parentSig[i] != 0.0){
                                sprintf(temp, "%i", parentSig[i]);
                            }else if ((i+2)%3 == 0 || (i+1)%3 == 0){
                                sprintf(temp, "%i", parentSig[i]);
                            }else
                                continue;
                        }
                        strncat(string, temp, 100-strlen(string));
                        if ((i+1)%3 == 0)
                            strncat(string, space, 100-strlen(space));
                    }
                    return string;
            }
}

extern const DSSSystem *DSCyclicalCaseSSystem(const DSCyclicalCase *cyclicalCase)
{
        return DSCaseSSystem(DSCyclicalCaseOriginalCase(cyclicalCase));
}

#if defined (__APPLE__) && defined (__MACH__)
#pragma mark Linear programming functions
#endif

extern const bool DSCyclicalCaseIsValid(const DSCyclicalCase *aSubcase, const bool strict)
{
        bool isValid = false;
        DSUInteger numberOfValidCases = 0, numberOfValidBlowingCases = 0;
        if (aSubcase == NULL) {
                DSError(M_DS_SUBCASE_NULL, A_DS_ERROR);
                goto bail;
        }
        numberOfValidCases = DSCyclicalCaseNumberOfValidSubcases(aSubcase);
        if (DSDesignSpaceUnstable(aSubcase->internalDesignspace) == true)
            numberOfValidBlowingCases = DSCyclicalCaseNumberOfValidBlowingSubcases(aSubcase);
        if (numberOfValidCases + numberOfValidBlowingCases > 0)
                isValid = true;
bail:
        return isValid;
}

extern const bool DSCyclicalCaseIsValidAtPoint(const DSCyclicalCase *aSubcase, const DSVariablePool * variablesToFix);

extern const bool DSCyclicalCaseIsValidAtSlice(const DSCyclicalCase *cyclicalCase, const DSVariablePool * lowerBounds, const DSVariablePool *upperBounds, const bool strict)
{
        bool isValid = false;
        DSUInteger numberValid;
        DSUInteger numberBlowingValid = 0;
        DSDesignSpace * ds;
        DSDictionary * validCases;
        if (cyclicalCase == NULL) {
                DSError(M_DS_SUBCASE_NULL, A_DS_ERROR);
                goto bail;
        }
        ds = cyclicalCase->internalDesignspace;
        if (ds == NULL) {
                DSError(M_DS_DESIGN_SPACE_NULL, A_DS_ERROR);
                goto bail;
        }
        numberValid = DSDesignSpaceNumberOfValidCases(ds);
        if (DSDesignSpaceUnstable(ds) == true)
            numberBlowingValid = DSDictionaryCount(ds->unstableCases);
        if (numberValid + numberBlowingValid == 0)
                goto bail;
        if (strict == true)
                validCases = DSDesignSpaceCalculateAllValidCasesForSlice(ds, lowerBounds, upperBounds);
        else
                validCases = DSDesignSpaceCalculateAllValidCasesForSliceNonStrict(ds, lowerBounds, upperBounds);
        if (validCases == NULL)
                goto bail;
        if (DSDictionaryCount(validCases) != 0)
                isValid = true;
        DSDictionaryFreeWithFunction(validCases, DSCaseFree);
bail:
        return isValid;
}

extern DSDictionary * DSCyclicalCaseVerticesForSlice(const DSCyclicalCase *aSubcase, const DSVariablePool * lowerBounds, const DSVariablePool *upperBounds, const DSUInteger numberOfVariables, const char ** variables);

#if defined (__APPLE__) && defined (__MACH__)
#pragma mark - Utility functions
#endif

extern void DSCyclicalCaseAddConstraints(DSCyclicalCase * cyclicalCase, const char ** strings, DSUInteger numberOfConstraints)
{
        DSDesignSpace * ds;
        if (cyclicalCase == NULL) {
                DSError(M_DS_CASE_NULL ": Cyclical Case is Null", A_DS_ERROR);
                goto bail;
        }
        ds = cyclicalCase->internalDesignspace;
                if (ds == NULL) {
                        DSError(M_DS_DESIGN_SPACE_NULL, A_DS_ERROR);
                        goto bail;
                }
                DSDesignSpaceAddConstraints(ds, strings, numberOfConstraints);
bail:
        return;
}

extern DSDictionary * DSCyclicalCaseCalculateAllValidSubcasesByResolvingCyclicalCases(DSCyclicalCase *cyclicalCase)
{
        DSDictionary * caseDictionary = NULL;
        char * subcaseString = NULL;
        if (cyclicalCase == NULL) {
                DSError(M_DS_CASE_NULL ": Cyclical Case is Null", A_DS_ERROR);
                goto bail;
        }
        caseDictionary = DSDesignSpaceCalculateAllValidCasesByResolvingCyclicalCases(cyclicalCase->internalDesignspace);
bail:
        if (subcaseString != NULL)
                DSSecureFree(subcaseString);
        return caseDictionary;
}

extern void DSDesignSpaceCalculateAllValidCasesByResolvingCyclicalCasesUnstable (DSDesignSpace *ds,
                                                                                 DSDictionary *caseDictionary)
{
    DSUInteger i, numberBlowingCases, ii, j;
    DSUnstableCase *uCase;
    DSCase *aCase;
    bool strict = true;
    DSVariablePool *lowerBounds, *upperBounds;
    
    
    // First determine total number of blowing subcases
    numberBlowingCases = DSDesignSpaceNumberOfValidBlowingCases(ds, true);
    if (numberBlowingCases == 0)
        goto bail;
    
    // Second loop to generate array of Cases validCasesBlowing.
    for (i = 0; i < DSDictionaryCount(ds->unstableCases); i++){
        uCase = DSDictionaryValueForName(ds->unstableCases, ds->unstableCases->names[i]);
        for (ii = 0; ii < DSDictionaryCount(uCase->ValidCases); ii++){
            
                aCase = DSCaseCopy(DSDictionaryValueForName(uCase->ValidCases, uCase->ValidCases->names[ii]));
            
                lowerBounds = DSVariablePoolCopy(aCase->Xi);
                upperBounds = DSVariablePoolCopy(aCase->Xi);
            
                for ( j = 0; j < DSVariablePoolNumberOfVariables(lowerBounds); j++){
                    DSVariablePoolSetValueForVariableWithName(lowerBounds,
                                                              DSVariablePoolVariableAtIndex(lowerBounds,j)->name,
                                                              1E-6);
                    DSVariablePoolSetValueForVariableWithName(upperBounds,
                                                              DSVariablePoolVariableAtIndex(lowerBounds,j)->name,
                                                              1E6);
                }
                if (DSCaseIsValidAtSlice(aCase, lowerBounds, upperBounds, strict) == true){
                        DSDictionaryAddValueWithName(caseDictionary, aCase->caseIdentifier, aCase);
                } else {
                    DSCaseFree(aCase);
                }
                DSVariablePoolFree(lowerBounds);
                DSVariablePoolFree(upperBounds);
        }
    }
    
bail:
    return;
    
}

extern void DSDesignSpaceCalculateAllValidCasesForSliceByResolvingCyclicalCasesUnstable (DSDesignSpace *ds,
                                                                                         DSDictionary *caseDictionary,
                                                                                         const DSVariablePool * lower,
                                                                                         const DSVariablePool * upper,
                                                                                         bool strict)
{
    DSUInteger i, numberBlowingCases, ii;
    DSUnstableCase *uCase;
    DSCase *aCase;
//    bool strict = false;
    
    
    // First determine total number of blowing subcases
    numberBlowingCases = DSDesignSpaceNumberOfValidBlowingCases(ds, strict);
    if (numberBlowingCases == 0)
        goto bail;
    
    // Second loop to generate array of Cases validCasesBlowing.
    for (i = 0; i < DSDictionaryCount(ds->unstableCases); i++){
        uCase = DSDictionaryValueForName(ds->unstableCases, ds->unstableCases->names[i]);
        for (ii = 0; ii < DSDictionaryCount(uCase->ValidCases); ii++){
            aCase = DSCaseCopy(DSDictionaryValueForName(uCase->ValidCases, uCase->ValidCases->names[ii]));
            if (DSCaseIsValidAtSlice(aCase, lower, upper, strict) == true){
                    DSDictionaryAddValueWithName(caseDictionary, aCase->caseIdentifier, aCase);
            } else {
                DSCaseFree(aCase);
            }
        }
    }
    
    
bail:
    return;
    
}

extern DSDictionary * DSCyclicalCaseCalculateAllValidSubcasesForSliceByResolvingCyclicalCases(DSCyclicalCase *cyclicalCase,
                                                                                              const DSVariablePool * lower,
                                                                                              const DSVariablePool * upper,
                                                                                              bool strict)
{
        DSDictionary * caseDictionary = NULL;
        if (cyclicalCase == NULL) {
                DSError(M_DS_CASE_NULL ": Cyclical Case is Null", A_DS_ERROR);
                goto bail;
        }
        if (lower == NULL || upper == NULL) {
                DSError(M_DS_VAR_NULL, A_DS_ERROR);
                goto bail;
        }
        caseDictionary = DSDesignSpaceCalculateAllValidCasesForSliceByResolvingCyclicalCases(cyclicalCase->internalDesignspace, lower, upper, strict);
bail:
        return caseDictionary;
}


extern DSDictionary * DSCyclicalCaseCalculateAllValidSubcases(const DSCyclicalCase * cyclicalCase)
{
        const DSDesignSpace * ds;
        DSDictionary * caseDictionary = NULL;
        DSUInteger i, numberValid = 0, numberValidSlice = 0;
        DSUInteger validCaseNumbers = 0;
        DSCase * aCase = NULL;
        if (cyclicalCase == NULL) {
                DSError(M_DS_CASE_NULL ": Cyclical Case is Null", A_DS_ERROR);
                goto bail;
        }
        numberValidSlice = 0;
        caseDictionary = DSDictionaryAlloc();
        ds = cyclicalCase->internalDesignspace;
        if (ds == NULL) {
                DSError(M_DS_DESIGN_SPACE_NULL, A_DS_ERROR);
                goto bail;
        }
        numberValid = DSDesignSpaceNumberOfValidCases(ds);
        if (numberValid == 0) {
                goto bail;
        }
        for (i = 0; i < numberValid; i++) {
                validCaseNumbers = atoi(ds->validCases->names[i]);
                aCase = DSDesignSpaceCaseWithCaseNumber(ds, validCaseNumbers);
                DSDictionaryAddValueWithName(caseDictionary, DSCaseIdentifier(aCase), aCase);
        }
bail:
        return caseDictionary;
}

extern DSDictionary * DSCyclicalCaseCalculateAllValidSubcasesForSlice(const DSCyclicalCase * cyclicalCase,
                                                                      const DSVariablePool *lower,
                                                                      const DSVariablePool *upper)
{
        const DSDesignSpace * ds;
        DSDictionary * caseDictionary = NULL;
        DSUInteger i, numberValid = 0, numberValidSlice = 0;
        DSUInteger validCaseNumbers = 0;
        DSCase * aCase = NULL;
        if (cyclicalCase == NULL) {
                DSError(M_DS_CASE_NULL ": Cyclical Case is Null", A_DS_ERROR);
                goto bail;
        }
        numberValidSlice = 0;
        caseDictionary = DSDictionaryAlloc();
        ds = cyclicalCase->internalDesignspace;
        if (ds == NULL) {
                DSError(M_DS_DESIGN_SPACE_NULL, A_DS_ERROR);
                goto bail;
        }
        numberValid = DSDesignSpaceNumberOfValidCases(ds);
        if (numberValid == 0) {
                goto bail;
        }
        for (i = 0; i < numberValid; i++) {
                validCaseNumbers = atoi(ds->validCases->names[i]);
                aCase = DSDesignSpaceCaseWithCaseNumber(ds, validCaseNumbers);
                if (DSCaseIsValidAtSlice(aCase, lower, upper, true) == true) {
                        DSDictionaryAddValueWithName(caseDictionary, DSCaseIdentifier(aCase), aCase);
                } else {
                        DSCaseFree(aCase);
                }
        }
bail:
        return caseDictionary;
}

extern DSDictionary * DSCyclicalCaseVerticesForSlice(const DSCyclicalCase *cyclicalCase,
                                                     const DSVariablePool * lowerBounds,
                                                     const DSVariablePool *upperBounds,
                                                     const DSUInteger numberOfVariables,
                                                     const char ** variables)
{
        const DSDesignSpace * ds = DSCyclicalCaseInternalDesignSpace(cyclicalCase);
        DSDictionary * caseDictionary = NULL;
        DSVertices * vertices = NULL;
        DSUInteger i, numberValid = 0;
        DSUInteger validCaseNumbers = 0;
        DSCase * aCase = NULL;
        if (ds == NULL) {
                DSError(M_DS_DESIGN_SPACE_NULL, A_DS_ERROR);
                goto bail;
        }
        caseDictionary = DSDictionaryAlloc();
        ds = cyclicalCase->internalDesignspace;
        if (ds == NULL) {
                DSError(M_DS_DESIGN_SPACE_NULL, A_DS_ERROR);
                goto bail;
        }
        numberValid = DSDesignSpaceNumberOfValidCases(ds);
        if (numberValid == 0) {
                goto bail;
        }
        for (i = 0; i < numberValid; i++) {
                validCaseNumbers = atoi(ds->validCases->names[i]);
                aCase = DSDesignSpaceCaseWithCaseNumber(ds, validCaseNumbers);
                if (DSCaseIsValidAtSlice(aCase, lowerBounds, upperBounds, true) == true) {
                        vertices = DSCaseVerticesForSlice(aCase, lowerBounds, upperBounds, numberOfVariables, variables);
                        DSDictionaryAddValueWithName(caseDictionary, DSCaseIdentifier(aCase), aCase);
                }
                DSCaseFree(aCase);
        }
bail:
        return caseDictionary;
}

extern DSDictionary * DSCyclicalCaseVerticesFor2DSlice(const DSCyclicalCase *cyclicalCase,
                                                       const DSVariablePool * lowerBounds,
                                                       const DSVariablePool *upperBounds,
                                                       const char * xVariable,
                                                       const char *yVariable)
{
        const DSDesignSpace * ds = DSCyclicalCaseInternalDesignSpace(cyclicalCase);
        DSDictionary * caseDictionary = NULL;
        DSVertices * vertices = NULL;
        DSUInteger i, numberValid = 0;
        DSUInteger validCaseNumbers = 0;
        DSCase * aCase = NULL;
        if (ds == NULL) {
                DSError(M_DS_DESIGN_SPACE_NULL, A_DS_ERROR);
                goto bail;
        }
        caseDictionary = DSDictionaryAlloc();
        ds = cyclicalCase->internalDesignspace;
        if (ds == NULL) {
                DSError(M_DS_DESIGN_SPACE_NULL, A_DS_ERROR);
                goto bail;
        }
        numberValid = DSDesignSpaceNumberOfValidCases(ds);
        if (numberValid == 0) {
                goto bail;
        }
        for (i = 0; i < numberValid; i++) {
                validCaseNumbers = atoi(ds->validCases->names[i]);
                aCase = DSDesignSpaceCaseWithCaseNumber(ds, validCaseNumbers);
                if (DSCaseIsValidAtSlice(aCase, lowerBounds, upperBounds, true) == true) {
                        vertices = DSCaseVerticesFor2DSlice(aCase, lowerBounds, upperBounds, xVariable, yVariable);
                        if (vertices != NULL)
                                DSDictionaryAddValueWithName(caseDictionary, DSCaseIdentifier(aCase), aCase);
                }
                DSCaseFree(aCase);
        }
bail:
        return caseDictionary;
}


#if defined(__APPLE__) && defined (__MACH__)
#pragma mark - Data Serialization
#endif


extern DSCyclicalCaseMessage * DSCyclicalCaseEncode(const DSCyclicalCase * aCase)
{
        DSCyclicalCaseMessage * message = NULL;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        message = DSSecureMalloc(sizeof(DSCyclicalCaseMessage));
        dscyclical_case_message__init(message);
        message->originalcase = DSCaseEncode(aCase->originalCase);
        message->internaldesignspace = DSDesignSpaceEncode(aCase->internalDesignspace);
bail:
        return message;
}

extern DSCyclicalCase * DSCyclicalCaseFromCyclicalCaseMessage(const DSCyclicalCaseMessage * message)
{
        DSCyclicalCase * aCase = NULL;
        if (message == NULL) {
                printf("message is NULL\n");
                goto bail;
        }
        aCase = DSSecureCalloc(sizeof(DSCyclicalCase), 1);
        aCase->originalCase = DSCaseFromCaseMessage(message->originalcase);
        aCase->internalDesignspace = DSDesignSpaceFromDesignSpaceMessage(message->internaldesignspace);
bail:
        return aCase;
}

extern DSCyclicalCase * DSCyclicalCaseDecode(size_t length, const void * buffer)
{
        DSCyclicalCase * aCase = NULL;
        DSCyclicalCaseMessage * message;
        message = dscyclical_case_message__unpack(NULL, length, buffer);
        aCase = DSCyclicalCaseFromCyclicalCaseMessage(message);
        dscyclical_case_message__free_unpacked(message, NULL);
bail:
        return aCase;
}


