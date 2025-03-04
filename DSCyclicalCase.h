/**
 * \file DSCyclicalCase.h
 * \brief Header file with functions for dealing with subcases.
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

#include "DSTypes.h"
#include "DSErrors.h"
#include "DSMemoryManager.h"
#include "DSExpression.h"
#include "DSMatrix.h"
#include "DSMatrixArray.h"
#include "DSSSystem.h"
#include "DSCase.h"
#include "DSGMASystem.h"
#include "DSDesignSpace.h"
#include "DSStack.h"

#ifndef __DS_SUBCASE__
#define __DS_SUBCASE__

/**
 *\addtogroup M_DS_Messages
 * Messages for DSCase related errors is M_DS_CASE_NULL.
 */
/*\{*/
#define M_DS_SUBCASE_NULL                  M_DS_NULL ": Subcase is NULL"
/*\}*/

#if defined (__APPLE__) && defined (__MACH__)
#pragma mark - Allocation, deallocation and initialization -
#endif

extern DSCyclicalCase * DSCyclicalCaseForCaseInDesignSpace(const DSDesignSpace * ds, const DSCase * aCase);
extern void DSCyclicalCaseFree(DSCyclicalCase * aSubcase);

extern void DSCyclicalCaseInit3dSignature(DSDesignSpace *collapsed, const DSDesignSpace *original, const DSCase *aCase);
extern void DSSubCaseGenerate3dSignature(const DSDesignSpace *ds, DSCase *aCase, DSuIntegerMatrix *three_digit, DSuIntegerMatrix *location);

#if defined (__APPLE__) && defined (__MACH__)
#pragma mark - Getter -
#endif

extern const DSVariablePool * DSCyclicalCaseXd(const DSCyclicalCase * cyclicalCase);
extern const DSVariablePool * DSCyclicalCaseXi(const DSCyclicalCase * cyclicalCase);
extern const DSVariablePool * DSCyclicalCaseMainCycleVariables(const DSCyclicalCase *cyclicalCase,
                                                               DSUInteger cycle);
extern const DSVariablePool * DSCyclicalCaseSecondaryCycleVariables(const DSCyclicalCase *cyclicalCase,
                                                                    DSUInteger cycle);
extern const DSDesignSpace * DSCyclicalCaseInternalDesignSpace(const DSCyclicalCase * cyclicalCase);
extern const DSCase * DSCyclicalCaseOriginalCase(const DSCyclicalCase * subcase);

extern const DSUInteger DSCyclicalCaseNumberOfValidSubcases(const DSCyclicalCase *cyclicalCase);
extern const DSUInteger DSCyclicalCaseNumberOfCycles(const DSCyclicalCase *cyclicalCase);
extern const DSUInteger DSCyclicalCaseNumberOfValidBlowingSubcases(const DSCyclicalCase *cyclicalCase);
extern const DSUInteger DSCyclicalCaseNumberOfSubcases(const DSCyclicalCase * cyclicalCase);
extern const DSUInteger DSCyclicalCaseDSDictionaryNumberOfValidSubCasesAndValidBlowingSubcases(const DSDictionary *cyclicalCases);

extern DSCase * DSCyclicalCaseSubcaseWithCaseNumber(const DSCyclicalCase * cyclicalCase, const DSUInteger subcaseNumber);
extern const DSCyclicalCase * DSCyclicalCaseCyclicalSubcaseWithCaseNumber(const DSCyclicalCase * cyclicalCase, const DSUInteger subcaseNumber);
//extern DSCase * DSCyclicalCaseSubcaseWithCaseSignature(const DSCyclicalCase * cyclicalCase, const DSUInteger * signature);

#if defined (__APPLE__) && defined (__MACH__)
#pragma mark - Getter -
#endif

extern const DSUInteger DSCyclicalCaseNumberOfEquations(const DSCyclicalCase *cyclicalCase);

extern DSExpression ** DSCyclicalCaseEquations(const DSCyclicalCase *cyclicalCase);

extern const DSUInteger DSCyclicalCaseNumberOfConditions(const DSCyclicalCase *cyclicalCase);
extern DSExpression ** DSCyclicalCaseConditions(const DSCyclicalCase *cyclicalCase);
extern DSExpression ** DSCyclicalCaseLogarithmicConditions(const DSCyclicalCase *cyclicalCase);

extern const DSUInteger DSCyclicalCaseNumberOfBoundaries(const DSCyclicalCase *cyclicalCase);
extern DSExpression ** DSCyclicalCaseBoundaries(const DSCyclicalCase *cyclicalCase);
extern DSExpression ** DSCyclicalCaseLogarithmicBoundaries(const DSCyclicalCase *cyclicalCase);

extern DSUInteger DSCyclicalCaseNumber(const DSCyclicalCase *cyclicalCase);
extern const char * DSCyclicalCaseIdentifier(const DSCyclicalCase *cyclicalCase);
extern const DSUInteger * DSCyclicalCaseSignature(const DSCyclicalCase *cyclicalCase);
extern char * DSCyclicalCaseSignatureToString(const DSCyclicalCase *cyclicalCase);
extern const DSSSystem *DSCyclicalCaseSSystem(const DSCyclicalCase *cyclicalCase);

//extern double DSCyclicalCaseLogarithmicGain(const DSCase *aCase, const char *XdName, const char *XiName);


#if defined (__APPLE__) && defined (__MACH__)
#pragma mark - Linear programming -
#endif

extern const bool DSCyclicalCaseIsValid(const DSCyclicalCase *aSubcase, const bool strict);
//extern const bool DSCyclicalCaseIsValidAtPoint(const DSCyclicalCase *aSubcase, const DSVariablePool * variablesToFix);
extern const bool DSCyclicalCaseIsValidAtSlice(const DSCyclicalCase *aSubcase, const DSVariablePool * lowerBounds, const DSVariablePool *upperBounds, const bool strict);

//extern const bool DSCyclicalCaseSubaseWithCaseNumberIsValid(const DSDesignSpace *ds, const DSUInteger caseNumber);
//extern const bool DSCyclicalCaseSubcaseWithCaseSignatureIsValid(const DSDesignSpace *ds, const DSUInteger * signature);
//extern const bool DSCyclicalCaseSubcaseWithCaseSignatureListIsValid(const DSDesignSpace *ds, const DSUInteger firstTerm, ...);


#if defined (__APPLE__) && defined (__MACH__)
#pragma mark - Utility -
#endif

extern DSDictionary * DSCyclicalCaseCalculateAllValidSubcasesByResolvingCyclicalCases(DSCyclicalCase *cyclicalCase);
extern void DSDesignSpaceCalculateAllValidCasesByResolvingCyclicalCasesUnstable (DSDesignSpace *ds,
                                                                                          DSDictionary *caseDictionary);
extern void DSDesignSpaceCalculateAllValidCasesForSliceByResolvingCyclicalCasesUnstable (DSDesignSpace *ds,
                                                                                         DSDictionary *caseDictionary,
                                                                                         const DSVariablePool * lower,
                                                                                         const DSVariablePool * upper,
                                                                                         bool strict);

extern DSDictionary * DSCyclicalCaseCalculateAllValidSubcasesForSliceByResolvingCyclicalCases(DSCyclicalCase *cyclicalCase,
                                                                                              const DSVariablePool * lower,
                                                                                              const DSVariablePool * upper,
                                                                                              bool strict);

extern DSDictionary * DSCyclicalCaseCalculateAllValidSubcases(const DSCyclicalCase * cyclicalCase);

extern DSDictionary * DSCyclicalCaseCalculateAllValidSubcasesForSlice(const DSCyclicalCase * cyclicalCase,
                                                                      const DSVariablePool *lower,
                                                                      const DSVariablePool *upper);

//extern DSDictionary * DSCyclicalCaseVerticesForSlice(const DSCyclicalCase *cyclicalCase,
//                                                     const DSVariablePool * lowerBounds,
//                                                     const DSVariablePool *upperBounds,
//                                                     const DSUInteger numberOfVariables,
//                                                     const char ** variables);

extern DSDictionary * DSCyclicalCaseVerticesFor2DSlice(const DSCyclicalCase *cyclicalCase,
                                                       const DSVariablePool * lowerBounds,
                                                       const DSVariablePool *upperBounds,
                                                       const char * xVariable,
                                                       const char *yVariable);

#if defined(__APPLE__) && defined (__MACH__)
#pragma mark - Data Serialization
#endif


extern DSCyclicalCaseMessage * DSCyclicalCaseEncode(const DSCyclicalCase * aCase);
extern DSCyclicalCase * DSCyclicalCaseFromCyclicalCaseMessage(const DSCyclicalCaseMessage * message);
extern DSCyclicalCase * DSCyclicalCaseDecode(size_t length, const void * buffer);

#endif
