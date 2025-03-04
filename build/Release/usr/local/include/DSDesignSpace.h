/**
 * \file DSDesignSpace.h
 * \brief Header file with functions for dealing with Design Spaces
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

#ifndef __DS_DESIGN_SPACE__
#define __DS_DESIGN_SPACE__

#include "DSTypes.h"
#include "DSErrors.h"
#include "DSDataSerialization.pb-c.h"

#ifdef __cplusplus
__BEGIN_DECLS
#endif

#define M_DS_DESIGN_SPACE_NULL              M_DS_NULL ": Design Space is NULL"

#define DS_DESIGN_SPACE_FLAG_SERIAL                      0x01
#define DS_DESIGN_SPACE_FLAG_CYCLICAL                    0x02
#define DS_DESIGN_SPACE_FLAG_RESOLVE_CO_DOMINANCE        0x04
#define DS_DESIGN_SPACE_FLAG_UNSTABLE                    0x08
#define DS_DESIGN_SPACE_FLAG_CONSERVATIONS               0x10
#define DS_DESIGN_SPACE_FLAG_CO_DOMINANCE_ADJUST_STOICHIOMETRY          0x20
#define DS_DESIGN_SPACE_FLAG_CO_DOMINANCE_SKIP_OVERLAPPING_PHENOTYPES   0x40
#define DS_DESIGN_SPACE_FLAG_MASS_BALANCES                              0x80

// these flags are contained in ds->modifierFlags2
#define DS_DESIGN_SPACE_FLAG_C0_DOMINANCE_ADJUST_BOUNDARIES             0x01

#if defined (__APPLE__) && defined (__MACH__)
#pragma mark - Allocation, deallocation and initialization
#endif

DSDesignSpace * DSDesignSpaceAlloc(void);
void DSDesignSpaceFree(DSDesignSpace * ds);

#if defined (__APPLE__) && defined (__MACH__)
#pragma mark - Factory -
#endif

extern DSDesignSpace * DSDesignSpaceByParsingStringList(const char * string, const DSVariablePool * const Xd_a, ...);
extern DSDesignSpace * DSDesignSpaceByParsingStrings(char * const * const strings, const DSVariablePool * const Xd_a, const DSUInteger numberOfEquations);
extern DSDesignSpace * DSDesignSpaceByParsingStringsWithXi(char * const * const strings, const DSVariablePool * const Xd_a, const DSVariablePool * const Xi, const DSUInteger numberOfEquations);

#if defined (__APPLE__) && defined (__MACH__)
#pragma mark - Setters -
#endif

extern void DSDesignSpaceSetGMA(DSDesignSpace * ds, DSGMASystem *gma);
extern void DSDesignSpaceAddConditions(DSDesignSpace *ds, const DSMatrix * Cd, const DSMatrix * Ci, const DSMatrix * delta);
extern void DSDesignSpaceSetSerial(DSDesignSpace *ds, bool serial);
extern void DSDesignSpaceSetCyclical(DSDesignSpace *ds, bool cyclical);
extern void DSDesignSpaceSetResolveCoDominance(DSDesignSpace *ds, bool Codominance);
extern void DSDesignSpaceSetAdjustCodominantStoichiometry(DSDesignSpace *ds, bool adjust);
extern void DSDesignSpaceSetSkipOverlappingCodominantPhenotypes(DSDesignSpace *ds, bool adjust);
extern void DSDesignSpaceSetShouldConsiderMassBalances(DSDesignSpace *ds, bool mass_balance);
extern void DSDesignSpaceSetAdjustCodominantBoundaries(DSDesignSpace *ds, bool adjust_boundaries);
extern void DSDesignSpaceSetUnstable(DSDesignSpace *ds, bool Unstable);
extern void DSDesignSpaceSetResolveConservations(DSDesignSpace *ds, bool Conservations);
extern void DSDesignSpaceSetNumberOfConservations(DSDesignSpace *ds, DSUInteger numberOfConservations);
extern void DSDesignSpaceSetNumberOfInheritedConservations(DSDesignSpace *collapsed, const DSDesignSpace *original);
extern void DSDesignSpaceInitializeMassBalances(DSDesignSpace *ds,
                                                const char ** fin_strings,
                                                const char ** fout_strings,
                                                const char ** signature_string,
                                                DSUInteger numberOfMassBalances,
                                                const DSVariablePool * metabolicBlocks,
                                                const char ** S_string, DSUInteger rows, DSUInteger columns,
                                                const char ** rxns);

#if defined (__APPLE__) && defined (__MACH__)
#pragma mark - Getters -
#endif

extern bool DSDesignSpaceSerial(const DSDesignSpace *ds);
extern bool DSDesignSpaceCyclical(const DSDesignSpace *ds);
extern bool DSDesignSpaceResolveCoDominance(const DSDesignSpace *ds);
extern bool DSDesignSpaceUnstable(const DSDesignSpace *ds);
extern bool DSDesignSpaceConserved(const DSDesignSpace *ds);
extern bool DSDesignSpaceShouldConsiderMassBalances(const DSDesignSpace * ds);
extern bool DSDesignSpaceAdjustCodominantStoichiometry(const DSDesignSpace *ds);
extern bool DSDesignSpaceSkipOverlappingCodominantPhenotypes(const DSDesignSpace *ds);
extern bool DSDesignSpaceShouldAdjustCodominantBoundaries(DSDesignSpace *ds);


extern const char * DSDesignSpaceFinAtIndex(const DSDesignSpace *ds, DSUInteger n);
extern const char * DSDesignSpaceFoutAtIndex(const DSDesignSpace *ds, DSUInteger n);
extern DSUInteger DSDesignSpaceNumberOfMetabolicBlocks(const DSDesignSpace *ds);

extern const DSVariablePool * DSDesignSpaceXi(const DSDesignSpace *ds);
extern const DSVariablePool * DSDesignSpaceXd(const DSDesignSpace *ds);

extern const DSUInteger DSDesignSpaceNumberOfEquations(const DSDesignSpace *ds);

extern DSExpression ** DSDesignSpaceEquations(const DSDesignSpace *ds);
extern const DSUInteger * DSDesignSpaceSignature(const DSDesignSpace *ds);
extern  char * DSDesignSpaceDominantSignature(const DSDesignSpace * ds,
                                                    const DSVariablePool * Xi,
                                                    const DSVariablePool * Xd);
extern char * DSDesignSpaceSignatureToString(const DSDesignSpace *ds);

extern const DSUInteger DSDesignSpaceNumberOfValidCases(const DSDesignSpace *ds);
extern const DSUInteger DSDesignSpaceNumberOfValidBlowingCases(const DSDesignSpace *ds, bool strict);
extern const DSUInteger DSDesignSpaceNumberOfCases(const DSDesignSpace *ds);
extern DSUInteger DSDesignSpaceNumberOfConservations(const DSDesignSpace *ds);
extern DSUInteger DSDesignSpaceNumberOfBoundaries(const DSDesignSpace *ds);

extern DSCase * DSDesignSpaceCaseWithCaseNumber(const DSDesignSpace * ds, const DSUInteger caseNumber);
extern DSCase * DSDesignSpaceCaseWithCaseIdentifier(const DSDesignSpace * ds, const char * identifer);
extern DSCase * DSDesignSpaceCaseWithCaseSignature(const DSDesignSpace * ds, const DSUInteger * signature);
//extern DSCase * DSDesignSpaceCaseWithCaseSignatureList(const DSDesignSpace *ds, const DSUInteger firstTerm, ...);

extern const bool DSDesignSpaceCaseWithCaseNumberIsValid(const DSDesignSpace *ds, const DSUInteger caseNumber);
extern const bool DSDesignSpaceCaseWithCaseSignatureIsValid(const DSDesignSpace *ds, const DSUInteger * signature);
//extern const bool DSDesignSpaceCaseWithCaseSignatureListIsValid(const DSDesignSpace *ds, const DSUInteger firstTerm, ...);

extern const DSGMASystem * DSDesignSpaceGMASystem(const DSDesignSpace * ds);
extern const DSDictionary * DSDesignSpaceCyclicalCaseDictionary(const DSDesignSpace *ds);
//extern DSDictionary * DSDesignSpaceCycleDictionaryForSignature(const DSDesignSpace * ds, const DSUInteger * signature);

extern const char * DSDesignSpaceCasePrefix(const DSDesignSpace * ds);

#if defined (__APPLE__) && defined (__MACH__)
#pragma mark - Utility -
#endif

extern void * DSDesignSpaceTermListForAllStrings(const char ** strings, const DSUInteger numberOfEquations);
extern void DSDesignSpaceAddConstraints(DSDesignSpace * ds, const char ** strings, DSUInteger numberOfConstraints);
extern void DSDesignSpacePrint(const DSDesignSpace * ds);

#if defined (__APPLE__) && defined (__MACH__)
#pragma mark Case and Case validity
#endif

extern DSUInteger * DSDesignSpaceCaseNumbersWithPrefix(const DSDesignSpace * ds, const DSUInteger sizeOfPrefix, const DSUInteger *prefix, DSUInteger * numberOfCases)
;

extern DSCase ** DSDesignSpaceCalculateCases(DSDesignSpace *ds, const DSUInteger numberOfCase, DSUInteger *cases);
extern DSCase ** DSDesignSpaceCalculateValidCasesByPrunning(DSDesignSpace *ds);
extern DSCase ** DSDesignSpaceCalculateAllValidCases(DSDesignSpace *ds);
extern DSDictionary * DSDesignSpaceCalculateAllValidCasesForSliceByResolvingCyclicalCases(DSDesignSpace *ds,
                                                                                          const DSVariablePool * lower,
                                                                                          const DSVariablePool * upper,
                                                                                          bool strict);
extern DSDictionary * DSDesignSpaceCalculateAllValidCasesByResolvingCyclicalCases(DSDesignSpace *ds);
extern DSDictionary * DSDesignSpaceCalculateAllValidCasesForSlice(DSDesignSpace *ds, const DSVariablePool *lower, const DSVariablePool *upper);
extern DSDictionary * DSDesignSpaceCalculateAllValidCasesForSliceNonStrict(DSDesignSpace *ds, const DSVariablePool *lower, const DSVariablePool *upper);

#if defined (__APPLE__) && defined (__MACH__)
#pragma mark Cyclical Cases and Cyclical Case validity
#endif

extern DSUInteger DSDesignSpaceNumberOfCyclicalCases(const DSDesignSpace * ds);
extern const DSCyclicalCase * DSDesignSpaceCyclicalCaseWithCaseNumber(const DSDesignSpace *ds, DSUInteger caseNumber);
extern const DSCyclicalCase * DSDesignSpaceCyclicalCaseWithCaseIdentifier(const DSDesignSpace * ds, const char * identifer);
extern void DSDesignSpaceCalculateCyclicalCase(DSDesignSpace *ds, DSCase * aCase);
extern void DSDesignSpaceCalculateCyclicalCases(DSDesignSpace *ds);

extern void DSDesignSpaceCalculateValidityOfCases(DSDesignSpace *ds);
extern DSDictionary * DSDesignSpaceCalculateValidityOfCaseSet(DSDesignSpace *ds, DSUInteger numberOfCases, DSCase ** cases);

#if defined (__APPLE__) && defined (__MACH__)
#pragma mark Unstable Cases
#endif

extern DSUInteger * DSBlowUpSignatureForCaseNumber(const DSUInteger caseNumber, const DSVariablePool *freeDependentVariable);
extern DSUnstableCase * DSCalculateUnstableCase(DSDesignSpace *ds, DSCase *aCase);
extern void DSDesignSpaceCalculateUnstableCases(DSDesignSpace *ds);
extern void DSDesignSpaceCalculateUnstableCase(DSDesignSpace *ds, DSCase *aCase);
extern const bool DSCaseIsCyclical(const DSCase *aCase);


#if defined(__APPLE__) && defined (__MACH__)
#pragma mark - Data Serialization
#endif


extern DSDesignSpaceMessage * DSDesignSpaceEncode(const DSDesignSpace * aCase);
extern DSExtensionDataMessage * DSExtensionDataEncode(const DSCycleExtensionData *extensionData);
extern DSVectorMessage * DSVectorEncode(const DSUInteger * data, DSUInteger n);
extern DSMassBalanceDataMessage * DSMassBalanceEncode(const DSMassBalanceData *data);
extern DSDesignSpace * DSDesignSpaceFromDesignSpaceMessage(const DSDesignSpaceMessage * message);
extern DSDesignSpace * DSDesignSpaceDecode(size_t length, const void * buffer);

#ifdef __cplusplus
__END_DECLS
#endif

#endif
