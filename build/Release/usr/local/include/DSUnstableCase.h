//
//  DSUnstableCase.h
//  designspace
//
//  Created by Miguel Angel Valderrama Gomez on 12/4/18.
//

#ifndef DSUnstableCase_h
#define DSUnstableCase_h

#include "DSTypes.h"
#include "DSErrors.h"
#include "DSDataSerialization.pb-c.h"
#include <glpk.h>


#define DSDSUnstable(x)                         ((x)->unstableCases)

extern void DSuCaseFree(DSUnstableCase *uCase);
extern bool DSuCaseIsValid(DSCase *aCase, bool strict);

extern void dsUnstableCaseGetAdditionalConstraintMatrices(DSUnstableCase *uCase,
                                                          const DSUInteger *bSignature);
extern DSUnstableCase * DSUnstableCaseAddBoundaryMatrices(DSCase *aCase, const DSDesignSpace *ds);
extern void DSUnstableCaseExpandConditionMatrices(DSUnstableCase *uCase);
extern void DSUnstableCaseGetKnife(DSUnstableCase *uCase, const DSMatrix *gAd);
extern void DSUnstableCaseGetXd_b(DSUnstableCase *uCase, DSMatrix *Ad);
extern void DSUnstableCaseIdentifyBlowingDependentVariables(DSUnstableCase *uCase);
extern void DSUnstableCasePrintLinearProblem(glp_prob *lp, const DSUnstableCase *uCase, const DSUInteger *bSignature);
extern void DSUnstableCaseDetermineBlowingBehavior(glp_prob *lp, DSUnstableCase *uCase, const DSUInteger *bSignature);
extern void DSUnstableCaseCreateBoundaryMatrices(DSUnstableCase *uCase);
extern void DSUnstableCaseCreateBoundaryMatrices_alt(DSUnstableCase *uCase);
extern void DSUnstableCaseCreateBoundaryMatrices_alt2(DSUnstableCase *uCase);



extern void DSUnstableCaseGetEqualityAndKnife(const DSMatrix *pInverse,
                                              const DSMatrixArray *gaussArray,
                                              DSUnstableCase *uCase);

extern DSMatrixArray * dsUnstableCaseLinearProblemAddEqualityConstraints(const DSMatrix *A,
                                                                         const DSMatrix *B,
                                                                         const DSUnstableCase *uCase);

extern DSMatrixArray * dsUnstableCaseLinearProblemAddKnifeEdgeConditions(const DSMatrix *A_,
                                                                         const DSMatrix *B_,
                                                                         const DSUnstableCase *uCase);

extern const DSMatrix *DSUnstableCaseGetSubSetPseudoInverse(const DSVariablePool *Xd_t, const DSVariablePool *Xd_e,
                                                            const DSVariablePool *Xd_b,
                                                            const DSMatrix *Ad);


extern const bool DSUnstableCaseConditionsAreValid(DSUnstableCase *uCase, const DSUInteger *bSignature);
extern glp_prob * dsUnstableCaseLinearProblemForMatrices(const DSMatrix *A,
                                                         const DSMatrix *B,
                                                         DSUnstableCase *uCase,
                                                         const DSUInteger *bSignature);

extern glp_prob * dsUnstableCaseLinearProblemForCaseValidity(const DSMatrix * U,
                                                             const DSMatrix *zeta,
                                                             DSUnstableCase *uCase,
                                                             const DSUInteger * bSignature);

extern const DSCase ** DSUnstableCaseListAllSubcases(const DSCase * aCase, const DSDesignSpace *ds);
extern DSUInteger DSUnstableCaseSubcasesCount(const DSCase * aCase, const DSDesignSpace *ds);


#endif /* DSUnstableCase_h */
