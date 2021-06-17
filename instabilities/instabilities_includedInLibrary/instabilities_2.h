//
//  designspacetest.c
//
//
//  Created by Miguel Valderrama on 10/23/18.
//
//

// function headers to be added:

// DSUnstableCase.h

// macros.
//#define DSDSUnstable(x)                         ((x)->unstableCases)
//
//
//void dsUnstableCaseGetAdditionalConstraintMatrices(DSUnstableCase *uCase, const DSUInteger *bSignature);                    //added
//void DSUnstableCaseExpandConditionMatrices(DSUnstableCase *uCase);                                                          //added
//void DSUnstableCaseGetXd_b(DSUnstableCase *uCase, DSMatrixArray *gaussArray);                                               //added
//void DSUnstableCaseIdentifyBlowingDependentVariables(DSUnstableCase *uCase);                                                //added
//void DSUnstableCasePrintLinearProblem(glp_prob *lp, const DSUnstableCase *uCase, const DSUInteger *bSignature);             //added
//void DSUnstableCaseDetermineBlowingBehavior(glp_prob *lp, DSUnstableCase *uCase, const DSUInteger *bSignature);             //added
//void DSUnstableCaseCreateBoundaryMatrices(DSUnstableCase *uCase);                                                           //added
//
//
//
//DSMatrixArray * dsUnstableCaseLinearProblemAddEqualityConstraints(const DSMatrix *A,const DSMatrix *B, const DSUnstableCase *uCase);    //added
//
//
//DSMatrixArray * dsUnstableCaseLinearProblemAddKnifeEdgeConditions(const DSMatrix *A_,
//                                                            const DSMatrix *B_, const DSUnstableCase *uCase);               //added
//
//const DSMatrix *DSUnstableCaseGetSubSetPseudoInverse(const DSVariablePool *Xd_t, const DSVariablePool *Xd_e,                //added
//                                                     const DSVariablePool *Xd_b,
//                                                     const DSMatrix *Ad);
//
//
//
//const bool DSUnstableCaseConditionsAreValid(DSUnstableCase *uCase, const DSUInteger *bSignature);                           //added
//
//glp_prob * dsUnstableCaseLinearProblemForMatrices(const DSMatrix *A, const DSMatrix *B, DSUnstableCase *uCase,              //added
//                                                  const DSUInteger *bSignature);
//
//glp_prob * dsUnstableCaseLinearProblemForCaseValidity(const DSMatrix * U, const DSMatrix *zeta, DSUnstableCase *uCase,      //added
//                                                      const DSUInteger * bSignature);
//
//
//
//// DSDesignSpace.h
//DSUInteger * DSBlowUpSignatureForCaseNumber(const DSUInteger caseNumber, const DSVariablePool *freeDependentVariable);      //added
//DSUnstableCase * DSCalculateUnstableCase(DSDesignSpace *ds, DSCase *aCase);                                                 //added
//void DSDesignSpaceCalculateUnstableCases(DSDesignSpace *ds);                                                                //added
//void DSDesignSpaceCalculateUnstableCase(DSDesignSpace *ds, DSCase *aCase);                                                  //added
//const bool DSCaseIsCyclical(const DSCase *aCase);                                                                           //added
//
//
////DSMatrix.h
//const DSMatrix * DSMatrixPseudoInverse(const DSMatrix *matrix);                                         // added
//DSUInteger * DSMatrixSwapRows(DSMatrix *matrix, DSUInteger rowA, DSUInteger rowB);                      // added
//DSMatrix * DSMatrixSwapRowsBack(DSMatrix *matrix, DSUInteger *swappingVector);                          // added
//DSMatrixArray * DSMatrixGaussElimination(const DSMatrix *Ad, const DSMatrix *Ai, const DSMatrix *B);    // added







// functions not to be added.
void DSGetVertices(DSUnstableCase *uCase);

DSVertices * DSUnstableCaseVerticesFor2DSlice(DSUnstableCase *uCase, const DSVariablePool * lowerBounds, const DSVariablePool *upperBounds, const char * xVariable, const char *yVariable);
void printValidityOfCases(DSDesignSpace * ds);
void printPool(const DSVariablePool *pool);
#define DSDSGMA(x)                              ((x)->gma)
void printDictionary(DSDictionary *dictionary);








