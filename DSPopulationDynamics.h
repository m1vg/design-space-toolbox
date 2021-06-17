//
//  DSPopulationDynamics.h
//  designspace
//
//  Created by Miguel Angel Valderrama Gomez on 8/12/20.
//

#ifndef PopulationDynamics_h
#define PopulationDynamics_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#include "DSCase.h"
#include "DSTypes.h"
#include "DSErrors.h"
#include "DSVariable.h"
#include "DSMemoryManager.h"
#include "DSMatrix.h"

#ifdef __cplusplus
__BEGIN_DECLS
#endif

#if defined(__APPLE__) && defined(__MACH__)
#pragma mark Utility functions
#endif

const extern double DSPopDynamicsMutationRateForTransition(const DSVariablePool *oPoint1,
                                                      const DSVariablePool *oPoint2,
                                                      const double l,
                                                      const double d,
                                                      const DSVariablePool * identity);

const DSVariablePool * DSPopDynamicsGetDifferences(const DSVariablePool *oPoint1,
                                             const DSVariablePool *oPoint2,
                                             const bool log_scale);

const DSVariablePool * DSPopDynamicsCalculateParametricMutationRates(const DSVariablePool *Differences,
                                                               const DSVariablePool *identity,
                                                               const double l,
                                                               const double d);

const double DSPopDynamicsCalculateMutationRate(const DSVariablePool * BiasFactor);


#ifdef __cplusplus
__END_DECLS
#endif

#endif /* PopulationDynamics_h */
