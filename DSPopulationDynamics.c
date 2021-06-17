//
//  PopulationDynamics.c
//  designspace
//
//  Created by Miguel Angel Valderrama Gomez on 8/12/20.
//

#include "DSPopulationDynamics.h"

#include <stdbool.h>
#include <string.h>
#include <math.h>

const extern double DSPopDynamicsMutationRateForTransition(const DSVariablePool *oPoint1,
                                                     const DSVariablePool *oPoint2,
                                                     const double l,
                                                     const double d,
                                                     const DSVariablePool * identity)
{
    
        /*This function calculates the phenotype-specific mutation rate for the transition from the operating point oPoint1 to oPoint2. Three major steps are involved:
         1. Determine differences between the two operating points
         2. Determine bias factor
         3. Calculate mutation rate as the product of parameters
         */
        
        double mutationRate = INFINITY;
        
        if (oPoint1 == NULL || oPoint2 == NULL || identity == NULL){
            DSError(M_DS_NULL, A_DS_ERROR);
            goto bail;
        }
        
        const DSVariablePool  * Differences = NULL;
        const DSVariablePool * BiasFactor = NULL;
        
        Differences = DSPopDynamicsGetDifferences(oPoint1, oPoint2, false);
        BiasFactor = DSPopDynamicsCalculateParametricMutationRates(Differences, identity, l, d);
        mutationRate = DSPopDynamicsCalculateMutationRate(BiasFactor);
        
        if (Differences != NULL)
            DSVariablePoolFree((DSVariablePool *)Differences);
        if (BiasFactor != NULL)
            DSVariablePoolFree((DSVariablePool *)BiasFactor);
        
bail:
        return mutationRate;
}


const DSVariablePool * DSPopDynamicsGetDifferences(const DSVariablePool *oPoint1,
                                             const DSVariablePool *oPoint2,
                                             const bool log_scale){
        DSVariablePool * Differences = NULL;
        double diff;
        const DSVariable ** names;
    
        
        if (oPoint1 == NULL || oPoint2 == NULL){
            DSError(M_DS_NULL, A_DS_ERROR);
            goto bail;
        }
    
        if (DSVariablePoolIsReadWriteAdd(oPoint1) == false)
            DSVariablePoolSetReadWriteAdd((DSVariablePool *) oPoint1);
        if (DSVariablePoolIsReadWriteAdd(oPoint2) == false)
            DSVariablePoolSetReadWriteAdd((DSVariablePool *) oPoint2);
        
        if (DSVariablePoolNumberOfVariables(oPoint1) != DSVariablePoolNumberOfVariables(oPoint2)){
            DSError(M_DS_WRONG "Number of variables does not agree", A_DS_ERROR);
            goto bail;
        }
        
        DSUInteger i;
        Differences = DSVariablePoolCopy(oPoint1);
        names = DSVariablePoolAllVariables(oPoint1);
    
        for (i=0; i<DSVariablePoolNumberOfVariables(oPoint1); i++){
                if (log_scale == false){
                    
                    diff = DSVariablePoolValueForVariableWithName(oPoint2, DSVariableName(names[i]))/DSVariablePoolValueForVariableWithName(oPoint1, DSVariableName(names[i]));
                    DSVariablePoolSetValueForVariableWithName(Differences,
                                                              DSVariableName(names[i]),
                                                              diff);
                }else{
                    diff = log10(DSVariablePoolValueForVariableWithName(oPoint2, DSVariableName(names[i]))) - log10(DSVariablePoolValueForVariableWithName(oPoint1, DSVariableName(names[i])));
                    
                    DSVariablePoolSetValueForVariableWithName(Differences,
                                                              DSVariableName(names[i]),
                                                              diff);
                }
        }
bail:
    return Differences;
}

const DSVariablePool * DSPopDynamicsCalculateParametricMutationRates(const DSVariablePool *Differences,
                                                               const DSVariablePool *identity,
                                                               const double l,
                                                               const double d){
    
        DSVariablePool * ParametricMutationRates = NULL;
        DSUInteger i, c;
        double bias_factor, difference;
        const DSVariable ** names;
        bool increasing = true;
    
    
        if (Differences == NULL || identity == NULL){
            DSError(M_DS_NULL, A_DS_ERROR);
            goto bail;
        }
    
        if (DSVariablePoolNumberOfVariables(Differences) != DSVariablePoolNumberOfVariables(identity)){
            DSError(M_DS_WRONG "Number of variables does not agree", A_DS_ERROR);
            goto bail;
        }
    
        /* This function calculates the bias factor that considers the magnitude of the difference and the direction of change. It uses the identity of the pool to decide how to calculate the bias. Each type of parameter has two options: increase or decrease.
         
             Case 1: K
             Case 2: b (beta)
             Case 3: amax (alpha_max)
             Case 4: amin (alpha_min)
         
         */
        names = DSVariablePoolAllVariables(Differences);
        ParametricMutationRates = DSVariablePoolCopy(Differences);
    
        
        for(i=0; i<DSVariablePoolNumberOfVariables(Differences); i++){
            difference = DSVariablePoolValueForVariableWithName(Differences,DSVariableName(names[i]));
            c = DSVariablePoolValueForVariableWithName(identity, DSVariableName(names[i]));
            
            if (difference > 1.0)
                increasing = true;
            else
                increasing = false; 
            if (d == 0.0)
                c = 6;
            
            difference = fabs(log10(difference));
                
                    switch (c) {
                            
                        case 1: // K
                            if (increasing == true){ // exp(-d/(lambda * delta))
                                bias_factor = exp( -difference/(l * d) );
                            }else{  // exp(-d/(lambda / delta^2))
                                bias_factor = exp(-difference / (l / (d )) );
                            }
                            DSVariablePoolSetValueForVariableWithName(ParametricMutationRates,
                                                                      DSVariableName(names[i]),
                                                                      bias_factor);
                            break;
                            
                        case 2: // beta
                            if (increasing == true){ // exp(-d/(lambda * delta))
                                bias_factor = exp(-difference/(l * d ));
                            }else{  // exp(-d/(lambda / delta))
                                bias_factor = exp(-difference/(l / d ));
                            }
                            DSVariablePoolSetValueForVariableWithName(ParametricMutationRates,
                                                                      DSVariableName(names[i]),
                                                                      bias_factor);
                            break;
                            
                            
                        case 3: // a_max
                            if (increasing == true){ // exp(-d/(lambda / delta))
                                bias_factor = exp(-difference/(l / d ));
                            }else{ // exp(-d/(lambda * delta))
                                bias_factor = exp(-difference/(l * d ));
                            }
                            DSVariablePoolSetValueForVariableWithName(ParametricMutationRates,
                                                                      DSVariableName(names[i]),
                                                                      bias_factor);
                            break;
                            
                        case 4: // a_min
                            if (increasing == true){ // exp(-d/(lambda * delta))
                                bias_factor = exp(-difference/(l * d ));
                                
                            }else{ // exp(-d/(lambda / delta))
                                bias_factor = exp(-difference/(l / d ));
                            }
                            DSVariablePoolSetValueForVariableWithName(ParametricMutationRates,
                                                                      DSVariableName(names[i]),
                                                                      bias_factor);
                            break;
                            
                        case 5: // a
                            if (increasing == true){ // exp(-d/(lambda / delta))
                                
                                bias_factor = exp(-difference/(l / d ));
                                
                            }else{ // exp(-d/(lambda * delta))
                                
                                bias_factor = exp(-difference/(l * d ));
                            }
                            DSVariablePoolSetValueForVariableWithName(ParametricMutationRates,
                                                                      DSVariableName(names[i]),
                                                                      bias_factor);
                            break;
                            
                        case 6: // a different kind
                            printf("Warning, function DSPopDynamicsCalculateParametricMutationRates is using case 6 \n");
                            bias_factor = exp(-difference/l);
                            DSVariablePoolSetValueForVariableWithName(ParametricMutationRates,
                                                                      DSVariableName(names[i]),
                                                                      bias_factor);
                            break;
                        case 7: // when d == 0
                            bias_factor = 1.0;
                            DSVariablePoolSetValueForVariableWithName(ParametricMutationRates,
                                                                      DSVariableName(names[i]),
                                                                      bias_factor);
                            
                            break;
                            
                        default:
                            break;
                    }
        }
    bail:
        return ParametricMutationRates;
}

const double DSPopDynamicsCalculateMutationRate(const DSVariablePool * BiasFactor){
    
        double mutationRate = 1.0;
        
        if (BiasFactor == NULL){
            DSError(M_DS_NULL, A_DS_ERROR);
            goto bail;
        }
    
        DSUInteger i;
        const DSVariable ** names;
        names = DSVariablePoolAllVariables(BiasFactor);
    
    
        for (i=0; i<DSVariablePoolNumberOfVariables(BiasFactor); i++){
            mutationRate = DSVariablePoolValueForVariableWithName(BiasFactor,
                                                                  DSVariableName(names[i]))*mutationRate;
        }
    
        bail:
            return mutationRate;
}

