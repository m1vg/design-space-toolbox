//
//  designspacetest.c
//
//
//  Created by Miguel Valderrama on 08/28/18.
//
//
#include <stdio.h>
#include <string.h>
#include <designspace/DSStd.h>
#include <gsl/gsl_matrix_double.h>


int main(int argc, const char ** argv) {
//        int i;
//        char * strings[3] = {"\0"};
//        strings[0] = strdup("x1. = a11 + 2*b31*x3 - b11*x1 - 2*b12*(x1^2)");
//        strings[1] = strdup("x2. = b12*(x1^2) - b23*x2 - b22*x2");
//        strings[2] = strdup("x3. = a31 + b23*x2 - b31*x3 - b33*x3");
//        printf("The equations are: %s\n %s\n %s\n", *strings, *(strings+1), *(strings+2) );
//        printf("\n aqui \n");


        char * strings[4] = {"\0"};
        strings[0] = strdup("x1. =  + x4*k7 + 2*a31 - x1*b11 - x1*k8 - 2*x2*b22");
        strings[1] = strdup("x2. =  + k1 - x4*k6");
        strings[2] = strdup("x3. =  + x2*b23 - x3*b31");
        strings[3] = strdup("x4.= x1*k8-x4*k7");
        printf("The equations are: %s\n %s\n %s\n %s\n", *strings, *(strings+1), *(strings+2), *(strings+3));

        DSDesignSpace * ds;
//        DSExpression ** expr = NULL;
//        DSGMASystem *gma = NULL;

        
        ds = DSDesignSpaceByParsingStrings(strings, NULL, 4);
        DSDesignSpaceCalculateCyclicalCases(ds);
        printf("Number of valid cases is: %i\n", DSDesignSpaceNumberOfValidCases(ds));
        printf("DSDesignSpaceNumberOfValidCases passed!\n");

        // Now let's generate a case.

        DSCase *aCase = NULL;
        DSUInteger * termSignature, caseNumber = 2; // old case number 27

        termSignature = DSCaseSignatureForCaseNumber(caseNumber, ds->gma);
        aCase = DSCaseWithTermsFromDesignSpace(ds, termSignature, DSDesignSpaceCasePrefix(ds));

        printf("The signature of the case is: %u \n", *aCase->signature);
        printf("The case number is: %u \n", aCase->caseNumber);
        printf("The case identifier is: %s \n", aCase->caseIdentifier);
        printf("The term signature used to generate the case was: %u %u %u %u %u %u %u %u \n", *(termSignature),
        *(termSignature+1), *(termSignature+2), *(termSignature+3), *(termSignature+4), *(termSignature+5), *(termSignature+6), *(termSignature+7) );

        // header for the print function for condition matrices Cd, Ci, U, delta, zeta.
        void printMatrix(DSMatrix *Matrix);
        void printPool(DSVariablePool *pool);
        void printMatrixArray(DSMatrixArray * Array);

        // let's take a look at the condition matrix Cd
        printf("The condition Matrix Cd is: \n");
        printMatrix(aCase->Cd);

        printf("The condition Matrix Ci is: \n");
        printMatrix(aCase->Ci);


        // Let's check the matrix A
        DSMatrix *nullspace = NULL, *A = NULL;
        A = DSSSystemA(aCase->ssys);
        printf("The Matrix A is: \n");
        printMatrix(A);

        // Now let's calculate the left null space of A
        nullspace = DSMatrixLeftNullspace(A);
        printf("The left null space of A is: \n");
        printMatrix(nullspace);

        // and the problematic rows are:
        DSMatrix *problematic = NULL;
        problematic = DSMatrixIdenticalRows(nullspace);
        printf("The problematic equations are: \n");
        printMatrix(problematic);

        //the problematic equations are:
        DSMatrix * problematicEquations = NULL;
        DSMatrix * dsSubcaseProblematicEquations(const DSCase * aCase);

        problematicEquations = dsSubcaseProblematicEquations(aCase);

        // Now let's calculate the problematic terms
        DSMatrixArray * problematicTerms = NULL;
        DSMatrixArray * dsSubcaseProblematicTerms(const DSCase *aCase, const DSMatrix * dependentEquations);

        problematicTerms = dsSubcaseProblematicTerms(aCase, problematicEquations);

        printf("The problematic equations are: \n");
        printMatrix(problematicEquations);

        printf("The Problematic Terms are: \n");
        printMatrixArray(problematicTerms);

        // Now let's calculate the coefficientArray
        DSMatrixArray * coefficientArray = NULL;
        DSMatrixArray * dsSubcaseCoefficientsOfInterest(const DSCase * aCase, const DSMatrixArray * problematicTerms);

        coefficientArray = dsSubcaseCoefficientsOfInterest(aCase, problematicTerms);
        printf("The Array of coefficients is: \n");
        printMatrixArray(coefficientArray);

        // now let's see how a subcase is built!

        DSDesignSpace * dsCyclicalCaseCollapsedSystem(const DSCase * aCase,
                                              const DSDesignSpace * original,
                                              DSMatrix * problematicEquations,
                                              const DSMatrixArray * coefficientArray);
        DSDesignSpace * subcase = NULL;
//        subcase = dsCyclicalCaseCollapsedSystem(aCase, ds, problematicEquations, coefficientArray);

        char ** systemEquations = NULL;
        DSCycleExtensionData * extensionData;

        DSCycleExtensionData * dsCycleExtensionDataInitForCyclicalCase(const DSCase * aCase,
                                                               const DSDesignSpace * original);

//        extensionData = dsCycleExtensionDataInitForCyclicalCase(aCase, ds);

        char ** dsCyclicalCaseEquationsSplitVariables(const DSCase * aCase,
                                                     const DSDesignSpace * original,
                                                     DSMatrix * problematicEquations,
                                                     const DSMatrixArray * coefficientArray,
                                                     DSCycleExtensionData * extensionData);



//        systemEquations = dsCyclicalCaseEquationsSplitVariables(aCase,
//                                                                ds,
//                                                                problematicEquations,
//                                                                coefficientArray,
//                                                                extensionData);

        DSExpression ** dsCyclicalCaseEquationsForCycle(const DSCase * aCase,
                                                       const DSDesignSpace * original,
                                                       const DSMatrixArray * coefficientArray,
                                                       const DSUInteger cycleNumber,
                                                       const DSUInteger primaryCycleVariable,
                                                       const DSUInteger numberSecondaryVariables,
                                                       const DSUInteger * secondaryCycleVariables,
                                                       const DSMatrix * LI,
                                                       const DSMatrix * Lc,
                                                       const DSMatrix * MBn,
                                                       const DSVariablePool * yn,
                                                       const DSVariablePool * yc);


char * dsCyclicalCaseEquationForFlux(const DSCase * aCase,
                                            const DSDesignSpace * original,
                                            const DSMatrixArray * coefficientArray,
                                            const DSUInteger variableIndex,
                                            const DSUInteger fluxIndex,
                                            const bool positiveFlux,
                                            const DSUInteger cycleNumber,
                                            const DSUInteger primaryVariable,
                                            const DSUInteger numberSecondaryVariables,
                                            const DSUInteger * secondaryCycleVariables,
                                            const DSMatrix * LI,
                                            const DSMatrix * Lc,
                                            const DSMatrix * MBn,
                                            const DSVariablePool * yc);


void dsCyclicalCaseSolutionOfPartitionedMatrices(const DSCase * aCase,
                                                        const DSUInteger numberOfSecondaryVariables,
                                                        const DSUInteger * secondaryVariables,
                                                        DSMatrix ** LI,
                                                        DSMatrix **Lc,
                                                        DSMatrix **MBn,
                                                        DSVariablePool ** yn,
                                                        DSVariablePool ** yc);


void dsCyclicalCasePartitionSolutionMatrices(const DSCase * aCase,
                                                    const DSUInteger numberOfSecondaryVariables,
                                                    const DSUInteger * secondaryVariables,
                                                    DSMatrix ** ADn,
                                                    DSMatrix ** ADc,
                                                    DSMatrix ** AIn,
                                                    DSMatrix ** Bn,
                                                    DSVariablePool ** yn,
                                                    DSVariablePool ** yc);


DSUInteger dsCyclicalCaseSecondaryCycleVariableIndicesForCycle(DSMatrix * problematicEquations,
                                                                      DSUInteger cycleNumber,
                                                                      DSUInteger * primaryVariables,
                                                                      DSUInteger ** cycleIndices);


DSUInteger dsCyclicalCaseAllSecondaryCycleVariables(const DSMatrix * problematicEquations,
                                                           const DSMatrixArray * coefficientArray,
                                                           const DSUInteger numberOfCycles,
                                                           const DSUInteger * primaryVariables,
                                                           DSUInteger ** cycleIndices,
                                                           double ** coefficients);


DSUInteger dsCyclicalCasePrimaryCycleVariableIndices(const DSCase * aCase,
                                                            DSMatrix * problematicEquations,
                                                            DSUInteger ** primaryVariables);


DSMatrix * dsCyclicalCaseExpandLcMatrix(const DSVariablePool * Xd,
                                               const DSMatrix * Lc,
                                               const DSVariablePool * yc);

DSStack * dsCyclicalCaseCreateAugmentedSystems(const DSCase * aCase,
                                                      const DSDesignSpace * original,
                                                      DSMatrix * problematicEquations,
                                                      const DSMatrixArray * problematicTerms,
                                                      const DSMatrixArray * coefficientArray);

DSUInteger dsCyclicalCaseNumberOfAugmentedSystems(const DSDesignSpace * original,
                                                         const DSMatrix * problematicEquations);
void dsCyclicalCaseAugmentedEquationsForCycleALT(char ** systemEquations,
                                                     const DSCase * aCase,
                                                     const DSDesignSpace * original,
                                                     const DSMatrix * problematicMatrix,
                                                     const DSMatrixArray * coefficientArray,
                                                     const DSUInteger cycleNumber,
                                                     const DSUInteger primaryCycleVariable,
                                                     const DSUInteger numberSecondaryVariables,
                                                     const DSUInteger * secondaryVariables,
                                                     const DSMatrix * LI,
                                                     const DSMatrix * Lc,
                                                     const DSMatrix * Mb,
                                                     const DSVariablePool * yn,
                                                     const DSVariablePool * yc);

                   DSDesignSpace * dsCyclicalCaseAugmentedSystemForSubdominantDecays(const DSCase * aCase,
                                                                         const DSDesignSpace * original,
                                                                         DSMatrix * problematicEquations,
                                                                         const DSMatrixArray * problematicTerms,
                                                                         const DSMatrixArray * coefficientArray,
                                                                         const DSUInteger *subdominantDecaySpecies,
                                                                         const DSUInteger *subdominantDecayTerm);


                   char ** dsCyclicalCaseOriginalEquationsWithEquilibriumConstraintsALT(const DSCase * aCase,
                                                                         const DSDesignSpace * original,
                                                                         const DSUInteger numberSecondaryVariables,
                                                                         const DSUInteger * secondaryVariables,
                                                                         const double * coefficientMultipliers,
                                                                         const DSMatrix * LI,
                                                                         const DSMatrix * Lc,
                                                                         const DSMatrix * Mb,
                                                                         const DSVariablePool * yn,
                                                                         const DSVariablePool * yc);

void dsAddConstraintsForSubdominantDecays(DSDesignSpace * subcase, const DSCase * aCase, const DSDesignSpace * original, const DSMatrix * problematicEquations, const DSMatrixArray * coefficientArray, const DSUInteger * subdominantDecays, const DSUInteger * subdominantDecayTerms);

DSDesignSpace * dsCyclicalCaseCreateUniqueAugmentedSystem(const DSCase *aCase, const DSGMASystem * modifiedGMA, const DSMatrix  * problematicEquations, const DSExpression ** augmentedEquations, const DSUInteger * subdominantDecays);

void dsCyclicalCaseEquilibriumEquationForVariableALT(DSUInteger index,
                                                         char ** systemEquations,
                                                         const DSCase * aCase,
                                                         const DSDesignSpace * original,
                                                         const DSMatrix * LI,
                                                         const DSMatrix * Lc,
                                                         const DSMatrix * Mb,
                                                         const DSVariablePool * yn,
                                                         const DSVariablePool * yc);


     subcase = dsCyclicalCaseCollapsedSystem(aCase, ds, problematicEquations, coefficientArray);

extensionData = dsCycleExtensionDataInitForCyclicalCase(aCase, ds);

systemEquations = dsCyclicalCaseEquationsSplitVariables(aCase,
                                                        ds,
                                                        problematicEquations,
                                                        coefficientArray,
                                                        extensionData);

       // now let's print the system equations.
       printf("System Equations are (from dsCyclicalCaseEquationsSplitVariables): \n");

       for (int y=0; y<4; y++)
       {


        for (int i=0; i < 100; i++)
        {
        printf("%c", *(*(systemEquations+y)+i ) );
        }
        printf("\n");

        }

        // Now the idea is to understand how these sytemEquations are being generated. In order to do that
        // we need to print the variables that are being generated in the function dsCyclicalCaseEquationsSplitVariables(...)

        DSMatrix * Mb = NULL, *LI = NULL, *Lc = NULL;
        DSMatrix * TLI = NULL, * TLc = NULL, * TMb = NULL;
        DSVariablePool * yn = NULL, *yc = NULL;
        char *name, ** systemEquations2 = NULL;
        const DSVariablePool *Xd;
        bool error=false;
        DSUInteger i, j, k, index, numberOfCycles, numberSecondaryVariables, *primaryVariables = NULL, *secondaryVariables = NULL;
        DSUInteger numberAllSecondaryVariables, * allSecondaryVariables = NULL;
        double  * coefficientMultipliers = NULL;


        numberOfCycles = dsCyclicalCasePrimaryCycleVariableIndices(aCase, problematicEquations, &primaryVariables);
        printf("The number of cycles is: %i \n", numberOfCycles);
        printf("The index of primary variable is: %u \n", (*primaryVariables));

        numberAllSecondaryVariables = dsCyclicalCaseAllSecondaryCycleVariables(problematicEquations,
                                                                               coefficientArray, numberOfCycles,
                                                                               primaryVariables,
                                                                               &allSecondaryVariables, &coefficientMultipliers);

        printf("The number of secondary variables is: %u \n", numberAllSecondaryVariables);
        printf("all secondary variables are: %u, %u \n", *(allSecondaryVariables + 0), *(allSecondaryVariables + 1) );
        printf("Coefficient multipliers: %f \n", *(coefficientMultipliers));

        i = 0;

        numberSecondaryVariables = dsCyclicalCaseSecondaryCycleVariableIndicesForCycle(problematicEquations, i, primaryVariables, &secondaryVariables);
        dsCyclicalCaseSolutionOfPartitionedMatrices(aCase, numberSecondaryVariables, secondaryVariables, &TLI, &TLc, &TMb, &yn, &yc);

        // let's print now the number of SecondaryVariables and the matrices TLI, TLc, TMb, yn and yc.

        printf("numberSecondaryVariables is: %u \n", numberSecondaryVariables);
        printf("secondaryVariables are: %u", *secondaryVariables );
        printf("The Matrix TLI is: \n");
        printMatrix(TLI);

        printf("The Matrix TLc is: \n");
        printMatrix(TLc);

        printf("The Matrix TMb is: \n");
        printMatrix(TMb);

        printf("The pool yn is: \n");
        printPool(yn);

        printf("The pool yc is: \n");
        printPool(yc);

        printf("The Matrix Ad of the case 27 is: \n");
        printMatrix(DSSSystemAd(aCase->ssys));

        printf("The Matrix Ai of the case 27 is: \n");
        printMatrix(DSSSystemAi(aCase->ssys));

        printf("The Matrix B of the case 27 is: \n");
        printMatrix(DSSSystemB(aCase->ssys));

        // Now let's take a look at the other matrices and pools

        DSMatrix * tempMatrix, *Ad, *Ai, *B, *ADn, *ADc, *AIn, *Bn, *Mn, *MBn;

        Ad = DSSSystemAd(aCase->ssys);
        Ai = DSSSystemAi(aCase->ssys);
        B = DSSSystemB(aCase->ssys);
        tempMatrix = DSMatrixSubMatrixIncludingRows(Ad, numberSecondaryVariables, secondaryVariables);
        ADc = DSMatrixSubMatrixExcludingColumns(tempMatrix, numberSecondaryVariables, secondaryVariables);
        ADn = DSMatrixSubMatrixIncludingColumns(tempMatrix, numberSecondaryVariables, secondaryVariables);
        AIn = DSMatrixSubMatrixIncludingRows(Ai, numberSecondaryVariables, secondaryVariables);
        Bn = DSMatrixSubMatrixIncludingRows(B, numberSecondaryVariables, secondaryVariables);


        printf("The Matrix temp is: \n");
        printMatrix(tempMatrix);

        printf("The Matrix ADc is: \n");
        printMatrix(ADc);

        printf("The Matrix ADn is: \n");
        printMatrix(ADn);

        printf("The Matrix AIn is: \n");
        printMatrix(AIn);

        printf("The Matrix Bn is: \n");
        printMatrix(Bn);

        printf("SecondaryVariables is: %u %u  \n", *(secondaryVariables), *(secondaryVariables+1) );

//        DSMatrix * Mb = NULL, *LI = NULL, *Lc = NULL;
        Mb = DSMatrixCalloc(numberAllSecondaryVariables, 1);
        LI = DSMatrixCalloc(numberAllSecondaryVariables, DSVariablePoolNumberOfVariables(DSCaseXi(aCase)));
        Lc = DSMatrixCalloc(numberAllSecondaryVariables, DSVariablePoolNumberOfVariables(DSCaseXd(aCase)));



        Mn = DSMatrixInverse(ADn);

        LI = DSMatrixByMultiplyingMatrix(Mn, AIn);
        Lc = DSMatrixByMultiplyingMatrix(Mn, ADc);
        MBn = DSMatrixByMultiplyingMatrix(Mn, Bn);




        printf("The Matrix LI is: \n");
        printMatrix(LI);

        printf("The Matrix Lc is: \n");
        printMatrix(Lc);

        printf("The Matrix MBn is: \n");
        printMatrix(MBn);



        DSMatrixMultiplyByScalar(LI, 1.);
        DSMatrixMultiplyByScalar(Lc, 1.);

        printf("The Matrix LI after is: \n");
        printMatrix(LI);

        printf("The Matrix Lc after is: \n");
        printMatrix(Lc);

        DSMatrix * Mb2 = NULL, *LI2 = NULL, *Lc2 = NULL;

        Mb2 = DSMatrixCalloc(numberAllSecondaryVariables, 1);
        LI2 = DSMatrixCalloc(numberAllSecondaryVariables, DSVariablePoolNumberOfVariables(DSCaseXi(aCase)));
        Lc2 = DSMatrixCalloc(numberAllSecondaryVariables, DSVariablePoolNumberOfVariables(DSCaseXd(aCase)));

                        for (j = 0; j < DSMatrixRows(LI); j++){
                        for (k = 0; k < numberAllSecondaryVariables; k++) {
                                if (secondaryVariables[j] == allSecondaryVariables[k]) {
                                        index = k;
                                        break;
                                }
                        }
                        DSMatrixSetDoubleValue(Mb2, index, 0, DSMatrixDoubleValue(Mb, j, 0));
                        for (k = 0; k < DSMatrixColumns(LI); k++) {
                                DSMatrixSetDoubleValue(LI2, index, k, DSMatrixDoubleValue(LI, j, k));
                        }
                        for (k = 0; k < DSMatrixColumns(TLc); k++) {
                                name = DSVariableName(DSVariablePoolVariableAtIndex(yc, k));
                                DSMatrixSetDoubleValue(Lc2, index, DSVariablePoolIndexOfVariableWithName(DSCaseXd(aCase), name),
                                                       DSMatrixDoubleValue(Lc, j, k));
                        }
                }

        printf("The Matrix LI2 after is: \n");
        printMatrix(LI2);

        printf("The Matrix Mb2 after is: \n");
        printMatrix(Mb2);

        printf("The Matrix Lc2 after is: \n");
        printMatrix(Lc2);

         yn = DSVariablePoolAlloc();
        yc = DSVariablePoolAlloc();
        for (i = 0; i < numberAllSecondaryVariables; i++) {
                DSVariablePoolAddVariableWithName(yn, DSVariableName(DSVariablePoolVariableAtIndex(DSCaseXd(aCase), allSecondaryVariables[i])));
        }
        for (i = 0; i < DSVariablePoolNumberOfVariables(DSCaseXd(aCase)); i++) {
                name = DSVariableName(DSVariablePoolVariableAtIndex(DSCaseXd(aCase), i));
                if (DSVariablePoolHasVariableWithName(yn, name) == false) {
                        DSVariablePoolAddVariableWithName(yc, name);
                }
        }

        char ** systemEquations3 = NULL;
        systemEquations3 = dsCyclicalCaseOriginalEquationsWithEquilibriumConstraintsALT(aCase, ds, numberAllSecondaryVariables, allSecondaryVariables, coefficientMultipliers, LI2, Lc2, Mb2, yn, yc);

        // now let's print the system equations.
       printf("\n aqui1 \n");
       printf("System Equations are (from dsCyclicalCaseOriginalEquationsWithEqConstraintsALT): \n");

       for (int y=0; y<4; y++)
       {
        for (int i=0; i < 100; i++)
        {
        printf("%c", *(*(systemEquations3+y)+i ) );
        }
        printf("\n");

        }



        for (i = 0; i < numberOfCycles; i++) {
        dsCyclicalCaseAugmentedEquationsForCycleALT(systemEquations3, aCase, ds, problematicEquations, coefficientArray, i, primaryVariables[i], numberAllSecondaryVariables, allSecondaryVariables, LI2, Lc2, Mb2, yn, yc);
        }

           // now let's print the system equations.
       printf("System Equations are (from dsCyclicalCaseAugmentedEquationsForCycleALT): \n");

       for (int y=0; y<4; y++)
       {


        for (int i=0; i < 100; i++)
        {
        printf("%c", *(*(systemEquations3+y)+i ) );
        }
        printf("\n");

        }






        printf("\n End Of File \n");




















        DSDesignSpaceFree(ds);

        return 0;
}









void printMatrix(DSMatrix *Matrix){

        for (int rows=0; rows < Matrix->rows; rows++)
                {
                    for(int columns=0; columns < Matrix->columns; columns++)
                    {
                    printf("%f     ", gsl_matrix_get(Matrix->mat, rows, columns));
                    }
                printf("\n");
        }


}

void printPool(DSVariablePool *pool){

        for(int j = 0; j < pool->numberOfVariables; j++)
        {
        printf("Name %s ; Value %f  \n", pool->variables[j]->name, pool->variables[j]->value );
        }
}


void printMatrixArray(DSMatrixArray * Array){

        for (int m=0; m < DSMatrixArrayNumberOfMatrices(Array); m++)
        {
        DSMatrix Matrix = *DSMatrixArrayMatrix(Array, m);
        printf("Printing Matrix Nr. %i \n", m);

        for (int rows=0; rows < Matrix.rows; rows++)
        {
            for(int columns=0; columns < Matrix.columns; columns++)
            {
            printf("%f     ", gsl_matrix_get(Matrix.mat, rows, columns));
            }
        printf("\n");
        }

}
}

char ** dsCyclicalCaseEquationsSplitVariables(const DSCase * aCase,
                                                     const DSDesignSpace * original,
                                                     DSMatrix * problematicEquations,
                                                     const DSMatrixArray * coefficientArray,
                                                     DSCycleExtensionData * extensionData)
{
        DSMatrix * Mb = NULL, *LI = NULL, *Lc = NULL;
        DSMatrix * TLI = NULL, * TLc = NULL, * TMb = NULL;
        DSVariablePool * yn = NULL, *yc = NULL;
        char *name, ** systemEquations = NULL;
        const DSVariablePool *Xd;
        bool error=false;
        DSUInteger i, j, k, index, numberOfCycles, numberSecondaryVariables, *primaryVariables = NULL, *secondaryVariables = NULL;
        DSUInteger numberAllSecondaryVariables, * allSecondaryVariables = NULL;
        double  * coefficientMultipliers = NULL;
        //        DSCycleExtensionData * newExtensionData;
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
        if (coefficientArray == NULL)
                goto bail;
        numberOfCycles = dsCyclicalCasePrimaryCycleVariableIndices(aCase, problematicEquations, &primaryVariables);
        numberAllSecondaryVariables = dsCyclicalCaseAllSecondaryCycleVariables(problematicEquations,
                                                                               coefficientArray, numberOfCycles,
                                                                               primaryVariables,
                                                                               &allSecondaryVariables, &coefficientMultipliers);
        if (primaryVariables == NULL) {
                goto bail;
        }

        Mb = DSMatrixCalloc(numberAllSecondaryVariables, 1);
        LI = DSMatrixCalloc(numberAllSecondaryVariables, DSVariablePoolNumberOfVariables(DSCaseXi(aCase)));
        Lc = DSMatrixCalloc(numberAllSecondaryVariables, DSVariablePoolNumberOfVariables(DSCaseXd(aCase)));
        for (i = 0; i < numberOfCycles; i++) {
                numberSecondaryVariables = dsCyclicalCaseSecondaryCycleVariableIndicesForCycle(problematicEquations, i, primaryVariables, &secondaryVariables);
                dsCyclicalCaseSolutionOfPartitionedMatrices(aCase, numberSecondaryVariables, secondaryVariables, &TLI, &TLc, &TMb, &yn, &yc);
                if (TLI == NULL) {
                        DSSecureFree(secondaryVariables);
                        goto bail;
                }
                for (j = 0; j < DSMatrixRows(TLI); j++){
                        for (k = 0; k < numberAllSecondaryVariables; k++) {
                                if (secondaryVariables[j] == allSecondaryVariables[k]) {
                                        index = k;
                                        break;
                                }
                        }
                        DSMatrixSetDoubleValue(Mb, index, 0, DSMatrixDoubleValue(TMb, j, 0));
                        for (k = 0; k < DSMatrixColumns(TLI); k++) {
                                DSMatrixSetDoubleValue(LI, index, k, DSMatrixDoubleValue(TLI, j, k));
                        }
                        for (k = 0; k < DSMatrixColumns(TLc); k++) {
                                name = DSVariableName(DSVariablePoolVariableAtIndex(yc, k));
                                DSMatrixSetDoubleValue(Lc, index, DSVariablePoolIndexOfVariableWithName(DSCaseXd(aCase), name),
                                                       DSMatrixDoubleValue(TLc, j, k));
                        }
                }
                DSMatrixFree(TLc);
                DSMatrixFree(TLI);
                DSMatrixFree(TMb);
                DSVariablePoolFree(yc);
                DSVariablePoolFree(yn);
                DSSecureFree(secondaryVariables);
        }
        yn = DSVariablePoolAlloc();
        yc = DSVariablePoolAlloc();
        for (i = 0; i < numberAllSecondaryVariables; i++) {
                DSVariablePoolAddVariableWithName(yn, DSVariableName(DSVariablePoolVariableAtIndex(DSCaseXd(aCase), allSecondaryVariables[i])));
        }
        for (i = 0; i < DSVariablePoolNumberOfVariables(DSCaseXd(aCase)); i++) {
                name = DSVariableName(DSVariablePoolVariableAtIndex(DSCaseXd(aCase), i));
                if (DSVariablePoolHasVariableWithName(yn, name) == false) {
                        DSVariablePoolAddVariableWithName(yc, name);
                }
        }

//        goto bail;
//        if (secondaryVariables == NULL && numberSecondaryVariables > 0) {
//                goto bail;
//        }
//        if (numberSecondaryVariables > 0) {
//                dsCyclicalCaseSolutionOfPartitionedMatrices(aCase, numberSecondaryVariables, secondaryVariables, &LI, &Lc, &Mb, &yn, &yc);
//                if (DSCaseNumber(aCase) == 49) {
//                        printf("49 LI:\n");
//                        DSMatrixPrint(LI);
//                }
//        }
        //        systemEquations = dsCyclicalCaseOriginalEquationsWithEquilibriumConstraints(aCase, original, numberSecondaryVariables, secondaryVariables, LI, Lc, Mb, yn, yc);
        //        systemEquations = dsCyclicalCaseOriginalCaseEquationsWithEquilibriumConstraints(aCase, numberSecondaryVariables, secondaryVariables, LI, Lc, Mb, yn, yc);
        //        extensionData = DSSecureCalloc(sizeof(extensionData), 1);
        systemEquations = dsCyclicalCaseOriginalEquationsWithEquilibriumConstraintsALT(aCase, original, numberAllSecondaryVariables, allSecondaryVariables, coefficientMultipliers, LI, Lc, Mb, yn, yc);
        Xd = DSGMASystemXd(DSDesignSpaceGMASystem(original));
        //        extensionData->numberCycles = numberOfCycles;
        //        extensionData->cycleVariables = DSSecureCalloc(sizeof(DSUInteger), numberOfCycles);
        //        extensionData->fluxEquations =  DSSecureCalloc(sizeof(DSExpression **), numberOfCycles);
        //        extensionData->fluxIndex =  DSSecureCalloc(sizeof(DSUInteger *), numberOfCycles);
        //        extensionData->numberOfFluxes = DSSecureCalloc(sizeof(DSUInteger), numberOfCycles);
        if (systemEquations == NULL) {
                goto bail;
        }
        for (i = 0; i < DSDesignSpaceNumberOfEquations(original); i++) {
                if (systemEquations[i] == NULL) {
                        error = true;
                }
        }
        if (error == true) {
                for (i = 0; i < DSDesignSpaceNumberOfEquations(original); i++) {
                        if (systemEquations[i])
                                DSSecureFree(systemEquations[i]);
                }
                DSSecureFree(systemEquations);
                goto bail;
        }
        for (i = 0; i < numberOfCycles; i++) {
                dsCyclicalCaseAugmentedEquationsForCycleALT(systemEquations, aCase, original, problematicEquations, coefficientArray, i, primaryVariables[i], numberAllSecondaryVariables, allSecondaryVariables, LI, Lc, Mb, yn, yc);
                //                printf("%s\n", systemEquations[primaryVariables[i]]);
                //                numberSecondaryVariables = dsCyclicalCaseSecondaryCycleVariableIndicesForCycle(problematicEquations, i, primaryVariables, &secondaryVariables);
                //                cycleEquations = dsCyclicalCaseEquationsForCycle(aCase,
                //                                                                 original,
                //                                                                 coefficientArray,
                //                                                                 i,
                //                                                                 primaryVariables[i],
                //                                                                 numberSecondaryVariables,
                //                                                                 secondaryVariables,
                //                                                                 LI,
                //                                                                 Lc,
                //                                                                 Mb,
                //                                                                 yn,
                //                                                                 yc);
                //                if (cycleEquations == NULL) {
                //                        DSSecureFree(secondaryVariables);
                //                        break;
                //                }
                //                dsExtensionDataForCycle(aCase, original, extensionData, i, primaryVariables[i], numberSecondaryVariables, secondaryVariables);
                //                for (j = 0; j < numberSecondaryVariables+1; j++) {
                //                        if (j == 0) {
                //                                index = primaryVariables[i];
                //                        } else {
                //                                index = secondaryVariables[j-1];
                //                        }
                //                        DSSecureFree(systemEquations[index]);
                //                        systemEquations[index] = DSExpressionAsString(cycleEquations[j]);
                //                        DSExpressionFree(cycleEquations[j]);
                //                }
                //                if (secondaryVariables != NULL) {
                //                        DSSecureFree(cycleEquations);
                //                        DSSecureFree(secondaryVariables);
                //                }
        }
        if (i != numberOfCycles) {
                for (i = 0; i < DSDesignSpaceNumberOfEquations(original); i++) {
                        DSSecureFree(systemEquations[i]);
                }
                DSSecureFree(systemEquations);
                systemEquations = NULL;
        }
bail:
        if (primaryVariables != NULL)
                DSSecureFree(primaryVariables);
        if (allSecondaryVariables != NULL)
                DSSecureFree(allSecondaryVariables);
        if (coefficientMultipliers != NULL)
                DSSecureFree(coefficientMultipliers);
        if (LI != NULL)
                DSMatrixFree(LI);
        if (Lc != NULL)
                DSMatrixFree(Lc);
        if (Mb != NULL)
                DSMatrixFree(Mb);
        if (yn != NULL)
                DSVariablePoolFree(yn);
        if (yc != NULL)
                DSVariablePoolFree(yc);
        return systemEquations;
}

DSUInteger dsCyclicalCaseNumberOfAugmentedSystems(const DSDesignSpace * original,
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

DSStack * dsCyclicalCaseCreateAugmentedSystems(const DSCase * aCase,
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
        DSDictionary * validity;
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
                        validity = DSDesignSpaceCalculateAllValidCasesByResolvingCyclicalCases(subcase);
                        if (validity == NULL) {
                                DSDesignSpaceFree(subcase);
                        } else {
                                DSStackPush(augmentedSystemsStack, subcase);
                        }
                        DSDictionaryFree(validity);
                }
        }
        DSSecureFree(decayEquations);
        DSSecureFree(decayTerms);
        DSSecureFree(numberOfequations);
bail:
        return augmentedSystemsStack;
}


DSMatrix * dsCyclicalCaseExpandLcMatrix(const DSVariablePool * Xd,
                                               const DSMatrix * Lc,
                                               const DSVariablePool * yc)
{
        DSUInteger i, j, index;
        DSMatrix * newLc = NULL;
        if (Xd == NULL || yc == NULL) {
                DSError(M_DS_VAR_NULL, A_DS_ERROR);
                goto bail;
        }
        if (Lc == NULL) {
                DSError(M_DS_MAT_NULL, A_DS_ERROR);
                goto bail;
        }
        newLc = DSMatrixCalloc(DSMatrixRows(Lc), DSVariablePoolNumberOfVariables(Xd));
        for (i = 0; i < DSMatrixRows(Lc); i++) {
                for (j = 0; j < DSMatrixColumns(Lc); j++) {
                        index = DSVariablePoolIndexOfVariableWithName(Xd, DSVariableName(DSVariablePoolVariableAtIndex(yc, j)));
                        DSMatrixSetDoubleValue(newLc, i, index, DSMatrixDoubleValue(Lc, i, j));
                }
        }
bail:
        return newLc;
}


DSUInteger dsCyclicalCasePrimaryCycleVariableIndices(const DSCase * aCase,
                                                            DSMatrix * problematicEquations,
                                                            DSUInteger ** primaryVariables)
{
        DSUInteger numberOfCycles = 0, * cycleIndices;
        DSUInteger i, j, k, index;
        DSMatrix * Ad, *temp, *nullspace;
        double value, matrixValue;
        DSUInteger max;
        if (primaryVariables == NULL) {
                DSError(M_DS_NULL ": Pointer to hold primary cycle variable indices cannot be null", A_DS_ERROR);
                goto bail;
        }
        *primaryVariables = NULL;
        if (problematicEquations == NULL) {
                goto bail;
        }
        numberOfCycles = DSMatrixColumns(problematicEquations);
        *primaryVariables = DSSecureCalloc(sizeof(DSUInteger), numberOfCycles);
        cycleIndices = DSSecureCalloc(sizeof(DSUInteger), DSMatrixRows(problematicEquations));
        Ad = DSSSystemAd(DSCaseSSys(aCase));
        for (i = 0; i < numberOfCycles; i++) {
                k = 0;
                for (j = 0; j < DSMatrixRows(problematicEquations); j++) {
                        if (DSMatrixDoubleValue(problematicEquations, j, i) == 0)
                                continue;
                        cycleIndices[k++] = j;
                }
                if (k == 0) {
                        DSSecureFree(*primaryVariables);
                        *primaryVariables = NULL;
                        break;
                }
                if (k == 1) {
                        (*primaryVariables)[i] = cycleIndices[0];
                        continue;
                }
                temp = DSMatrixSubMatrixIncludingRowsAndColumns(Ad, k, k, cycleIndices, cycleIndices);
                nullspace = DSMatrixRightNullspace(temp);
                printf("showing nullspace: \n ");
                printMatrix(nullspace);
                if (nullspace == NULL) {
                        max = 0;
                } else {
                        value = NAN;
                        max = 0;
                        index = 0;
                        for (j = 0; j < DSMatrixRows(nullspace); j++) {
                                matrixValue = fabs(DSMatrixDoubleValue(nullspace, j, 0));
                                if (matrixValue < 1e-14)
                                        continue;
                                if (isnan(value)) {
                                        max = j;
                                        value = matrixValue;
                                        continue;

                                }
                                if (matrixValue < value) {
                                        value = matrixValue;
                                        max = j;
                                        printf("entre!");
                                }
                                index++;
                                if (index >= k) {
                                        break;
                                }
                        }
                }
                (*primaryVariables)[i] = cycleIndices[max];
                if (nullspace != NULL)
                        DSMatrixFree(nullspace);
                DSMatrixFree(temp);

        }
        DSSecureFree(cycleIndices);
        DSMatrixFree(Ad);

bail:
        return numberOfCycles;
}

DSUInteger dsCyclicalCaseAllSecondaryCycleVariables(const DSMatrix * problematicEquations,
                                                           const DSMatrixArray * coefficientArray,
                                                           const DSUInteger numberOfCycles,
                                                           const DSUInteger * primaryVariables,
                                                           DSUInteger ** cycleIndices,
                                                           double ** coefficients)
{
        DSUInteger numberSecondaryVariables = 0;
        DSUInteger i, j, k, index, * coefficientIndices = NULL;
        if (cycleIndices == NULL) {
                DSError(M_DS_NULL ": Pointer to hold primary cycle variable indices cannot be null", A_DS_ERROR);
                goto bail;
        }
        *cycleIndices = NULL;
        *coefficients = NULL;
        if (problematicEquations == NULL) {
                goto bail;
        }
        coefficientIndices = DSSecureCalloc(sizeof(DSUInteger), DSMatrixArrayNumberOfMatrices(coefficientArray));
        for (i = 0; i < DSMatrixColumns(problematicEquations); i++) {
                k = 0;
                for (j = 0; j < DSMatrixRows(problematicEquations); j++) {
                        if (DSMatrixDoubleValue(problematicEquations, j, i) == 0)
                                continue;
                        if (primaryVariables[i] == j) {
                                coefficientIndices[i] = k;
                                continue;
                        }
                        numberSecondaryVariables++;
                        k++;
                }
        }
        if (numberSecondaryVariables == 0) {
                goto bail;
        }
        *cycleIndices = DSSecureCalloc(sizeof(DSUInteger), numberSecondaryVariables);
        *coefficients = DSSecureCalloc(sizeof(double), numberSecondaryVariables);
        index = 0;
        for (i = 0; i < DSMatrixColumns(problematicEquations); i++) {
                k = 0;
                for (j = 0; j < DSMatrixRows(problematicEquations); j++) {
                        if (DSMatrixDoubleValue(problematicEquations, j, i) == 0)
                                continue;
                        if (primaryVariables[i] == j)
                                continue;
                        (*coefficients)[index] = DSMatrixArrayDoubleWithIndices(coefficientArray, i, k++, 0)/DSMatrixArrayDoubleWithIndices(coefficientArray, i, coefficientIndices[i], 0);
                        (*cycleIndices)[index++] = j;
                }
        }
bail:
        if (coefficientIndices != NULL)
                DSSecureFree(coefficientIndices);
        return numberSecondaryVariables;
}

DSUInteger dsCyclicalCaseSecondaryCycleVariableIndicesForCycle(DSMatrix * problematicEquations,
                                                                      DSUInteger cycleNumber,
                                                                      DSUInteger * primaryVariables,
                                                                      DSUInteger ** cycleIndices)
{
        DSUInteger numberOfCycleVariables = 0;
        DSUInteger i, j;
        if (cycleIndices == NULL) {
                DSError(M_DS_NULL ": Pointer to hold primary cycle variable indices cannot be null", A_DS_ERROR);
                goto bail;
        }
        *cycleIndices = NULL;
        if (problematicEquations == NULL) {
                goto bail;
        }
        for (i = 0; i < DSMatrixRows(problematicEquations); i++) {
                if (DSMatrixDoubleValue(problematicEquations, i, cycleNumber) == 0)
                        continue;
                numberOfCycleVariables++;
        }
        if (numberOfCycleVariables == 0) {
                goto bail;
        }
        numberOfCycleVariables -= 1;
        if (numberOfCycleVariables == 0) {
                goto bail;
        }
        *cycleIndices = DSSecureCalloc(sizeof(DSUInteger), numberOfCycleVariables);
        j = 0;
        for (i = 0; i < DSMatrixRows(problematicEquations); i++) {
                if (DSMatrixDoubleValue(problematicEquations, i, cycleNumber) == 0)
                        continue;
                if (i != primaryVariables[cycleNumber]) {
                        (*cycleIndices)[j] = i;
                        j++;
                }
                if (j == numberOfCycleVariables)
                        break;
        }
bail:
        return numberOfCycleVariables;
}

void dsCyclicalCasePartitionSolutionMatrices(const DSCase * aCase,
                                                    const DSUInteger numberOfSecondaryVariables,
                                                    const DSUInteger * secondaryVariables,
                                                    DSMatrix ** ADn,
                                                    DSMatrix ** ADc,
                                                    DSMatrix ** AIn,
                                                    DSMatrix ** Bn,
                                                    DSVariablePool ** yn,
                                                    DSVariablePool ** yc)
{
        DSMatrix * tempMatrix, *Ad, *Ai, *B;
        const DSSSystem * ssystem;
        DSUInteger i;
        char * name;
        if (ADn == NULL || ADc == NULL || AIn == NULL || Bn == NULL) {
                DSError(M_DS_NULL ": Matrix pointers to hold partitioned matrices cannot be null", A_DS_ERROR);
                goto bail;
        }
        *ADn = NULL;
        *ADc = NULL;
        *AIn = NULL;
        *Bn = NULL;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (secondaryVariables == NULL) {
                goto bail;
        }
        ssystem = DSCaseSSys(aCase);
        Ad = DSSSystemAd(ssystem);
        Ai = DSSSystemAi(ssystem);
        B = DSSSystemB(ssystem);
        tempMatrix = DSMatrixSubMatrixIncludingRows(Ad, numberOfSecondaryVariables, secondaryVariables);
        *ADc = DSMatrixSubMatrixExcludingColumns(tempMatrix, numberOfSecondaryVariables, secondaryVariables);
        *ADn = DSMatrixSubMatrixIncludingColumns(tempMatrix, numberOfSecondaryVariables, secondaryVariables);
        DSMatrixFree(tempMatrix);
        *AIn = DSMatrixSubMatrixIncludingRows(Ai, numberOfSecondaryVariables, secondaryVariables);
        *Bn = DSMatrixSubMatrixIncludingRows(B, numberOfSecondaryVariables, secondaryVariables);
        *yn = DSVariablePoolAlloc();
        *yc = DSVariablePoolAlloc();
        DSMatrixFree(Ad);
        DSMatrixFree(Ai);
        DSMatrixFree(B);
        for (i = 0; i < numberOfSecondaryVariables; i++) {
                name = DSVariableName(DSVariablePoolVariableAtIndex(DSSSystemXd(ssystem), secondaryVariables[i]));
                DSVariablePoolAddVariableWithName(*yn, name);
        }
        for (i = 0; i < DSVariablePoolNumberOfVariables(DSSSystemXd(ssystem)); i++) {
                name = DSVariableName(DSVariablePoolVariableAtIndex(DSSSystemXd(ssystem), i));
                if (DSVariablePoolHasVariableWithName(*yn, name) == false) {
                        DSVariablePoolAddVariableWithName(*yc, name);
                }
        }
bail:
        return;
}

void dsCyclicalCaseSolutionOfPartitionedMatrices(const DSCase * aCase,
                                                        const DSUInteger numberOfSecondaryVariables,
                                                        const DSUInteger * secondaryVariables,
                                                        DSMatrix ** LI,
                                                        DSMatrix **Lc,
                                                        DSMatrix **MBn,
                                                        DSVariablePool ** yn,
                                                        DSVariablePool ** yc)
{
        DSMatrix *ADn = NULL, * AIn = NULL, * ADc = NULL, * Bn = NULL, * Mn = NULL;
        if (LI == NULL || Lc == NULL || MBn == NULL) {
                DSError(M_DS_NULL ": Matrix pointers to hold partitioned matrices cannot be null", A_DS_ERROR);
                goto bail;
        }
        *LI = NULL;
        *Lc = NULL;
        *MBn = NULL;
        if (aCase == NULL) {
                DSError(M_DS_CASE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (secondaryVariables == NULL) {
                goto bail;
        }
        dsCyclicalCasePartitionSolutionMatrices(aCase, numberOfSecondaryVariables, secondaryVariables, &ADn, &ADc, &AIn, &Bn, yn, yc);
        if (ADn == NULL || ADc == NULL || AIn == NULL || Bn == NULL) {
                goto bail;
        }
        Mn = DSMatrixInverse(ADn);
        if (Mn == NULL) {
                goto bail;
        }
        *LI = DSMatrixByMultiplyingMatrix(Mn, AIn);
        *Lc = DSMatrixByMultiplyingMatrix(Mn, ADc);
        *MBn = DSMatrixByMultiplyingMatrix(Mn, Bn);
        DSMatrixMultiplyByScalar(*LI, 1.);
        DSMatrixMultiplyByScalar(*Lc, 1.);
bail:
        if (ADn != NULL)
                DSMatrixFree(ADn);
        if (ADc != NULL)
                DSMatrixFree(ADc);
        if (AIn != NULL)
                DSMatrixFree(AIn);
        if (Bn != NULL)
                DSMatrixFree(Bn);
        if (Mn != NULL)
                DSMatrixFree(Mn);
        return;
}

char * dsCyclicalCaseEquationForFlux(const DSCase * aCase,
                                            const DSDesignSpace * original,
                                            const DSMatrixArray * coefficientArray,
                                            const DSUInteger variableIndex,
                                            const DSUInteger fluxIndex,
                                            const bool positiveFlux,
                                            const DSUInteger cycleNumber,
                                            const DSUInteger primaryVariable,
                                            const DSUInteger numberSecondaryVariables,
                                            const DSUInteger * secondaryCycleVariables,
                                            const DSMatrix * LI,
                                            const DSMatrix * Lc,
                                            const DSMatrix * MBn,
                                            const DSVariablePool * yc)
{
        char * fluxEquationString = NULL;
        DSExpression * fluxEquation = NULL;
        const DSMatrix * A;
        DSMatrix *Kd, *LKi, *LKd, *Kc, *Ks, *Ki, *temp;
        const DSGMASystem * gma = NULL;
        DSUInteger i, j, index;
        double denominator = 0, numerator = 0;
        DSExpression * (*termFunction)(const DSGMASystem *, const DSUInteger, DSUInteger);
        // Checks...
        gma = DSDesignSpaceGMASystem(original);
        denominator = DSMatrixArrayDoubleWithIndices(coefficientArray, cycleNumber, 0, 0);
        if (variableIndex == primaryVariable) {
                numerator = denominator;
        } else {
                for (i = 0; i < numberSecondaryVariables; i++) {
                        if (secondaryCycleVariables[i] == variableIndex) {
                                numerator = DSMatrixArrayDoubleWithIndices(coefficientArray, cycleNumber, i+1, 0);
                        }
                }
        }
        if (numerator == 0) {
                goto bail;
        }
        if (positiveFlux == true) {
                A = DSGMASystemAlpha(gma);
                Kd = DSMatrixArrayMatrix(DSGMASystemGd(gma), variableIndex);
                Kd = DSMatrixSubMatrixIncludingRowList(Kd, 1, fluxIndex);
                Ki = DSMatrixArrayMatrix(DSGMASystemGi(gma), variableIndex);
                Ki = DSMatrixSubMatrixIncludingRowList(Ki, 1, fluxIndex);
                termFunction = DSGMASystemPositiveTermForEquations;
        } else {
                A = DSGMASystemBeta(gma);
                Kd = DSMatrixArrayMatrix(DSGMASystemGd(gma), variableIndex);
                Kd = DSMatrixSubMatrixIncludingRowList(Kd, 1, fluxIndex);
                Ki = DSMatrixArrayMatrix(DSGMASystemGi(gma), variableIndex);
                Ki = DSMatrixSubMatrixIncludingRowList(Ki, 1, fluxIndex);
                termFunction = DSGMASystemNegativeTermForEquations;
        }
        if (numberSecondaryVariables > 0) {
                Ks = DSMatrixSubMatrixIncludingColumns(Kd, numberSecondaryVariables, secondaryCycleVariables);
                LKi = DSMatrixByMultiplyingMatrix(Ks, LI);
                Kc = DSMatrixCalloc(DSMatrixRows(Kd), DSMatrixColumns(Kd));
                temp = dsCyclicalCaseExpandLcMatrix(DSCaseXd(aCase), Lc, yc);
                for (i = 0; i < DSMatrixRows(Kc); i++) {
                        for (j = 0; j < DSMatrixColumns(Lc); j++) {
                                index = DSVariablePoolIndexOfVariableWithName(DSCaseXd(aCase), DSVariableName(DSVariablePoolVariableAtIndex(yc, j)));
                                DSMatrixSetDoubleValue(Kc, i, index, DSMatrixDoubleValue(Kd, i, index));
                        }
                }
                LKd = DSMatrixByMultiplyingMatrix(Ks, temp);
                DSMatrixFree(temp);
                temp = DSMatrixByMultiplyingMatrix(Ks, MBn);
                DSMatrixAddByMatrix(LKd, Kc);
                DSMatrixAddByMatrix(Ki, LKi);
                DSMatrixSetDoubleValue(temp, 0, 0, DSMatrixDoubleValue(temp, 0, 0)+log10(DSMatrixDoubleValue(A, variableIndex, fluxIndex)));
                fluxEquation = DSExpressionFromPowerlawInMatrixForm(0, LKd, DSCaseXd(aCase), Ki, DSCaseXi(aCase), temp);
                fluxEquationString = DSSecureCalloc(sizeof(char), 1000);
                asprintf(&fluxEquationString, "%i*%s*%lf", (positiveFlux ? 1 : -1), DSExpressionAsString(fluxEquation), numerator/denominator);
                DSExpressionFree(fluxEquation);
        }// else {
                fluxEquation = termFunction(gma, variableIndex, fluxIndex);
                fluxEquationString = DSSecureCalloc(sizeof(char), 1000);
                asprintf(&fluxEquationString, "%s*%lf", DSExpressionAsString(fluxEquation), numerator/denominator);
                DSExpressionFree(fluxEquation);
//        }
bail:
        return fluxEquationString;
}

DSExpression ** dsCyclicalCaseEquationsForCycle(const DSCase * aCase,
                                                       const DSDesignSpace * original,
                                                       const DSMatrixArray * coefficientArray,
                                                       const DSUInteger cycleNumber,
                                                       const DSUInteger primaryCycleVariable,
                                                       const DSUInteger numberSecondaryVariables,
                                                       const DSUInteger * secondaryCycleVariables,
                                                       const DSMatrix * LI,
                                                       const DSMatrix * Lc,
                                                       const DSMatrix * MBn,
                                                       const DSVariablePool * yn,
                                                       const DSVariablePool * yc)
{
        DSUInteger i, j, index, numberOfX, pcount = 0, ncount = 0;
        const DSSSystem * ssys;
        char * string = NULL, *name, *flux;
        double value;
        DSExpression ** cycleEquations = NULL;

//        DSMatrix * LI;
//        DSMatrix * Lc;
//        DSMatrix * MBn;
//        DSVariablePool * yn;
//        DSVariablePool * yc;
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
        if (secondaryCycleVariables == NULL && numberSecondaryVariables > 0) {
                DSError(M_DS_NULL ": Array of secondary cycle variables is null", A_DS_ERROR);
                goto bail;
        }
        if (numberSecondaryVariables > 0) {
//                dsCyclicalCaseSolutionOfPartitionedMatrices(aCase, numberSecondaryVariables, secondaryCycleVariables, &LI, &Lc, &MBn, &yn, &yc);
                if (LI == NULL || Lc == NULL || MBn == NULL) {
                        goto bail;
                }
        }
        ssys = DSCaseSSystem(aCase);
        string = DSSecureCalloc(sizeof(char *), 1000);
        // asprintf introduces memory error....
        asprintf(&string, "%s. = ", DSVariableName(DSVariablePoolVariableAtIndex(DSSSystemXd(ssys), primaryCycleVariable)));
        for (i = 0; i < numberSecondaryVariables+1; i++) {
                if (i == 0) {
                        index = primaryCycleVariable;
                } else {
                        index = secondaryCycleVariables[i-1];
                }
                for (j = 1; j <= DSDesignSpaceSignature(original)[2*index]; j++) {
                        if (j == DSCaseSignature(aCase)[2*index]) {
                                continue;
                        }
                        flux = dsCyclicalCaseEquationForFlux(aCase, original,
                                                             coefficientArray,
                                                             index,
                                                             j-1,
                                                             true,
                                                             cycleNumber,
                                                             primaryCycleVariable,
                                                             numberSecondaryVariables,
                                                             secondaryCycleVariables,
                                                             LI, Lc, MBn, yc);
                        if (flux != NULL) {
                                asprintf(&string, "%s + %s", string, flux);
                                DSSecureFree(flux);
                                pcount++;
                        }
                }
                for (j = 1; j <= DSDesignSpaceSignature(original)[2*index+1]; j++) {
                        if (j == DSCaseSignature(aCase)[2*index+1]) {
                                continue;
                        }
                        flux = dsCyclicalCaseEquationForFlux(aCase, original,
                                                             coefficientArray,
                                                             index,
                                                             j-1,
                                                             false,
                                                             cycleNumber,
                                                             primaryCycleVariable,
                                                             numberSecondaryVariables,
                                                             secondaryCycleVariables,
                                                             LI, Lc, MBn, yc);
                        if (flux != NULL) {
                                asprintf(&string, "%s + %s", string, flux);
                                DSSecureFree(flux);
                                ncount++;
                        }
                }
        }
        printf("%s\n", string);
        if (pcount == 0 || ncount == 0) {
                goto bail;
        }
        cycleEquations = DSSecureCalloc(sizeof(DSExpression *), numberSecondaryVariables+1);
        cycleEquations[0] = DSExpressionByParsingString(string);
        for (i = 0; i < numberSecondaryVariables; i++) {
                name = DSVariableName(DSVariablePoolVariableAtIndex(DSSSystemXd(ssys), secondaryCycleVariables[i]));
                index = DSVariablePoolIndexOfVariableWithName(yn, name);
                asprintf(&string, "%s = 10^%lf", name, DSMatrixDoubleValue(MBn, index, 0));
                numberOfX = DSVariablePoolNumberOfVariables(DSSSystemXi(ssys));
                for (j = 0; j < numberOfX; j++) {
                        name = DSVariableName(DSVariablePoolVariableAtIndex(DSSSystemXi(ssys), j));
                        value = -DSMatrixDoubleValue(LI, index, j);
                        if (value == 0.0)
                                continue;
                        if (value == 1.0)
                                asprintf(&string, "%s*%s", string, name);
                        else
                                asprintf(&string, "%s*%s^%lf", string, name, value);
                }
                numberOfX = DSVariablePoolNumberOfVariables(yc);
                for (j = 0; j < numberOfX; j++) {
                        name = DSVariableName(DSVariablePoolVariableAtIndex(yc, j));
                        value = -DSMatrixDoubleValue(Lc, index, j);
                        if (value == 0.0)
                                continue;
                        if (value == 1.0)
                                asprintf(&string, "%s*%s", string, name);
                        else
                                asprintf(&string, "%s*%s^%lf", string, name, value);
                }
                cycleEquations[i+1] = DSExpressionByParsingString(string);
        }

bail:
        if (string != NULL)
                DSSecureFree(string);
        return cycleEquations;
}

void dsCyclicalCaseAugmentedEquationsForCycleALT(char ** systemEquations,
                                                     const DSCase * aCase,
                                                     const DSDesignSpace * original,
                                                     const DSMatrix * problematicMatrix,
                                                     const DSMatrixArray * coefficientArray,
                                                     const DSUInteger cycleNumber,
                                                     const DSUInteger primaryCycleVariable,
                                                     const DSUInteger numberSecondaryVariables,
                                                     const DSUInteger * secondaryVariables,
                                                     const DSMatrix * LI,
                                                     const DSMatrix * Lc,
                                                     const DSMatrix * Mb,
                                                     const DSVariablePool * yn,
                                                     const DSVariablePool * yc)
{
        DSExpression *fluxEquation;
        char * string, *name;
        DSUInteger i, j, k, l, index, eta;
        const DSUInteger *signature;
        DSMatrix * C, *Kd, *Ki;
        DSMatrix * Ks, *Kn, *LKi, *LKd, *temp;
        const DSGMASystem * gma;
        double denominator = 0, numerator = 0;

        if (original == NULL) {
                DSError(M_DS_DESIGN_SPACE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (problematicMatrix == NULL) {
                DSError(M_DS_MAT_NULL, A_DS_ERROR);
                goto bail;
        }
        l = 0;
        for (i = 0; i < DSDesignSpaceNumberOfEquations(original); i++) {
                if (DSMatrixDoubleValue(problematicMatrix, i, cycleNumber) == 0.0f) {
                        continue;
                }
                if (primaryCycleVariable == i) {
                        denominator = DSMatrixArrayDoubleWithIndices(coefficientArray, cycleNumber, l, 0);
                        break;
                }
                l++;
        }
        signature = DSDesignSpaceSignature(original);
        gma = DSDesignSpaceGMASystem(original);
        string = DSSecureCalloc(sizeof(char), 1000);
        DSSecureFree(systemEquations[primaryCycleVariable]);
        systemEquations[primaryCycleVariable] = string;
        name = DSVariableName(DSVariablePoolVariableAtIndex(DSGMASystemXd(gma), primaryCycleVariable));
        sprintf(string, "%s. = ", name);
        if (string != systemEquations[primaryCycleVariable]) {
                DSSecureFree(systemEquations[primaryCycleVariable]);
                systemEquations[primaryCycleVariable] = string;
        }
        l = 0;
        for (i = 0; i < 2*DSDesignSpaceNumberOfEquations(original); i++) {
                if (DSMatrixDoubleValue(problematicMatrix, i/2, cycleNumber) == 0.0f) {
                        continue;
                }
                index = i/2;
                if (i % 2 == 0) {
                        if (l >= DSMatrixRows(DSMatrixArrayMatrix(coefficientArray, cycleNumber))) // This Check was added to avoid an error that was appearing.
                                break;
                        C = DSMatrixSubMatrixIncludingRowList(DSGMASystemAlpha(gma), 1, i/2);
                        Kd = DSMatrixCopy(DSMatrixArrayMatrix(DSGMASystemGd(gma), i/2));
                        Ki = DSMatrixCopy(DSMatrixArrayMatrix(DSGMASystemGi(gma), i/2));
                        if (index != primaryCycleVariable) {
                                numerator = DSMatrixArrayDoubleWithIndices(coefficientArray, cycleNumber, l++, 0);
                                dsCyclicalCaseEquilibriumEquationForVariableALT(index, systemEquations, aCase, original, LI, Lc, Mb, yn, yc);
                        } else {
                                numerator = denominator;
                                l++;
                        }
                } else {
                        C = DSMatrixSubMatrixIncludingRowList(DSGMASystemBeta(gma), 1, i/2);
                        Kd = DSMatrixCopy(DSMatrixArrayMatrix(DSGMASystemHd(gma), i/2));
                        Ki = DSMatrixCopy(DSMatrixArrayMatrix(DSGMASystemHi(gma), i/2));

                }
                temp = DSMatrixTranspose(C);
                DSMatrixFree(C);
                C = temp;
                if (numberSecondaryVariables > 0) { // block commented out by Miguel
                        Ks = DSMatrixSubMatrixIncludingColumns(Kd, numberSecondaryVariables, secondaryVariables);
                        Kn = DSMatrixCopy(Kd);
                        for (j = 0; j < DSMatrixRows(Kd); j++) {
                                for (k = 0; k < numberSecondaryVariables; k++) {
                                        index = secondaryVariables[k];
                                        DSMatrixSetDoubleValue(Kn, j, index, 0.0);
                                }
                        }
                        LKi = DSMatrixByMultiplyingMatrix(Ks, LI);
                        LKd = DSMatrixByMultiplyingMatrix(Ks, Lc);
                        DSMatrixMultiplyByScalar(LKd, -1.);
                        DSMatrixMultiplyByScalar(LKi, -1.);
                        DSMatrixAddByMatrix(LKd, Kn);
//                        DSMatrixFree(Kd);
//                        Kd = LKd;
//                        DSMatrixAddByMatrix(Ki, LKi);
                        temp = DSMatrixByMultiplyingMatrix(Ks, Mb);
                        for (j = 0; j < DSMatrixRows(temp); j++) {
                                DSMatrixSetDoubleValue(temp, j, 0, numerator/denominator*DSMatrixDoubleValue(C, j, 0)*pow(10, DSMatrixDoubleValue(temp, j, 0)));
                        }
                        DSMatrixFree(C);
                        DSMatrixFree(LKi);
                        C = temp;
                        DSMatrixFree(Kn);
                        DSMatrixFree(Ks);
                }
                for (j = 0; j < signature[i]; j++) {
                        if (j + 1 == DSCaseSignature(aCase)[i])
                                continue;
                        fluxEquation = DSExpressionFromPowerlawInMatrixForm(j, Kd, DSGMASystemXd(gma), Ki, DSGMASystemXi(gma), C);
                        name = DSExpressionAsString(fluxEquation);
                        if (i % 2 == 0) {
                                sprintf(string, "%s + %s", systemEquations[primaryCycleVariable], name);
                        } else {
                                if ( primaryCycleVariable < i && i < (primaryCycleVariable+1)*2 ){
                                sprintf(string, "%s - %s", systemEquations[primaryCycleVariable], name);
                            }
                        }
                        if (string != systemEquations[primaryCycleVariable]) {
                                DSSecureFree(systemEquations[primaryCycleVariable]);
                                systemEquations[primaryCycleVariable] = string;
                        }
                        DSExpressionFree(fluxEquation);
//                        if (DSCaseNumber(aCase) == 126 && DSDesignSpaceCyclical(original) == false) {
//                                printf("%i,%i: %s\n", i, j, name);
//                        }
                        DSSecureFree(name);
                }
                DSMatrixFree(Kd);
                DSMatrixFree(Ki);
                DSMatrixFree(C);
        }

bail:
        return;
}

DSDesignSpace * dsCyclicalCaseAugmentedSystemForSubdominantDecays(const DSCase * aCase,
                                                                         const DSDesignSpace * original,
                                                                         DSMatrix * problematicEquations,
                                                                         const DSMatrixArray * problematicTerms,
                                                                         const DSMatrixArray * coefficientArray,
                                                                         const DSUInteger *subdominantDecaySpecies,
                                                                         const DSUInteger *subdominantDecayTerm)
{
        DSDesignSpace * augmentedSystem = NULL, *modifierDesignSpace = NULL;
        DSGMASystem * gma = NULL;
        DSExpression ** augmentedEquations = NULL, ** dsEquations, ** caseEquations;
        DSUInteger i, j, k, l,*signature;
        double subCoefficient, value;
        char ** subcaseEquations;
        DSUInteger positiveTerms;
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
                positiveTerms = 0;
                l = 0;
                for (j = 0; j < DSMatrixRows(problematicEquations); j++) {
                        if (DSMatrixDoubleValue(problematicEquations, j, i) == 0)
                                continue;
                        positiveTerms += DSDesignSpaceSignature(original)[2*j]-1;
                        if (j == subdominantDecaySpecies[i]) {
                                subCoefficient = DSMatrixArrayDoubleWithIndices(coefficientArray, i, l, 0);
                        }
                        l++;
                }
                l = 0;
                for (j = 0; j < DSMatrixRows(problematicEquations); j++) {
                        if (DSMatrixDoubleValue(problematicEquations, j, i) == 0)
                                continue;
                        value = DSMatrixArrayDoubleWithIndices(coefficientArray, i, l, 0);
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
                if (positiveTerms == 0)
                        goto bail;
                j = subdominantDecaySpecies[i];
                for (k = 0; k < DSMatrixColumns(DSGMASystemBeta(gma)); k++) {
                        if (k+1 == aCase->signature[2*j+1]) {
                                DSMatrixSetDoubleValue((DSMatrix *)DSGMASystemBeta(gma), j, k, 0.0f);
                        } else {
                                DSMatrixSetDoubleValue((DSMatrix *)DSGMASystemBeta(gma), j, k,
                                                       DSMatrixDoubleValue(DSGMASystemBeta(gma), j, k));
                        }
                }
                augmentedEquations[i] = DSExpressionAddExpressions(augmentedEquations[i], DSGMASystemNegativeTermsForEquations(gma, j));
        }
        augmentedSystem = dsCyclicalCaseCreateUniqueAugmentedSystem(aCase,
                                                                    gma,
                                                                    problematicEquations,
                                                                    (const DSExpression **)augmentedEquations,
                                                                    subdominantDecaySpecies);
        DSDesignSpaceAddConditions(augmentedSystem, DSCaseCd(aCase), DSCaseCi(aCase), DSCaseDelta(aCase));
        dsAddConstraintsForSubdominantDecays(augmentedSystem,
                                             aCase,
                                             original,
                                             problematicEquations,
                                             coefficientArray,
                                             subdominantDecaySpecies,
                                             subdominantDecayTerm);
        DSDesignSpaceSetSerial(augmentedSystem, true);
        DSDesignSpaceSetCyclical(augmentedSystem, true);
        dsEquations = DSDesignSpaceEquations(original);
        for (i = 0; i < DSGMASystemNumberOfEquations(gma); i++) {
                subcaseEquations[i] = DSExpressionAsString(dsEquations[i]);
                signature[2*i] = DSCaseSignature(aCase)[2*i];
                signature[2*i+1] = DSCaseSignature(aCase)[2*i+1];
                DSExpressionFree(dsEquations[i]);
        }
        DSSecureFree(dsEquations);
        dsEquations = DSDesignSpaceEquations(augmentedSystem);
        caseEquations = DSCaseEquations(aCase);
        for (i = 0; i < DSMatrixColumns(problematicEquations); i++) {
                for (j = 0; j < DSMatrixRows(problematicEquations); j++) {
                        if (DSMatrixDoubleValue(problematicEquations, j, i) == 0)
                                continue;
                        signature[2*j] = 0;
                        signature[2*j+1] = 0;
                        DSSecureFree(subcaseEquations[j]);
                        subcaseEquations[j] = DSExpressionAsString(caseEquations[j]);
                }
                j = subdominantDecaySpecies[i];
                DSSecureFree(subcaseEquations[j]);
                subcaseEquations[j] = DSExpressionAsString(dsEquations[j]);
        }
        for (i = 0; i < DSDesignSpaceNumberOfEquations(augmentedSystem); i++) {
                DSExpressionFree(dsEquations[i]);
                DSExpressionFree(caseEquations[i]);
        }
        DSSecureFree(dsEquations);
        DSSecureFree(caseEquations);
        modifierDesignSpace = DSDesignSpaceByParsingStringsWithXi(subcaseEquations,
                                                                  DSGMASystemXd_a(DSDesignSpaceGMASystem(original)),
                                                                  DSGMASystemXi(DSDesignSpaceGMASystem(original)),
                                                                  DSDesignSpaceNumberOfEquations(original));
        modifierDesignSpace->Cd = DSMatrixCopy(augmentedSystem->Cd);
        modifierDesignSpace->Ci = DSMatrixCopy(augmentedSystem->Ci);
        modifierDesignSpace->delta = DSMatrixCopy(augmentedSystem->delta);
        DSCyclicalCaseDesignSpaceCalculateSubCyclicalCases(augmentedSystem, modifierDesignSpace, signature);
        for (i = 0; i < DSDesignSpaceNumberOfEquations(original); i++) {
                DSSecureFree(subcaseEquations[i]);
        }
        DSSecureFree(subcaseEquations);
        DSSecureFree(signature);
bail:
        if (augmentedEquations != NULL) {
                for (i = 0; i < DSMatrixColumns(problematicEquations); i++) {
                        if (augmentedEquations[i] != NULL)
                                DSExpressionFree(augmentedEquations[i]);
                }
                DSSecureFree(augmentedEquations);
        }
        if (gma != NULL) {
                DSGMASystemFree(gma);
        }
        if (modifierDesignSpace != NULL) {
                DSDesignSpaceFree(modifierDesignSpace);
        }
        return augmentedSystem;
}

char ** dsCyclicalCaseOriginalEquationsWithEquilibriumConstraintsALT(const DSCase * aCase,
                                                                         const DSDesignSpace * original,
                                                                         const DSUInteger numberSecondaryVariables,
                                                                         const DSUInteger * secondaryVariables,
                                                                         const double * coefficientMultipliers,
                                                                         const DSMatrix * LI,
                                                                         const DSMatrix * Lc,
                                                                         const DSMatrix * Mb,
                                                                         const DSVariablePool * yn,
                                                                         const DSVariablePool * yc)
{
        DSExpression ** equations, **caseEquations, *fluxEquation;
        char ** systemEquations = NULL, * string, *name;
        DSUInteger i, j, k, index;
        const DSUInteger *signature;
        DSMatrix * C, *Kd, *Ki;
        DSMatrix * Ks, *Kn, *LKi, *LKd, *temp;
        const DSGMASystem * gma;
        if (original == NULL) {
                DSError(M_DS_DESIGN_SPACE_NULL, A_DS_ERROR);
                goto bail;
        }
        gma = DSDesignSpaceGMASystem(original);
        if (numberSecondaryVariables == 0) {
                systemEquations = DSSecureCalloc(sizeof(char *), DSDesignSpaceNumberOfEquations(original));
                equations = DSDesignSpaceEquations(original);
                caseEquations = DSCaseEquations(aCase);
                for (i = 0; i < DSDesignSpaceNumberOfEquations(original); i++) {
                        name = DSVariableName(DSVariablePoolVariableAtIndex(DSGMASystemXd(gma), i));
                        if (DSVariablePoolHasVariableWithName(DSGMASystemXd_t(gma), name) == true) {
                                systemEquations[i] = DSExpressionAsString(equations[i]);
                        } else {
                                systemEquations[i] = DSExpressionAsString(caseEquations[i]);
                        }
                        DSExpressionFree(equations[i]);
                        DSExpressionFree(caseEquations[i]);
                }
                DSSecureFree(equations);
                DSSecureFree(caseEquations);
                goto bail;
        }
        if (LI == NULL || Lc == NULL || Mb == NULL) {
                goto bail;
        }
        if (yn == NULL || yc == NULL) {
                goto bail;
        }
        systemEquations = DSSecureCalloc(sizeof(char *), DSDesignSpaceNumberOfEquations(original));
        signature = DSDesignSpaceSignature(original);
        equations = DSDesignSpaceEquations(original);
        caseEquations = DSCaseEquations(aCase);
        for (i = 0; i < 2*DSDesignSpaceNumberOfEquations(original); i++) {
                name = DSVariableName(DSVariablePoolVariableAtIndex(DSGMASystemXd(gma), i/2));
                if (DSVariablePoolHasVariableWithName(DSGMASystemXd_t(gma), name) == false) {
                        if (i % 2 == 0) {
                                //                                string = DSSecureCalloc(sizeof(char), 1000);
                                systemEquations[i/2] = DSExpressionAsString(caseEquations[i/2]);//DSExpressionAsString(caseEquations[i/2]);
                        }
                        continue;
                }
                if (i % 2 == 0) {
                        string = DSSecureCalloc(sizeof(char), 1000);
                        systemEquations[i/2] = string;
                        sprintf(string, "%s. = ", name);
                        if (string != systemEquations[i/2]) {
                                DSSecureFree(systemEquations[i/2]);
                                systemEquations[i/2] = string;
                        }
                        C = DSMatrixSubMatrixIncludingRowList(DSGMASystemAlpha(gma), 1, i/2);
                        Kd = DSMatrixCopy(DSMatrixArrayMatrix(DSGMASystemGd(gma), i/2));
                        Ki = DSMatrixCopy(DSMatrixArrayMatrix(DSGMASystemGi(gma), i/2));
                } else {
                        C = DSMatrixSubMatrixIncludingRowList(DSGMASystemBeta(gma), 1, i/2);
                        Kd = DSMatrixCopy(DSMatrixArrayMatrix(DSGMASystemHd(gma), i/2));
                        Ki = DSMatrixCopy(DSMatrixArrayMatrix(DSGMASystemHi(gma), i/2));
                }
                temp = DSMatrixTranspose(C);
                DSMatrixFree(C);
                C = temp;
//                if (numberSecondaryVariables > 0) {
//                        Ks = DSMatrixSubMatrixIncludingColumns(Kd, numberSecondaryVariables, secondaryVariables);
//                        Kn = DSMatrixCopy(Kd);
//                        for (j = 0; j < DSMatrixRows(Kd); j++) {
//                                for (k = 0; k < numberSecondaryVariables; k++) {
//                                        index = secondaryVariables[k];
//                                        DSMatrixSetDoubleValue(Kn, j, index, 0.0);
//                                }
//                        }
//                        LKi = DSMatrixByMultiplyingMatrix(Ks, LI);
//                        LKd = DSMatrixByMultiplyingMatrix(Ks, Lc);
//                        DSMatrixMultiplyByScalar(LKd, -1.);
//                        DSMatrixMultiplyByScalar(LKi, -1.);
//                        DSMatrixAddByMatrix(LKd, Kn);
//                        DSMatrixAddByMatrix(Ki, LKi);
//                        temp = DSMatrixByMultiplyingMatrix(Ks, Mb);
//                        for (j = 0; j < DSMatrixRows(temp); j++) {
//                                DSMatrixSetDoubleValue(C, j, 0, DSMatrixDoubleValue(C, j, 0)*pow(10, DSMatrixDoubleValue(temp, j, 0)));
//                        }
//                        DSMatrixFree(Kd);
//                        Kd = LKd;
//                        DSMatrixFree(temp);
//                        DSMatrixFree(LKi);
//                        DSMatrixFree(Kn);
//                        DSMatrixFree(Ks);
//                }
                for (j = 0; j < numberSecondaryVariables; j++) {
                        if (i/2 == secondaryVariables[j]) {
                                DSMatrixMultiplyByScalar(C, coefficientMultipliers[j]);
                                break;
                        }
                }
                for (j = 0; j < signature[i]; j++) {
                        fluxEquation = DSExpressionFromPowerlawInMatrixForm(j, Kd, DSGMASystemXd(gma), Ki, DSGMASystemXi(gma), C);


                        name = DSExpressionAsString(fluxEquation);
                        //                      let's print matrices Kd, variablepool Xd
//                        printf("\nThe value of j is: %u \n", j);
//                        printf("\nThe matrix Kd is: \n");
//                        printMatrix(Kd);
//                        printf("\nThe matrix Ki is: \n");
//                        printMatrix(Ki);
//                        printf("\nThe matrix C is: \n");
//                        printMatrix(C);
//
//                        printf("\nThe corresponding fluxEquation in string form is\n");
//                        printf(name);
//                        printf("\n");



                        if (i % 2 == 0) {
                                sprintf(string, "%s + %s", systemEquations[i/2], name);
                        } else {
                                sprintf(string, "%s - %s", systemEquations[i/2], name);
                        }
                        if (string != systemEquations[i/2]) {
                                DSSecureFree(systemEquations[i/2]);
                                systemEquations[i/2] = string;
                        }
                        DSExpressionFree(fluxEquation);
                        DSSecureFree(name);
                }
                DSMatrixFree(Kd);
                DSMatrixFree(Ki);
                DSMatrixFree(C);
        }
        for (i = 0; i < DSDesignSpaceNumberOfEquations(original); i++) {
//                if (DSCaseNumber(aCase) == 126 && DSDesignSpaceCyclical(original) == false) {
//                        printf("%s\n", systemEquations[i]);
//                }
                DSExpressionFree(equations[i]);
                DSExpressionFree(caseEquations[i]);
        }
        DSSecureFree(equations);
        DSSecureFree(caseEquations);
bail:
        return systemEquations;
}


void dsAddConstraintsForSubdominantDecays(DSDesignSpace * subcase, const DSCase * aCase, const DSDesignSpace * original, const DSMatrix * problematicEquations, const DSMatrixArray * coefficientArray, const DSUInteger * subdominantDecays, const DSUInteger * subdominantDecayTerms)
{
        const DSGMASystem * gma = NULL;
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
        gma = DSDesignSpaceGMASystem(original);
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
                                value = log10((subCoefficient/coefficient)*DSMatrixDoubleValue(DSGMASystemBeta(gma), subdominantDecays[i], subdominantDecayTerms[i])
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
        return;
}

DSDesignSpace * dsCyclicalCaseCreateUniqueAugmentedSystem(const DSCase *aCase, const DSGMASystem * modifiedGMA, const DSMatrix  * problematicEquations, const DSExpression ** augmentedEquations, const DSUInteger * subdominantDecays)
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
                equations[i] = NULL;
                asprintf(&equations[i], "%s = %s", temp_lhs, temp_rhs);
                DSExpressionFree(eqLHS);
                DSSecureFree(temp_rhs);
                DSSecureFree(temp_lhs);
        }
        ds = DSDesignSpaceByParsingStringsWithXi(equations, Xda, DSGMASystemXi(modifiedGMA), DSCaseNumberOfEquations(aCase));
        for (i = 0; i < DSCaseNumberOfEquations(aCase); i++) {
                DSSecureFree(equations[i]);
                DSExpressionFree(caseEquations[i]);
        }
        DSVariablePoolFree(Xda);
        DSSecureFree(equations);
        DSSecureFree(caseEquations);
bail:
        return ds;
}

void dsCyclicalCaseEquilibriumEquationForVariableALT(DSUInteger index,
                                                         char ** systemEquations,
                                                         const DSCase * aCase,
                                                         const DSDesignSpace * original,
                                                         const DSMatrix * LI,
                                                         const DSMatrix * Lc,
                                                         const DSMatrix * Mb,
                                                         const DSVariablePool * yn,
                                                         const DSVariablePool * yc)
{
        DSExpression *fluxEquation, ** caseEquations;
        char * string, *name;
        DSUInteger i;
        DSMatrix * C;
        if (original == NULL) {
                DSError(M_DS_DESIGN_SPACE_NULL, A_DS_ERROR);
                goto bail;
        }
        if (LI == NULL || Lc == NULL || Mb == NULL) {
                DSError(M_DS_MAT_NULL, A_DS_ERROR);
                goto bail;
        }
        if (yn == NULL || yc == NULL) {
                DSError(M_DS_VAR_NULL, A_DS_ERROR);
                goto bail;
        }
//        C = DSMatrixCopy(Mb);
//        for (i = 0; i < DSMatrixRows(C); i++) {
//                DSMatrixSetDoubleValue(C, i, 0, pow(10, DSMatrixDoubleValue(C, i, 0)));
//        }
//        name = DSVariableName(DSVariablePoolVariableAtIndex(DSCaseXd(aCase), index));
//        i = DSVariablePoolIndexOfVariableWithName(yn, name);
//        Lc = DSMatrixCopy(Lc);
//        LI = DSMatrixCopy(LI);
//        DSMatrixMultiplyByScalar((DSMatrix *)Lc, -1.0f);
//        DSMatrixMultiplyByScalar((DSMatrix *)LI, -1.0f);
//        fluxEquation = DSExpressionFromPowerlawInMatrixForm(i, Lc, DSCaseXd(aCase), LI, DSCaseXi(aCase), C);
//        string = DSExpressionAsString(fluxEquation);
//        DSExpressionFree(fluxEquation);
//        DSSecureFree(systemEquations[index]);
//        asprintf(&(systemEquations[index]), "%s = %s", name, string);
//        DSSecureFree(string);
//        DSMatrixFree((DSMatrix *)Lc);
//        DSMatrixFree((DSMatrix *)LI);
//        DSMatrixFree((DSMatrix *)C);


          caseEquations = DSCaseEquations(aCase);
          systemEquations[index] = DSExpressionAsString(caseEquations[index]);
          DSSecureFree(caseEquations);

bail:
        return;
}
