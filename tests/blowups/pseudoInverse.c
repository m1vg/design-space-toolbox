#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <designspace/DSStd.h>

int main() {

	DSInteger n_row = 2;
	DSInteger n_col = 2;
	double mat[2][2] =  {{ 0, 2},		// Example Rank-deficient matrix
		             	 { 0, -1}};
	DSUInteger i = 0;
	DSUInteger j = 0;

	DSMatrix * gA = DSMatrixAlloc (n_row, n_col);
	for (i = 0; i < n_row; i++)
		for (j = 0; j < n_col; j++)
	   		DSMatrixSetDoubleValue (gA, i, j, mat[i][j]);


    DSMatrix *U, *V, *S, *invS, *trU, *trV, *pA, *VinvS;
    DSMatrixArray *array;

    array = DSMatrixSVD(gA);
    // unpack matrices
    S = DSMatrixArrayMatrix(array, 0);
    U = DSMatrixArrayMatrix(array, 1);
    V = DSMatrixArrayMatrix(array, 2);

    printf("The matrix S is: \n");
    DSMatrixPrint(S);

    printf("The matrix U is: \n");
    DSMatrixPrint(U);

    printf("The matrix V is: \n");
    DSMatrixPrint(V);

    // calculate inverse of S, invS
    invS = DSMatrixCalloc(n_row, n_row);
    for (i = 0; i < n_row; i++){
        if (DSMatrixDoubleValue(S, 0, i) != 0){
            DSMatrixSetDoubleValue(invS, i, i, 1/DSMatrixDoubleValue(S, 0, i));
        }
    }

    // Transpose V
    trV = DSMatrixTranspose(V);

    // Transpose U
    trU = DSMatrixTranspose(U);

    // Calculate pseudo inverse of matrix A
    VinvS = DSMatrixByMultiplyingMatrix(trV, invS);
    pA = DSMatrixByMultiplyingMatrix(VinvS, trU);

    //print Result
    printf("The matrix gA is: \n");
    DSMatrixPrint(gA);
    printf("The matrix pA is: \n");
    DSMatrixPrint(pA);

	return 0;
}