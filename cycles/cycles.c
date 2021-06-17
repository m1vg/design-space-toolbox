
#include <stdio.h>
#include <string.h>
#include <designspace/DSStd.h>

int main(int argc, const char ** argv) {
    int i;
    char * strings[10] = {"\0"};
    strings[0] = strdup("X1. = a1*S*(K5^-2)*(X5^2)*(D5^-1) - b1*(K1^-2)*(X1^2)*(D1^-1)");
    strings[1] = strdup("X2. = a2*(r1^-1)*(D21^-1) + a2*(K21^-2)*(X1^2)*(D21^-1) - b2*(K2^-1)*(X2^1)*(D23^-1)*(D22^-1) - b2*(r3^-1)*(K23^-n)*(X3^2)*(K2^-1)*(X2^1)*(D23^-1)*(D22^-1)");
    strings[2] = strdup("X3. = a31*X2 + 2*a32*X4 - 2*a4*X3 - b3*X3");
    strings[3] = strdup("X4. = a4*X3 - a32*X4 - b4*X4");
    strings[4] = strdup("0  = 1 + (K1^-2)*(X1^2) -D1 ");
    strings[5] = strdup("0 = 1 + (K21^-2)*(X1^2) - D21 ");
    strings[6] = strdup("0 = 1 + (K2^-1)*(X2^1) - D22");
    strings[7] = strdup("0 = 1 + (K23^-2)*(X3^2) - D23 ");
    strings[8] = strdup("0 = 1 + (K5^-2)*(X5^2) - D5");
    strings[9] = strdup("0 = X1 + X5 - XT");
    DSDesignSpace * ds;
    DSVariablePool * Xd_a = NULL, *Xi = NULL;
    
    printf("initializing program \n");
    
    char * xd[6] = {"\0"};
    xd[0] = strdup("D1");
    xd[1] = strdup("D21");
    xd[2] = strdup("D22");
    xd[3] = strdup("D23");
    xd[4] = strdup("D5");
    xd[5] = strdup("XT");

    Xd_a = DSVariablePoolAlloc();
    DSVariablePoolAddVariableWithName(Xd_a, xd[0]);
    DSVariablePoolAddVariableWithName(Xd_a, xd[1]);
    DSVariablePoolAddVariableWithName(Xd_a, xd[2]);
    DSVariablePoolAddVariableWithName(Xd_a, xd[3]);
    DSVariablePoolAddVariableWithName(Xd_a, xd[4]);
    DSVariablePoolAddVariableWithName(Xd_a, xd[5]);
    
    char * xi[18] = {"\0"};
    xi[0] = strdup("K1");
    xi[1] = strdup("K2");
    xi[2] = strdup("K21");
    xi[3] = strdup("K23");
    xi[4] = strdup("K5");
    xi[5] = strdup("S");
    xi[6] = strdup("X5");
    xi[7] = strdup("a1");
    xi[8] = strdup("a2");
    xi[9] = strdup("a31");
    xi[10] = strdup("a32");
    xi[11] = strdup("a4");
    xi[12] = strdup("b1");
    xi[13] = strdup("b2");
    xi[14] = strdup("b3");
    xi[15] = strdup("b4");
    xi[16] = strdup("r1");
    xi[17] = strdup("r3");
    
    Xi = DSVariablePoolAlloc();
    for (i = 0; i < 18 ; i++ )
        DSVariablePoolAddVariableWithName(Xi, xi[i]);


    printf("The variable pool Xd_a is: \n");
    DSVariablePoolPrint(Xd_a);
    
    printf("The variable pool Xi is: \n");
    DSVariablePoolPrint(Xi);
    
    ds = DSDesignSpaceByParsingStringsWithXi(strings, Xd_a, Xi,  10);
    printf("design space constructed! \n");
    
    DSDesignSpaceSetSerial(ds, true);
    printf("setting serial \n");
    
    DSDesignSpaceCalculateCyclicalCases(ds);
    printf("Number of valid cases is: %i\n", DSDesignSpaceNumberOfValidCases(ds));
    printf("DSDesignSpaceNumberOfValidCases passed!\n");
    
    DSDesignSpaceFree(ds);
    printf("DSDesignSpaceFree passed!\n");
    
    return 0;
}
