//
//  designspacetest.c
//
//
//  Created by Jason Lomnitz on 5/13/13.
//
//
#include <stdio.h>
#include <string.h>
#include <designspace/DSStd.h>

int main(int argc, const char ** argv) {
        int i;
        char * strings[3] = {"\0"};
        strings[0] = strdup("x1. = a11 + 2*b31*x3 - b11*x1 - 2*b12*(x1^2)");
        strings[1] = strdup("x2. = b12*(x1^2) - b23*x2 - b22*x2");
        strings[2] = strdup("x3. = a31 + b23*x2 - b31*x3 - b33*x3");
        DSDesignSpace * ds;
        DSExpression ** expr = NULL;
        
        ds = DSDesignSpaceByParsingStrings(strings, NULL, 3);
        DSDesignSpaceCalculateCyclicalCases(ds);
        printf("Number of valid cases is: %i\n", DSDesignSpaceNumberOfValidCases(ds));
        printf("DSDesignSpaceNumberOfValidCases passed!\n");
        
        DSDesignSpaceFree(ds);
        printf("DSDesignSpaceFree passed!\n");
        
        return 0;
}