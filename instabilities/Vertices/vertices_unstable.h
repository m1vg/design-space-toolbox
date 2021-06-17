//
//  designspacetest.c
//
//
//  Created by Miguel Valderrama on 10/23/18.
//
//

void printPool(const DSVariablePool *pool);
#define DSDSGMA(x)                              ((x)->gma)
void printDictionary(DSDictionary *dictionary);


void setRemainingParameterValues(DSVariablePool *lowerBounds, DSVariablePool *upperBounds);

void loopValid(DSCase *aCase);
