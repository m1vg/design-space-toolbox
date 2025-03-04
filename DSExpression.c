/**
 * \file DSExpression.c
 * \brief Implementation file with functions for dealing with mathematical
 *        expressions.
 *
 * \details The DSExpression object is used internally to parse mathematical
 *          expressions and to evaluate these expressions using a variable pool.
 *          The DSExpression object supports scalar adition, multiplication, 
 *          powers and real valued functions of real variables.
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
 *
 * \todo Add options to register custom functions.
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "DSErrors.h"
#include "DSMemoryManager.h"
#include "DSVariable.h"
#include "DSExpression.h"
#include "DSExpressionTokenizer.h"
#include "DSMatrix.h"

#define DS_EXPRESSION_CONSTANT_BRANCH_INDEX     0
#define DS_EXPRESSION_STRING_INIT_LENGTH        10000 // Default value is 1000

static void dsExpressionToStringAdditionOperator(const DSExpression *current, char ** string, DSUInteger *length);


/**
 * \brief Allocates a node in the DSExpression tree that holds a constant
 *        double value.
 *
 * \details The function allocates a node in the expression tree and initializes
 *          it such that it is internally recognized as a constant double value,
 *          which is set upon calling this function.
 *
 * \param value A const double with the value assigned at the node.
 * \return A DSExpression pointer to the node initialized as a constant.
 *
 * \warning This function should not be called directly by the user and is 
 *          therefore not exposed; however, it is not a static function as it is
 *          called by the parsing tools in the DSExpressionGrammar files.
 */
extern DSExpression * dsExpressionAllocWithConstant(const double value)
{
        DSExpression *newNode = NULL;
        newNode = DSSecureCalloc(1, sizeof(DSExpression));
        DSExpressionSetConstant(newNode, value);
        newNode->type = DS_EXPRESSION_TYPE_CONSTANT;
        return newNode;
}

/**
 * \brief Allocates a node in the DSExpression tree that holds an operator node.
 *
 * \details The function allocates a node in the expression tree and initializes
 *          it such that it is internally recognized as a operator node with a,
 *          specified operator type: a '+', '*' or '^'.
 *          
 *
 * \param value A const char with the operator code.
 * \return A DSExpression pointer to the node initialized as an operator.
 * 
 * \warning This function should not be called directly by the user and is 
 *          therefore not exposed; however, it is not a static function as it is
 *          called by the parsing tools in the DSExpressionGrammar files.
 */
extern DSExpression * dsExpressionAllocWithOperator(const char op_code)
{
        DSExpression *newNode = NULL;
        switch (op_code) {
                case '=':
                        newNode = DSSecureCalloc(1, sizeof(DSExpression));
                        DSExpressionSetOperator(newNode, op_code);
                        break;
                case '<':
                        newNode = DSSecureCalloc(1, sizeof(DSExpression));
                        DSExpressionSetOperator(newNode, op_code);
                        break;
                case '>':
                        newNode = DSSecureCalloc(1, sizeof(DSExpression));
                        DSExpressionSetOperator(newNode, op_code);
                        break;
                case '.':
                        newNode = DSSecureCalloc(1, sizeof(DSExpression));
                        DSExpressionSetOperator(newNode, '.');
                        break;
                case '+':
                        newNode = DSSecureCalloc(1, sizeof(DSExpression));
                        DSExpressionSetOperator(newNode, '+');
                        /* First branch reserved for constants */
                        newNode->branches = DSSecureMalloc(sizeof(DSExpression *));
                        newNode->branches[DS_EXPRESSION_CONSTANT_BRANCH_INDEX] = dsExpressionAllocWithConstant(0.0);
                        newNode->numberOfBranches = 1;
                        break;
                case '*':
                        newNode = DSSecureCalloc(1, sizeof(DSExpression));
                        DSExpressionSetOperator(newNode, '*');
                        /* First branch reserved for constants */
                        newNode->branches = DSSecureMalloc(sizeof(DSExpression *));
                        newNode->branches[DS_EXPRESSION_CONSTANT_BRANCH_INDEX] = dsExpressionAllocWithConstant(1.0);
                        newNode->numberOfBranches = 1;
                        break;
                case '^':
                        newNode = DSSecureCalloc(1, sizeof(DSExpression));
                        DSExpressionSetOperator(newNode, op_code);
                        break;
                case '-':
                        DSError(M_DS_WRONG ": DSExpression does not internally use '-' operators", A_DS_ERROR);
                        break;
                case '/':
                        DSError(M_DS_WRONG ": DSExpression does not internally use '/' operators", A_DS_ERROR);
                        break;
                default:
                        DSError(M_DS_WRONG ": DSExpression found unrecognized operator.", A_DS_ERROR);
                        break;
        }
        return  newNode;
}

extern DSExpression * dsExpressionAllocWithVariableName(const char * name)
{
        DSExpression *newNode = NULL;
        if (name == NULL) {
                DSError(M_DS_WRONG ": name of variable is NULL", A_DS_ERROR);
                goto bail;
        }
        if (strlen(name) == 0) {
                DSError(M_DS_WRONG ": name of variable is empty", A_DS_ERROR);
                goto bail;
        }
        newNode = DSSecureCalloc(1, sizeof(DSExpression));
        DSExpressionSetVariable(newNode, strdup(name));
bail:
        return newNode;       
}

extern void DSExpressionFree(DSExpression *root)
{
        DSUInteger i;
        if (root == NULL) {
                DSError(M_DS_NULL ": Expression to free is NULL", A_DS_ERROR);
                goto bail;
        }
        for (i = 0; i < DSExpressionNumberOfBranches(root); i++) {
                DSExpressionFree(DSExpressionBranchAtIndex(root, i));
        }
        if (root->branches != NULL)
                DSSecureFree(root->branches);
        switch (DSExpressionType(root)) {
                case DS_EXPRESSION_TYPE_FUNCTION:
                case DS_EXPRESSION_TYPE_VARIABLE:
                        DSSecureFree(DSExpressionVariable(root));
                        break;
                default:
                        break;
        }
        DSSecureFree(root);
bail:
        return;
}

extern DSExpression * DSExpressionCopy(const DSExpression * expression)
{
        DSExpression * root = NULL;
        char * string = NULL;
        if (expression == NULL) {
                DSError(M_DS_NULL ": Expression to copy is NULL", A_DS_ERROR);
                goto bail;
        }
        string = DSExpressionAsString(expression);
        root = DSExpressionByParsingString(string);
        DSSecureFree(string);
bail:
        return root;
}

#if defined(__APPLE__) && defined (__MACH__)
#pragma mark - Factory functions
#endif

extern DSExpression * DSExpressionByParsingString(const char *string)
{
        DSExpression *root = NULL;
        void *parser = NULL;
        struct expression_token *tokens, *current;
        parse_expression_s parsed;
        if (string == NULL) {
                DSError(M_DS_WRONG ": String to parse is NULL", A_DS_ERROR);
                goto bail;
        }
        if (strlen(string) == 0) {
                DSError(M_DS_WRONG ": String to parse is empty", A_DS_WARN);
                goto bail;
        }
        tokens = DSExpressionTokenizeString(string);
        if (tokens == NULL) {
                DSError(M_DS_PARSE ": Token stream is NULL", A_DS_ERROR);
                goto bail;
        }
        parser = DSExpressionParserAlloc(DSSecureMalloc);
        current = tokens;
        parsed.wasSuccesful = true;
        while (current != NULL) {
                if (DSExpressionTokenType(current) == DS_EXPRESSION_TOKEN_START) {
                        current = DSExpressionTokenNext(current);
                        continue;
                }
                DSExpressionParser(parser, 
                                   DSExpressionTokenType(current), 
                                   (DSExpression *)current,
                                   &parsed);
                current = DSExpressionTokenNext(current);
        }
        DSExpressionParser(parser, 
                           0, 
                           NULL,
                           &parsed);
        DSExpressionParserFree(parser, DSSecureFree);
        DSExpressionTokenFree(tokens);
        if (parsed.wasSuccesful == true)
                root = parsed.root;
        else
                DSExpressionFree(parsed.root);
bail:
        return root;
}

extern DSExpression * DSExpressionAddExpressions(DSExpression *lvalue, DSExpression *rvalue)
{
        DSExpression * newRoot = NULL;
        if (lvalue == NULL && rvalue == NULL) {
                DSError(M_DS_NULL ": Expression is NULL", A_DS_ERROR);
                goto bail;
        }
        if (lvalue == NULL) {
                newRoot = rvalue;
                goto bail;
        }
        if (rvalue == NULL) {
                newRoot = lvalue;
                goto bail;
        }
        if (DSExpressionType(lvalue) == DS_EXPRESSION_TYPE_OPERATOR) {
                if (DSExpressionOperator(lvalue) == '=') {
                        newRoot = dsExpressionAllocWithOperator('=');
                        DSExpressionAddBranch(newRoot, DSExpressionAddExpressions(DSExpressionBranchAtIndex(lvalue, 0), rvalue));
                        DSExpressionAddBranch(newRoot, DSExpressionAddExpressions(DSExpressionBranchAtIndex(lvalue, 1), rvalue));
                        goto bail;
                }
                if (DSExpressionOperator(lvalue) == '<') {
                        newRoot = dsExpressionAllocWithOperator('<');
                        DSExpressionAddBranch(newRoot, DSExpressionAddExpressions(DSExpressionBranchAtIndex(lvalue, 0), rvalue));
                        DSExpressionAddBranch(newRoot, DSExpressionAddExpressions(DSExpressionBranchAtIndex(lvalue, 1), rvalue));
                        goto bail;
                }
                if (DSExpressionOperator(lvalue) == '>') {
                        newRoot = dsExpressionAllocWithOperator('>');
                        DSExpressionAddBranch(newRoot, DSExpressionAddExpressions(DSExpressionBranchAtIndex(lvalue, 0), rvalue));
                        DSExpressionAddBranch(newRoot, DSExpressionAddExpressions(DSExpressionBranchAtIndex(lvalue, 1), rvalue));
                        goto bail;
                }
        }
        newRoot = dsExpressionAllocWithOperator('+');
        DSExpressionAddBranch(newRoot, lvalue);
        DSExpressionAddBranch(newRoot, rvalue);
bail:
        return newRoot;
}

extern DSExpression * DSExpressionSubstractExpressions(DSExpression *lvalue, DSExpression *rvalue)
{
        DSExpression * newRoot = NULL, *temp;
        DSUInteger i;
        if (lvalue == NULL && rvalue == NULL) {
                DSError(M_DS_NULL ": Expression is NULL", A_DS_ERROR);
                goto bail;
        }
        if (lvalue == NULL) {
                newRoot = dsExpressionAllocWithOperator('*');
                DSExpressionAddBranch(newRoot, dsExpressionAllocWithConstant(-1.0));
                DSExpressionAddBranch(newRoot, rvalue);
                goto bail;
        }
        if (rvalue == NULL) {
                newRoot = lvalue;
                goto bail;
        }
        newRoot = dsExpressionAllocWithOperator('+');
        DSExpressionAddBranch(newRoot, lvalue);
        if (DSExpressionType(rvalue) == DS_EXPRESSION_TYPE_OPERATOR) {
                if (DSExpressionOperator(rvalue) == '+') {
                        for (i = 0; i < DSExpressionNumberOfBranches(rvalue); i++) {
                                temp = dsExpressionAllocWithOperator('*');
                                DSExpressionAddBranch(temp, dsExpressionAllocWithConstant(-1.0));
                                DSExpressionAddBranch(temp, DSExpressionBranchAtIndex(rvalue, i));
                                DSExpressionAddBranch(newRoot, temp);
                                rvalue->branches[i] = NULL;
                        }
                        rvalue->numberOfBranches = 0;
                        DSExpressionFree(rvalue);
                        goto bail;
                }
        }
        temp = dsExpressionAllocWithOperator('*');
        DSExpressionAddBranch(temp, dsExpressionAllocWithConstant(-1.0));
        DSExpressionAddBranch(temp, rvalue);
        DSExpressionAddBranch(newRoot, temp);
bail:
        return newRoot;
}

extern DSExpression * DSExpressionMultiplyExpressionByConstant(DSExpression *expression, double constant)
{
        DSExpression * newRoot = NULL, *temp;
        if (expression == NULL) {
                goto bail;
        }
        newRoot = dsExpressionAllocWithOperator('*');
        temp = dsExpressionAllocWithConstant(constant);
        DSExpressionAddBranch(newRoot, expression);
        DSExpressionAddBranch(newRoot, temp);
bail:
        return newRoot;
}

static DSExpression * dsExpressionCompressConstantVariableNode(const DSExpression * current, const DSVariablePool * assumedConstant)
{
        DSUInteger i;
        DSExpression * compressed = NULL;
        if (current == NULL) {
                goto bail;
        }
        if (DSExpressionType(current) == DS_EXPRESSION_TYPE_VARIABLE) {
                if (DSVariablePoolHasVariableWithName(assumedConstant, DSExpressionVariable(current)) == true) {
                        compressed = dsExpressionAllocWithConstant(DSVariablePoolValueForVariableWithName(assumedConstant, DSExpressionVariable(current)));
                        goto bail;
                } else {
                        compressed = DSExpressionCopy(current);
                        goto bail;
                }
        }
        if (DSExpressionType(current) == DS_EXPRESSION_TYPE_FUNCTION ||
            DSExpressionType(current) == DS_EXPRESSION_TYPE_CONSTANT) {
                compressed = DSExpressionCopy(current);
                goto bail;
        }
        if (DSExpressionType(current) != DS_EXPRESSION_TYPE_OPERATOR) {
                DSError(M_DS_WRONG "Expression Node is Undefined", A_DS_ERROR);
                goto bail;
        }
        compressed = dsExpressionAllocWithOperator(DSExpressionOperator(current));
        for (i = 0; i < DSExpressionNumberOfBranches(current); i++) {
                DSExpressionAddBranch(compressed,
                                      dsExpressionCompressConstantVariableNode(DSExpressionBranchAtIndex(current, i),
                                                                               assumedConstant));
        }
bail:
        return compressed;
}

extern DSExpression * DSExpressionByCompressingConstantVariables(const DSExpression *expression, const DSVariablePool * assumedConstant)
{
        DSExpression * newExpression = NULL;
        newExpression = dsExpressionCompressConstantVariableNode(expression, assumedConstant);
bail:
        return newExpression;
}

#if defined(__APPLE__) && defined (__MACH__)
#pragma mark - Branch adding functions
#endif

static const bool dsExpressionOperatorBranchIsZero(DSExpression * expression)
{
        bool isZero = false;
        if (expression == NULL) {
                DSError(M_DS_NULL ": Branch being added is NULL", A_DS_ERROR);
                goto bail;
        }
        if (DSExpressionType(expression) != DS_EXPRESSION_TYPE_OPERATOR) {
                goto bail;
        }
        switch (DSExpressionOperator(expression)) {
                case '+':
                        if (DSExpressionNumberOfBranches(expression) == 0) {
                                isZero = true;
                                break;
                        } if (DSExpressionNumberOfBranches(expression) == 1) {
                                if (DSExpressionConstant(DSExpressionBranchAtIndex(expression, 0)) == 0.0)
                                        isZero = true;
                                break;
                        }
                        break;
                case '*':
                        if (DSExpressionNumberOfBranches(expression) >= 1) {
                                if (DSExpressionConstant(DSExpressionBranchAtIndex(expression, 0)) == 0.0)
                                        isZero = true;
                                break;
                        }
                        break;
                default:
                        break;
        }
bail:
        return isZero;
}

static void DSExpressionAddNonConstantBranch(DSExpression *expression, DSExpression *branch)
{
        if (expression == NULL) {
                DSError(M_DS_NULL ": Expression root is NULL", A_DS_ERROR);
                goto bail;
        }
        if (branch == NULL) {
                DSError(M_DS_NULL ": Branch being added is NULL", A_DS_ERROR);
                goto bail;
        }
        if (DSExpressionType(expression) != DS_EXPRESSION_TYPE_OPERATOR) {
                DSError(M_DS_WRONG ": Expression root is not an operator", A_DS_ERROR);
                goto bail;
        }
        if (dsExpressionOperatorBranchIsZero(branch) == true) {
                DSExpressionFree(branch);
                goto bail;
        }
        if (expression->branches == NULL)
                expression->branches = DSSecureMalloc(sizeof(DSExpression *)*(DSExpressionNumberOfBranches(expression)+1));
        else
                expression->branches = DSSecureRealloc(expression->branches, 
                                                       sizeof(DSExpression *)*(DSExpressionNumberOfBranches(expression)+1));
        expression->branches[DSExpressionNumberOfBranches(expression)] = branch;
        expression->numberOfBranches++;
bail:
        return;
}

static void DSExpressionAddConstantBranch(DSExpression *expression, DSExpression *branch)
{
        DSExpression * newBranch;
        double constant;
        if (expression == NULL) {
                DSError(M_DS_NULL ": Expression root is NULL", A_DS_ERROR);
                goto bail;
        }
        if (branch == NULL) {
                DSError(M_DS_NULL ": Branch being added is NULL", A_DS_ERROR);
                goto bail;
        }
        if (DSExpressionType(expression) != DS_EXPRESSION_TYPE_OPERATOR) {
                DSError(M_DS_WRONG ": Expression root is not an operator", A_DS_ERROR);
                goto bail;
        }
        if (DSExpressionType(branch) != DS_EXPRESSION_TYPE_CONSTANT) {
                DSError(M_DS_WRONG ": branch expression is not a constant", A_DS_ERROR);
                goto bail;
        }
        constant = DSExpressionConstant(branch);
        switch (DSExpressionOperator(expression)) {
                case '+':
                        newBranch = DSExpressionBranchAtIndex(expression, 0);
                        if (newBranch == NULL) {
                                DSError(M_DS_NULL ": Constant branch is null" , A_DS_ERROR);
                                break;
                        }
                        DSExpressionSetConstant(newBranch,
                                                DSExpressionConstant(newBranch)+constant);
                        DSExpressionFree(branch);
                        break;
                case '*':
                        newBranch = DSExpressionBranchAtIndex(expression, 0);
                        if (newBranch == NULL) {
                                DSError(M_DS_NULL ": Constant branch is null" , A_DS_ERROR);
                                break;
                        }
                        DSExpressionSetConstant(newBranch,
                                                DSExpressionConstant(newBranch)*constant);
                        DSExpressionFree(branch);
                        break;
                case '>':
                case '<':
                case '=':
                        DSExpressionAddNonConstantBranch(expression, branch);
                        break;
                case '.':
                        DSExpressionAddNonConstantBranch(expression, branch);
                        break;
                case '^':
                        DSExpressionAddNonConstantBranch(expression, branch);
                        break;
                default:
                        break;
        }
bail:
        return;
}

static void dsExpressionAddBranchToFunction(DSExpression *expression, DSExpression *branch)
{
        if (expression == NULL) {
                DSError(M_DS_NULL ": Expression root is NULL", A_DS_ERROR);
                goto bail;
        }
        if (branch == NULL) {
                DSError(M_DS_NULL ": Branch being added is NULL", A_DS_ERROR);
                goto bail;
        }
        expression->branches = DSSecureMalloc(sizeof(DSExpression *)*(DSExpressionNumberOfBranches(expression)+1));
        expression->branches[DSExpressionNumberOfBranches(expression)] = branch;
        expression->numberOfBranches++;
        expression->type = DS_EXPRESSION_TYPE_FUNCTION;
bail:
        return;
}
extern void DSExpressionAddBranch(DSExpression *expression, DSExpression *branch)
{
        DSUInteger i;
        char expressionOperator = '?';
        if (expression == NULL) {
                DSError(M_DS_NULL ": Expression root is NULL", A_DS_ERROR);
                goto bail;
        }
        if (branch == NULL) {
                DSError(M_DS_NULL ": Branch being added is NULL", A_DS_ERROR);
                goto bail;
        }
        if (DSExpressionType(expression) == DS_EXPRESSION_TYPE_VARIABLE) {
                dsExpressionAddBranchToFunction(expression, branch);
                goto bail;
        }
        if (DSExpressionType(expression) != DS_EXPRESSION_TYPE_OPERATOR) {
                DSError(M_DS_WRONG ": Adding branch to non-operator expression", A_DS_ERROR);
                goto bail;
        }
        if (DSExpressionType(branch) == DS_EXPRESSION_TYPE_UNDEFINED) {
                DSError(M_DS_WRONG ": branch expression type is undefined", A_DS_ERROR);
                goto bail;
        }
        if (DSExpressionType(branch)  == DS_EXPRESSION_TYPE_VARIABLE || 
            DSExpressionType(branch) == DS_EXPRESSION_TYPE_FUNCTION) {
                DSExpressionAddNonConstantBranch(expression, branch);
                goto bail;
        }
        if (DSExpressionType(branch) == DS_EXPRESSION_TYPE_CONSTANT) {
                DSExpressionAddConstantBranch(expression, branch);
                goto bail;
        }
        if (DSExpressionNumberOfBranches(branch) < 2 && DSExpressionOperator(branch) != '.') {
                if (DSExpressionNumberOfBranches(branch) == 1) {
                        DSExpressionAddBranch(expression, DSExpressionBranchAtIndex(branch, 0));
                        branch->numberOfBranches = 0;
                        DSExpressionFree(branch);
                } else {
                        DSError(M_DS_WRONG ": branch has insufficient branches", A_DS_ERROR);
                }
                goto bail;
        }
        expressionOperator = DSExpressionOperator(branch);
        switch (DSExpressionOperator(expression)) {
                case '+':
                        if (expressionOperator == '+') {
                                for (i = 0; i < DSExpressionNumberOfBranches(branch); i++)
                                        DSExpressionAddBranch(expression, DSExpressionBranchAtIndex(branch, i));
                                branch->numberOfBranches = 0;
                                DSExpressionFree(branch);
                        } else {
                                DSExpressionAddNonConstantBranch(expression, branch);
                        }
                        break;
                case '*':
                        if (expressionOperator == '*') {
                                for (i = 0; i < DSExpressionNumberOfBranches(branch); i++)
                                        DSExpressionAddBranch(expression, DSExpressionBranchAtIndex(branch, i));
                                branch->numberOfBranches = 0;
                                DSExpressionFree(branch);
                        } else {
                                DSExpressionAddNonConstantBranch(expression, branch);
                        }
                        break;
                case '=':
                        DSExpressionAddNonConstantBranch(expression, branch);
                        break;
                case '<':
                        DSExpressionAddNonConstantBranch(expression, branch);
                        break;
                case '>':
                        DSExpressionAddNonConstantBranch(expression, branch);
                        break;
                case '.':
                        DSExpressionAddNonConstantBranch(expression, branch);
                        break;
                case '^':
                        DSExpressionAddNonConstantBranch(expression, branch);
                        break;
                default:
                        DSError(M_DS_WRONG ": Operator for expression root is undefined", A_DS_ERROR);
                        break;
        }
bail:
        return;
}

#if defined(__APPLE__) && defined (__MACH__)
#pragma mark - Expression properties
#endif


#if defined(__APPLE__) && defined (__MACH__)
#pragma mark - Utility functions
#endif

#define ds_function_index_log    0
#define ds_function_index_ln     1
#define ds_function_index_log10  2
#define ds_function_index_cos    3
#define ds_function_index_sin    4
#define ds_function_index_abs    5
#define ds_function_index_sign   6
#define ds_function_index_sqrt   7
#define ds_function_index_real   8
#define ds_function_index_imag   9




static double dsExpressionEvaluateMathematicalFunction(const DSExpression *function, const DSVariablePool * pool)
{
        double value = 0, eval = NAN;
        double complex complex_value;
        int functionIndex;
        DSVariablePool *functionNames = NULL;
        if (function == NULL) {
                DSError(M_DS_NULL ": Expression node is null", A_DS_ERROR);
                goto bail;
        }
        if (DSExpressionType(function) != DS_EXPRESSION_TYPE_FUNCTION) {
                DSError(M_DS_WRONG ": Expression node must be a function", A_DS_ERROR);
                goto bail;
        }
        functionNames = DSVariablePoolByParsingString("log : 1, ln : 1, log10 : 1, cos : 1, sin : 1, abs : 1, sign : 1, sqrt : 1, real : 1, imag : 1");
        if (DSVariablePoolHasVariableWithName(functionNames, DSExpressionVariable(function)) == false) {
                DSError(M_DS_WRONG ": Function name not recognized", A_DS_ERROR);
                goto bail;
        }
        functionIndex =DSVariablePoolIndexOfVariableWithName(functionNames, DSExpressionVariable(function));
        if (functionIndex == ds_function_index_real || functionIndex == ds_function_index_imag) {
                complex_value = DSExpressionEvaluateComplexWithVariablePool(DSExpressionBranchAtIndex(function, 0), pool);
        } else {
                value = DSExpressionEvaluateWithVariablePool(DSExpressionBranchAtIndex(function, 0), pool);
        }
        switch (functionIndex) {
                case ds_function_index_ln:
                        eval = log(value);
                        break;
                case ds_function_index_log:
                case ds_function_index_log10:
                        eval = log10(value);
                        break;
                case ds_function_index_cos:
                        eval = cos(value);
                        break;
                case ds_function_index_sin:
                        eval = sin(value);
                        break;
                case ds_function_index_abs:
                        eval = fabs(value);
                        break;
                case ds_function_index_sign:
                        if (value > 0.0f)
                                eval = 1.0f;
                        else if (value < 0.0f)
                                eval = -1.0f;
                        else
                                eval = 0.0f;
                        break;
                case ds_function_index_sqrt:
                        eval = sqrt(value);
                        break;
                case ds_function_index_real:
                        eval = creal(complex_value);
                        break;
                case ds_function_index_imag:
                        eval = cimag(complex_value);
                        break;
                default:
                        break;
        }
bail:
        if (functionNames != NULL)
                DSVariablePoolFree(functionNames);
        return eval;
}

extern double DSExpressionEvaluateWithVariablePool(const DSExpression *expression, const DSVariablePool *pool)
{
        double value = NAN;
        DSVariable *variable = NULL;
        DSUInteger i;
        if (expression == NULL) {
                DSError(M_DS_NULL ": Expression is NULL", A_DS_ERROR);
                goto bail;
        }
        switch (DSExpressionType(expression)) {
                case DS_EXPRESSION_TYPE_VARIABLE:
                        if (pool != NULL) {
                                if (DSVariablePoolHasVariableWithName(pool, DSExpressionVariable(expression)) == true) {
                                        variable = DSVariablePoolVariableWithName(pool, DSExpressionVariable(expression));
                                        value = DSVariableValue(variable);
                                } else {
                                        DSError(M_DS_WRONG ": Variable pool does not have variable.", A_DS_ERROR);
                                        printf("[%s] Not Found.\n", DSExpressionVariable(expression));
                                }
                        } else {
                                DSError(M_DS_VAR_NULL, A_DS_ERROR);
                        }
                        
                        break;
                case DS_EXPRESSION_TYPE_CONSTANT:
                        value = DSExpressionConstant(expression);
                        break;
                case DS_EXPRESSION_TYPE_FUNCTION:
                        value=dsExpressionEvaluateMathematicalFunction(expression, pool);
                        break;
                case DS_EXPRESSION_TYPE_OPERATOR:
                        switch (DSExpressionOperator(expression)) {
                                case '+':
                                        value = 0;
                                        for (i = 0; i < DSExpressionNumberOfBranches(expression); i++)
                                                value += DSExpressionEvaluateWithVariablePool(DSExpressionBranchAtIndex(expression, i), pool);
                                        break;
                                case '*':
                                        value = 1;
                                        for (i = 0; i < DSExpressionNumberOfBranches(expression); i++)
                                                value *= DSExpressionEvaluateWithVariablePool(DSExpressionBranchAtIndex(expression, i), pool);
                                        break;
                                case '^':
                                        value = pow(DSExpressionEvaluateWithVariablePool(DSExpressionBranchAtIndex(expression, 0), pool),
                                                    DSExpressionEvaluateWithVariablePool(DSExpressionBranchAtIndex(expression, 1), pool));
                                        break;
                                default:
                                        DSError(M_DS_WRONG "Operators cannot be evaluated as a function", A_DS_WARN);
                                        value = NAN;
                                        goto bail;
                                        break;
                        }
                        break;
                default:
                        break;
        }
bail:
        return value;
}


static double complex dsExpressionEvaluateMathematicalFunctionComplex(const DSExpression *function, const DSVariablePool * pool)
{
        double complex value, eval = NAN;
        DSVariablePool *functionNames = NULL;
        if (function == NULL) {
                DSError(M_DS_NULL ": Expression node is null", A_DS_ERROR);
                goto bail;
        }
        if (DSExpressionType(function) != DS_EXPRESSION_TYPE_FUNCTION) {
                DSError(M_DS_WRONG ": Expression node must be a function", A_DS_ERROR);
                goto bail;
        }
        functionNames = DSVariablePoolByParsingString("log : 1, ln : 1, log10 : 1, cos : 1, sin : 1, abs : 1, sign : 1, sqrt : 1, real : 1, imag : 1");
        if (DSVariablePoolHasVariableWithName(functionNames, DSExpressionVariable(function)) == false) {
                DSError(M_DS_WRONG ": Function name not recognized", A_DS_ERROR);
                goto bail;
        }
        value = DSExpressionEvaluateComplexWithVariablePool(DSExpressionBranchAtIndex(function, 0), pool);
        switch (DSVariablePoolIndexOfVariableWithName(functionNames, DSExpressionVariable(function))) {
                case ds_function_index_ln:
                        eval = clog(value);
                        break;
                case ds_function_index_log:
                case ds_function_index_log10:
                        eval = log10(creal(value));
                        if (cimag(value) == 0.0) {
                                eval = log10(creal(value));
                        } else {
                                DSError(M_DS_NOT_IMPL ": Using log10 of real part.", A_DS_WARN);
                        }
                        break;
                case ds_function_index_cos:
                        eval = ccos(value);
                        break;
                case ds_function_index_sin:
                        eval = csin(value);
                        break;
                case ds_function_index_abs:
                        eval = cabs(value);
                        break;
                case ds_function_index_sign:
                        if (creal(value) > 0.0f)
                                eval = 1.0f;
                        else if (creal(value) < 0.0f)
                                eval = -1.0f;
                        else if (cimag(value) > 0.0f)
                                eval = 1.0f;
                        else if (cimag(value) < 0.0f)
                                eval = -1.0f;
                        else
                                eval = 0.0f;
                        break;
                case ds_function_index_sqrt:
                        eval = csqrt(value);
                        break;
                case ds_function_index_real:
                        eval = creal(value);
                        break;
                case ds_function_index_imag:
                        eval = cimag(value);
                        break;
                default:
                        break;
        }
bail:
        if (functionNames != NULL)
                DSVariablePoolFree(functionNames);
        return eval;
}

extern double complex DSExpressionEvaluateComplexWithVariablePool(const DSExpression *expression, const DSVariablePool *pool)
{
        double complex value = NAN;
        DSVariable *variable = NULL;
        DSUInteger i;
        if (expression == NULL) {
                DSError(M_DS_NULL ": Expression is NULL", A_DS_ERROR);
                goto bail;
        }
        switch (DSExpressionType(expression)) {
                case DS_EXPRESSION_TYPE_VARIABLE:
                        if (strcmp(DSExpressionVariable(expression), DSExpressionImaginaryNumber) == 0) {
                                value = I;
                        } else if (pool != NULL) {
                                if (DSVariablePoolHasVariableWithName(pool, DSExpressionVariable(expression)) == true) {
                                        variable = DSVariablePoolVariableWithName(pool, DSExpressionVariable(expression));
                                        value = DSVariableValue(variable);
                                } else {
                                        DSError(M_DS_WRONG ": Variable pool does not have variable.", A_DS_ERROR);
                                        printf("[%s] Not Found.\n", DSExpressionVariable(expression));
                                }
                        } else {
                                DSError(M_DS_VAR_NULL, A_DS_ERROR);
                        }
                        
                        break;
                case DS_EXPRESSION_TYPE_CONSTANT:
                        value = DSExpressionConstant(expression);
                        break;
                case DS_EXPRESSION_TYPE_FUNCTION:
                        value=dsExpressionEvaluateMathematicalFunctionComplex(expression, pool);
                        break;
                case DS_EXPRESSION_TYPE_OPERATOR:
                        switch (DSExpressionOperator(expression)) {
                                case '+':
                                        value = 0;
                                        for (i = 0; i < DSExpressionNumberOfBranches(expression); i++)
                                                value += DSExpressionEvaluateComplexWithVariablePool(DSExpressionBranchAtIndex(expression, i), pool);
                                        break;
                                case '*':
                                        value = 1;
                                        for (i = 0; i < DSExpressionNumberOfBranches(expression); i++)
                                                value *= DSExpressionEvaluateComplexWithVariablePool(DSExpressionBranchAtIndex(expression, i), pool);
                                        break;
                                case '^':
                                        value = cpow(DSExpressionEvaluateComplexWithVariablePool(DSExpressionBranchAtIndex(expression, 0), pool),
                                                    DSExpressionEvaluateComplexWithVariablePool(DSExpressionBranchAtIndex(expression, 1), pool));
                                        break;
                                default:
                                        DSError(M_DS_WRONG "Operators cannot be evaluated as a function", A_DS_WARN);
                                        value = NAN;
                                        goto bail;
                                        break;
                        }
                        break;
                default:
                        break;
        }
bail:
        return value;
}


extern DSExpression * DSExpressionEquationLHSExpression(const DSExpression *expression)
{
        DSExpression * lhs = NULL;
        if (expression == NULL) {
                DSError(M_DS_NULL ": Expression is NULL", A_DS_ERROR);
                goto bail;
        }
        if (DSExpressionType(expression) != DS_EXPRESSION_TYPE_OPERATOR) {
                DSError(M_DS_WRONG ": Expression is not an equation", A_DS_ERROR);
                goto bail;
        }
        if (DSExpressionOperator(expression) != '=' && DSExpressionOperator(expression) != '<' && DSExpressionOperator(expression) != '>' ) {
                DSError(M_DS_WRONG ": Expression is not an equation", A_DS_ERROR);
                goto bail;
        }
        if (DSExpressionNumberOfBranches(expression) < 2) {
                DSError(M_DS_WRONG ": Equation does not have a right hand side and left hand side", A_DS_ERROR);
                goto bail;                
        }
        switch (DSExpressionType(expression)) {
                case DS_EXPRESSION_TYPE_VARIABLE:
                        break;
                case DS_EXPRESSION_TYPE_CONSTANT:
                        break;
                case DS_EXPRESSION_TYPE_FUNCTION:
                        break;
                case DS_EXPRESSION_TYPE_OPERATOR:
                        switch (DSExpressionOperator(expression)) {
                                case '>':
                                case '<':
                                case '=':
                                        lhs = DSExpressionCopy(DSExpressionBranchAtIndex(expression, 0));
                                        break;
                                default:
                                        break;
                        }
                        break;
                default:
                        break;
        }
bail:
        return lhs;
}

extern DSExpression * DSExpressionEquationRHSExpression(const DSExpression *expression)
{
        DSExpression * rhs = NULL;
        if (expression == NULL) {
                DSError(M_DS_NULL ": Expression is NULL", A_DS_ERROR);
                goto bail;
        }
        if (DSExpressionType(expression) != DS_EXPRESSION_TYPE_OPERATOR) {
                DSError(M_DS_WRONG ": Expression is not an equation", A_DS_ERROR);
                goto bail;
        }
        if (DSExpressionOperator(expression) != '=' && DSExpressionOperator(expression) != '<' && DSExpressionOperator(expression) != '>') {
                DSError(M_DS_WRONG ": Expression is not an equation", A_DS_ERROR);
                goto bail;
        }
        switch (DSExpressionType(expression)) {
                case DS_EXPRESSION_TYPE_VARIABLE:
                        break;
                case DS_EXPRESSION_TYPE_CONSTANT:
                        break;
                case DS_EXPRESSION_TYPE_FUNCTION:
                        break;
                case DS_EXPRESSION_TYPE_OPERATOR:
                        switch (DSExpressionOperator(expression)) {
                                case '>':
                                case '<':
                                case '=':
                                        rhs = DSExpressionCopy(DSExpressionBranchAtIndex(expression, 1));
                                        break;
                                default:
                                        break;
                        }
                        break;
                default:
                        break;
        }
bail:
        return rhs;
}

static void dsExpressionVariablesInExpressionInternal(const DSExpression * current, DSVariablePool * pool)
{
        DSUInteger i;
        if (current == NULL) {
                DSError(M_DS_NULL ": Expression is NULL", A_DS_ERROR);
                goto bail;
        }
        switch (DSExpressionType(current)) {
                case DS_EXPRESSION_TYPE_VARIABLE:
                        if (strcmp(DSExpressionVariable(current), DSExpressionImaginaryNumber) == 0)
                                break;
                        if (DSVariablePoolHasVariableWithName(pool, DSExpressionVariable(current)) == false)
                                DSVariablePoolAddVariableWithName(pool, DSExpressionVariable(current));
                        break;
                case DS_EXPRESSION_TYPE_CONSTANT:
                        break;
                case DS_EXPRESSION_TYPE_FUNCTION:
                        dsExpressionVariablesInExpressionInternal(DSExpressionBranchAtIndex(current, 0), pool);
                        break;
                case DS_EXPRESSION_TYPE_OPERATOR:
                        for (i = 0; i < DSExpressionNumberOfBranches(current); i++)
                                dsExpressionVariablesInExpressionInternal(DSExpressionBranchAtIndex(current, i), pool);
                        break;
                default:
                        break;
        }
bail:
        return;
}

extern DSVariablePool * DSExpressionVariablesInExpression(const DSExpression * expression)
{
        DSVariablePool * variables = NULL;
        if (expression == NULL) {
                DSError(M_DS_NULL ": Expression is NULL", A_DS_ERROR);
                goto bail;
        }
        variables = DSVariablePoolAlloc();
        dsExpressionVariablesInExpressionInternal(expression, variables);
bail:
        return variables;
}

#if defined(__APPLE__) && defined (__MACH__)
#pragma mark - Utility functions
#endif

static bool operatorIsLowerPrecedence(char op1, char op2)
{
        bool isLower = false;
        char * precedence = ".^*+<>=";
        DSUInteger i, index1, index2;
        for (i = 0; i < strlen(precedence); i++) {
                if (precedence[i] == op1)
                        index1 = i;
                if (precedence[i] == op2)
                        index2 = i;
        }
        if (index2 > index1)
                isLower = true;
        return isLower;
}

static DSUInteger dsExpressionConstantNumberOfDecimals(double constant)
{
        DSUInteger count = 0;
        constant = fabs(constant);
        constant -= floor(constant);
        while (constant > 1e-14 && count < 16) {
                count++;
                constant*=10;
                constant -= floor(constant);
        }
        return count;
}

static void dsExpressionToStringInternal(const DSExpression *current, char ** string, DSUInteger *length)
{
        DSUInteger i;
        DSExpression * branch, * constantBranch;
        double constant;
        char temp[100] = {'\0'};
        if (current == NULL) {
                DSError(M_DS_NULL ": Node to print is nil", A_DS_ERROR);
                goto bail;
        }
        switch (DSExpressionType(current)) {
                case DS_EXPRESSION_TYPE_CONSTANT:
                        sprintf(temp, "%.*lf", dsExpressionConstantNumberOfDecimals(DSExpressionConstant(current)), DSExpressionConstant(current));
                        break;
                case DS_EXPRESSION_TYPE_VARIABLE:
                        sprintf(temp, "%s", DSExpressionVariable(current));
                        break;
                case DS_EXPRESSION_TYPE_OPERATOR:
                        constantBranch = DSExpressionBranchAtIndex(current, 0);
                        if (constantBranch == NULL) {
                                DSError(M_DS_NULL ": Constant branch is NULL", A_DS_ERROR);
                                break;
                        }
                        constant = DSExpressionConstant(constantBranch);
                        for (i=0; i < DSExpressionNumberOfBranches(current); i++) {
                                if (i == 0 && DSExpressionOperator(current) == '+' && constant == 0.0)
                                        continue;
                                if (i == 0 && DSExpressionOperator(current) == '*' && constant == 1.0)
                                        continue;
                                if (i == 0 && DSExpressionOperator(current) == '*' && constant == 0.0) {
                                        strncat(*string, "0", *length-strlen(*string));
                                        break;
                                }
                                if (i == 0 && DSExpressionOperator(current) == '*' && constant == -1.0) {
                                        strncat(*string, "-", *length-strlen(*string));
                                        continue;
                                }
                                branch = DSExpressionBranchAtIndex(current, i);
                                if (DSExpressionType(branch) == DS_EXPRESSION_TYPE_OPERATOR &&
                                    operatorIsLowerPrecedence(DSExpressionOperator(current), DSExpressionOperator(branch)))
                                        strncat(*string, "(", *length-strlen(*string));
                                dsExpressionToStringInternal(branch, string, length);
                                if (DSExpressionType(branch) == DS_EXPRESSION_TYPE_OPERATOR &&
                                    operatorIsLowerPrecedence(DSExpressionOperator(current), DSExpressionOperator(branch)))
                                        strncat(*string, ")", *length-strlen(*string));
                                if (i < DSExpressionNumberOfBranches(current)-1 || DSExpressionOperator(current) == '.') {
                                        if (DSExpressionOperator(current) == '+' &&
                                            DSExpressionType(DSExpressionBranchAtIndex(current, i+1)) == DS_EXPRESSION_TYPE_OPERATOR) {
                                                if (DSExpressionConstant(DSExpressionBranchAtIndex(DSExpressionBranchAtIndex(current, i+1), 0)) < 0) {
                                                        continue;
                                                }
                                        }
                                        sprintf(temp, "%c", DSExpressionOperator(current));
                                        strncat(*string,  temp, *length-strlen(*string));
                                        temp[0] = '\0';
                                }
                        }
                        break;
                case DS_EXPRESSION_TYPE_FUNCTION:
                        sprintf(temp, "%s(", DSExpressionVariable(current));
                        if (strlen(*string)+strlen(temp) >= *length) {
                                length += DS_EXPRESSION_STRING_INIT_LENGTH;
                                *string = DSSecureRealloc(string, sizeof(char)**length);
                        }
                        strncat(*string, temp, *length-strlen(*string));
                        dsExpressionToStringInternal(DSExpressionBranchAtIndex(current, 0), string, length);
                        strncat(*string, ")", *length-strlen(*string));
                        temp[0] = '\0';
                        break;
                default:
                        break;
        }
        if (strlen(*string)+strlen(temp) >= *length) {
                length += DS_EXPRESSION_STRING_INIT_LENGTH;
                *string = DSSecureRealloc(string, sizeof(char)**length);
        }
        strncat(*string,  temp, *length-strlen(*string));
bail:
        return;
}

extern char * DSExpressionAsString(const DSExpression *expression)
{
        DSUInteger length = DS_EXPRESSION_STRING_INIT_LENGTH;
        char * string = NULL;
        if (expression == NULL) {
                DSError(M_DS_NULL ": Node to print is null", A_DS_ERROR);
                goto bail;
        }
        string = DSSecureCalloc(sizeof(char), length);
        dsExpressionToStringInternal(expression, &string, &length);
bail:
        return string;
}

extern char * DSExpressionBranchAtIndexAsString(const DSExpression *expression, DSUInteger index){
        
    char * string = NULL;
    
    if (expression == NULL) {
            DSError(M_DS_NULL ": Node to print is null", A_DS_ERROR);
            goto bail;
    }
    
    if (index > DSExpressionNumberOfBranches(expression))
        goto bail;
    
    DSExpression *expression_n = NULL;
    DSUInteger length = DS_EXPRESSION_STRING_INIT_LENGTH, i;
    
    string = DSSecureCalloc(sizeof(char), length);
    expression_n = DSExpressionBranchAtIndex(expression, index);
    dsExpressionToStringInternal(expression_n, &string, &length);
    
    
bail:
    return string;
}

extern DSUInteger DSExpressionNumberOfTerms(const DSExpression *expression){
    
    DSUInteger numberOfTerms = 0;
    if (expression == NULL)
        goto bail;
    numberOfTerms = DSExpressionNumberOfBranches(expression);
bail:
    return numberOfTerms;
}

static void expressionToLatexStringInternal(const DSExpression *current, char ** string, DSUInteger *length, const DSDictionary * substitutionDict)
{
        DSUInteger i;
        DSExpression * branch, *constantBranch;
        double constant;
        char * name, * altName, * open, *close;
        char temp[100] = {'\0'};
        if (current == NULL) {
                DSError(M_DS_NULL ": Node to print is nil", A_DS_ERROR);
                goto bail;
        }
        switch (DSExpressionType(current)) {
                case DS_EXPRESSION_TYPE_CONSTANT:
                        sprintf(temp, "%.*lf", dsExpressionConstantNumberOfDecimals(DSExpressionConstant(current)), DSExpressionConstant(current));
                        break;
                case DS_EXPRESSION_TYPE_VARIABLE:
                        name = DSExpressionVariable(current);
                        altName = (char *)DSDictionaryValueForName(substitutionDict, name);
                        if (altName != NULL) {
                                name = altName;
                        }
                        sprintf(temp, "%s ", name);
                        break;
                case DS_EXPRESSION_TYPE_OPERATOR:
                        constantBranch = DSExpressionBranchAtIndex(current, 0);
                        if (constantBranch == NULL) {
                                DSError(M_DS_NULL ": Constant branch is NULL", A_DS_ERROR);
                                break;
                        }
                        constant = DSExpressionConstant(constantBranch);
                        for (i=0; i < DSExpressionNumberOfBranches(current); i++) {
                                if (i == 0 && DSExpressionOperator(current) == '+' && constant == 0.0)
                                        continue;
                                if (i == 0 && DSExpressionOperator(current) == '*' && constant == 1.0)
                                        continue;
                                if (i == 0 && DSExpressionOperator(current) == '*' && constant == 0.0) {
                                        strncat(*string, "0", *length-strlen(*string));
                                        break;
                                }
                                if (i == 0 && DSExpressionOperator(current) == '*' && constant == -1.0) {
                                        strncat(*string, "-", *length-strlen(*string));
                                        continue;
                                }
                                if (DSExpressionOperator(current) == '.') {
                                        sprintf(temp, "\\dot{");
                                        strncat(*string,  temp, *length-strlen(*string));
                                        temp[0] = '\0';
                                }
                                branch = DSExpressionBranchAtIndex(current, i);
                                if (DSExpressionType(branch) == DS_EXPRESSION_TYPE_OPERATOR &&
                                    operatorIsLowerPrecedence(DSExpressionOperator(current), DSExpressionOperator(branch)))
                                        strncat(*string, "(", *length-strlen(*string));
                                expressionToLatexStringInternal(branch, string, length, substitutionDict);
                                if (DSExpressionType(branch) == DS_EXPRESSION_TYPE_OPERATOR &&
                                    operatorIsLowerPrecedence(DSExpressionOperator(current), DSExpressionOperator(branch)))
                                        strncat(*string, ")", *length-strlen(*string));
                                if (DSExpressionOperator(current) == '.') {
                                        sprintf(temp, "}");
                                        strncat(*string,  temp, *length-strlen(*string));
                                        temp[0] = '\0';
                                }
                                if (i < DSExpressionNumberOfBranches(current)-1) {
                                        if (DSExpressionOperator(current) == '+' ||
                                            DSExpressionOperator(current) == '=' ||
                                            DSExpressionOperator(current) == '<' ||
                                            DSExpressionOperator(current) == '>') {
                                                if (DSExpressionOperator(current) == '+' &&
                                                    DSExpressionType(DSExpressionBranchAtIndex(current, i+1)) == DS_EXPRESSION_TYPE_OPERATOR) {
                                                        if (DSExpressionConstant(DSExpressionBranchAtIndex(DSExpressionBranchAtIndex(current, i+1), 0)) < 0) {
                                                                continue;
                                                        }
                                                }
                                                sprintf(temp, " %c ", DSExpressionOperator(current));
                                                strncat(*string,  temp, *length-strlen(*string));
                                                temp[0] = '\0';
                                        } else if (DSExpressionOperator(current) == '^') {
                                                sprintf(temp, "^{");
                                                strncat(*string,  temp, *length-strlen(*string));
                                                temp[0] = '\0';
                                        }
                                } else {
                                        if (DSExpressionOperator(current) == '^') {
                                                sprintf(temp, "}");
                                                strncat(*string,  temp, *length-strlen(*string));
                                                temp[0] = '\0';
                                        }
                                }
                        }
                        break;
                case DS_EXPRESSION_TYPE_FUNCTION:
                        name = DSExpressionVariable(current);
                        open = "(";
                        close = ")";
                        if (strcmp("log", name) == 0) {
                                name = "\\log";
                        } else if (strcmp("log10", name) == 0) {
                                name = "\\log_{10}";
                        } else  if (strcmp("ln", name) == 0) {
                                name = "\\ln";
                        } else  if (strcmp("sin", name) == 0) {
                                name = "\\sin";
                        } else  if (strcmp("cos", name) == 0) {
                                name = "\\cos";
                        } else if (strcmp("sqrt", name) == 0) {
                                name = "\\sqrt";
                                open = "{";
                                close = "}";
                        } else if (strcmp("real", name) == 0) {
                                name = "\\Re";
                        } else if (strcmp("imag", name) == 0) {
                                name = "\\Im";
                        }
                        sprintf(temp, "%s%s", name, open);
                        if (strlen(*string)+strlen(temp) >= *length) {
                                length += DS_EXPRESSION_STRING_INIT_LENGTH;
                                *string = DSSecureRealloc(string, sizeof(char)**length);
                        }
                        strncat(*string, temp, *length-strlen(*string));
                        expressionToLatexStringInternal(DSExpressionBranchAtIndex(current, 0), string, length, substitutionDict);
                        strncat(*string, close, *length-strlen(*string));
                        temp[0] = '\0';
                        break;
                default:
                        break;
        }
        if (strlen(*string)+strlen(temp) >= *length) {
                length += DS_EXPRESSION_STRING_INIT_LENGTH;
                *string = DSSecureRealloc(string, sizeof(char)**length);
        }
        strncat(*string,  temp, *length-strlen(*string));
bail:
        return;
}

extern char * DSExpressionAsLatexString(const DSExpression *expression, const DSDictionary * substitutionDict)
{
        DSUInteger length = DS_EXPRESSION_STRING_INIT_LENGTH;
        char * string = NULL;
        if (expression == NULL) {
                DSError(M_DS_NULL ": Node to print is nil", A_DS_ERROR);
                goto bail;
        }
        if (substitutionDict == NULL) {
                DSError(M_DS_DICTIONARY_NULL, A_DS_ERROR);
                goto bail;
        }
        string = DSSecureCalloc(sizeof(char), length);
        expressionToLatexStringInternal(expression, &string, &length, substitutionDict);
bail:
        return string;
}

extern void DSExpressionPrint(const DSExpression *expression)
{
        char *string = NULL;
        if (expression == NULL) {
                DSError(M_DS_NULL ": Node to print is null", A_DS_ERROR);
                goto bail;
        }
        string = DSExpressionAsString(expression);
        if (string == NULL)
                goto bail;
        if (strlen(string) != 0) {
                if (DSPrintf == NULL) {
                        printf("%s\n", string);
                } else {
                        DSPrintf("%s\n", string);
                }
        } else {
                if (DSPrintf == NULL) {
                        printf("0\n");
                } else {
                        DSPrintf("0\n");
                }
        }
        DSSecureFree(string);
bail:
        return;
}

extern bool DSExpressionIsEqualToExpression(const DSExpression * lhs, const DSExpression *rhs)
{
        bool areEqual = true;
        DSUInteger i;
        if (lhs == rhs) {
                goto bail;
        }
        if (lhs == NULL || rhs == NULL) {
                areEqual = false;
                goto bail;
        }
        if (DSExpressionType(lhs) != DSExpressionType(rhs)) {
                areEqual = false;
                goto bail;
        }
        switch (DSExpressionType(lhs)) {
                case DS_EXPRESSION_TYPE_CONSTANT:
                        if (DSExpressionConstant(lhs) != DSExpressionConstant(rhs))
                                areEqual = false;
                        break;
                case DS_EXPRESSION_TYPE_VARIABLE:
                        if (strcmp(DSExpressionVariable(lhs), DSExpressionVariable(rhs)) != 0)
                                areEqual = false;
                        break;
                case DS_EXPRESSION_TYPE_OPERATOR:
                        if (DSExpressionOperator(lhs) != DSExpressionOperator(rhs)) {
                                areEqual = false;
                                break;
                        }
                        if (DSExpressionNumberOfBranches(lhs) != DSExpressionNumberOfBranches(rhs)) {
                                areEqual = false;
                                break;
                        }
                        for (i = 0; i < DSExpressionNumberOfBranches(lhs); i++) {
                                
                                if (DSExpressionIsEqualToExpression(DSExpressionBranchAtIndex(lhs, i), DSExpressionBranchAtIndex(rhs, i)) == false) {
                                        areEqual = false;
                                        break;
                                }
                        }
                        break;
                case DS_EXPRESSION_TYPE_FUNCTION:
                        if (strcmp(DSExpressionVariable(lhs), DSExpressionVariable(lhs)) != 0) {
                                areEqual = false;
                                break;
                        }
                        areEqual = DSExpressionIsEqualToExpression(DSExpressionBranchAtIndex(lhs, 0), DSExpressionBranchAtIndex(rhs, 0));
                        break;
                default:
                        break;
        }
bail:
        return areEqual;
}

extern bool DSStringIsEqualToStringIgnoringOrder(const char *lhs_string, const char *rhs_string){
    
    DSExpression *ratio=NULL;
    bool areEqual = false;
    char ratio_string[1000];
    DSVariablePool *Xi=NULL;
    
    if (lhs_string == NULL || rhs_string == NULL) {
            areEqual = false;
            goto bail;
    }
    
    if (strcmp(lhs_string, rhs_string) == 0) {
            areEqual = true;
            goto bail;
    }
    
    sprintf(ratio_string, "%s/(%s)", lhs_string, rhs_string);
    ratio = DSExpressionByParsingStringSimplifyExpressionAsExpression(ratio_string);
    Xi = DSExpressionVariablesInExpression(ratio);
    
    if (DSVariablePoolNumberOfVariables(Xi) == 0)
        areEqual = true;
    else
        areEqual = false;
    
    if (ratio != NULL)
        DSSecureFree(ratio);
    if (Xi != NULL)
        DSVariablePoolFree(Xi);
    
bail:
    return areEqual;
    
}


extern bool DSExpressionIsEqualToExpressionIgnoringOrder(const DSExpression * lhs, const DSExpression *rhs){
    
    bool areEqual = false;
    char * lhs_string = NULL, *rhs_string = NULL;
    
    if (lhs == rhs) {
            areEqual = true;
            goto bail;
    }
    if (lhs == NULL || rhs == NULL) {
            areEqual = false;
            goto bail;
    }
    
    lhs_string = DSExpressionAsString(lhs);
    rhs_string = DSExpressionAsString(rhs);
    
    areEqual = DSStringIsEqualToStringIgnoringOrder(lhs_string, rhs_string);
    
    if (lhs_string !=NULL)
        DSSecureFree(lhs_string);
    if (rhs_string !=NULL)
        DSSecureFree(rhs_string);

    
bail:
    return areEqual;
    
}

extern DSExpression * DSExpressionByReplacingSubExpression(const DSExpression * expression, const DSExpression * target, const DSExpression * substitute)
{
        DSExpression * newExpression = NULL;
        DSExpression * new, * old;
        DSUInteger i;
        DSVariablePool * temp = DSVariablePoolAlloc();
        if (expression == NULL || target == NULL) {
                DSError(M_DS_NULL, A_DS_ERROR);
                goto bail;
        }
        switch (DSExpressionType(expression)) {
                case DS_EXPRESSION_TYPE_CONSTANT:
                case DS_EXPRESSION_TYPE_VARIABLE:
                        if (DSExpressionIsEqualToExpression(expression, target) == true) {
                                newExpression = DSExpressionCopy(substitute);
                        } else {
                                newExpression = DSExpressionCopy(expression);
                        }
                        break;
                case DS_EXPRESSION_TYPE_OPERATOR:
                case DS_EXPRESSION_TYPE_FUNCTION:
                        if (DSExpressionIsEqualToExpression(expression, target) == true) {
                                newExpression = DSExpressionCopy(substitute);
                                break;
                        }
                        newExpression = DSExpressionCopy(expression);
                        for (i = 0; i < DSExpressionNumberOfBranches(expression); i++) {
                                old = DSExpressionByReplacingSubExpression(DSExpressionBranchAtIndex(expression, i), target, substitute);
                                temp = DSExpressionVariablesInExpression(old);
                                if (DSVariablePoolNumberOfVariables(temp) == 0) {
                                        new = dsExpressionAllocWithConstant(DSExpressionEvaluateComplexWithVariablePool(old, temp));
                                        DSExpressionFree(old);
                                } else {
                                        new = old;
                                }
//                                new = dsExpressionCompressConstantVariableNode(old, temp);
                                if (DSExpressionIsEqualToExpression(DSExpressionBranchAtIndex(expression, i), new) == false) {
                                        DSExpressionFree(DSExpressionBranchAtIndex(newExpression, i));
                                        newExpression->branches[i] = new;
                                } else {
                                        DSExpressionFree(new);
                                }
                        }
                        break;
                default:
                        break;
        }
bail:
        DSVariablePoolFree(temp);
        return newExpression;
}

extern void DSExpressionRecastBranch(DSExpression *** array, DSExpression * expression, DSUInteger index, const char * prefix, DSUInteger *count) {
        DSUInteger i;
        DSExpression *new, *subs, *temp;
        char name[100];
        sprintf(name, "%s%i", prefix, *count);
        subs = dsExpressionAllocWithVariableName(name);
        for (i = index; i < *count; i++) {
                temp = DSExpressionByReplacingSubExpression((*array)[i], expression, subs);
                DSExpressionFree((*array)[i]);
                (*array)[i] = temp;
        }
        DSSecureFree(subs);
        *count += 1;
        *array = DSSecureRealloc(*array, sizeof(DSExpression **)**count);
        sprintf(name, "%s=replace", name);
        new = DSExpressionByParsingString(name);
        subs = dsExpressionAllocWithVariableName("replace");
        (*array)[*count-1] = DSExpressionByReplacingSubExpression(new, subs, expression);
        DSExpressionFree(new);
        DSExpressionFree(subs);
bail:
        return;
}

extern bool DSExpressionRecastExpression(DSExpression *** array, DSExpression * expression, DSUInteger index, const char * prefix, DSUInteger * count) {
        DSUInteger i;
        DSExpression *branch, *temp;
        bool changed = false;
        if (array == NULL) {
                DSError(M_DS_NULL, A_DS_ERROR);
                goto bail;
        }
        if (DSExpressionType(expression) != DS_EXPRESSION_TYPE_OPERATOR) {
                goto bail;
        }
        for (i = 0; i < DSExpressionNumberOfBranches(expression); i++) {
                if (changed == true)
                        break;
                branch = DSExpressionBranchAtIndex(expression, i);
                if (DSExpressionType(branch) != DS_EXPRESSION_TYPE_OPERATOR)
                        continue;
                if (operatorIsLowerPrecedence(DSExpressionOperator(expression), DSExpressionOperator(branch)) == true) {
                        temp = DSExpressionCopy(branch);
                        DSExpressionRecastBranch(array, temp, index, prefix, count);
                        DSSecureFree(temp);
                        changed = true;
                        break;
                }
                changed = DSExpressionRecastExpression(array, branch, index, prefix, count);
        }
bail:
        return changed;
}


extern DSExpression ** DSExpressionRecastSystemEquations(const DSExpression ** expressionArray, DSUInteger * numberOfEquations, const char * prefix) {
        DSExpression ** array = NULL, *temp;
        DSUInteger i;
        bool changed;
        if (expressionArray == NULL) {
                DSError(M_DS_NULL, A_DS_ERROR);
                goto bail;
        }
        if (*numberOfEquations == 0) {
                DSError(M_DS_WRONG, A_DS_ERROR);
                goto bail;
        }
        array = DSSecureCalloc(sizeof(DSExpression *), *numberOfEquations);
        for (i = 0; i < *numberOfEquations; i++) {
                array[i] = DSExpressionCopy(expressionArray[i]);
        }
        for (i = 0; i < *numberOfEquations; i++) {
                temp = NULL;
                changed = true;
                while (changed) {
                        changed = DSExpressionRecastExpression(&array, array[i], i, prefix, numberOfEquations);
                }
        }
bail:
        return array;
}



extern DSExpression * DSExpressionFromPowerlawInMatrixForm(const DSUInteger row, const DSMatrix * Kd, const DSVariablePool * Xd, const DSMatrix * Ki, const DSVariablePool *Xi, const DSMatrix * C)
{
        DSUInteger i;
        DSExpression * expression = NULL;
        char * name, *string = NULL, *temp;
        double value;
        string = DSSecureCalloc(sizeof(char), 100);
        temp = string;
        asprintf(&string, "%lf", DSMatrixDoubleValue(C, row, 0));
        if (temp != string) {
                DSSecureFree(temp);
                temp = string;
        }
        for (i = 0; i < DSVariablePoolNumberOfVariables(Xd); i++) {
                value = DSMatrixDoubleValue(Kd, row, i);
                if (value == 0.0)
                        continue;
                name = DSVariableName(DSVariablePoolVariableAtIndex(Xd, i));
                if (value == 1.) {
                        asprintf(&string, "%s*%s", string, name);
                } else {
                        asprintf(&string, "%s*%s^%lf", string, name, value);
                }
                if (temp != string) {
                        DSSecureFree(temp);
                        temp = string;
                }
        }
        for (i = 0; i < DSVariablePoolNumberOfVariables(Xi); i++) {
                value = DSMatrixDoubleValue(Ki, row, i);
                if (value == 0.0f)
                        continue;
                name = DSVariableName(DSVariablePoolVariableAtIndex(Xi, i));
                if (value == 1.) {
                        asprintf(&string, "%s*%s", string, name);
                } else {
                        asprintf(&string, "%s*%s^%lf", string, name, value);
                }
                if (temp != string) {
                        DSSecureFree(temp);
                        temp = string;
                }
        }
        expression = DSExpressionByParsingString(string);
        DSSecureFree(string);
bail:
        return expression;
}

static void processbranches(DSExpression * expression,
                            DSMatrix * Ki,
                            DSMatrix *C,
                            DSVariablePool * Xi,
                            double previous_exponent){
    
    if (expression == NULL || Ki == NULL || C == NULL || Xi == NULL)
        goto bail;
    
    DSUInteger i, ii, index;
    DSExpression *expression_i = NULL;
    DSVariablePool *tempPool = NULL, *Xd = NULL;
    double value, previous_exponent_i = 1.0;
    DSMatrix *Kd = NULL;
    char **ptr;

    Xd = DSVariablePoolAlloc();
    
//    printf("using processbranches. Received %s in form of matrices. Merging %s with the exponent %f \n",
//           DSExpressionAsString(DSExpressionFromPowerlawInMatrixForm(0, Kd, Xd, Ki, Xi, C)),
//           DSExpressionAsString(expression),
//           previous_exponent);
//    printf("The expression has %u branches \n", DSExpressionNumberOfBranches(expression));
//    printf("The operator is: %c \n", DSExpressionOperator(expression));

    if (DSExpressionOperator(expression) == '^'){
        previous_exponent_i = strtod(DSExpressionAsString(DSExpressionBranchAtIndex(expression, 1)),
                                   ptr);
        processbranches(DSExpressionBranchAtIndex(expression, 0), Ki, C, Xi, previous_exponent*previous_exponent_i);
        goto bail;
    }
    
    
    if (DSExpressionNumberOfBranches(expression) == 0 && DSExpressionOperator(expression) == '?'){
        tempPool = DSExpressionVariablesInExpression(expression);
//        DSMatrixSetDoubleValue(C, 0, 0, 1);
        if (DSVariablePoolNumberOfVariables(tempPool) == 1){
//            printf("process branches for expressions withoutbranches and 1 variable! the constant is: %f \n", expression->node.constant);
            index = DSVariablePoolIndexOfVariableWithName(Xi,
                                                          DSVariablePoolAllVariableNames(tempPool)[0]);
            value = DSMatrixDoubleValue(Ki, 0, index) + 1.0*previous_exponent;
            DSMatrixSetDoubleValue(Ki, 0, index, value);
        }
        if (tempPool != NULL){
            DSVariablePoolFree(tempPool);
            tempPool = NULL;
        }
        
    }
    
    
    for(i=0; i<DSExpressionNumberOfBranches(expression); i++){
        
            expression_i = DSExpressionBranchAtIndex(expression, i);
            tempPool = DSExpressionVariablesInExpression(expression_i);
        
//            printf("Branch %u: %s \n", i, DSExpressionAsString(expression_i));
            if (DSVariablePoolNumberOfVariables(tempPool) == 0){
//                    printf("branch %u has 0 variables \n", i);
                    if (DSExpressionOperator(expression) == '*')
                        DSMatrixSetDoubleValue(C, 0, 0, expression_i->node.constant);
                    if (tempPool != NULL){
                        DSVariablePoolFree(tempPool);
                        tempPool = NULL;
                    }
                    continue;
            }
        
            if (DSExpressionNumberOfBranches(expression_i) == 0){
                if (DSVariablePoolNumberOfVariables(tempPool) > 1){
                    previous_exponent_i = strtod(DSExpressionAsString(DSExpressionBranchAtIndex(expression_i, 1)),
                                               ptr);
                    processbranches(DSExpressionBranchAtIndex(expression_i, 0), Ki, C, Xi, previous_exponent_i);
                }else{
                    index = DSVariablePoolIndexOfVariableWithName(Xi,
                                                                  DSVariablePoolAllVariableNames(tempPool)[0]);
                    value = DSMatrixDoubleValue(Ki, 0, index) + 1.0*previous_exponent;
                    DSMatrixSetDoubleValue(Ki, 0, index, value);
                }
            }else{
                processbranches(expression_i, Ki, C, Xi, previous_exponent);
            }
            if (tempPool != NULL)
                DSVariablePoolFree(tempPool);
    }
//test1:
//    printf("using processbranches. Generated %s \n",
//           DSExpressionAsString(DSExpressionFromPowerlawInMatrixForm(0, Kd, Xd, Ki, Xi, C)));
    
bail:
    return;
    
}


extern DSExpression * DSExpressionSimplifyExpressionAsExpression(const DSExpression *expression){
    
    if (expression == NULL)
        goto bail;
    
//    printf("**************************************** start ****************************************  \n");
//    printf("Reporting from DSExpressionSimplifyExpressionAsExpression. The expression I got is: %s \n ", DSExpressionAsString(expression));
    
    DSMatrix * Ki = NULL, *Kd =NULL, *C = NULL;
    DSVariablePool *Xi = NULL, *Xd = NULL, *tempPool = NULL, *tempPool_ii = NULL;
    DSUInteger length = DS_EXPRESSION_STRING_INIT_LENGTH, i, index, ii;
    DSExpression *expression_i, *expression_ii, *expression_aux;
    double value, previous_exponent = 1.0;
    char **ptr;
    
    Xd = DSVariablePoolAlloc();
    Xi = DSExpressionVariablesInExpression(expression);
    if (DSVariablePoolNumberOfVariables(Xi) !=0)
        Ki = DSMatrixCalloc(1, DSVariablePoolNumberOfVariables(Xi));
    else
        Ki = DSMatrixCalloc(1, 1);
    C = DSMatrixCalloc(1, 1);
    
//    printf("The expression has %u branches \n", DSExpressionNumberOfBranches(expression));
//    printf("The operator is: %c \n", DSExpressionOperator(expression));
    
    if (DSExpressionNumberOfBranches(expression) == 0){
        if (DSVariablePoolNumberOfVariables(Xi) != 0)
            DSMatrixSetDoubleValue(C, 0, 0, 1.0);
        else
            DSMatrixSetDoubleValue(C, 0, 0, expression->node.constant);
        DSMatrixSetDoubleValue(Ki, 0, 0, 1.0);
    }

    if (DSExpressionOperator(expression) == '^'){
        previous_exponent = strtod(DSExpressionAsString(DSExpressionBranchAtIndex(expression, 1)),
                                   ptr);
        processbranches(DSExpressionBranchAtIndex(expression, 0), Ki, C, Xi, previous_exponent);
        goto bail0;
    }
    
    for(i=0; i<DSExpressionNumberOfBranches(expression); i++){
        
            expression_i = DSExpressionBranchAtIndex(expression, i);
            tempPool = DSExpressionVariablesInExpression(expression_i);
//            printf("Branch %u: %s \n", i, DSExpressionAsString(expression_i));
        
            if (DSVariablePoolNumberOfVariables(tempPool) == 0){
                    
                if (DSExpressionOperator(expression) == '*'){
//                    printf("branch %u has 0 variables. The constant is %f \n", i, expression_i->node.constant);
                    DSMatrixSetDoubleValue(C, 0, 0, expression_i->node.constant);
                }else{
                    printf("Generate an error! \n");
                }
                    if (tempPool != NULL){
                        DSVariablePoolFree(tempPool);
                        tempPool = NULL;
                    }
                    continue;
            }
            if (DSExpressionNumberOfBranches(expression_i) == 0){
                if (DSVariablePoolNumberOfVariables(tempPool) > 1){
//                    printf("branch %u has more than 1 variables \n", i);
                    previous_exponent = strtod(DSExpressionAsString(DSExpressionBranchAtIndex(expression_i, 1)),
                                               ptr);
                    processbranches(DSExpressionBranchAtIndex(expression_i, 0), Ki, C, Xi, previous_exponent);
                } else {
//                    printf("branch %u has 1 variable \n", i);
                    index = DSVariablePoolIndexOfVariableWithName(Xi,
                                                                  DSVariablePoolAllVariableNames(tempPool)[0]);
                    value = DSMatrixDoubleValue(Ki, 0, index) + 1.0;
                    DSMatrixSetDoubleValue(Ki, 0, index, value);
                }
            }else{
                    processbranches(expression_i, Ki, C, Xi, previous_exponent);
                }
        
            if (tempPool != NULL)
                DSVariablePoolFree(tempPool);
    }
bail0:
    expression_aux = DSExpressionFromPowerlawInMatrixForm(0, Kd, Xd, Ki, Xi, C);
    
//    printf("Reporting from DSExpressionSimplifyExpressionAsExpression. The expression I generated is: %s \n ", DSExpressionAsString(expression_aux));
//    printf("**************************************** end ****************************************  \n");
    
    if (Xd != NULL)
        DSVariablePoolFree(Xd);
    if (Xi != NULL)
        DSVariablePoolFree(Xi);
    if (Ki != NULL)
        DSMatrixFree(Ki);
    if (C != NULL)
        DSMatrixFree(C);
        
bail:
    return expression_aux;
    
}

extern char * DSExpressionSimplifyExpressionAsString(const DSExpression *expression){
    
    if (expression == NULL)
        goto bail;
    
    char * string = NULL;
    DSExpression *simplified_expression = NULL;
    
    simplified_expression = DSExpressionSimplifyExpressionAsExpression(expression);
    string = DSExpressionAsString(simplified_expression);
    
    if (simplified_expression != NULL)
        DSExpressionFree(simplified_expression);
    
bail:
    return string;
    
}

extern DSExpression * DSExpressionByParsingStringSimplifyExpressionAsExpression(const char * string){
    
    if (string == NULL)
        goto bail;
    
    DSExpression * simplified_expression = NULL, *string_as_expression=NULL;
    
    string_as_expression = DSExpressionByParsingString(string);
    simplified_expression = DSExpressionSimplifyExpressionAsExpression(string_as_expression);
    
    if (string_as_expression != NULL)
        DSExpressionFree(string_as_expression);
            
bail:
    return simplified_expression;
}

extern const char * DSExpressionByParsingStringSimplifyExpressionAsString(const char * string){
    
    if (string == NULL)
        goto bail;
    
    DSExpression * simplified_expression = NULL, *string_as_expression=NULL;
    char * simplified_string = NULL;
    
    string_as_expression = DSExpressionByParsingString(string);
    simplified_expression = DSExpressionSimplifyExpressionAsExpression(string_as_expression);
    simplified_string = DSExpressionAsString(simplified_expression);
    
    if (string_as_expression != NULL)
        DSExpressionFree(string_as_expression);
    if (simplified_expression != NULL)
        DSExpressionFree(simplified_expression);
            
bail:
    return simplified_string;
    
    
}


extern DSExpression * DSExpressionFromLogPowerlawInMatrixForm(const DSUInteger row, const DSMatrix * Kd, const DSVariablePool * Xd, const DSMatrix * Ki, const DSVariablePool *Xi, const DSMatrix * C)
{
        DSUInteger i;
        DSExpression * expression = NULL;
        char * name, *string = NULL, *temp;
        string = DSSecureCalloc(sizeof(char), 100);
        temp = string;
        asprintf(&string, "%lf", log10(DSMatrixDoubleValue(C, row, 0)));
        if (temp != string) {
                DSSecureFree(temp);
                temp = string;
        }
        for (i = 0; i < DSVariablePoolNumberOfVariables(Xd); i++) {
                if (DSMatrixDoubleValue(Kd, row, i) == 0.0f)
                        continue;
                name = DSVariableName(DSVariablePoolVariableAtIndex(Xd, i));
                asprintf(&string, "%s+%lf*%s", string, DSMatrixDoubleValue(Kd, row, i), name);
                if (temp != string) {
                        DSSecureFree(temp);
                        temp = string;
                }
        }
        for (i = 0; i < DSVariablePoolNumberOfVariables(Xi); i++) {
                if (DSMatrixDoubleValue(Ki, row, i) == 0.0f)
                        continue;
                name = DSVariableName(DSVariablePoolVariableAtIndex(Xi, i));
                asprintf(&string, "%s+%lf*%s", string, DSMatrixDoubleValue(Ki, row, i), name);
                if (temp != string) {
                        DSSecureFree(temp);
                        temp = string;
                }
        }
        expression = DSExpressionByParsingString(string);
        DSSecureFree(string);
bail:
        return expression;
}


