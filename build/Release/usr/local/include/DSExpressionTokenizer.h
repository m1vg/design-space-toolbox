/**
 * \file DSTokenizer.h
 * \brief
 *
 * \details
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

#include "DSTypes.h"
#include "DSErrors.h"
#include "DSMemoryManager.h"

#ifndef __DS_EXPRESSION_TOKENIZER__
#define __DS_EXPRESSION_TOKENIZER__

#include "DSExpressionGrammar.h"

#define DS_EXPRESSION_TOKEN_START      0                 //!< Token indicating the start of a tokenization.
#define DS_EXPRESSION_TOKEN_ID         TOKEN_EXPRESSION_ID  //!< Token indicating a variable identifier.
#define DS_EXPRESSION_TOKEN_VALUE      TOKEN_EXPRESSION_VALUE       //!< Token indicating a numerical value.
#define DS_EXPRESSION_TOKEN_EQUALS     TOKEN_EXPRESSION_EQUALS
#define DS_EXPRESSION_TOKEN_LT         TOKEN_EXPRESSION_LT
#define DS_EXPRESSION_TOKEN_MT         TOKEN_EXPRESSION_MT
#define DS_EXPRESSION_TOKEN_PLUS       TOKEN_EXPRESSION_PLUS
#define DS_EXPRESSION_TOKEN_MINUS      TOKEN_EXPRESSION_MINUS
#define DS_EXPRESSION_TOKEN_TIMES      TOKEN_EXPRESSION_TIMES
#define DS_EXPRESSION_TOKEN_DIVIDE     TOKEN_EXPRESSION_DIVIDE
#define DS_EXPRESSION_TOKEN_PRIME      TOKEN_EXPRESSION_PRIME
#define DS_EXPRESSION_TOKEN_NOT        TOKEN_EXPRESSION_NOT
#define DS_EXPRESSION_TOKEN_POWER      TOKEN_EXPRESSION_POWER
#define DS_EXPRESSION_TOKEN_LPAREN     TOKEN_EXPRESSION_LPAREN
#define DS_EXPRESSION_TOKEN_RPAREN     TOKEN_EXPRESSION_RPAREN




#define DSExpressionTokenNext(x)         ((x)->next)
#define DSExpressionTokenData(x)         ((x)->data)
#define DSExpressionTokenType(x)         ((x)->type)

#define DSExpressionTokenSetNext(x, y)   ((x)->next = (y))
#define DSExpressionTokenSetData(x, y)   ((x)->data = (y))
#define DSExpressionTokenSetType(x, y)   ((x)->type = (y))

/**
 * \brief A data structure representing a token used when parsing strings for
 * variable pools.
 *
 * \details This structures follows the convention used with the struct
 * variable_token and struct matrix_token, representing an ordered list
 * of tokens, as found by the tokenizers generated by the lex program.
 *
 * \see DSExpressionTokenizer()
 */
struct expression_token
{
        int type;                       //!< The current token code.
        union {
                char * name;            //!< Used for storing the name of a variable.
                double value;           //!< Used for storing the value of a constant.
        } data;                         //!< Union for holding either the name of a variable, or the value of a constant.
        struct expression_token *next;  //!< A pointer to the next token in the list.
};

/**
 * \brief Structure used when parsing a mathematical expression.
 *
 * \details This structure is used to parse a mathematical expression, it
 * holds (1) the root of the abstract syntax tree and a flag indicating if
 * any syntax errors were found.
 *
 */
typedef struct {
        DSExpression *root; //!< The pointer to the DSExpression representing the root of the syntax tree.
        bool wasSuccesful;  //!< Indicates if the parsing was succesful.
} parse_expression_s;

extern void DSExpressionAddBranch(DSExpression *expression, DSExpression *branch);

extern struct expression_token * DSExpressionTokenAlloc();
extern void DSExpressionTokenFree(struct expression_token *root);

extern void DSExpressionTokenSetString(struct expression_token *root, char *string);
extern void DSExpressionTokenSetDouble(struct expression_token *root, double value);

extern char * DSExpressionTokenString(struct expression_token *root);
extern double DSExpressionTokenDouble(struct expression_token *root);

extern struct expression_token * DSExpressionTokenizeString(const char *string);

#endif
