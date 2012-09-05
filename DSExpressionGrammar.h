/* This file was automatically generated.  Do not edit! */
#define DSExpressionParserTOKENTYPE DSExpression *
#define DSExpressionParserARG_PDECL ,void *parsed
void DSExpressionParser(void *yyp,int yymajor,DSExpressionParserTOKENTYPE yyminor DSExpressionParserARG_PDECL);
#if defined(YYTRACKMAXSTACKDEPTH)
int DSExpressionParserStackPeak(void *p);
#endif
void DSExpressionParserFree(void *p,void(*freeProc)(void *));
void *DSExpressionParserAlloc(void *(*mallocProc)(size_t));
#if !defined(NDEBUG)
void DSExpressionParserTrace(FILE *TraceFILE,char *zTracePrompt);
#endif
#define DSExpressionParserARG_STORE yypParser->parsed = parsed
#define DSExpressionParserARG_FETCH void *parsed = yypParser->parsed
#define DSExpressionParserARG_SDECL void *parsed;
#define TOKEN_EXPRESSION_RPAREN                         10
#define TOKEN_EXPRESSION_LPAREN                          9
#define TOKEN_EXPRESSION_VALUE                           8
#define TOKEN_EXPRESSION_POWER                           7
#define TOKEN_EXPRESSION_NOT                             6
#define TOKEN_EXPRESSION_TIMES                           5
#define TOKEN_EXPRESSION_DIVIDE                          4
#define TOKEN_EXPRESSION_MINUS                           3
#define TOKEN_EXPRESSION_PLUS                            2
#define TOKEN_EXPRESSION_ID                              1
#define INTERFACE 0
DSExpression *dsExpressionAllocWithVariableName(const char *name);
DSExpression *dsExpressionAllocWithConstant(const double value);
DSExpression *dsExpressionAllocWithOperator(const char op_code);
