%{
#include <math.h>
#include <stdio.h>
#define YYSTYPE double
int yylex (void);
void yyerror (const char *);
%}

%token LPAREN RPAREN ID NUM STR

%%
sexpr: atom                  {printf("matched sexpr\n");}
    | list
    ;
list: LPAREN members RPAREN  {printf("matched list\n");}
    | LPAREN RPAREN          {printf("matched empty list\n");}
    ;
members: sexpr               {printf("members 1\n");}
    | sexpr members          {printf("members 2\n");}
    ;
atom: ID                     {printf("ID\n");}
    | NUM                    {printf("NUM\n");}
    | STR                    {printf("STR\n");}
    ;
%%
