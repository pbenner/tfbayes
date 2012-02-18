%{
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <stdio.h>
#include <phylotree.hh>

#define YYSTYPE pt_node_t *

int yylex (void);
void yyerror (const char *);

pt_root_t* root;
extern char yytext[];

%}

%token LPAREN RPAREN ROOT NODE LEAF ID NUM STR

%%
root: LPAREN ROOT node node RPAREN { root = new pt_root_t(-1, $3, $4, ""); printf("now at root\n"); }
    | LPAREN LEAF NUM RPAREN       { root = new pt_root_t( 1, NULL, NULL, ""); }
    ;
node: LPAREN NODE node node RPAREN { $$ = new pt_node_t(-1, 0.0, $3, $4, ""); printf("now at node\n"); }
    | LPAREN LEAF NUM RPAREN       { $$ = new pt_leaf_t(-1, 0.0, "burp"); printf("now at leaf\n"); }
%%
