/* Copyright (C) 2013 Philipp Benner
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

// general options
////////////////////////////////////////////////////////////////////////////////

%pure-parser
%locations
%defines
%error-verbose
%lex-param   {context_t* context}
%parse-param {context_t* context}

// lexer definitions
////////////////////////////////////////////////////////////////////////////////

%code requires {
#include <tfbayes/config/partition-code-requires.hh>
}

%code {
#include <tfbayes/config/partition-code.hh>
}

// token definitions
////////////////////////////////////////////////////////////////////////////////

%token COLON COMMA DASH EXCLAMATION INTEGER LPAREN NEWLINE RPAREN LCPAREN RCPAREN SEMICOLON STRING

// grammar
////////////////////////////////////////////////////////////////////////////////
%%
partition_list:
      partition_list NEWLINE partition
      {
        if ($3) {
          context->partition_list.push_back($3->top().partition);
          delete($3);
        }
      }
    | partition
      {
        if ($1) {
          context->partition_list.push_back($1->top().partition);
          delete($1);
        }
      }
    ;
partition:
      partition COMMA subset
      {
        $$ = $1;
        $$->top().partition.push_back($3->top().subset);
        delete($3);
      }
    | subset
      {
        node_t node;
        node.partition.push_back($1->top().subset);
        $$ = new std::stack<node_t>();
        $$->push(node);
        delete($1);
      }
    | DASH
      {
        $$ = new std::stack<node_t>();
        $$->push(node_t());
      }
    | // skip empty lines
      { $$ = NULL; }
    ;
subset:
      string COLON integer COLON LCPAREN range_list RCPAREN
      {
        std::string baseline_tag = $1->top().string;
	size_t length            = $3->top().integer;
	dpm_subset_t subset      = model_id_t({baseline_tag, length});
        while (!$6->empty()) {
	  $6->top().range.length() = length;
          subset.insert($6->top().range);
          $6->pop();
        }
	delete($1);
	delete($3);
	delete($6);
	node_t node;
	node.subset = subset;
        $$ = new std::stack<node_t>();
        $$->push(node);
      }
    ;
range_list:
      range_list COMMA range
      {
        $$ = $1;
        $$->push($3->top());
	delete($3);
      }
    |
      range
      {
        $$ = $1;
      }
    ;
range:
      LPAREN integer COMMA integer RPAREN EXCLAMATION
      {
        size_t i = $2->top().integer;
        size_t j = $4->top().integer;
	delete($2);
	delete($4);
        node_t node;
        node.range = range_t(index_t(i, j), 0, true);
        $$ = new std::stack<node_t>();
        $$->push(node);
      }
    |
      LPAREN integer COMMA integer RPAREN
      {
        size_t i = $2->top().integer;
        size_t j = $4->top().integer;
	delete($2);
	delete($4);
        node_t node;
        node.range = range_t(index_t(i, j), 0, false);
        $$ = new std::stack<node_t>();
        $$->push(node);
      }
string:
      STRING
      {
        node_t node;
        node.string = std::string(yyget_text(context->scanner));
        $$ = new std::stack<node_t>();
        $$->push(node);
      }
integer:
      INTEGER
      {
        node_t node;
        node.integer = std::atoi(yyget_text(context->scanner));
        $$ = new std::stack<node_t>();
        $$->push(node);
      }
%%
