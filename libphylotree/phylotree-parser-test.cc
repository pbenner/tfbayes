/* Copyright (C) 2012 Philipp Benner
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <stdio.h>

#include <iostream>

#include <phylotree.hh>
#include <phylotree-parsetree.hh>

using namespace std;

extern char* yytext;
extern size_t line_count;
extern pt_parsetree_t* root;

int yyparse(void);

int yyerror(const char *msg) {
        printf("%s at line %ld near `%s'\n", msg, line_count, yytext);
        return 0;
}

int main(void) {

        yyparse();
        cout << root
             << endl;

        pt_node_t* pt_root = root->convert();

        delete(root);
        cout << pt_root << endl;

        return 0.0;
}
