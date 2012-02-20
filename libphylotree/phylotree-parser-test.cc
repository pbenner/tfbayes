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
#include <phylotree-parser.hh>

using namespace std;

int main(void) {

        yyparse();

        pt_node_t* pt_root = pt_parsetree->convert();

        cout << "Parsetree:" << endl
             << pt_parsetree << endl;

        cout << "Phylogenetic tree:" << endl
             << pt_root              << endl;

        pt_parsetree->destroy();
        pt_root->destroy();

        return 0.0;
}
