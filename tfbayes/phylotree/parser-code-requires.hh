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

#ifndef __TFBAYES_PHYLOTREE_PARSER_CODE_HH__
#define __TFBAYES_PHYLOTREE_PARSER_CODE_HH__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>

#include <tfbayes/phylotree/parsetree.hh>

#define YYSTYPE pt_parsetree_t *
#define YYINITDEPTH 100000
#define YYMAXDEPTH  100000

typedef void* yyscan_t;
typedef struct {
	pt_parsetree_t* pt_parsetree;
	yyscan_t scanner;
} context_t;

#endif /* __TFBAYES_PHYLOTREE_PARSER_CODE_HH__ */
