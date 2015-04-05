/* Copyright (C) 2015 Philipp Benner
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

#include <iostream>
#include <cstdlib>

#include <tfbayes/config/partition-lexer.hh>

using namespace std;

int
main(int argc, char *argv[])
{
        if (argc != 2) {
                cerr << "Usage: partition-parser-test INPUT_FILE"
                     << endl;
                exit(EXIT_FAILURE);
        }
        string filename(argv[1]);

        dpm_partition_list_t partition_list =
                parse_partition_list(filename);

        for (size_t i = 0; i < partition_list.size(); i++) {
                cout << partition_list[i]
                     << endl;
        }

        parse_partition_list_str("baseline-default:11:{(0, 334), (0, 230), (0, 230)!}, iohf");

        return 0;
}
