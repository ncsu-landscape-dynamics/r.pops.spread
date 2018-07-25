/*
 * SOD model
 *
 * Copyright (C) 2015-2017 by the authors.
 *
 * Authors: Zexi Chen (zchen22 ncsu edu)
 *          Vaclav Petras (wenzeslaus gmail com)
 *          Anna Petrasova (kratochanna gmail com)
 *
 * The code contained herein is licensed under the GNU General Public
 * License. You may obtain a copy of the GNU General Public License
 * Version 2 or later at the following locations:
 *
 * http://www.opensource.org/licenses/gpl-license.html
 * http://www.gnu.org/copyleft/gpl.html
 */

#include "popss/date.h"

extern "C" {
#include <grass/gis.h>
}

#include <map>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <string>

using std::string;
using std::cout;
using std::cerr;
using std::endl;

#ifdef SOD_TEST

int main(int argc, char *argv[])
{
    G_gisinit(argv[0]);

    struct GModule *module = G_define_module();

    G_add_keyword("test");
    G_add_keyword("date");
    G_add_keyword("raster");
    module->description = "Test of date computations for SOD model";

    if (G_parser(argc, argv))
        exit(EXIT_FAILURE);

    Date start(2018, 1, 1);
    Date end(2020, 12, 31);

    while (start < end) {
        cout << start << endl;
        if (start.isLastMonthOfYear()) {
            cout << "End of year" << endl;
        }
        start.increasedByMonth();
    }

    return 0;
}

#endif  // SOD_TEST
