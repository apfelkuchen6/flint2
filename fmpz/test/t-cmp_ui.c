/*
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("cmp_ui....");
    fflush(stdout);

    

    /* Compare with fmpz_cmp */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        ulong n;
        int lhs, rhs;

        fmpz_init(a);
        fmpz_init(b);

        fmpz_randtest(a, state, 200);

        n = n_randtest(state);
        fmpz_set_ui(b, n);

        lhs = fmpz_cmp(a, b);
        rhs = fmpz_cmp_ui(a, n);

        result = (lhs < 0) ? (rhs < 0) : ((lhs > 0) ? (rhs > 0) : (rhs == 0));
        if (result == 0)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), fmpz_print(a), flint_printf("\n");
            flint_printf("b = "), fmpz_print(b), flint_printf("\n");
            flint_printf("n = %wu\n", n);
            flint_printf("cmp(a, b) = %d\n", fmpz_cmp(a, b));
            flint_printf("cmp_ui(a, n) = %d\n", fmpz_cmp_ui(a, n));
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
