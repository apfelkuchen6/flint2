/*
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpf_mat.h"
#include "ulong_extras.h"

int
main(void)
{
    int i;
    FLINT_TEST_INIT(state);

    flint_printf("qr....");
    fflush(stdout);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        mpf_t dot, tmp;
        int j, k, l;
        mpf_mat_t A, Q, R, B;
        flint_bitcnt_t prec;

        slong m, n;

        m = n_randint(state, 10);
        n = n_randint(state, 10);
        prec = n_randint(state, 200) + 3;

        mpf_mat_init(A, m, n, prec);
        mpf_mat_init(Q, m, n, prec);
        mpf_mat_init(R, n, n, prec);
        mpf_mat_init(B, m, n, prec);
        mpf_init2(dot, prec);
        mpf_init2(tmp, prec);

        mpf_mat_randtest(A, state, prec);

        mpf_mat_qr(Q, R, A);

        mpf_mat_mul(B, Q, R);

        flint_mpf_set_ui(tmp, 1);
        mpf_div_2exp(tmp, tmp, prec - 2);
        for (j = 0; j < m; j++)
        {
            for (k = 0; k < n; k++)
            {
                mpf_sub(dot, mpf_mat_entry(A, j, k), mpf_mat_entry(B, j, k));
                mpf_abs(dot, dot);
                if (mpf_cmp(dot, tmp) > 0)
                {
                    flint_printf("FAIL:\n");
                    flint_printf("A:\n");
                    mpf_mat_print(A);
                    flint_printf("Q:\n");
                    mpf_mat_print(Q);
                    flint_printf("R:\n");
                    mpf_mat_print(R);
                    flint_printf("B:\n");
                    mpf_mat_print(B);
                    mpf_out_str(stdout, 10, 0, dot);
                    flint_printf("\n");
                    flint_printf("%d\n", prec);
                    fflush(stdout);
                    flint_abort();
                }
            }
        }

        for (j = 0; j < n; j++)
        {
            mpf_t norm;
            mpf_init2(norm, prec);
            for (l = 0; l < m; l++)
            {
                mpf_mul(tmp, mpf_mat_entry(Q, l, j), mpf_mat_entry(Q, l, j));
                mpf_add(norm, norm, tmp);
            }

            mpf_sub_ui(dot, norm, 1);
            mpf_abs(dot, dot);
            flint_mpf_set_ui(tmp, 1);
            mpf_div_2exp(tmp, tmp, prec - 3);
            if (flint_mpf_cmp_ui(norm, 0) != 0 && mpf_cmp(dot, tmp) > 0)
            {
                flint_printf("FAIL:\n");
                flint_printf("Q:\n");
                mpf_mat_print(Q);
                mpf_out_str(stdout, 10, 0, norm);
                flint_printf("\n");
                flint_printf("%d\n", j);
                fflush(stdout);
                flint_abort();
            }
            mpf_clear(norm);
            for (k = j + 1; k < n; k++)
            {
                flint_mpf_set_ui(dot, 0);
                for (l = 0; l < m; l++)
                {
                    mpf_mul(tmp, mpf_mat_entry(Q, l, j),
                            mpf_mat_entry(Q, l, k));
                    mpf_add(dot, dot, tmp);
                }

                flint_mpf_set_ui(tmp, 1);
                mpf_div_2exp(tmp, tmp, prec);
                mpf_abs(dot, dot);
                if (mpf_cmp(dot, tmp) > 0)
                {
                    flint_printf("FAIL:\n");
                    flint_printf("Q:\n");
                    mpf_mat_print(Q);
                    mpf_out_str(stdout, 10, 0, dot);
                    flint_printf("\n");
                    flint_printf("%d %d\n", j, k);
                    fflush(stdout);
                    flint_abort();
                }
            }
        }

        mpf_mat_clear(A);
        mpf_mat_clear(Q);
        mpf_mat_clear(R);
        mpf_mat_clear(B);
        mpf_clears(dot, tmp, NULL);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
