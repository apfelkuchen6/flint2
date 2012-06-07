/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2011 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "nmod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;
    flint_randinit(state);

    printf("xgcd_euclidean....");
    fflush(stdout);

    /* 
       Compare with result from gcd and check a*s + b*t = g
    */
    for (i = 0; i < 10000; i++)
    {
        nmod_poly_t a, b, c, g1, s, t, g2;

        mp_limb_t n;
        do n = n_randtest_not_zero(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);
        nmod_poly_init(g1, n);
        nmod_poly_init(g2, n);
        nmod_poly_init(s, n);
        nmod_poly_init(t, n);
        
        nmod_poly_randtest(a, state, n_randtest(state) % 200);
        nmod_poly_randtest(b, state, n_randtest(state) % 200);
        nmod_poly_randtest(c, state, n_randtest(state) % 200);
        
        nmod_poly_mul(a, a, c);
        nmod_poly_mul(b, b, c);

        nmod_poly_gcd_euclidean(g1, a, b);
        nmod_poly_xgcd_euclidean(g2, s, t, a, b);
        
        nmod_poly_mul(s, s, a);
        nmod_poly_mul(t, t, b);
        nmod_poly_add(s, s, t);

        result = (nmod_poly_equal(g1, g2) && nmod_poly_equal(s, g1) 
            && (g1->length == 0 || g1->coeffs[g1->length - 1] == 1));
        if (!result)
        {
            printf("FAIL:\n");
            nmod_poly_print(a), printf("\n\n");
            nmod_poly_print(b), printf("\n\n");
            nmod_poly_print(c), printf("\n\n");
            nmod_poly_print(g1), printf("\n\n");
            nmod_poly_print(g2), printf("\n\n");
            nmod_poly_print(s), printf("\n\n");
            nmod_poly_print(t), printf("\n\n");
            printf("n = %ld\n", n);
            abort();
        }
        
        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
        nmod_poly_clear(g1);
        nmod_poly_clear(g2);
        nmod_poly_clear(s);
        nmod_poly_clear(t);
    }

    /* Check aliasing of a and g */
    for (i = 0; i < 200; i++)
    {
        nmod_poly_t a, b, g, s, t;

        mp_limb_t n;
        do n = n_randtest(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(g, n);
        nmod_poly_init(s, n);
        nmod_poly_init(t, n);
        nmod_poly_randtest(a, state, n_randtest(state) % 200);
        nmod_poly_randtest(b, state, n_randtest(state) % 200);
        
        nmod_poly_xgcd_euclidean(g, s, t, a, b);
        nmod_poly_xgcd_euclidean(a, s, t, a, b);

        result = (nmod_poly_equal(a, g));
        if (!result)
        {
            printf("FAIL:\n");
            nmod_poly_print(a), printf("\n\n");
            nmod_poly_print(b), printf("\n\n");
            nmod_poly_print(g), printf("\n\n");
            printf("n = %ld\n", n);
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(g);
        nmod_poly_clear(s);
        nmod_poly_clear(t);
    }

    /* Check aliasing of b and g */
    for (i = 0; i < 200; i++)
    {
        nmod_poly_t a, b, g, s, t;

        mp_limb_t n;
        do n = n_randtest(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(g, n);
        nmod_poly_init(s, n);
        nmod_poly_init(t, n);
        nmod_poly_randtest(a, state, n_randtest(state) % 200);
        nmod_poly_randtest(b, state, n_randtest(state) % 200);
       
        nmod_poly_xgcd_euclidean(g, s, t, a, b);
        nmod_poly_xgcd_euclidean(b, s, t, a, b);

        result = (nmod_poly_equal(b, g));
        if (!result)
        {
            printf("FAIL:\n");
            nmod_poly_print(a), printf("\n\n");
            nmod_poly_print(b), printf("\n\n");
            nmod_poly_print(g), printf("\n\n");
            printf("n = %ld\n", n);
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(g);
        nmod_poly_clear(s);
        nmod_poly_clear(t);
    }

    /* Check aliasing of s and a */
    for (i = 0; i < 200; i++)
    {
        nmod_poly_t a, b, g, s, t;

        mp_limb_t n;
        do n = n_randtest(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(g, n);
        nmod_poly_init(s, n);
        nmod_poly_init(t, n);
        nmod_poly_randtest(a, state, n_randtest(state) % 200);
        nmod_poly_randtest(b, state, n_randtest(state) % 200);
       
        nmod_poly_xgcd_euclidean(g, s, t, a, b);
        nmod_poly_xgcd_euclidean(g, a, t, a, b);

        result = (nmod_poly_equal(s, a));
        if (!result)
        {
            printf("FAIL:\n");
            nmod_poly_print(a), printf("\n\n");
            nmod_poly_print(s), printf("\n\n");
            printf("n = %ld\n", n);
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(g);
        nmod_poly_clear(s);
        nmod_poly_clear(t);
    }

    /* Check aliasing of s and b */
    for (i = 0; i < 200; i++)
    {
        nmod_poly_t a, b, g, s, t;

        mp_limb_t n;
        do n = n_randtest(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(g, n);
        nmod_poly_init(s, n);
        nmod_poly_init(t, n);
        nmod_poly_randtest(a, state, n_randtest(state) % 200);
        nmod_poly_randtest(b, state, n_randtest(state) % 200);
       
        nmod_poly_xgcd_euclidean(g, s, t, a, b);
        nmod_poly_xgcd_euclidean(g, b, t, a, b);

        result = (nmod_poly_equal(s, b));
        if (!result)
        {
            printf("FAIL:\n");
            nmod_poly_print(b), printf("\n\n");
            nmod_poly_print(s), printf("\n\n");
            printf("n = %ld\n", n);
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(g);
        nmod_poly_clear(s);
        nmod_poly_clear(t);
    }

    /* Check aliasing of t and a */
    for (i = 0; i < 200; i++)
    {
        nmod_poly_t a, b, g, s, t;

        mp_limb_t n;
        do n = n_randtest(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(g, n);
        nmod_poly_init(s, n);
        nmod_poly_init(t, n);
        nmod_poly_randtest(a, state, n_randtest(state) % 200);
        nmod_poly_randtest(b, state, n_randtest(state) % 200);
       
        nmod_poly_xgcd_euclidean(g, s, t, a, b);
        nmod_poly_xgcd_euclidean(g, s, a, a, b);

        result = (nmod_poly_equal(t, a));
        if (!result)
        {
            printf("FAIL:\n");
            nmod_poly_print(a), printf("\n\n");
            nmod_poly_print(t), printf("\n\n");
            printf("n = %ld\n", n);
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(g);
        nmod_poly_clear(s);
        nmod_poly_clear(t);
    }

    /* Check aliasing of t and b */
    for (i = 0; i < 200; i++)
    {
        nmod_poly_t a, b, g, s, t;

        mp_limb_t n;
        do n = n_randtest(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(g, n);
        nmod_poly_init(s, n);
        nmod_poly_init(t, n);
        nmod_poly_randtest(a, state, n_randtest(state) % 200);
        nmod_poly_randtest(b, state, n_randtest(state) % 200);
       
        nmod_poly_xgcd_euclidean(g, s, t, a, b);
        nmod_poly_xgcd_euclidean(g, s, b, a, b);

        result = (nmod_poly_equal(t, b));
        if (!result)
        {
            printf("FAIL:\n");
            nmod_poly_print(b), printf("\n\n");
            nmod_poly_print(t), printf("\n\n");
            printf("n = %ld\n", n);
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(g);
        nmod_poly_clear(s);
        nmod_poly_clear(t);
    }

    flint_randclear(state);

    printf("PASS\n");
    return 0;
}
