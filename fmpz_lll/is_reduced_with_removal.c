/*
    Copyright (C) 2009, 2010 William Hart
    Copyright (C) 2009, 2010 Andy Novocin
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_lll.h"

int
fmpz_lll_is_reduced_with_removal(const fmpz_mat_t B, const fmpz_lll_t fl,
                                 const fmpz_t gs_B, int newd, flint_bitcnt_t prec)
{
    int res;
    fmpz_mat_t BB;


    if (gs_B == NULL)
        return fmpz_lll_is_reduced(B, fl, prec);

    /*
        The is_reduced checkers below could only return 1 when B has full row
        rank. Therefore, the definition of fmpz_lll_is_reduced needs to include
        a stripping out of *initial* zero rows (and possibly the corresponding
        columns in case rt == GRAM) followed by a call to one of these checkers.
    */
    if (fl->rt == Z_BASIS)
    {
        _fmpz_mat_read_only_window_init_strip_initial_zero_rows(BB, B);
    }
    else {
        _fmpz_mat_read_only_window_init_strip_initial_zero_rows_and_corresponding_cols(BB, B);
    }

    if (fmpz_lll_is_reduced_d_with_removal(BB, fl, gs_B, newd))
    {
        res = 1;
    }
    else if (fmpz_lll_is_reduced_mpfr_with_removal(BB, fl, gs_B, newd, prec))
    {
        res = 1;
    }
    else
    {
        if (fl->rt == Z_BASIS)
            res = fmpz_mat_is_reduced_with_removal(BB, fl->delta, fl->eta, gs_B, newd);
        else
            res = fmpz_mat_is_reduced_gram_with_removal(BB, fl->delta, fl->eta, gs_B, newd);
    }

    _fmpz_mat_read_only_window_clear(BB);
    return res;
}
