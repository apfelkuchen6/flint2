/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"


/* A = D - B*C, coefficients of D are clobbered */
static void _fq_nmod_mpoly_mulsub(fq_nmod_mpoly_t A,
                   mp_limb_t * Dcoeff, const ulong * Dexp, slong Dlen,
             const mp_limb_t * Bcoeffs, const ulong * Bexp, slong Blen,
             const mp_limb_t * Ccoeffs, const ulong * Cexp, slong Clen,
                         flint_bitcnt_t bits, slong N, const ulong * cmpmask,
                                                     const fq_nmod_ctx_t fqctx)
{
    slong d = fq_nmod_ctx_degree(fqctx);
    slong i, j;
    slong next_loc;
    slong heap_len = 2; /* heap zero index unused */
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    slong Di;
    slong Alen;
    mp_limb_t * Acoeffs = A->coeffs;
    ulong * Aexps = A->exps;
    ulong * exp, * exps;
    ulong ** exp_list;
    slong exp_next;
    slong * hind;
    mp_limb_t * t;
    int lazy_size = _n_fq_dot_lazy_size(Blen, fqctx);
    TMP_INIT;

    FLINT_ASSERT(Blen > 0);
    FLINT_ASSERT(Clen > 0);

    TMP_START;

    next_loc = Blen + 4; /* something bigger than heap can ever be */
    heap = (mpoly_heap_s *) TMP_ALLOC((Blen + 1)*sizeof(mpoly_heap_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(Blen*sizeof(mpoly_heap_t));
    store = store_base = (slong *) TMP_ALLOC(2*Blen*sizeof(slong));
    exps = (ulong *) TMP_ALLOC(Blen*N*sizeof(ulong));
    exp_list = (ulong **) TMP_ALLOC(Blen*sizeof(ulong *));
    for (i = 0; i < Blen; i++)
        exp_list[i] = exps + i*N;

    /* space for heap indices */
    hind = (slong *) TMP_ALLOC(Blen*sizeof(slong));
    for (i = 0; i < Blen; i++)
        hind[i] = 1;

    t = (mp_limb_t *) TMP_ALLOC(6*d*sizeof(mp_limb_t));

    /* start with no heap nodes and no exponent vectors in use */
    exp_next = 0;

    /* put (0, 0, Bexp[0] + Cexp[0]) on heap */
    x = chain + 0;
    x->i = 0;
    x->j = 0;
    x->next = NULL;

    heap[1].next = x;
    heap[1].exp = exp_list[exp_next++];

    if (bits <= FLINT_BITS)
        mpoly_monomial_add(heap[1].exp, Bexp + N*0, Cexp + N*0, N);
    else
        mpoly_monomial_add_mp(heap[1].exp, Bexp + N*0, Cexp + N*0, N);

    hind[0] = 2*1 + 0;

    Alen = 0;
    Di = 0;
    while (heap_len > 1)
    {
        exp = heap[1].exp;

        while (Di < Dlen && mpoly_monomial_gt(Dexp + N*Di, exp, N, cmpmask))
        {
            _fq_nmod_mpoly_fit_length(&Acoeffs, &A->coeffs_alloc, d,
                                      &Aexps, &A->exps_alloc, N, Alen + 1);

            mpoly_monomial_set(Aexps + N*Alen, Dexp + N*Di, N);
            _n_fq_set(Acoeffs + d*Alen, Dcoeff + d*Di, d);
            Alen++;
            Di++;
        }

        _fq_nmod_mpoly_fit_length(&Acoeffs, &A->coeffs_alloc, d,
                                  &Aexps, &A->exps_alloc, N, Alen + 1);

        mpoly_monomial_set(Aexps + N*Alen, exp, N);

        _nmod_vec_zero(t, 6*d);

        switch (lazy_size)
        {
        case 1:
            do {
                exp_list[--exp_next] = heap[1].exp;
                x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
                do {
                    hind[x->i] |= WORD(1);
                    *store++ = x->i;
                    *store++ = x->j;
                    _n_fq_madd2_lazy1(t, Bcoeffs + d*x->i, Ccoeffs + d*x->j, d);
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));
            _n_fq_reduce2_lazy1(t, d, fqctx->mod);
            break;

        case 2:
            do {
                exp_list[--exp_next] = heap[1].exp;
                x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
                do {
                    hind[x->i] |= WORD(1);
                    *store++ = x->i;
                    *store++ = x->j;
                    _n_fq_madd2_lazy2(t, Bcoeffs + d*x->i, Ccoeffs + d*x->j, d);
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));
            _n_fq_reduce2_lazy2(t, d, fqctx->mod);
            break;

        case 3:
            do {
                exp_list[--exp_next] = heap[1].exp;
                x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
                do {
                    hind[x->i] |= WORD(1);
                    *store++ = x->i;
                    *store++ = x->j;
                    _n_fq_madd2_lazy3(t, Bcoeffs + d*x->i, Ccoeffs + d*x->j, d);
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));
            _n_fq_reduce2_lazy3(t, d, fqctx->mod);
            break;

        default:
            do {
                exp_list[--exp_next] = heap[1].exp;
                x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
                do {
                    hind[x->i] |= WORD(1);
                    *store++ = x->i;
                    *store++ = x->j;
                    _n_fq_madd2(t, Bcoeffs + d*x->i,
                                   Ccoeffs + d*x->j, fqctx, t + 2*d);
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));
            break;
        }

        _n_fq_reduce2(Acoeffs + d*Alen, t, fqctx, t + 2*d);

        if (Di < Dlen && mpoly_monomial_equal(Dexp + N*Di, exp, N))
        {
            _n_fq_sub(Acoeffs + d*Alen, Dcoeff + d*Di, Acoeffs + d*Alen, d, fqctx->mod);
            Di++;
        }
        else
        {
            _n_fq_neg(Acoeffs + d*Alen, Acoeffs + d*Alen, d, fqctx->mod);            
        }

        Alen += !_n_fq_is_zero(Acoeffs + d*Alen, d);

        /* process nodes taken from the heap */
        while (store > store_base)
        {
            j = *--store;
            i = *--store;

            /* should we go right? */
            if (  (i + 1 < Blen)
                && (hind[i + 1] == 2*j + 1)
               )
            {
                x = chain + i + 1;
                x->i = i + 1;
                x->j = j;
                x->next = NULL;

                hind[x->i] = 2*(x->j+1) + 0;

                if (bits <= FLINT_BITS)
                    mpoly_monomial_add(exp_list[exp_next], Bexp + N*x->i,
                                                           Cexp + N*x->j, N);
                else
                    mpoly_monomial_add_mp(exp_list[exp_next], Bexp + N*x->i,
                                                              Cexp + N*x->j, N);

                exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
            }

            /* should we go up? */
            if (  (j + 1 < Clen)
               && ((hind[i] & 1) == 1)
               && ((i == 0) || (hind[i - 1] >= 2*(j + 2) + 1))
               )
            {
                x = chain + i;
                x->i = i;
                x->j = j + 1;
                x->next = NULL;

                hind[x->i] = 2*(x->j+1) + 0;

                if (bits <= FLINT_BITS)
                    mpoly_monomial_add(exp_list[exp_next], Bexp + N*x->i,
                                                           Cexp + N*x->j, N);
                else
                    mpoly_monomial_add_mp(exp_list[exp_next], Bexp + N*x->i,
                                                              Cexp + N*x->j, N);

                exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
            }
        }
    }

    FLINT_ASSERT(Di <= Dlen);
    if (Di < Dlen)
    {
        _fq_nmod_mpoly_fit_length(&Acoeffs, &A->coeffs_alloc, d,
                                  &Aexps, &A->exps_alloc, N, Alen + Dlen - Di);

        _nmod_vec_set(Acoeffs + d*Alen, Dcoeff + d*Di, d*(Dlen - Di));
        mpoly_copy_monomials(Aexps + N*Alen, Dexp + N*Di, Dlen - Di, N);
        Alen += Dlen - Di;
    }

    A->coeffs = Acoeffs;
    A->exps = Aexps;
    A->length = Alen;

    TMP_END;
}


/*
    Q = A/B return 1 if division was exact.
    nmainvars is the number of main vars.
*/
int fq_nmod_mpolyuu_divides(
    fq_nmod_mpolyu_t Q,
    const fq_nmod_mpolyu_t A,
    const fq_nmod_mpolyu_t B,
    slong nmainvars,
    const fq_nmod_mpoly_ctx_t ctx)
{
    flint_bitcnt_t bits = A->bits;
    fq_nmod_mpoly_struct * Bcoeff = B->coeffs;
    slong Blen = B->length;
    ulong * Bexp = B->exps;
    fq_nmod_mpoly_struct * Acoeff = A->coeffs;
    slong Alen = A->length;
    ulong * Aexp = A->exps;
    fq_nmod_mpoly_struct * a, * b, * q;
    slong N;
    ulong * cmpmask;    /* cmp mask for lesser variables */
    fq_nmod_mpoly_t T, S;
    int success;
    ulong maskhi = 0;   /* main variables are in lex */
    int lt_divides;
    slong i, j, s;
    slong next_loc, heap_len;
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    slong * hind;
    ulong mask, exp, maxexp = Aexp[Alen - 1];
    TMP_INIT;

    FLINT_ASSERT(Alen > 0);
    FLINT_ASSERT(Blen > 0);
    FLINT_ASSERT(bits == A->bits);
    FLINT_ASSERT(bits == B->bits);
    FLINT_ASSERT(bits == Q->bits);

    TMP_START;

    N = mpoly_words_per_exp(bits, ctx->minfo);
    cmpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, bits, ctx->minfo);

    /* alloc array of heap nodes which can be chained together */
    next_loc = Blen + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap1_s *) TMP_ALLOC((Blen + 1)*sizeof(mpoly_heap1_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(Blen*sizeof(mpoly_heap_t));
    store = store_base = (slong *) TMP_ALLOC(2*Blen*sizeof(slong));

    /* space for flagged heap indices */
    hind = (slong *) TMP_ALLOC(Blen*sizeof(slong));
    for (i = 0; i < B->length; i++)
        hind[i] = 1;

    /* mask with high bit set in each field of main exponent vector */
    mask = 0;
    for (i = 0; i < nmainvars; i++)
        mask = (mask << (FLINT_BITS/nmainvars)) + (UWORD(1) << (FLINT_BITS/nmainvars - 1));

    /* s is the number of terms * (latest quotient) we should put into heap */
    s = Blen;

    /* insert (-1, 0, exp2[0]) into heap */
    heap_len = 2;
    x = chain + 0;
    x->i = -WORD(1);
    x->j = 0;
    x->next = NULL;
    HEAP_ASSIGN(heap[1], Aexp[0], x);

    Q->length = 0;

    fq_nmod_mpoly_init3(T, 16, bits, ctx);
    fq_nmod_mpoly_init3(S, 16, bits, ctx);

    while (heap_len > 1)
    {
        exp = heap[1].exp;

        if (mpoly_monomial_overflows1(exp, mask))
        {
            goto not_exact_division;
        }

        fq_nmod_mpolyu_fit_length(Q, Q->length + 1, ctx);
        lt_divides = mpoly_monomial_divides1(Q->exps + Q->length, exp, Bexp[0], mask);

        T->length = 0;

        do
        {
            x = _mpoly_heap_pop1(heap, &heap_len, maskhi);
            do
            {
                *store++ = x->i;
                *store++ = x->j;
                if (x->i != -WORD(1))
                    hind[x->i] |= WORD(1);

                if (x->i == -WORD(1))
                {
                    a = Acoeff + x->j;
                    fq_nmod_mpoly_fit_length(S, T->length + a->length, ctx);
                    S->length = _fq_nmod_mpoly_add(
                                    S->coeffs, S->exps,
                                    T->coeffs, T->exps, T->length,
                                    a->coeffs, a->exps, a->length,
                                                      N, cmpmask, ctx->fqctx);
                }
                else
                {
                    b = Bcoeff + x->i;
                    q = Q->coeffs + x->j;
                    _fq_nmod_mpoly_mulsub(S,
                                    T->coeffs, T->exps, T->length,
                                    b->coeffs, b->exps, b->length,
                                    q->coeffs, q->exps, q->length,
                                                bits, N, cmpmask, ctx->fqctx);
                }
                fq_nmod_mpoly_swap(S, T, ctx);

            } while ((x = x->next) != NULL);
        } while (heap_len > 1 && heap[1].exp == exp);

        /* process nodes taken from the heap */
        while (store > store_base)
        {
            j = *--store;
            i = *--store;

            if (i == -WORD(1))
            {
                /* take next dividend term */
                if (j + 1 < Alen)
                {
                    x = chain + 0;
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;
                    _mpoly_heap_insert1(heap, Aexp[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
                }
            }
            else
            {
                /* should we go right? */
                if (  (i + 1 < Blen)
                   && (hind[i + 1] == 2*j + 1)
                   )
                {
                    x = chain + i + 1;
                    x->i = i + 1;
                    x->j = j;
                    x->next = NULL;
                    hind[x->i] = 2*(x->j + 1) + 0;
                    _mpoly_heap_insert1(heap, Bexp[x->i] + Q->exps[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
                }
                /* should we go up? */
                if (j + 1 == Q->length)
                {
                    s++;
                } else if (  ((hind[i] & 1) == 1)
                          && ((i == 1) || (hind[i - 1] >= 2*(j + 2) + 1))
                          )
                {
                    x = chain + i;
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;
                    hind[x->i] = 2*(x->j + 1) + 0;
                    _mpoly_heap_insert1(heap, Bexp[x->i] + Q->exps[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
                }
            }
        }

        if (T->length == 0)
        {
            continue;
        }

        if (mpoly_monomials_overflow_test(T->exps, T->length, bits, ctx->minfo))
        {
            goto not_exact_division;
        }

        q = Q->coeffs + Q->length;
        FLINT_ASSERT(q->bits == bits);
        b = Bcoeff + 0;
        if (!_fq_nmod_mpoly_divides_monagan_pearce(q,
                            T->coeffs, T->exps, T->length,
                            b->coeffs, b->exps, b->length,
                                              bits, N, cmpmask, ctx->fqctx))
        {
            goto not_exact_division;
        }

        if (!lt_divides || (exp^maskhi) < (maxexp^maskhi))
        {
            goto not_exact_division;
        }

        /* put newly generated quotient term back into the heap if necessary */
        if (s > 1)
        {
            i = 1;
            x = chain + i;
            x->i = i;
            x->j = Q->length;
            x->next = NULL;
            hind[x->i] = 2*(x->j + 1) + 0;
            _mpoly_heap_insert1(heap, Bexp[x->i] + Q->exps[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
        }
        s = 1;
        Q->length++;
    }

    success = 1;

cleanup:

    fq_nmod_mpoly_clear(T, ctx);
    fq_nmod_mpoly_clear(S, ctx);

    TMP_END;

    return success;

not_exact_division:

    success = 0;
    Q->length = 0;
    goto cleanup;
}
