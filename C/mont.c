
#include <assert.h>
#include <string.h>

#include "params.h"
#include "uint.h"
#include "fp.h"
#include "mont.h"

void proj_cswap(proj *P, proj *Q, bool c)
{
    fp_cswap(&P->x,&Q->x,c);
    fp_cswap(&P->z,&Q->z,c);
}

void xDBLADD(proj *R, proj *S, proj const *P, proj const *Q, proj const *PQ, proj const *A)
{
    fp a, b, c, d;

    fp_add3(&a, &Q->x, &Q->z);
    fp_sub3(&b, &Q->x, &Q->z);
    fp_add3(&c, &P->x, &P->z);
    fp_sub3(&d, &P->x, &P->z);
    fp_sq2(&R->x, &c);
    fp_sq2(&S->x, &d);
    fp_mul2(&c, &b);
    fp_mul2(&d, &a);
    fp_sub3(&b, &R->x, &S->x);
    fp_add3(&a, &A->z, &A->z); /* multiplication by 2 */
    fp_mul3(&R->z, &a, &S->x);
    fp_add3(&S->x, &A->x, &a);
    fp_add2(&R->z, &R->z); /* multiplication by 2 */
    fp_mul2(&R->x, &R->z);
    fp_mul2(&S->x, &b);
    fp_sub3(&S->z, &c, &d);
    fp_add2(&R->z, &S->x);
    fp_add3(&S->x, &c, &d);
    fp_mul2(&R->z, &b);
    fp_sq2(&d, &S->z);
    fp_sq2(&b, &S->x);
    fp_mul3(&S->x, &PQ->z, &b);
    fp_mul3(&S->z, &PQ->x, &d);
}

void xDBL(proj *Q, proj const *A, proj const *P)
{
    fp a, b, c;
    fp_add3(&a, &P->x, &P->z);
    fp_sq1(&a);
    fp_sub3(&b, &P->x, &P->z);
    fp_sq1(&b);
    fp_sub3(&c, &a, &b);
    fp_add2(&b, &b); fp_add2(&b, &b); /* multiplication by 4 */
    fp_mul2(&b, &A->z);
    fp_mul3(&Q->x, &a, &b);
    fp_add3(&a, &A->z, &A->z); /* multiplication by 2 */
    fp_add2(&a, &A->x);
    fp_mul2(&a, &c);
    fp_add2(&a, &b);
    fp_mul3(&Q->z, &a, &c);
}

void xADD(proj *S, proj const *P, proj const *Q, proj const *PQ)
{
    fp a, b, c, d;
    fp_add3(&a, &P->x, &P->z);
    fp_sub3(&b, &P->x, &P->z);
    fp_add3(&c, &Q->x, &Q->z);
    fp_sub3(&d, &Q->x, &Q->z);
    fp_mul2(&a, &d);
    fp_mul2(&b, &c);
    fp_add3(&c, &a, &b);
    fp_sub3(&d, &a, &b);
    fp_sq1(&c);
    fp_sq1(&d);
    fp_mul3(&S->x, &PQ->z, &c);
    fp_mul3(&S->z, &PQ->x, &d);
}

/* Montgomery ladder. */
/* P must not be the unique point of order 2. */
/* not constant-time! */
void xMUL(proj *Q, proj const *A, proj const *P, uint const *k)
{
    proj R = *P;
    const proj Pcopy = *P; /* in case Q = P */

    Q->x = fp_1;
    Q->z = fp_0;

    unsigned long i = 64 * LIMBS;
    while (--i && !uint_bit(k, i));

    do {

        bool bit = uint_bit(k, i);

        if (bit) { proj T = *Q; *Q = R; R = T; } /* not constant-time */
        //fp_cswap(&Q->x, &R.x, bit);
        //fp_cswap(&Q->z, &R.z, bit);

        xDBLADD(Q, &R, Q, &R, &Pcopy, A);

        if (bit) { proj T = *Q; *Q = R; R = T; } /* not constant-time */
        //fp_cswap(&Q->x, &R.x, bit);
        //fp_cswap(&Q->z, &R.z, bit);

    } while (i--);
}

/* computes the isogeny with kernel point K of order k */
/* returns the new curve coefficient A and the image of two points P1 and P2 */
/* (obviously) not constant time in k */
void xISOG(proj *A, proj *P1, proj *P2, proj const *K, uint64_t k)
{
    assert (k >= 3);
    assert (k % 2 == 1);

    fp tmp0, tmp1;
    fp T[4] = {K->z, K->x, K->x, K->z};
    proj Q1, Q2;

    fp_mul3(&Q1.x,  &P1->x, &K->x);
    fp_mul3(&tmp0, &P1->z, &K->z);
    fp_sub2(&Q1.x,  &tmp0);

    fp_mul3(&Q1.z,  &P1->x, &K->z);
    fp_mul3(&tmp0, &P1->z, &K->x);
    fp_sub2(&Q1.z,  &tmp0);

    fp_mul3(&Q2.x,  &P2->x, &K->x);
    fp_mul3(&tmp0, &P2->z, &K->z);
    fp_sub2(&Q2.x,  &tmp0);

    fp_mul3(&Q2.z,  &P2->x, &K->z);
    fp_mul3(&tmp0, &P2->z, &K->x);
    fp_sub2(&Q2.z,  &tmp0);

    proj M[3] = {*K};
    xDBL(&M[1], A, K);

    for (uint64_t i = 1; i < k / 2; ++i) {

        if (i >= 2)
            xADD(&M[i % 3], &M[(i - 1) % 3], K, &M[(i - 2) % 3]);

        fp_mul3(&tmp0, &M[i % 3].x, &T[0]);
        fp_mul3(&tmp1, &M[i % 3].z, &T[1]);
        fp_add3(&T[0], &tmp0, &tmp1);

        fp_mul2(&T[1], &M[i % 3].x);

        fp_mul3(&tmp0, &M[i % 3].z, &T[2]);
        fp_mul3(&tmp1, &M[i % 3].x, &T[3]);
        fp_add3(&T[2], &tmp0, &tmp1);

        fp_mul2(&T[3], &M[i % 3].z);


        fp_mul3(&tmp0, &P1->x, &M[i % 3].x);
        fp_mul3(&tmp1, &P1->z, &M[i % 3].z);
        fp_sub2(&tmp0, &tmp1);
        fp_mul2(&Q1.x,  &tmp0);

        fp_mul3(&tmp0, &P1->x, &M[i % 3].z);
        fp_mul3(&tmp1, &P1->z, &M[i % 3].x);
        fp_sub2(&tmp0, &tmp1);
        fp_mul2(&Q1.z,  &tmp0);
    
        fp_mul3(&tmp0, &P2->x, &M[i % 3].x);
        fp_mul3(&tmp1, &P2->z, &M[i % 3].z);
        fp_sub2(&tmp0, &tmp1);
        fp_mul2(&Q2.x,  &tmp0);

        fp_mul3(&tmp0, &P2->x, &M[i % 3].z);
        fp_mul3(&tmp1, &P2->z, &M[i % 3].x);
        fp_sub2(&tmp0, &tmp1);
        fp_mul2(&Q2.z,  &tmp0);
    }

    fp_mul2(&T[0], &T[1]);
    fp_add2(&T[0], &T[0]); /* multiplication by 2 */

    fp_sq1(&T[1]);

    fp_mul2(&T[2], &T[3]);
    fp_add2(&T[2], &T[2]); /* multiplication by 2 */

    fp_sq1(&T[3]);

    /* Ax := T[1] * T[3] * Ax - 3 * Az * (T[1] * T[2] - T[0] * T[3]) */
    fp_mul3(&tmp0, &T[1], &T[2]);
    fp_mul3(&tmp1, &T[0], &T[3]);
    fp_sub2(&tmp0, &tmp1);
    fp_mul2(&tmp0, &A->z);
    fp_add3(&tmp1, &tmp0, &tmp0); fp_add2(&tmp0, &tmp1); /* multiplication by 3 */

    fp_mul3(&tmp1, &T[1], &T[3]);
    fp_mul2(&tmp1, &A->x);

    fp_sub3(&A->x, &tmp1, &tmp0);

    /* Az := Az * T[3]^2 */
    fp_sq1(&T[3]);
    fp_mul2(&A->z, &T[3]);

    /* X := X * Xim^2, Z := Z * Zim^2 */
    fp_sq1(&Q1.x);
    fp_sq1(&Q1.z);
    fp_mul2(&P1->x, &Q1.x);
    fp_mul2(&P1->z, &Q1.z);

    fp_sq1(&Q2.x);
    fp_sq1(&Q2.z);
    fp_mul2(&P2->x, &Q2.x);
    fp_mul2(&P2->z, &Q2.z);
}

/* Evaluates a Montgomery 4-isogeny at the point P */
/* We assume that A is given in affine form (i.e. *A.z = fp_1) */
void eval_4_isog(proj const *A, proj *P){

    assert (!memcmp(&A->z, &fp_1, sizeof(fp)));

    fp tmp0, tmp1, tmp2;            /* create temporary storage of intermediate results */

    fp_mul3(&tmp0, &P->x, &P->z);   /* tmp0 = XZ */
    fp_sub2(&P->x, &P->z);          /* X = X - Z */
    fp_sq1(&P->x);                  /* X = (X-Z)^2 */
    P->z = tmp0;                    /* Z = XZ */

    fp_sub3(&tmp2, &A->x, &fp_2);   /* tmp2 = A - 2 */
    fp_mul2(&tmp2, &P->x);          /* tmp2 = (A - 2)X */
    fp_mul2(&tmp2, &P->z);          /* tmp2 = (A - 2)XZ */
    fp_add3(&tmp1, &A->x, &fp_2);   /* tmp1 = A + 2 */
    fp_add3(&tmp0, &P->z, &P->z);   /* tmp0 = 2Z */
    fp_add2(&tmp0, &tmp0);          /* tmp0 = 4Z */
    fp_add2(&tmp0, &P->x);          /* tmp0 = X + 4Z */
    fp_mul2(&tmp1, &P->z);          /* tmp1 = (A + 2) Z */
    fp_add2(&tmp1, &P->x);          /* tmp1 = X + (A + 2) Z */
    fp_mul3(&P->x, &tmp0, &tmp1);   /* X = (X + 4Z)(X + (A + 2)Z) */
    P->z = tmp2;                    /* Z = (A - 2)XZ */
}
