#include <string.h>
#include <assert.h>

#include "uint.h"
#include "fp.h"
#include "mont.h"
#include "csidh.h"
#include "rng.h"

#include <stdio.h>

void csidh_private(private_key *priv) /* Secret key generation */
{
    memset(&priv->e, 0, sizeof(priv->e));
    for (size_t i = 0; i < NUM_PRIMES; ) {
        int8_t buf[64];
        randombytes(buf, sizeof(buf));
        for (size_t j = 0; j < sizeof(buf); ++j) {
            if (buf[j] <= MAX_EXPONENT && buf[j] >= -MAX_EXPONENT) {
                priv->e[i / 2] |= (buf[j] & 0xf) << i % 2 * 4;
                if (++i >= NUM_PRIMES)
                    break;
            }
        }
    }
}

void action(public_key *out, public_key const *in, private_key const *priv)
{
    uint cofactor;
    proj kernel_point;
    int8_t e[NUM_PRIMES];       /* Exponent vector */
    proj A = {in->A, fp_1};     /* Initial Montgomery parameter A */
    proj Aright = A;
    proj Aleft = A;
    proj right_kernel = in->P;  /* Point walking to the right */
    proj left_kernel = in->Q;   /* Point walking to the left */
    proj P = right_kernel;
    proj Q = left_kernel;

    for (size_t i = 0; i < NUM_PRIMES; ++i) {
        /* Unpacking private key */
        int8_t t = (int8_t) (priv->e[i / 2] << i % 2 * 4) >> 4;
        e[i] = t;                               

    }
    
    /* Assert Montgomery parameter is affine */
    assert(!memcmp(&A.z, &fp_1, sizeof(fp)));               

    for (int8_t k = -MAX_EXPONENT; k < MAX_EXPONENT; ++k) {

        proj right_kernel = P;  /* Point walking to the right */
        proj left_kernel = Q;   /* Point walking to the left */

        for (size_t i = 0; i < NUM_PRIMES; ++i) {

            cofactor = uint_1;

            for (size_t j = (i+1); j < NUM_PRIMES; ++j) {
                /* Constructing the current cofactor */
                uint_mul3_64(&cofactor, &cofactor, primes[j]);
                /* This is just the previous one divided by primes[i] */
            }

            uint ell_i;
            uint_set(&ell_i, primes[i]);
            bool c = (e[i] <= k);

            proj_cswap(&Aright, &Aleft, c);
            proj_cswap(&right_kernel, &left_kernel, c);
            proj_cswap(&P, &Q, c);

            xMUL(&kernel_point, &Aright, &right_kernel, &cofactor);
            /* kernel_point = cofactor * right_kernel; this has order primes[i] */
            xISOG(&Aright, &Q, &right_kernel, &kernel_point, primes[i]);
            /* We computed the isogeny with kernel_point, and push through Q and right_kernel */
            xMUL(&left_kernel, &Aleft, &left_kernel, &ell_i);
            /* Now both left_kernel and right_kernel have order equal to cofactor */

            proj_cswap(&Aright, &Aleft, c);
            proj_cswap(&right_kernel, &left_kernel, c);
            proj_cswap(&P, &Q, c);
        }

        fp_inv(&Aleft.z);
        fp_mul2(&Aleft.x, &Aleft.z);
        Aleft.z = fp_1;

        fp_inv(&Aright.z);
        fp_mul2(&Aright.x, &Aright.z);
        Aright.z = fp_1;

        eval_4_isog(&Aleft, &P);
        Aleft = Aright;
    }

    out->A = Aright.x;
    out->P = P;
    out->Q = Q;
}

bool csidh(public_key *out, public_key const *in, private_key const *priv)
{
    action(out, in, priv);
    return true;
}