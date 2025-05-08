#ifndef MONT_H
#define MONT_H

#include "params.h"

void proj_cswap(proj *P, proj *Q, bool c);
void xDBL(proj *Q, proj const *A, proj const *P);
void xADD(proj *S, proj const *P, proj const *Q, proj const *PQ);
void xDBLADD(proj *R, proj *S, proj const *P, proj const *Q, proj const *PQ, proj const *A);
void xMUL(proj *Q, proj const *A, proj const *P, uint const *k);
void xISOG(proj *A, proj *P1, proj *P2, proj const *K, uint64_t k);
void eval_4_isog(proj const *A, proj *P);

#endif
