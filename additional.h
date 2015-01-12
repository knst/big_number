/*
 * File:   additional.h
 * Author: knst
 *
 * Created on 27 Май 2010 г., 3:50
 */

#ifndef _ADDITIONAL_H
#define _ADDITIONAL_H

#include "BN.h"
#include "BNsign.h"
#include <vector>


BN gcdEuclidean(BN x,BN y);
BN gcdBinary(BN x,BN y);
BN gcdLehmer(BN x,BN y);
BN gcdInverseEuclidean(BN a,BN mod);
BN gcdInverseEuclideanBinary(BN x,BN mod);
BN gcdExtendedEuclideanBinary(BN xx, BN yy);
vector <BN> multi_inverse(const vector <BN> &x, const BN &mod);
BN Garner(vector <BN> m, vector <BN> v);
BN CTO(vector <BN> m, vector <BN> v);
BN expSimutaneousMul(vector <BN> G, vector <BN> exp, BN mod);
vector <BN> expSimutaneousMulPrecomputation(vector <BN> g, vector <BN> exp, BN mod);

//fixed-base
vector <BN> expFixedBaseWindowPrecomputation(BN g, int t, BN mod);
BN expFixedBaseWindow(const vector <BN> & g, BN exp, BN mod);
BN expFixedBaseEuclidean(const vector<BN> & g, BN exp, BN mod);

//exponent-recording
vector <int> karyStringReplacementRepresentation(BN exp, int k);
BN expkaryStringReplacement(const BN &g, const vector <int> & exp, const BN & mod, int k);

#endif /* _ADDITIONAL_H */

