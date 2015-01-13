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

#endif /* _ADDITIONAL_H */

