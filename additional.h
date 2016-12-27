#ifndef _ADDITIONAL_H
#define _ADDITIONAL_H

#include "BN.h"
#include "BNsign.h"

#include <vector>

BN gcdEuclidean(BN x,BN y);
BN gcdBinary(BN x,BN y);
BN gcdLehmer(BN x,BN y);
BN gcdInverseEuclidean(const BN& a, const BN& mod);
BN gcdInverseEuclideanBinary(BN x,BN mod);
BN gcdExtendedEuclideanBinary(BN xx, BN yy);
std::vector <BN> multi_inverse(const std::vector <BN>& x, const BN &mod);
BN Garner(const std::vector <BN>& m, const std::vector <BN>& v);
BN CTO(const std::vector <BN>& m, const std::vector <BN>& v);

#endif /* _ADDITIONAL_H */

