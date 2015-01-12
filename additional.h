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


class SDbn {
private:
    int lenght;
    char *d;
public:
    int Lenght();
    SDbn(const SDbn&);
    SDbn(const BN&);
    SDbn operator = (const SDbn&);
    char operator [] (int);
};

class Fbc {
private:
    BN ** G;
    BN m;
    int h;
    int v;
    int b;
    int a;
    int BinaryToInt (int t, int *e);
public:
    Fbc(BN g, BN mod, int bitcount);
    ~Fbc();
    BN fixed_base_comb(BN exp);
};



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
BN expSignDigitRightToLeft(BN g, BN exponent, BN mod);
vector <int> karyStringReplacementRepresentation(BN exp, int k);
BN expkaryStringReplacement(const BN &g, const vector <int> & exp, const BN & mod, int k);

#endif /* _ADDITIONAL_H */

