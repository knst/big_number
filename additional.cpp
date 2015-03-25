#include "additional.h"

#include "BN.h"
#include "BNsign.h"

#include <algorithm>
#include <iostream>
#include <vector>

using namespace std;

BN gcdEuclidean(BN a,BN b)
{
    while(!b.is0()) {
        a = a % b;
        a.swap(b);
    }
    return a;
}

BN gcdInverseEuclidean(const BN& a, const BN& mod)
{
    BN startMod(mod);

    BN A(a % mod);
    BN B(mod);

    if(A.is0())
        return move(BN::bn0());

    BNsign x0;
    BNsign x1(BN::bn0());
    BNsign x2(BN::bn1());
    BN Div;
    BN Mod;
    while (true) {
        B.divmod(A, Div, Mod);
        if (Mod.is0())
            break;
        x0 = move(x1);
        x1 = move(x2);

        x2 = x0 - x1 * (BNsign)Div;

        B = move(A);
        A = move(Mod);
    }

    if(A != BN::bn1())
        return move(BN::bn0());
    if(x2.sign && !x2.value.is0())
        return move(startMod -= x2.value);
    return move(x2.value);
}

BN gcdBinary(BN a,BN b)
{
    if(a.is0())
        return b;

    if(b.is0())
        return a;

    size_t acount = a.countzeroright();
    size_t bcount = b.countzeroright();

    a >>= acount;
    b >>= bcount;

    // TODO why swap(a, b) is faster than a.swap(b) ?
    if(a < b)
        swap(a, b);

    while(!a.is0()) {
        if(a >= b) {
            a -= b;
            a >>= a.countzeroright();
        } else {
            b -= a;
            b >>= b.countzeroright();
        }
    }
    return move(b <<= min(acount,bcount));
}

BN gcdInverseEuclideanBinary(BN xx, BN mod)
{
    BN& yy = mod;
    BN g(BN::bn1());
    size_t xcount=xx.countzeroright();
    size_t ycount=yy.countzeroright();
    if (xcount && ycount)
        return BN::bn0();

    xx >>= min(xcount,ycount);
    yy >>= min(xcount,ycount);

    BN u=xx;
    BN v=yy;
    BNsign x=xx;
    BNsign y=yy;

    BNsign a(BN::bn1());
    BNsign b(BN::bn0());
    BNsign c(BN::bn0());
    BNsign d(BN::bn1());

    do {
        size_t uZeros = u.countzeroright();
        size_t vZeros = v.countzeroright();
        for (size_t i = 0; i < uZeros; ++i) {
            if(!a.value.isEven() || !b.value.isEven()) {
                a += y;
                b -= x;
            }
            a.value >>= 1;
            b.value >>= 1;
        }
        for (size_t i = 0; i < vZeros; ++i) {
            if(!c.value.isEven() || !d.value.isEven()) {
                c += y;
                d -= x;
            }
            c.value >>= 1;
            d.value >>= 1;
        }
        u >>= uZeros;
        v >>= vZeros;
        if(u >= v) {
            u -= v;
            a -= c;
            b -= d;
        } else {
            v -= u;
            c -= a;
            d -= b;
        }
    }
    while(!u.is0());

    if(v != BN::bn1())
        return move(BN::bn0());

    if(c.sign && !c.value.is0())
        return move(mod -= c.value % mod);
    return move(c.value % mod);
}

vector <BN> multi_inverse(const vector <BN> &x, const BN &mod)
{
    int count=x.size();
    vector <BN> a;
    vector <BN> x_inverse;

    a.reserve(count);
    x_inverse.reserve(count);

    a.push_back(x[0]);

    for( vector<BN>::const_iterator xi=x.begin()+1;
        xi!=x.end();
        xi++)
        a.push_back(a.back() * (*xi) % mod);

    BN curr = gcdInverseEuclideanBinary(a[0],mod);
    for(int i=x.size()-1;i>0;i--)
    {
        x_inverse.push_back(a[i-1] * curr % mod);    //index = i
        curr = x[i] * curr % mod;
        //x_inverse.push_back( (a[i-1] * curr).reduction_special(mod));
        //curr = (x[i] * curr).reduction_special(mod);
    }
    x_inverse.push_back(curr);                //index = 0
    reverse(x_inverse.begin(),x_inverse.end());
    return x_inverse;
}

BN Garner(const std::vector <BN>& m, const std::vector <BN>& v)
{
    if(m.size()!=v.size())
        throw invalid_argument("Garner: Sizes of arrays are different");
    int t = m.size();
    vector<BN> C(t);
    BN u;
    for(int i=1;i<t;i++)
    {
        C[i] = BN::bn1();
        for(int j=0; j<i; j++)
        {
            u = gcdInverseEuclideanBinary(m[j],m[i]);
            C[i] = u * C[i] % m[i];
        }
    }
    u = v[0];
    BN x = u;
    BN Pm(BN::bn1());
    for(int i=1;i<t;i++)
    {
        BN xmi = x%m[i];
        BN vix = (v[i] >= xmi ? v[i] - xmi : m[i] + v[i] - xmi);
        //u = (v[i] - x) * C[i] % m[i];
        u = vix * C[i] % m[i];
        Pm = Pm * m[i-1];
        x += u * Pm;
    }
    return x;
}

BN CTO(const std::vector <BN>& m, const std::vector <BN>& v)
{
    if(m.size()!=v.size())
        throw invalid_argument("CTO: Sizes of arrays are different");

    int t = m.size();
    BN M = m[0];
    for (int i=1;i<t;i++)
        M = M * m[i];
    BN x(BN::bn0());
    for (int i=0;i<t;i++)
    {
        BN Mi = M/m[i];
        x += v[i] * (M/m[i]) * gcdInverseEuclideanBinary(M/m[i],m[i]) % M;
    }
    return x % M;
}
