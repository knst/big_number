#include "additional.h"

#include "BN.h"
#include "BNsign.h"

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <vector>

using namespace std;

bt2s Abs(bt2s x) {
    return x<0?-x:x;
}

BN gcdEuclidean(BN a,BN b)
{
    while(!b.is0()) {
        a = a % b;
        a.swap(b);
    }
    return a;
}

BN gcdInverseEuclidean(BN a, BN mod)
{
    BN start_mod = mod;
    a = a % mod;

    if(a.is0())
        return BN::bn0();

    BNsign x0;
    BNsign x1(BN::bn0());
    BNsign x2(BN::bn1());
    while (true) {
        BN Div;
        BN Mod;
        mod.divmod(a, Div, Mod);
        if (Mod.is0())
            break;
        x0 = move(x1);
        x1 = move(x2);

        x2 = x0 - x1 * (BNsign)Div;

        mod = move(a);
        a = move(Mod);
    }

    if(a != BN::bn1())
        return BN::bn0();
    if(x2.value.is0())
        x2.sign = false;
    if(!x2.sign)
        return x2.value;
    else
        return start_mod - x2.value;
}

BN gcdBinary(BN a,BN b)
{
    if(a.is0())
        return b;

    if(b.is0())
        return a;

    int acount = a.countzeroright();
    int bcount = b.countzeroright();

    a = a >> acount;
    b = b >> bcount;

    // TODO why swap(a, b) is faster than a.swap(b) ?
    if(a < b)
        swap(a, b);

    while(!a.is0()) {
        if(a >= b) {
            a.subappr(b);
            a = a >> a.countzeroright();
        } else {
            b.subappr(a);
            b = b >> b.countzeroright();
        }
    }
    return b << min(acount,bcount);
}

BN gcdExtendedEuclideanBinary(BN xx, BN yy)
{
    BN g(BN::bn1());
    int xcount=xx.countzeroright();
    int ycount=yy.countzeroright();
    g=g<<min(xcount,ycount);
    xx=xx>>min(xcount,ycount);
    yy=yy>>min(xcount,ycount);

    BN u=xx;
    BN v=yy;
    BNsign x=xx;
    BNsign y=yy;

    BNsign a(BN::bn1());
    BNsign b(BN::bn0());
    BNsign c(BN::bn0());
    BNsign d(BN::bn1());

    do {
        while(u.isEven())
        {
            u=u>>1;
            if(a.value.isEven() && b.value.isEven()) {
                a.value=a.value>>1;
                b.value=b.value>>1;
            } else {
                a=(a+y);
                a.value=a.value>>1;
                b=(b-x);
                b.value=b.value>>1;
            }
        }
        while(v.isEven())
        {
            v=v>>1;
            if(c.value.isEven()&&d.value.isEven()) {
                c.value=c.value>>1;
                d.value=d.value>>1;
            } else {
                c=(c+y);
                c.value=c.value>>1;
                d=(d-x);
                d.value=d.value>>1;
            }
        }
        if(u>=v) {
            u=u-v;
            a=a-c;
            b=b-d;
        } else {
            v=v-u;
            c=c-a;
            d=d-b;
        }
    }
    while(!u.is0());
    //u.is0() == true;
    BNsign A=c;
    BNsign B=d;
    A.PrintSign();
    B.PrintSign();
    (g*v).PrintDec();
    return g*v;
    // A*x + B*y == v
}

BN gcdInverseEuclideanBinary(BN xx, BN mod)
{
    if(mod.is0())
        return BN::bn0();
    if (xx == BN::bn1())
        return xx;

    BN g(BN::bn1());
    int xcount=xx.countzeroright();
    int ycount=mod.countzeroright();
    if(min(xcount,ycount))
        return BN::bn0();
    //min(xcount,ycount) = 0!
//    g=g<<min(xcount,ycount);
//    xx=xx>>min(xcount,ycount);
//    mod=mod>>min(xcount,ycount);

    BN u=xx;
    BN v=mod;
    BNsign x=xx;
    BNsign y=mod;

    BNsign a(BN::bn1());
    BNsign b(BN::bn0());
    BNsign c(BN::bn0());
    BNsign d(BN::bn1());
    if(mod.isEven())
        do
        {
            int ucount=u.countzeroright();
            u=u>>ucount;
            while(ucount)
            {
                int abcount=min(a.value.countzeroright(),b.value.countzeroright());
                if(abcount>ucount)
                {
                    a.value=a.value>>ucount;
                    b.value=b.value>>ucount;
                    ucount=0;
                }
                else if(abcount)
                {
                    a.value=a.value>>abcount;
                    b.value=b.value>>abcount;
                    ucount-=abcount;
                }
                else
                {
                    a=(a+y);
                    a.value=a.value>>1;
                    b=(b-x);
                    b.value=b.value>>1;
                    ucount--;
                }
            }
            int vcount=v.countzeroright();
            v=v>>vcount;
            while(vcount)
            {
                int cdcount=min(c.value.countzeroright(),d.value.countzeroright());
                if(cdcount>vcount) {
                    c.value=c.value>>vcount;
                    d.value=d.value>>vcount;
                    vcount=0;
                } else if(cdcount) {
                    c.value=c.value>>cdcount;
                    d.value=d.value>>cdcount;
                    vcount-=cdcount;
                } else {
                    c=(c+y);
                    c.value=c.value>>1;
                    d=(d-x);
                    d.value=d.value>>1;
                    vcount--;
                }
            }
            if(u>=v) {
                u=u-v;
                a=a-c;
                b=b-d;
            } else {
                v=v-u;
                c=c-a;
                d=d-b;
            }
        }
        while(!u.is0());
    else
        do
        {
            int ucount=u.countzeroright();
            u=u>>ucount;
            while(ucount)
            {
                int acount = a.value.countzeroright();
                if(acount) {
                    a.value=a.value>>min(acount,ucount);
                    ucount-=min(acount,ucount);
                } else {
                    a=(a+y);
                    a.value=a.value>>1;
                    ucount--;
                }
            }
            int vcount=v.countzeroright();
            v=v>>vcount;
            while(vcount) {
                int ccount = c.value.countzeroright();
                if(ccount) {
                    c.value=c.value>>min(ccount,vcount);
                    vcount-=min(ccount,vcount);
                } else {
                    c=(c+y);
                    c.value=c.value>>1;
                    vcount--;
                }
            }
            if(u>=v) {
                u=u-v;
                a=a-c;
            } else {
                v=v-u;
                c=c-a;
            }
        }
        while(!u.is0());

    if(v!=BN::bn1())
        return BN::bn0();

    if(c.value.is0())
        c.sign=false;
    if(c.sign)
        return mod-(c.value)%mod;
    else
        return c.value%mod;
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

//#define DEBUG_LEHMER
BN gcdLehmer(BN x,BN y)
{
    if(x<y)
        x.swap(y);

    while(y.digitCount()>1)
    {
        bt2s xp = x[x.digitCount()-1];
        bt2s yp = y[y.digitCount()-1];
        bt2s A = 1;
        bt2s B = 0;
        bt2s C = 0;
        bt2s D = 1;
        if(x.digitCount()==y.digitCount())
            while( (yp+C)!=0 && (yp+D)!=0 )
            {
                bt2s q = (xp+A)/(yp+C);
                bt2s qp= (xp+B)/(yp+D);
                if(q!=qp)
                    break;
                bt2s t;
                t = A - q*C;    A = C;    C = t;
                t = B - q*D;    B = D;    D = t;
                t = xp- q*yp;   xp= yp;   yp= t;
            }
        if(B == 0) {
            x = x % y;
            x.swap(y);
        } else {
            BN T = (B<0 ? x.mulbase(Abs(A)) - y.mulbase(Abs(B)) : y.mulbase(Abs(B)) - x.mulbase(Abs(A)));
            BN U = (D<0 ? x.mulbase(Abs(C)) - y.mulbase(Abs(D)) : y.mulbase(Abs(D)) - x.mulbase(Abs(C)));
            x=T;
            y=U;
        }
    }
    return gcdBinary(x,y);
}

BN Garner(const std::vector <BN>& m, const std::vector <BN>& v)
{
    if(m.size()!=v.size())
        throw "Garner: size's M and V are different!";
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
        x = x + u * Pm;
    }
    return x;
}

BN CTO(const std::vector <BN>& m, const std::vector <BN>& v)
{
    if(m.size()!=v.size())
        throw "CTO: size's M and V are different!";

    int t = m.size();
    BN M = m[0];
    for (int i=1;i<t;i++)
        M = M * m[i];
    BN x(BN::bn0());
    for (int i=0;i<t;i++)
    {
        BN Mi = M/m[i];
        x = x + v[i] * (M/m[i]) * gcdInverseEuclideanBinary(M/m[i],m[i]) % M;
    }
    return x%M;
}
