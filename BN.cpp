/*
 * File:   BN.cpp
 * Author: knst
 *
 * Created on 27 Май 2010 г., 3:16
 */

#include <cstdio>
#include <cstdlib>
#include <exception>
#include <limits>
#include <stdexcept>
#include <iostream>
#include <time.h>

#include "BN.h"
#include "additional.h"

using namespace std;

// maximal size for cache-optimazed multiplication is:
// n < bt4_max * / (bt - 1) * bt_max^3
constexpr bt4 MaximalSizeForFastMul = numeric_limits<bt4>::max() / bmax / bmax / bmax * (bmax - 1) - 1;

void BN::InitMemory(int type)
{
    switch(type) {
    case 1: {
        ba.resize(rbc);
        ba.assign(sizeof(bmax), bmax);
        break;
    }
    case -1:
        for(size_t i = 0; i < rbc; ++i)
            ba[i] = rand() % bsize;
        Norm();
        break;

    default:
        throw "Unknow type (BN::BN)";
    }
}

void BN::Norm()
{
    for(rbc = ba.size() - 1; rbc > 0 && ba[rbc] == 0; rbc--);
        rbc++;
}

BN::BN()
: rbc(1)
, ba(rbc + 1)
{
}

BN::BN(uint64_t basecount, int type)
: rbc(1)
, ba(basecount + 1)
{
    if (type == 1 || type == -1) {
        rbc = basecount;
        InitMemory(type);
    } else if (type != 0)
        throw invalid_argument("BN constructor: invalid type " + to_string(type));
}

BN::BN(uint64_t x)
: rbc((sizeof(uint64_t) + bz - 1) / bz)
, ba(rbc + 1)
{
    for(size_t i = 0; i < rbc; i++, x >>= bz8)
        ba[i] = static_cast<bt>(x);
    Norm();
}

BN::BN(const BN& bn)
: rbc(bn.rbc)
, ba(bn.ba)
{
}

BN::BN(BN&& bn)
: rbc(bn.rbc)
, ba(move(bn.ba))
{
}

BN::BN(const vector<bt>& _ba, size_t _rbc)
: rbc(_rbc)
, ba(_ba)
{
    if (!_rbc)
        Norm();
}

BN::BN(const BN& bn, size_t start, size_t count)
: rbc(count ? count : bn.rbc - start)
, ba(rbc + 1)
{
    size_t last = min(count, bn.rbc - start);
    for(size_t i = 0; i < last; i++)
        ba[i] = bn.ba[i + start];
    Norm();
}

BN::BN(const string &str,const int &status)
{
    if(status == 1)
    {
        BN bn(1,0);
        for (auto i : str) {
            if (i < '0' || i > '9')
                continue;
            bn = bn.mulbase(10) + BN(i - '0');
        }
        rbc = bn.rbc;
        ba = move(bn.ba);
        return;
    }

    size_t length = str.size();
    ba.resize((length + bz * 2 - 1) / (bz * 2));

    size_t index = ba.size() - 1;
    size_t shift = bz * 2 - (ba.size() * bz * 2 - length);

    for (size_t i = 0; i < length; ++i) {
        bt d;
        if (str[i] >= 'A' && str[i] <= 'F')
            d = str[i] - 'A' + 10;
        else if (str[i] >= 'a' && str[i] <= 'f')
            d = str[i] - 'a' + 10;
        else if (str[i] >= '0' && str[i] <= '9')
            d = str[i] - '0';
        else
            throw invalid_argument(string("BN constructor: invalid char '") + str[i] + "' in HEX string");
        ba[index] = (ba[index] << 4) | d;
        shift--;
        if (shift == 0) {
            shift = bz * 2;
            index--;
        }
    }
    Norm();
}

void BN::swap(BN& bn) {
    ba.swap(bn.ba);
    std::swap(rbc, bn.rbc);
}

BN & BN::operator = (const BN& bn)
{
    if (this == &bn)
        return *this;
    ba = bn.ba;
    rbc = bn.rbc;

    return *this;
}

BN & BN::operator = (BN&& bn)
{
    if (this == &bn)
        return *this;
    ba = move(bn.ba);
    rbc = bn.rbc;

    return *this;
}

const BN BN::operator + (const BN&bn)const {
    const BN& a = rbc > bn.rbc ? *this : bn;
    const BN& b = rbc > bn.rbc ? bn : *this;

    // TODO: replace to rbc + 1 in future, some error need It
    vector<bt> result(a.rbc + 2);

    bt over = 0;
    for(size_t pos = 0; pos < b.rbc; pos++) {
        result[pos] = a.ba[pos] + b.ba[pos];
        bt over2 = (result[pos] < b.ba[pos]);

        result[pos] = result[pos] + over;
        bt over3 = (result[pos] < over);
        over = over2 + over3;
    }

    for(size_t pos = b.rbc; pos < a.rbc; pos++) {
        result[pos] = a.ba[pos] + over;
        over = (result[pos] < over);
    }

    result[a.rbc] = over;
    BN resultBn(move(result));
    resultBn.Norm();
    return resultBn;
}


BN & BN::operator ++()
{
    size_t index = 0;
    do {
        ++ba[index];
        ++index;
    } while (index < rbc && ba[index - 1] == 0);
    if (index == rbc) {
        if (ba.size() == rbc + 1)
            ba.push_back(1);
        else
            ba[rbc] = 1;
        ++rbc;
    }
    return *this;
}

const BN BN::operator - (const BN& bn)const
{
    return this->sub(bn);
}

BN & BN::operator --()
{
    size_t index = 0;
    do {
        --ba[index];
        ++index;
    } while (index < rbc && ba[index - 1] == bmax);

    // Normalization:
    if (index == rbc)
        while (rbc && ba[rbc-1] == 0)
            --rbc;
    return *this;
}

BN BN::mulbt(size_t t) const
{
    if(t == 0)
        return *this;

    BN res(rbc + t, 0);
    for(size_t i = 0; i < rbc; ++i)
        res.ba[i + t] = ba[i];
    res.rbc = rbc + t;
    res.ba.resize(res.rbc);
    return res;
}

BN BN::divbt(size_t t) const
{
    if(t == 0)
        return *this;

    if(t >= rbc)
        return (BN) 0;

    BN res(rbc - t, 0);
    res.ba.assign(ba.begin() + t, ba.end());
    res.rbc = rbc - t;
    return res;
}

BN BN::modbt(size_t t) const
{
    if(t == 0)
        return BN(1,0);

    if(t >= rbc)
            return *this;

    BN res(move(vector<bt>(ba.begin(), ba.begin() + t)));
    res.Norm();
    res.ba.resize(res.rbc);
    return res;
}

const BN BN::mulbase(const bt &multiplier)const
{
    BN result(rbc+1, 0);
    bt2 curr = 0;
    for(size_t i = 0; i < rbc; i++, curr>>=bz8)
        result.ba[i] = curr += static_cast<bt2>(ba[i]) * multiplier;
    if (curr) {
        result.ba[rbc] = curr;
        result.rbc = rbc + 1;
    } else
        result.rbc = rbc;
    result.ba.resize(result.rbc);
    return result;
}

BN& BN::mulbaseappr(const bt &multiplier)
{
    bt2 curr = 0;
    ba.resize(rbc + 1);
    for (size_t i =0; i < rbc; ++i, curr >>= bz8)
        ba[i] = curr += ba[i] * multiplier;
    if (curr) {
        ba[rbc] = curr;
        ++rbc;
    }
    ba.resize(rbc);
    return *this;
}

const BN BN::operator * (const BN&bn)const {
    if (bn.rbc < MaximalSizeForFastMul && rbc < MaximalSizeForFastMul)
        return fast_mul(bn);

    // Else classical O(n*n) multiplication.
    if(bn.rbc == 1)
        return mulbase(bn.ba[0]);

    if (rbc == 1)
        return bn.mulbase(ba[0]);

    // Tested: b * a is faster than a * b
    const BN& b = rbc > bn.rbc ? *this : bn;
    const BN& a = rbc > bn.rbc ? bn : *this;

    BN result(a.rbc + b.rbc + 1, 0);
    for (size_t i = 0; i < b.rbc; ++i) {
        bt2 curr = 0;
        bt2 x = b.ba[i];
        for (size_t j = 0; j < a.rbc; ++j) {
            curr = (curr >> bz8) + result.ba[i + j] + x * a.ba[j];
            result.ba[i + j] = curr;
        }
        result.ba[i + a.rbc] = curr >> bz8;
    }
    result.Norm();
    return result;
}

const BN BN::fast_mul (const BN& bn) const {
    size_t n = rbc;
    size_t m = bn.rbc;

    BN result(n + m, 0);

    bt4 t = 0;
    for(size_t s = 0; s < m + n - 1; s++) {

        size_t end_index = min(n - 1, s);
        size_t start_index = s >= m ? s - m + 1 : 0;
        for(size_t i = start_index, j = s - start_index; i <= end_index; i++, --j)
            t += static_cast<bt2>(ba[i]) * bn.ba[j];


        result.ba[s] = t;
        t = t >> bz8;
    }

    result.ba[m + n - 1] = t;
    result.Norm();
    return move(result);
}

BN BN::karatsuba_add(size_t start, size_t count) const {
    BN result(count + 1, 0);

    bt2 res = 0;
    for(size_t pos = 0; pos < count; pos++) {
        res = res + (bt2) ba[start + pos] + (bt2) ba[start + count + pos];
        result.ba[pos] = res;
        res >>= bz8;
    }

    result.ba[count] = res;
    result.rbc = count + 1;
    return result;
}

BN BN::add_appr (const BN&bn, size_t mul_bt) {
    size_t result_len = max(rbc, bn.rbc + mul_bt);
    if(result_len < ba.size()) {
        bt2 res = 0;
        size_t pos = 0;
        size_t m = min(rbc-mul_bt, bn.rbc);
        for(; pos < m; pos++) {
            res += static_cast<bt2>(ba[pos + mul_bt]) + bn.ba[pos];
            ba[pos + mul_bt] = res;
            res >>= bz8;
        }

        for(; pos < bn.rbc; pos++) {
            res += bn.ba[pos];
            ba[pos + mul_bt] = res;
            res >>= bz8;
        }

        while(res) {
            res += ba[pos + mul_bt];
            ba[pos + mul_bt] = res;
            res >>= bz8;
            pos ++;
        }

        Norm();
        return *this;
    }
    *this = (*this) + bn.mulbt(mul_bt);
    return *this;
}


BN BN::karatsubaRecursive(BN & bn, size_t start, size_t len) {
    size_t n = len / 2;
    if (n < karacuba_const) {
        BN result(len + len + 2, 0);
        bt4 t = 0;
        for(size_t s = 0; s < len + len - 1; s++) {
            size_t end_index = min(len - 1, s);
            size_t start_index = s >= len ? s - len + 1 : 0;
            for(size_t i = start_index; i <= end_index; i++) {
                t += static_cast<bt2>(ba[i + start]) * bn.ba[s - i + start];
            }
            result.ba[s] = t;
            t = t >> bz8;
        }

        result.ba[len + len - 1] = t;
        result.Norm();
        return result;
    }

    BN& U = *this;
    BN& V = bn;

    // A = u1.caracuba(v1);
    BN A(U.karatsubaRecursive(V, start + n, n + len % 2));
    // B = u0.caracuba(v0);
    BN B(U.karatsubaRecursive(V, start, n));

    //  C = (u0 + u1).caracuba(v0 + v1)
    BN u01(U.karatsuba_add(start, n));
    BN v01(V.karatsuba_add(start, n));
    BN C(u01.karatsubaRecursive(v01, 0, n+1));

//    res = B + A.mulbt(2*n);
    BN res(B.ba, A.rbc + 2 * n);
    res.ba.resize(A.rbc + 2 + 2 * n);

    for(size_t i = 0; i < A.rbc; ++i)
        res.ba[i + 2 * n] = A.ba[i];

    return res + C.mulbt(n) - (A + B).mulbt(n);
}

const BN BN::karatsuba(const BN& bn)const {
    size_t x = rbc;
    size_t y = bn.rbc;
    size_t len = max(x, y);
    if(min(x, y) < karacuba_const)
        return this->fast_mul(bn);

    BN U(*this);
    BN V(bn);
    U.ba.resize(len + 1);
    V.ba.resize(len + 1);
    return U.karatsubaRecursive(V, 0, len);
}

const BN BN::karatsuba_old(const BN& bn)const {
    size_t x = rbc;
    size_t y = bn.rbc;

    size_t M = max(x,y);
    size_t n = (M + 1) / 2;
    if(min(x,y) < karacuba_const)
        return move(this->fast_mul(bn));

    const BN& U = *this;
    const BN& V = bn;

    BN u0(U, 0, n);
    BN v0(V, 0, n);

    BN u1(U, n, n);
    BN v1(V, n, n);

    BN A = u1.karatsuba_old(v1);
    BN B = u0.karatsuba_old(v0);
    BN C = (u0 + u1).karatsuba_old(v0+v1);

    return move(A.mulbt(2*n) + (C-A-B).mulbt(n) + B);
}

const BN BN::divbase(const bt& diviser) const
{
    if(diviser == 0)
        throw "Div by 0";
    BN result(rbc, 0);
    bt2 curr=0;
    for (size_t i = rbc - 1; i < rbc; --i) {
        curr <<= bz8;
        curr += ba[i];
        result.ba[i] = curr / diviser;
        curr %= diviser;
    }
    result.Norm();
    return result;
}

BN& BN::divbaseappr(const bt &diviser)
{
    if(diviser == 0)
            throw "Div by 0";
    bt2 curr = 0;
    for(size_t i = rbc - 1; i < rbc; --i) {
        curr <<= bz8;
        curr += ba[i];
        ba[i] = curr / diviser;
        curr %= diviser;
    }
    Norm();
    return *this;
}

const BN BN::modbase(const bt &diviser)const
{
    if(diviser == 0)
        throw "Div by 0";
    BN result(1, 0);
    bt2 curr=0;
    for (size_t i = rbc - 1; i < rbc; --i) {
        curr <<= bz8;
        curr += ba[i];
        curr %= diviser;
    }
    result.ba[0]=curr;
    result.rbc = 1;
    return result;
}

BN& BN::modbaseappr(const bt &diviser)
{
    if(diviser == 0)
        throw "Div by 0";
    bt2 curr = 0;
    for(size_t i = rbc - 1; i < rbc; --i) {
        curr <<= bz8;
        curr += ba[i];
        ba[i] = curr / diviser;
        curr %= diviser;
    }
    ba.resize(1);
    ba[0] = curr;
    rbc = 1;
    return *this;
}

void BN::subappr(const BN& bn)
{
    bool flag = 0;
    size_t pos = 0;

    for (; pos < bn.rbc; ++pos) {
        bt2s res = static_cast<bt2s>(ba[pos]) - bn.ba[pos] - flag;
        ba[pos] = static_cast<bt>(res);
        flag = (res < 0);
    }

    if (flag) {
        while (!ba[pos]--)
            ++pos;
    }

    Norm();
}

BN BN::sub(const BN& bn)const
{
    BN result(rbc, 0);

    bool flag = 0;
    size_t pos = 0;

    for (; pos < bn.rbc && pos < rbc; ++pos) {
        bt2s res = static_cast<bt2s>(ba[pos]) - bn.ba[pos] - flag;
        result.ba[pos] = static_cast<bt>(res);
        flag = (res < 0);
    }

    for (; flag && pos < rbc; ++pos) {
        result.ba[pos] = ba[pos] - 1;
        flag = (result.ba[pos] > ba[pos]);
    }

    for(;pos < rbc; pos++)
        result.ba[pos] = ba[pos];

    result.Norm();
    return result;
}

void BN::divmod(const BN& bn, BN& div, BN& mod) const
{
    if(bn.is0())
        throw "Div by 0";

    if(bn.rbc == 1) {
        div = move(this -> divbase(bn.ba[0]));
        mod = move(this -> modbase(bn.ba[0]));
        return;
    }

    if(*this < bn) {
        div = 0;
        mod = *this;
        return;
    }

    bt d = bsize  / (bn.ba[bn.rbc-1] + 1);

    BN delimoe(d == 1 ? *this : this->mulbase(d));
    BN delitel(d == 1 ? bn : bn.mulbase(d));

    delimoe.ba.resize(delimoe.ba.size() + 2);
    delitel.ba.resize(delitel.ba.size() + 1);

    size_t n = delitel.rbc;
    size_t m = delimoe.rbc - delitel.rbc + 1;

    div.ba.resize(m + 2);

    vector<bt> temp(n + 1);
    for (size_t j = m; j <= m; --j) {
        bt2 q = (delimoe.ba[j + n] * bsize + delimoe.ba[j + n - 1]) / delitel.ba[n-1];
        bt2 r = (delimoe.ba[j + n] * bsize + delimoe.ba[j + n - 1]) % delitel.ba[n-1];

        bool doAdjust = true;
        if (q == bsize || q * delitel.ba[n-2] > bsize * r + delimoe.ba[j + n - 2]) {
            --q;
            r += delitel.ba[n-1];
            if (q == bsize || (r < bsize && q * delitel.ba[n-2] > bsize * r + delimoe.ba[j + n - 2])) {
                --q;
                doAdjust = false;
            }
        }

        if (!q)
            continue;

        // Calculation temp = delitel * q
        bt2 x = 0;
        for(size_t i = 0; i < n; ++i) {
            x = (x >> bz8) + q * delitel.ba[i];
            temp[i] = x;
        }
        temp[n] = x >> bz8;

        bool qDecremented = false;
        if (doAdjust) {
            // temp1 * b^j <= delimoe ?

            size_t index = n;
            while (index <= n && temp[index] == delimoe.ba[index + j])
                --index;
            if (index <= n && temp[index] > delimoe.ba[index + j]) {
                --q;
                qDecremented = true;
            }
        }

        // delimoe = delimoe - temp * b^(i-n)
        {
            bt2s res = 0;
            size_t pos = 0;

            if (qDecremented)
                for (; pos <= n; ++pos) {
                    res = (res >> 8) + delimoe.ba[j + pos] + delitel.ba[pos] - temp[pos];
                    delimoe.ba[j + pos] = static_cast<bt>(res);
                }
            else
                for (; pos <= n; ++pos) {
                    res = (res >> 8) + delimoe.ba[j + pos] - temp[pos];
                    delimoe.ba[j + pos] = static_cast<bt>(res);
                }

            if (res) {
                while (!delimoe.ba[j + pos]--)
                    ++pos;
            }
        }

        div.ba[j] = q;
    }
    div.Norm();
    delimoe.Norm();
    if (d != 1)
        mod = move(delimoe.divbase(d));
    else
        mod = move(delimoe);
}

const BN BN::operator / (const BN&bn)const
{
    BN div;
    BN mod;
    divmod(bn, div, mod);
    return move(div);
}

const BN BN::operator % (const BN& bn)const
{
    BN div;
    BN mod;
    divmod(bn, div, mod);
    return move(mod);
}

const BN BN::operator >> (int shift) const {
    if(shift == 0)
        return *this;
    if(shift < 0)
        return (*this) >> (-shift);

    size_t baseshift = shift / bz8;
    size_t realshift = shift - baseshift * bz8;

    if (realshift == 0)
        return divbt(baseshift);

    if (baseshift >= rbc)
        return BN(1, 0);

    BN result(rbc - baseshift, 0);
    for (size_t i = 0; i < rbc - baseshift; ++i) {
        result.ba[i] =
            (ba[i + baseshift] >> realshift) |
            (ba[i + baseshift + 1] << (bz8 - realshift));
    }
    result.Norm();
    return result;
}

const BN BN::operator << (int shift) const {
    if(shift == 0)
        return *this;

    if (shift < 0)
        return *this >> (-shift);

    size_t baseshift = shift / bz8;
    size_t realshift = shift - baseshift * bz8;

    if (realshift == 0)
        return mulbt(baseshift);

    BN result(rbc + baseshift + 1, 0);
    result.ba[baseshift] = ba[0] << realshift;
    for(size_t i = 1; i <= rbc; i++) {
        result.ba[i + baseshift] =
            (ba[i - 1] >> (bz8 - realshift)) |
            (ba[i] << realshift);
    }
    result.Norm();
    return result;
}

bool BN::operator < (const BN&bn)const
{
    // TODO: replace to compare ba and bn.ba
    if(rbc > bn.rbc)
        return false;
    if(rbc < bn.rbc)
        return true;
    for(size_t i = rbc - 1; i < rbc; i--)
    {
        if(ba[i] > bn.ba[i])
            return false;
        if(ba[i] < bn.ba[i])
            return true;
    }
    return false;
}

bool BN::operator <= (const BN&bn)const
{
    return !(*this > bn);
}

bool BN::operator > (const BN&bn)const
{
    if(rbc > bn.rbc)
        return true;
    if(rbc < bn.rbc)
        return false;
    for(size_t i = rbc - 1; i < rbc; i--)
    {
        if(ba[i] > bn.ba[i])
            return true;
        if(ba[i] < bn.ba[i])
            return false;
    }
    return false;
}

bool BN::operator >= (const BN&bn)const
{
    return !(*this < bn);
}

bool BN::operator ==(const BN&bn)const
{
    if (rbc != bn.rbc)
        return false;
    for (size_t i = 0; i < rbc; ++i)
        if (ba[i] != bn.ba[i])
            return false;
    return true;
}

bool BN::operator !=(const BN&bn)const
{
    return !(*this == bn);
}

bt BN::operator [](size_t index)const
{
    return ba[index];
}

size_t BN::digitCount() const
{
    return rbc;
}

size_t BN::bitCount() const
{
    size_t x = 0;
    bt value = ba[rbc-1];
    while (value) {
        ++x;
        value >>= 1;
    }
    return (rbc - 1) * bz8 + x;
}

BN BN::reduction_barrett_precomputation() const {
    return ( (BN)1 ).mulbt(2*rbc) / *this;
}

BN BN::reduction_barrett(const BN& mod, const BN& mu) const {
    // m = m[k-1]...m[1]m[0], rbc = k
    size_t k = mod.rbc;
    if(k*2 < this->rbc) {
        return (*this)%mod;
    }

    BN x = *this;
    BN q1 = x.divbt(k-1);
    BN q2 = q1*mu;
    BN q3 = q2.divbt(k+1);
    BN r1 = this->modbt(k+1);
    BN r2 = (q3 * mod).modbt(k+1);
    BN r;
    if(r1 >= r2)
            r = r1 - r2;
        else
            r = ((BN) 1).mulbt(k+1) + r1 - r2;
    while(r >= mod)
        r = r - mod;
    return r;
}



BN BN::reduction_special(const BN &mod) const {
    size_t t = mod.rbc;

    // Bt = b^t;
    BN Bt(t+1,0);
    Bt.ba[t] = (bt) 1;
    Bt.rbc = t+1;

    BN c = Bt - mod;

    BN q = this -> divbt(t);
    BN r = *this - q.mulbt(t);

    BN r_sum = r;

    while(!q.is0()) {
        BN Q = (q*c).divbt(t);
        BN R = q*c - Q.mulbt(t);
        r_sum = r_sum + R;
        q = Q;
        r = R;
    }
    while(r_sum>=mod)
        r_sum=r_sum-mod;
    return r_sum;
}

BN BN::Pow(uint64_t power) const
{
    BN res(1);
    if (power == 0)
        return res;
    BN t = *this;
    do {
        if (power & 1)
            res = res * t;
        power >>= 1;
        if (power)
            t = t.Qrt();
    } while (power);
    return res;
}

BN BN::PowMod(uint64_t power, const BN& mod) const
{
    BN res(1);
    if (power == 0)
        return res;
    BN t = *this % mod;
    do {
        if (power & 1)
            res = res * t % mod;
        power >>= 1;
        if (power)
            t = t.Qrt() % mod;
    } while (power);
    return res;
}

BN BN::PowMod(const BN& power, const BN& mod) const {
    if(power.is0())
        return BN(1);

    BN Res(1);
    BN t = *this % mod;

    size_t len = power.bitCount();
    bt mask = 1;
    const bt *curr = &*power.ba.begin();
    for(size_t i = 0; i < len; i++) {
        if(!mask) {
            mask = 1;
            ++curr;
        }
        if((*curr) & mask)
            Res = Res * t % mod;
        if (i + 1 != len)
            t = t.Qrt() % mod;
        mask <<= 1;
    }
    return Res;
}

BN BN::PowModBarrett(const BN& power, const BN& mod) const {

    if(power.is0())
        return (BN) 1;

    BN mu = mod.reduction_barrett_precomputation();
    BN Res(1);
    BN t = (*this) % mod;

    int len = power.bitCount();
    bt mask = 1;
    const bt *curr = &*power.ba.begin();
    for(int i = 0; i < len; i++) {
        if(!mask) {
            mask = 1;
            ++curr;
        }
        if( (*curr) & mask)
            Res = (Res*t).reduction_barrett(mod, mu);

        if (i + 1 != len)
            t = t.Qrt().reduction_barrett(mod, mu);
        mask <<= 1;
    }
    return Res;
}


BN BN::expRightToLeft(const BN& exponent, const BN& mod)const {

    if(exponent.is0())
        return (BN) 1;

    BN A = 1;
    BN S = (*this) % mod;

    int exponent_len = exponent.bitCount();
    bt exponent_mask = (bt)  1;
    const bt *exponent_current_base = &*exponent.ba.begin();

    for(int i=0;i<exponent_len;i++) {
        if(!exponent_mask) {
            exponent_mask = (bt) 1;
            ++exponent_current_base;
        }

        if( (*exponent_current_base) & exponent_mask) {
            A = A * S % mod;
        }

        if (i + 1 != exponent_len)
            S = S.Qrt() % mod;
        exponent_mask <<= 1;
    }
    return A;
}

BN BN::expLeftToRight(const BN& exponent, const BN& mod) const {
    if(exponent.is0())
        return (BN) 1;

    BN A = 1;
    BN g = *this % mod;

    int exponent_len = exponent.bitCount();
    bt exponent_mask = (bt) 1;
    int start_shift_exponent_mask = (exponent_len - 1 ) % bz8;
    exponent_mask <<= start_shift_exponent_mask;
    const bt * exponent_current_base = &*exponent.ba.begin() + (exponent.rbc - 1);
    for(int i = 0; i < exponent_len; i++) {
        if(!exponent_mask) {
            exponent_mask = (bt) 1 << (bz8 - 1);
            --exponent_current_base;
            //printf("(%x)",*exponent_current_base);
        }

        if( (*exponent_current_base) & exponent_mask) {
            A = A.Qrt() % mod * g % mod;
        }
        else {
            A = A.Qrt() % mod;
        }

        exponent_mask >>= 1;
    }
    return A;

}

vector <BN> BN::expLeftToRightK_aryPrecomputation(const BN& mod) const {
    BN g = *this % mod;
    vector <BN> garr(bsize);
    garr[0] = BN(1);
    for(bt2 i = 1; i < bsize; i++) {
        garr[i] = garr[i-1] * g % mod;
    }
    return garr;
}

BN BN::expLeftToRightK_ary(const BN& exponent, const BN& mod, const vector<BN>& g) const {
    if(exponent.is0())
        return (BN) 1;

    BN A = 1;
    for(int i = exponent.rbc - 1; i >= 0; i--) {
        for(int k = 0; k < (int)bz8; k++)
            A = A.Qrt() % mod;
        A = A * g[exponent.ba[i]] % mod;
    }
    return A;
}

vector <BN> BN::expLeftToRightK_aryVarPrecomputation(const BN& mod, int K) const {
    int Kmax = (1 << K);
    BN g = *this % mod;
    vector <BN> garr(Kmax);
    garr[0] = BN(1);
    for(int i = 1; i < Kmax; i++) {
        garr[i] = garr[i-1] * g % mod;
    }
    return garr;
}

BN BN::expLeftToRightK_aryVar(BN exponent, BN mod, vector <BN> g, int K) const {
    if(exponent.is0())
        return (BN) 1;

    BN A = 1;

    int x = K;
    for(int i = exponent.rbc * bz8 - 1; i >= K; i -= K) {
        x = i;
        for(int k = 0; k < K; k++)
            A = A.Qrt() % mod;
        int curr = 0;
        for(int k = 0; k < K; k++) {
            curr <<= 1;
            curr |= exponent.bitI(i-k);
        }
        A = A * g[curr] % mod;
    }

    uint32_t curr = 0;
    for(int i = x - K; i >= 0; i--) {
        A = A.Qrt();
        curr <<= 1;
        curr |= exponent.bitI(i);
    }
    A = A * g[curr] % mod;

    return A;
}

vector <BN> BN::expLeftToRightK_aryModifPrecomputation(BN mod) const {
    BN g = *this % mod;

    vector <BN> garr(bsize);

    garr[0] = (BN) 1;
    garr[1] = g;
    garr[2] = g.Qrt() % mod;
    for(bt2 i = 1; i < bsize/2; i++)
        garr[2*i+1] = garr[2*i-1] * garr[2] % mod;
    return garr;
}

BN BN::expLeftToRightK_aryMod(BN exponent, BN mod, vector <BN> g) const {
    if(exponent.is0())
        return (BN) 1;

    BN A = (BN) 1;
    for(int i = exponent.rbc - 1; i >= 0; i--) {
        bt ei = exponent.ba[i];

        int hi = 0;
        if(ei != 0) {
            while(! (ei & 1)) {
                ei >>= 1;
                hi++;
            }
        }

        for(int k = 0; k < (int)bz8 - hi; k++)
            A = A.Qrt() % mod;
        A = A * g[ei] % mod;
        for(int k = 0; k < hi; k++)
            A = A.Qrt() % mod;
    }
    return A;

}

vector <BN> BN::expSlidingWindowPrecomputation(BN mod, int k) const {
    int k_pow = 2 << (k-1);
    vector <BN> garr (k_pow);
    BN g = *this % mod;
    garr[0] = (BN) 1;
    garr[1] = g;
    garr[2] = g.Qrt() % mod;
    for(int i = 1; i < k_pow/2; i++)
        garr[2*i+1] = garr[2*i-1] * garr[2] % mod;
    return garr;
}

BN BN::expSlidingWindow(BN exponent, BN mod, vector <BN> g, int k) const {
    BN A = (BN) 1;
    int i = exponent.bitCount() - 1;
    while (i >= 0) {
        if(exponent.bitI(i) == 0) {
            A = A.Qrt() % mod;
            i--;
            continue;
        }
        int l = max(i - k + 1, 0);
        while(exponent.bitI(l) == 0)
            l++;

        int gx = 0;
        for(int j = i; j >= l; j--)
            gx = (gx << 1) | exponent.bitI(j);
        for(int j = 0; j < i - l + 1; j++)
            A = A.Qrt() % mod;
        A = A * g[gx] % mod;
        i = l - 1;
    }
    return A;
}

vector <BN> BN::expBest_SlidePrecomp(BN mod) const {
    vector <BN> garr (bsize);
    BN mu = mod.reduction_barrett_precomputation();
    BN g = this -> reduction_barrett(mod, mu);
    garr[0] = (BN) 1;
    garr[1] = g;
    garr[2] = g.Qrt().reduction_barrett(mod,mu);
    for(bt2 i = 1; i < bsize/2; i++)
        garr[2*i+1] = (garr[2*i-1] * garr[2]).reduction_barrett(mod,mu);
    return garr;

}



BN BN::expBest_Slide(BN exponent, BN mod, vector <BN> g) const {
    BN A = (BN) 1;
    BN mu = mod.reduction_barrett_precomputation();
    int i = exponent.bitCount() - 1;
    int k = bz8;
    while (i >= 0) {
        if(exponent.bitI(i) == 0) {
            A = A.Qrt().reduction_barrett(mod, mu);
            i--;
            continue;
        }
        int l = max(i - k + 1, 0);
        while(exponent.bitI(l) == 0)
            l++;

        int gx = 0;
        for(int j = i; j >= l; j--)
            gx = (gx << 1) | exponent.bitI(j);
        for(int j = 0; j < i - l + 1; j++)
            A = A.Qrt().reduction_barrett(mod, mu);
        A = (A * (g[gx])).reduction_barrett(mod, mu);
        i = l - 1;
    }
    return A;

}


bt2 inverse(bt2 a, bt2 mod) {

    bt2 start_mod = mod;
    a = a % mod;
    if(a == 0)
        return 0;
    bt2s x0;
    bt2s x1 = 0;
    bt2s x2 = 1;
    while(mod % a) {
        bt2 q = mod/a;
        x0 = x1;
        x1 = x2;

        x2 = x0 - x1 * (bt2s) q;
        bt2 new_mod = a;
        a = mod % a;
        mod = new_mod;
    }

    if(a !=  1)
                return 0;
        if(x2 < 0)
                return start_mod + x2;
        return x2;
}


BN BN::Sqrt()const
{
    if(is0())
        return BN(1,0);
    BN x((this->rbc+1)/2,0);
    x.ba[(this->rbc+1)/2]=1;
    x.rbc=(this->rbc+1)/2+1;
    BN x0;
    do {
        x0=x;
        x=( (*this)/x+x )>>1;
    }
    while(x0>x);
    return x0;
}

BN BN::fastQrt() const
{
    size_t n = rbc;

    BN result(n + n, 0);

    bt4 t = 0;
    for(size_t s = 0; s < n + n - 1; s++) {

        size_t start_index = s >= n ? s - n + 1 : 0;
        size_t end_index = min(n - 1, s);
        while (start_index < end_index) {
            t += static_cast<bt4>(2 * static_cast<bt2>(ba[start_index]) * ba[end_index]);
            ++start_index;
            --end_index;
        }
        if (start_index == end_index)
            t += static_cast<bt2>(ba[start_index]) * ba[end_index];

        result.ba[s] = t;
        t = t >> bz8;
    }

    result.ba[n + n - 1] = t;
    result.Norm();
    return move(result);
}

BN BN::Qrt() const
{
    if (rbc < MaximalSizeForFastMul)
        return fastQrt();

    BN res(2 * rbc + 1, 0);
    for (size_t i = 0; i < rbc; ++i) {
        bt4 cuv = res.ba[2 * i] + static_cast<bt2>(ba[i]) * ba[i];
        res.ba[2 * i] = cuv;
        for (size_t j = i + 1; j < rbc; ++j) {
            cuv = static_cast<bt2>(res.ba[i + j]) +
                (static_cast<bt4>((static_cast<bt2>(ba[i]) * ba[j])) << 1) +
                (cuv >> bz8);
            res.ba[i + j] = cuv;
        }
        cuv = res.ba[i + rbc] + (cuv >> bz8);
        res.ba[i + rbc] = cuv;
        res.ba[i + rbc + 1] += (cuv >> bz8);
    }
    res.Norm();
    return res;
}

int BN::countzeroright()const
{
    if(ba[0] & 1)
            return 0;
    size_t count = 0;
    while(count < rbc && !ba[count])
            count++;
    if (count == rbc)
        return bz8 * (count - 1);

    bt last = ba[count];
    size_t d = 0;
    while(!(last & 1)) {
        ++d;
        last >>= 1;
    }
    return count * bz8 + d;
}

bool BN::bitI(size_t index)const {
    if(index >= bz8 * rbc)
        return false;

    bt mask = 1;
    mask <<= (index % bz8);

    if (ba[index / bz8] & mask)
        return true;
    return false;
}

BN::operator uint64_t () const {
    uint64_t result = 0;
    for(size_t i = min(rbc, sizeof(uint64_t)/bz); i; i--)
        result = (result << bz8) | ba[i-1];
    return result;
}

void BN::PrintHex(bool newstr)const
{
        for(int i=rbc-1;i>=0;i--)
                printf("%0*x",(int)bz*2,ba[i]);
        if(newstr)
                printf("\n");
}

void BN::PrintDec(bool newstr)const
{
        BN bn=*this;
        int slen=rbc*bz*4+1;                // длина 10-ного числа не более чем в 2 раза больше 16-ричного
        char *s=new char [slen];
        int count=0;
        if(bn.is0())
                s[count++]='0';
        while(!bn.is0())
        {
                s[count++]=(bn.modbase(10)).ba[0]+'0';
                bn.divbaseappr(10);
        }
        for(int i=count-1;i>=0;i--)
                printf("%c",s[i]);
        if(newstr)
                printf("\n");
        delete []s;
}

bool BN::is0() const
{
    if (rbc > 1 || ba[0])
        return false;
    return true;
}

bool BN::isEven() const
{
    if(ba[0] & 1)
        return false;
    return true;
}

