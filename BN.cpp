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

void BN::InitMemory(int type)
{
    switch(type) {
    case 1: {
        bt zero = 0;
        for(size_t i = 0; i < rbc; ++i)
            ba[i] = ~zero;
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

int BN::Norm()
{
    for(rbc = ba.size() - 1; rbc > 0 && ba[rbc] == 0; rbc--);
        rbc++;
    return rbc;
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
        throw std::invalid_argument("BN constructor: invalid type " + to_string(type));
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

BN::BN(const BN& bn, size_t start, size_t count) {
    if(!count)
        count = bn.rbc - start;

    rbc = count;
    ba.resize(rbc + 1);
    for(size_t i = 0; i < count && i + start < bn.rbc; i++)
        ba[i] = bn.ba[i + start];

    size_t bbstart = max<size_t>(0, bn.rbc - start);
    for(size_t i = bbstart; i < ba.size(); i++)
        ba[i] = 0;
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

        int length = str.size();
        rbc=(length+bz*2-1)/(bz*2);
        ba.resize(rbc + 1);
        InitMemory(0);
        bt z;
        for(int i=0;i<length;i++)
        {
                if(str[i]>='A'&&str[i]<='F')
                {
                        int x=(length-i-1)/(bz*2);
                        ba[x]<<=4;
                        z=str[i]-'A'+10;
                        ba[x]|=z;
                }
                if(str[i]>='a'&&str[i]<='f')
                {
                        int x=(length-i-1)/(bz*2);
                        ba[x]<<=4;
                        z=str[i]-'a'+10;
                        ba[x]|=z;
                }
                if(str[i]>='0'&&str[i]<='9')
                {
                        int x=(length-i-1)/(bz*2);
                        ba[x]<<=4;
                        z=str[i]-'0';
                        ba[x]|=z;
                }
        }
        Norm();
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
    size_t result_len = max(rbc, bn.rbc);
    BN result(result_len + 1, 0);

    bt2 res = 0;
    size_t pos = 0;
    size_t m = min(rbc, bn.rbc);
    for(; pos < m; pos++) {
        res += static_cast<bt2>(ba[pos]) + bn.ba[pos];
        result.ba[pos] = res;
        res >>= bz8;
    }

    for(; pos < rbc; pos++) {
        res += ba[pos];
        result.ba[pos] = res;
        res >>= bz8;
    }

    for(; pos < bn.rbc; pos++) {
        res += bn.ba[pos];
        result.ba[pos] = res;
        res >>= bz8;
    }

    result.ba[pos] = res;
    result.Norm();
    return result;
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

    BN res(t, 0);
    res.ba.assign(ba.begin(), ba.begin() + t);
    res.Norm();
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
    return result;
}

BN& BN::mulbaseappr(const bt &multiplier)
{
    bt2 curr = 0;
    if (rbc + 1 > ba.size()) {
        vector<bt> ba_new(rbc + 2);
        for (size_t i =0; i < rbc; ++i, curr >>= bz8)
            ba_new[i] = curr += ba[i] * multiplier;
        ba = ba_new;
    } else {
        for (size_t i =0; i < rbc; ++i, curr >>= bz8)
            ba[i] = curr += ba[i] * multiplier;
    }
    if (curr) {
        ba[rbc] = curr;
        ++rbc;
    }
    return *this;
}

const BN BN::operator * (const BN&bn)const {
    // maximal size for cache-optimazed multiplication is:
    // n < bt4_max * / (bt - 1) * bt_max^3

    constexpr bt4 maximalLen = numeric_limits<bt4>::max() / bmax / bmax / bmax * (bmax - 1) - 1;

    // If fast mul is possible, do It.
    if (bn.rbc < maximalLen && rbc < maximalLen)
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

    BN result(n + m + 1, 0);

    bt4 t = 0;
    for(size_t s = 0; s < m + n; s++) {

        size_t end_index = min(n, s);
        size_t start_index = s > m ? s - m + 1 : 0;
        for(size_t i = start_index; i <= end_index; i++) {
            t += static_cast<bt2>(ba[i]) * bn.ba[s-i];
        }

        result.ba[s] = t;
        t = t >> bz8;
    }

    result.ba[m+n] = t;
    result.Norm();
    return move(result);
}

BN BN::karatsuba_add(const BN & bn, int start_1, int count_1, int start_2, int count_2) const {
    const BN & bn1 = *this;
    const BN & bn2 = bn;

    size_t result_len = max(count_1, count_2);
    BN result(result_len + 1, 0);

    size_t max_1 = min<size_t>(count_1, bn1.rbc - start_1);
    size_t max_2 = min<size_t>(count_2, bn2.rbc - start_2);
    size_t m_min = min<size_t>(max_1, max_2);

    bt2 res = 0;
    size_t pos = 0;
    for(; pos < m_min; pos++) {
        res = res + (bt2) bn1.ba[start_1 + pos] + (bt2) bn2.ba[start_2 + pos];
        result.ba[pos] = res;
        res >>= bz8;
    }
    for(; pos < max_1; pos++) {
        res = res + (bt2) bn1.ba[start_1 + pos];
        result.ba[pos] = res;
        res >>= bz8;
    }

    for(; pos < max_2;pos++) {
        res = res + (bt2) bn2.ba[start_2 + pos];
        result.ba[pos] = res;
        res >>= bz8;
    }

    result.ba[pos] = res;
    for(pos++; pos <= result_len; pos++)
        result.ba[pos] = 0;

    result.Norm();
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


BN BN::karatsubaRecursive(const BN & bn, size_t start, size_t len) const {
    size_t n = len / 2;
    if (len / 2 < karacuba_const) {
        BN U(*this, start, len);
        BN V(bn, start, len);
        return U.fast_mul(V);
    }

    const BN & U = *this;
    const BN & V = bn;

    BN A = U.karatsubaRecursive(V, start + n, n);   // A = u1.caracuba(v1);
    BN B = U.karatsubaRecursive(V, start, n);       // B = u0.caracuba(v0);
    BN C = U.karatsuba_add(U, start, n, start + n, n).
            karatsuba
            (V.karatsuba_add(V, start, n, start + n, n));   //  (u0 + u1).caracuba(v0 + v1)

//    res = A.mulbt(2*n);
    BN res(B.ba, A.rbc + 2 * n);
    res.ba.resize(A.ba.size() + 2 * n);

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

    const BN & U = *this;
    const BN & V = bn;
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
    if (index >= rbc)
        throw logic_error("Error in BN::operator[]. Index is too large");
    return ba[index];
}

size_t BN::basecount() const
{
        return rbc;
}

size_t BN::bitcount() const
{
    size_t x = 0;
    bt value = ba[rbc-1];
    while (value) {
        ++x;
        value >>= 1;
    }
    return (rbc - 1) * bz8 + x;
}

BN BN::transformationMontgomery(const BN & mod, bt m1) const {
    int k = mod.rbc;
    BN y(rbc + k, 0);                     // выделяем память с запасом
    y = *this;
    for(int i = 0; i < k; i++) {
        bt u = y.ba[i] * m1;
        y.add_appr(mod.mulbase(u), i);
    }
    y = y.divbt(k);
    if(y >= mod)
        return y - mod;
    return y;
}

BN BN::mulMontgomery(const BN& bn, const BN& mod, bt m1) const {
//    вернуть проверки, если сделать функцию не приватной
//    if(gcdBinary(mod, (BN)bsize) != (BN) 1)
//        throw "montgomery: gcd != 1\n";
//    if(*this >= mod || bn >= mod)
//        throw "montgomery: *this >= mod || bn >= mod\n";

    BN A = 0;
    size_t n = mod.basecount();
    for(size_t i = 0; i < n && i < this->rbc; i++) {
        bt u = (A.ba[0] + this->ba[i] * bn[0]) * m1;
        A = (A + bn.mulbase(this->ba[i]) + mod.mulbase(u)).divbt(1);

    }
    for(size_t i = this->rbc; i < n; i++) {
        bt u = A.ba[0] * m1;
        A = (A + mod.mulbase(u)).divbt(1);
    }

    if(A >= mod)
        return A - mod;
    return A;
}

BN BN::reduction_montgomery(const BN& mod, bt m1, const BN& T) const {
    // TODO: not use, not testing
    if(gcdBinary(mod, (BN)bsize) != (BN) 1)
        throw "montgomery: gcd != 1\n";
    int n = mod.bitcount();
    if(T >= mod.mulbt(n))                       //mod * R; R = b^n;
        throw "montgomery: T >= mod * R\n";
    BN A = T;
    BN u(n,0);
    for(int i = 0; i < n; i++) {
        u.ba[i] = (A[i] * m1);
        A = A + mod.mulbase(u[i]).mulbt(i);
    }
    A = A.divbt(n);
    if(A >= mod)
        A = A - mod;
    return A;
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
        BN Res(1);
        if(power==(int64_t)0)
                return Res;
        BN t=(*this);
        while(power)
        {
                if(power&(int64_t)1)
                        Res=Res*t;
                t=t*t;
                power>>=(int64_t)1;
        }
        return Res;
}

BN BN::PowMod(uint64_t power, const BN& mod) const
{
        BN Res(1);
        if(power==(int64_t)0)
                return Res;
        BN t=(*this)%mod;
        while(power)
        {
                if(power&(int64_t)1)
                        Res=(Res*t)%mod;
                t=(t*t)%mod;
                power>>=(int64_t)1;
        }
        return Res%mod;
}

BN BN::PowMod(const BN& power, const BN& mod)const {

    if(power.is0())
        return (BN) 1;

    BN Res(1);
    BN t = (*this) % mod;

    int len = power.bitcount();
    bt mask = 1;
    const bt *curr = &*power.ba.begin();
    for(int i = 0; i < len; i++) {
        if(!mask) {
            mask = 1;
            ++curr;
        }
        if((*curr) & mask)
            Res = Res * t % mod;

        t=t.Qrt() % mod;
        mask <<= 1;
    }
    return Res % mod;
}

BN BN::PowModBarrett(const BN& power, const BN& mod) const {

    if(power.is0())
        return (BN) 1;

    BN mu = mod.reduction_barrett_precomputation();
    BN Res(1);
    BN t = (*this) % mod;

    int len = power.bitcount();
    bt mask = 1;
    const bt *curr = &*power.ba.begin();
    for(int i = 0; i < len; i++) {
        if(!mask) {
            mask = 1;
            ++curr;
        }
        if( (*curr) & mask)
            Res = (Res*t).reduction_barrett(mod, mu);

        t = t.Qrt().reduction_barrett(mod, mu);
        mask <<= 1;
    }
    return Res.reduction_barrett(mod, mu);
}


BN BN::expRightToLeft(const BN& exponent, const BN& mod)const {

    if(exponent.is0())
        return (BN) 1;

    BN A = 1;
    BN S = (*this) % mod;

    int exponent_len = exponent.bitcount();
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


        S = S.Qrt() % mod;
        exponent_mask <<= 1;
    }
    return A % mod;
}

BN BN::expLeftToRight(const BN& exponent, const BN& mod) const {
    if(exponent.is0())
        return (BN) 1;

    BN A = 1;
    BN g = *this % mod;

    int exponent_len = exponent.bitcount();
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
    for(int i = 1; i < bsize; i++) {
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
    for(int i = 1; i < bsize/2; i++)
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
    int i = exponent.bitcount() - 1;
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
    for(int i = 1; i < bsize/2; i++)
        garr[2*i+1] = (garr[2*i-1] * garr[2]).reduction_barrett(mod,mu);
    return garr;

}



BN BN::expBest_Slide(BN exponent, BN mod, vector <BN> g) const {
    BN A = (BN) 1;
    BN mu = mod.reduction_barrett_precomputation();
    int i = exponent.bitcount() - 1;
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


BN BN::expMontgomery(BN exponent, BN mod) const {
    if(gcdBinary(mod, (BN)bsize) != (BN) 1) {
        cout<<"%";
        return PowMod(exponent, mod);
        throw "expMontgomery: gcdBinary(mod, b) != 1\n";
    }
    bt mod1 = bsize - inverse(mod[0], bsize);

    BN R = ((BN) 1).mulbt(mod.basecount());
    BN x = *this % mod;
    BN x1 = x.mulMontgomery(R.Qrt() % mod, mod, mod1);

    BN A = R % mod;

    for(int i = exponent.bitcount(); i >= 0; i--) {
        A = A.Qrt().transformationMontgomery(mod, mod1);
        if(exponent.bitI(i))
            A = A.mulMontgomery(x1, mod, mod1);
    }
    A = A.transformationMontgomery(mod, mod1);
    return A;
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

BN BN::Qrt()const
{
        BN res(2*rbc+1,0);
        for(size_t i = 0; i < rbc; ++i) {
                bt4 cuv=res.ba[2*i]+((bt2)ba[i])*ba[i];
                res.ba[2*i]=cuv;
                for (size_t j = i + 1; j < rbc; ++j) {
                        cuv=(bt4)res.ba[i+j]+((bt4)((bt2)ba[i]*ba[j])<<1)+(cuv>>bz8);
                        res.ba[i+j]=cuv;
                }
                //(*((bt2*)(res.ba+i+rbc)))+=cuv>>bz8;
                cuv=(res.ba[i+rbc+1]<<bz8)+res.ba[i+rbc]+(cuv>>bz8);
                res.ba[i+rbc]=cuv;
                res.ba[i+rbc+1]=(cuv>>bz8);
        }
        res.Norm();
        return res;
}

int BN::countzeroright()const
{
        if(ba[0] & (bt)1)
                return 0;
        int count=0;
        while(!ba[count])
                count++;
        bt mask=1;
        int d=0;
        while(!(mask&ba[count]))
        {
                mask<<=1;
                d++;
        }
        return count*bz8+d;
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

