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

constexpr size_t karatsubaMinimalSize = 50;

void BN::InitMemory(int type)
{
    switch(type) {
    case 1: {
        ba.assign(sizeof(bmax), bmax);
        break;
    }
    case -1:
        for(size_t i = 0; i < ba.size(); ++i)
            ba[i] = rand() % bsize;
        Norm();
        break;

    default:
        throw "Unknow type (BN::BN)";
    }
}

void BN::Norm()
{
    while (ba.size() > 1 && ba.back() == 0)
        ba.pop_back();
}

BN::BN()
: ba(1)
{
}

BN::BN(uint64_t basecount, int type)
: ba(basecount)
{
    if (type == 1 || type == -1) {
        InitMemory(type);
        Norm();
    } else if (type != 0)
        throw invalid_argument("BN constructor: invalid type " + to_string(type));
}

BN::BN(uint64_t x)
: ba((sizeof(uint64_t) + bz - 1) / bz)
{
    for(size_t i = 0; i < ba.size(); i++, x >>= bz8)
        ba[i] = static_cast<bt>(x);
    Norm();
}

BN::BN(const BN& bn)
: ba(bn.ba)
{
}

BN::BN(BN&& bn)
: ba(move(bn.ba))
{
}

BN::BN(const vector<bt>& _ba, size_t _rbc)
: ba(_ba)
{
    if (_rbc)
        ba.resize(_rbc);
    else
        Norm();
}

BN::BN(const BN& bn, size_t start, size_t count)
: ba(count ? count : bn.ba.size() - start + 1)
{
    size_t last = min(count, bn.ba.size() - start);
    if (bn.ba.size() < start)
        last = 0;
    for(size_t i = 0; i < last; i++)
        ba[i] = bn.ba[i + start];
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
}

BN & BN::operator = (const BN& bn)
{
    if (this == &bn)
        return *this;
    ba = bn.ba;

    return *this;
}

BN & BN::operator = (BN&& bn)
{
    if (this == &bn)
        return *this;
    ba = move(bn.ba);

    return *this;
}

const BN BN::operator + (const BN&bn)const {
    const BN& a = ba.size() > bn.ba.size() ? *this : bn;
    const BN& b = ba.size() > bn.ba.size() ? bn : *this;

    vector<bt> result(a.ba.size() + 1);

    bt over = 0;
    for(size_t pos = 0; pos < b.ba.size(); pos++) {
        result[pos] = a.ba[pos] + b.ba[pos];
        bt over2 = (result[pos] < b.ba[pos]);

        result[pos] = result[pos] + over;
        bt over3 = (result[pos] < over);
        over = over2 + over3;
    }

    for(size_t pos = b.ba.size(); pos < a.ba.size(); pos++) {
        result[pos] = a.ba[pos] + over;
        over = (result[pos] < over);
    }

    result.back() = over;
    BN resultBn(move(result));
    resultBn.Norm();
    return resultBn;
}


BN & BN::operator ++()
{
    for (auto& i : ba) {
        ++i;
        if (i != 0)
            return *this;
    }
    ba.push_back(1);
    return *this;
}

const BN BN::operator - (const BN& bn) const
{
    BN result(ba.size(), 0);

    bool flag = 0;
    size_t pos = 0;

    for (; pos < bn.ba.size() && pos < ba.size(); ++pos) {
        bt2s res = static_cast<bt2s>(ba[pos]) - bn.ba[pos] - flag;
        result.ba[pos] = static_cast<bt>(res);
        flag = (res < 0);
    }

    for (; flag && pos < ba.size(); ++pos) {
        result.ba[pos] = ba[pos] - 1;
        flag = (result.ba[pos] > ba[pos]);
    }

    for(;pos < ba.size(); pos++)
        result.ba[pos] = ba[pos];

    result.Norm();
    return move(result);
}

BN & BN::operator --()
{
    for (auto& i : ba) {
        if (i) {
            --i;
            return *this;
        }
        --i;
    }
    Norm();

    return *this;
}

BN BN::mulbt(size_t t) const
{
    if(t == 0)
        return *this;

    BN res(ba.size() + t, 0);
    for(size_t i = 0; i < ba.size(); ++i)
        res.ba[i + t] = ba[i];
    res.Norm();
    return res;
}

BN BN::divbt(size_t t) const
{
    if(t == 0)
        return *this;

    if(t >= ba.size())
        return BN::bn0();

    BN res(ba.size() - t, 0);
    res.ba.assign(ba.begin() + t, ba.end());
    return res;
}

BN BN::modbt(size_t t) const
{
    if(t == 0)
        return BN::bn0();

    if(t >= ba.size())
        return *this;

    BN res(move(vector<bt>(ba.begin(), ba.begin() + t)));
    res.Norm();
    return res;
}

const BN BN::mulbase(const bt &multiplier)const
{
    BN result(ba.size() + 1, 0);
    bt2 curr = 0;
    for(size_t i = 0; i < ba.size(); i++, curr>>=bz8)
        result.ba[i] = curr += static_cast<bt2>(ba[i]) * multiplier;
    if (curr) {
        result.ba[ba.size()] = curr;
    } else
    result.Norm();
    return result;
}

BN& BN::mulbaseappr(const bt &multiplier)
{
    bt2 curr = 0;
    for(auto& i : ba) {
        i = curr += static_cast<bt2>(i) * multiplier;
        curr >>= bz8;
    }
    if (curr)
        ba.push_back(curr);
    return *this;
}

const BN BN::operator * (const BN& bn)const {
    if (min(bn.ba.size(), ba.size()) > karatsubaMinimalSize)
        return karatsubaMultiplication(bn);

    if (max(bn.ba.size(), ba.size()) < MaximalSizeForFastMul)
        return fastMultiplication(bn);

    return classicMultiplication(bn);
}

const BN BN::classicMultiplication(const BN& bn) const {
    // classical O(n*n) multiplication.
    if(bn.ba.size() == 1)
        return mulbase(bn.ba.front());

    if (ba.size() == 1)
        return bn.mulbase(ba.front());

    // Tested: b * a is faster than a * b
    const BN& b = ba.size() > bn.ba.size() ? *this : bn;
    const BN& a = ba.size() > bn.ba.size() ? bn : *this;

    BN result(a.ba.size() + b.ba.size(), 0);
    for (size_t i = 0; i < b.ba.size(); ++i) {
        bt2 curr = 0;
        bt2 x = b.ba[i];
        for (size_t j = 0; j < a.ba.size(); ++j) {
            curr = (curr >> bz8) + result.ba[i + j] + x * a.ba[j];
            result.ba[i + j] = curr;
        }
        result.ba[i + a.ba.size()] = curr >> bz8;
    }
    result.Norm();
    return result;
}

const BN BN::fastMultiplication (const BN& bn) const {
    size_t n = ba.size();
    size_t m = bn.ba.size();

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

    result.ba.back() = t;
    result.Norm();
    return move(result);
}

vector<bt> karatsubaSum(
    const vector<bt>& a,
    const vector<bt>& b
) {
    size_t n = a.size();
    size_t m = b.size();
    if (n != m)
        throw "wtf";

    vector<bt> result(n + 1);

    bt2 sum = 0;
    for(size_t pos = 0; pos < n; ++pos) {
        sum += a[pos] + b[pos];
        result[pos] = sum;
        sum >>= bz8;
    }

    result.back() = sum;
    return result;
}

vector<bt> karatsubaSum2(const vector<bt>& u, size_t start, size_t n, size_t m) {
    vector<bt> result(m + 1);

    bt2 sum = 0;
    for (size_t pos = 0; pos < n; ++pos) {
        sum += u[start + pos] + u[start + pos + n];
        result[pos] = sum;
        sum >>= bz8;
    }
    if (n != m) {
        sum += u[start + n + n];
        result[n] = sum;
        sum >>= bz8;
    }

    result[m] = sum;
    return move(result);
}


vector<bt> karatsubaRecursive(
    const vector<bt>& U,
    const vector<bt>& V,
    size_t start,
    size_t count
) {
    size_t len = count;
    size_t n = len / 2;
    size_t m = len - n;
    if (n < karatsubaMinimalSize) {
        vector<bt> a(U.begin() + start, U.begin() + start + count);
        vector<bt> b(V.begin() + start, V.begin() + start + count);
        BN A(a);
        BN B(b);
        return A.fastMultiplication(B).raw();
    }

    const vector<bt>& A = karatsubaRecursive(U, V, start + n, m);
    const vector<bt>& B = karatsubaRecursive(U, V, start, n);

    const vector<bt>& u01 = karatsubaSum2(U, start, n, m);
    const vector<bt>& v01 = karatsubaSum2(V, start, n, m);
    const vector<bt>& C = karatsubaRecursive(u01, v01, 0, m + 1);

    BN Ab(A);
    BN Bb(B);
    BN Cb(C);

    BN result = Bb + Ab.mulbt(n + n) + Cb.mulbt(n) - (Ab + Bb).mulbt(n);
    return result.raw();
}

const BN BN::karatsubaMultiplication(const BN& bn) const {
    size_t x = ba.size();
    size_t y = bn.ba.size();
    size_t len = max(x, y);
    if(min(x, y) < karatsubaMinimalSize)
        return this->fastMultiplication(bn);

    vector<bt> U = ba;
    vector<bt> V = bn.ba;
    U.resize(len);
    V.resize(len);
    return karatsubaRecursive(U, V, 0, len);
}

const BN BN::karatsuba_old(const BN& bn)const {
    size_t x = ba.size();
    size_t y = bn.ba.size();

    size_t M = max(x,y);
    size_t n = (M + 1) / 2;
    if(min(x,y) < karatsubaMinimalSize)
        return move(this->fastMultiplication(bn));

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
    BN result(ba.size(), 0);
    bt2 curr=0;
    for (size_t i = ba.size() - 1; i < ba.size(); --i) {
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
    for(size_t i = ba.size() - 1; i < ba.size(); --i) {
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
    bt2 curr=0;
    for (size_t i = ba.size()- 1; i < ba.size(); --i) {
        curr <<= bz8;
        curr += ba[i];
        curr %= diviser;
    }

    BN result = BN::bn0();
    result.ba[0]=curr;

    return result;
}

BN& BN::modbaseappr(const bt &diviser)
{
    if(diviser == 0)
        throw "Div by 0";
    bt2 curr = 0;
    for(size_t i = ba.size() - 1; i < ba.size(); --i) {
        curr <<= bz8;
        curr += ba[i];
        ba[i] = curr / diviser;
        curr %= diviser;
    }
    ba.resize(1);
    ba[0] = curr;
    return *this;
}

void BN::divmod(const BN& bn, BN& div, BN& mod) const
{
    if(bn.is0())
        throw "Div by 0";

    if(bn.ba.size() == 1) {
        div = move(this -> divbase(bn.ba.front()));
        mod = move(this -> modbase(bn.ba.front()));
        return;
    }

    if(*this < bn) {
        div = BN::bn0();
        mod = *this;
        return;
    }

    bt d = bsize  / (bn.ba.back() + 1);

    BN delimoe(d == 1 ? *this : this->mulbase(d));
    BN delitel(d == 1 ? bn : bn.mulbase(d));

    size_t n = delitel.ba.size();
    size_t m = delimoe.ba.size() - n + 1;

    delimoe.ba.resize(delimoe.ba.size() + 2);
    delitel.ba.resize(delitel.ba.size() + 1);

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

    if (baseshift >= ba.size())
        return BN::bn0();

    BN result(ba.size() - baseshift, 0);
    for (size_t i = 0; i < ba.size() - baseshift - 1; ++i) {
        result.ba[i] =
            (ba[i + baseshift] >> realshift) |
            (ba[i + baseshift + 1] << (bz8 - realshift));
    }
    result.ba.back() = ba.back() >> realshift;

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

    BN result(ba.size() + baseshift + 1, 0);
    result.ba[baseshift] = ba[0] << realshift;
    for(size_t i = 1; i < ba.size(); i++) {
        result.ba[i + baseshift] =
            (ba[i - 1] >> (bz8 - realshift)) |
            (ba[i] << realshift);
    }
    result.ba.back() = ba.back() >> (bz8 - realshift);
    result.Norm();
    return result;
}

bool BN::operator < (const BN& bn) const
{
    // TODO: maybe can be replaced to compare ba and bn.ba
    if(ba.size() > bn.ba.size())
        return false;
    if(ba.size() < bn.ba.size())
        return true;
    for(size_t i = ba.size() - 1; i < ba.size(); i--)
    {
        if(ba[i] > bn.ba[i])
            return false;
        if(ba[i] < bn.ba[i])
            return true;
    }
    return false;
}

bool BN::operator <= (const BN&bn) const
{
    return !(*this > bn);
}

bool BN::operator > (const BN& bn) const
{
    return bn < *this;
}

bool BN::operator >= (const BN& bn) const
{
    return !(*this < bn);
}

bool BN::operator ==(const BN& bn) const
{
    return ba == bn.ba;
}

bool BN::operator !=(const BN& bn) const
{
    return !(*this == bn);
}

bt BN::operator [](size_t index) const
{
    return ba[index];
}

size_t BN::digitCount() const
{
    return ba.size();
}

size_t BN::bitCount() const
{
    size_t x = 0;
    bt value = ba.back();
    while (value) {
        ++x;
        value >>= 1;
    }
    return (ba.size() - 1) * bz8 + x;
}

BN BN::reduction_barrett_precomputation() const {
    return BN::bn1().mulbt(2*ba.size()) / *this;
}

BN BN::reduction_barrett(const BN& mod, const BN& mu) const {
    // m = m[k-1]...m[1]m[0], rbc = k
    size_t k = mod.ba.size();
    if(k*2 < ba.size()) {
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
            r = BN::bn1().mulbt(k+1) + r1 - r2;
    while(r >= mod)
        r = r - mod;
    return r;
}



BN BN::reduction_special(const BN &mod) const {
    size_t t = mod.ba.size();

    // Bt = b^t;
    BN Bt(t+1,0);
    Bt.ba[t] = (bt) 1;

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
    if (power == 0)
        return BN::bn1();

    BN res(BN::bn1());
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
    if (power == 0)
        return BN::bn1();

    BN res(BN::bn1());
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
        return BN::bn1();

    BN res(BN::bn1());
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
            res = res * t % mod;
        if (i + 1 != len)
            t = t.Qrt() % mod;
        mask <<= 1;
    }
    return res;
}

BN BN::PowModBarrett(const BN& power, const BN& mod) const {
    if(power.is0())
        return BN::bn1();


    BN mu = mod.reduction_barrett_precomputation();
    BN res(BN::bn1());
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
            res = (res*t).reduction_barrett(mod, mu);

        if (i + 1 != len)
            t = t.Qrt().reduction_barrett(mod, mu);
        mask <<= 1;
    }
    return res;
}


BN BN::expRightToLeft(const BN& exponent, const BN& mod)const {

    if(exponent.is0())
        return BN::bn1();

    BN A(BN::bn1());
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
        return BN::bn1();

    BN A(BN::bn1());
    BN g = *this % mod;

    int exponent_len = exponent.bitCount();
    bt exponent_mask = 1;
    int start_shift_exponent_mask = (exponent_len - 1 ) % bz8;
    exponent_mask <<= start_shift_exponent_mask;
    const bt * exponent_current_base = &*exponent.ba.begin() + (exponent.ba.size() - 1);
    for(int i = 0; i < exponent_len; i++) {
        if(!exponent_mask) {
            exponent_mask = (bt) 1 << (bz8 - 1);
            --exponent_current_base;
            //printf("(%x)",*exponent_current_base);
        }

        if( (*exponent_current_base) & exponent_mask) {
            A = A.Qrt() % mod * g % mod;
        } else {
            A = A.Qrt() % mod;
        }

        exponent_mask >>= 1;
    }
    return A;

}

vector <BN> BN::expLeftToRightK_aryPrecomputation(const BN& mod) const {
    BN g = *this % mod;
    vector <BN> garr(bsize);
    garr[0] = BN::bn1();
    for(bt2 i = 1; i < bsize; i++) {
        garr[i] = garr[i-1] * g % mod;
    }
    return garr;
}

BN BN::expLeftToRightK_ary(const BN& exponent, const BN& mod, const vector<BN>& g) const {
    if(exponent.is0())
        return BN::bn1();

    BN A(BN::bn1());
    for(int i = exponent.ba.size() - 1; i >= 0; i--) {
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
    garr[0] = BN::bn1();
    for(int i = 1; i < Kmax; i++) {
        garr[i] = garr[i-1] * g % mod;
    }
    return garr;
}

BN BN::expLeftToRightK_aryVar(BN exponent, BN mod, vector <BN> g, int K) const {
    if(exponent.is0())
        return BN::bn1();

    BN A(BN::bn1());

    int x = K;
    for(int i = exponent.ba.size() * bz8 - 1; i >= K; i -= K) {
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

    vector<BN> garr(bsize);

    garr[0] = BN::bn1();
    garr[1] = g;
    garr[2] = g.Qrt() % mod;
    for(bt2 i = 1; i < bsize/2; i++)
        garr[2*i+1] = garr[2*i-1] * garr[2] % mod;
    return garr;
}

BN BN::expLeftToRightK_aryMod(BN exponent, BN mod, vector <BN> g) const {
    if(exponent.is0())
        return BN::bn1();

    BN A(BN::bn1());
    for(int i = exponent.ba.size() - 1; i >= 0; i--) {
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
    garr[0] = BN::bn1();
    garr[1] = g;
    garr[2] = g.Qrt() % mod;
    for(int i = 1; i < k_pow/2; i++)
        garr[2*i+1] = garr[2*i-1] * garr[2] % mod;
    return garr;
}

BN BN::expSlidingWindow(BN exponent, BN mod, vector <BN> g, int k) const {
    BN A(BN::bn1());
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
    garr[0] = BN::bn1();
    garr[1] = g;
    garr[2] = g.Qrt().reduction_barrett(mod,mu);
    for(bt2 i = 1; i < bsize/2; i++)
        garr[2*i+1] = (garr[2*i-1] * garr[2]).reduction_barrett(mod,mu);
    return garr;

}



BN BN::expBest_Slide(BN exponent, BN mod, vector <BN> g) const {
    BN A(BN::bn1());
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

BN BN::Sqrt() const
{
    if(is0())
        return BN::bn1();

    size_t rbc2 = (ba.size() + 1) / 2 + 1;
    BN x(rbc2, 0);
    x.ba.back() = 1;

    BN x0;
    do {
        x0 = x;
        x = ((*this)/x + x) >> 1;
    } while(x0 > x);
    return x0;
}

BN BN::fastQrt() const
{
    size_t n = ba.size();

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
    if (ba.size() < MaximalSizeForFastMul)
        return fastQrt();

    BN res(2 * ba.size() + 1, 0);
    for (size_t i = 0; i < ba.size(); ++i) {
        bt4 cuv = res.ba[2 * i] + static_cast<bt2>(ba[i]) * ba[i];
        res.ba[2 * i] = cuv;
        for (size_t j = i + 1; j < ba.size(); ++j) {
            cuv = static_cast<bt2>(res.ba[i + j]) +
                (static_cast<bt4>((static_cast<bt2>(ba[i]) * ba[j])) << 1) +
                (cuv >> bz8);
            res.ba[i + j] = cuv;
        }
        cuv = res.ba[i + ba.size()] + (cuv >> bz8);
        res.ba[i + ba.size()] = cuv;
        res.ba[i + ba.size() + 1] += (cuv >> bz8);
    }
    res.Norm();
    return res;
}

int BN::countzeroright()const
{
    if (ba[0] & 1)
        return 0;
    if (is0())
        return 0;

    size_t count = 0;
    while(!ba[count])
        ++count;

    bt last = ba[count];
    size_t result = count * bz8;
    while(!(last & 1)) {
        ++result;
        last >>= 1;
    }
    return result;
}

bool BN::bitI(size_t index)const {
    if(index >= bz8 * ba.size())
        return false;

    bt mask = 1;
    mask <<= (index % bz8);

    if (ba[index / bz8] & mask)
        return true;
    return false;
}

BN::operator uint64_t () const {
    uint64_t result = 0;
    for(size_t i = min(ba.size(), sizeof(uint64_t)/bz); i; i--)
        result = (result << bz8) | ba[i-1];
    return result;
}

void BN::PrintHex(bool newstr)const
{
        for(int i=ba.size() -1;i>=0;i--)
                printf("%0*x",(int)bz*2,ba[i]);
        if(newstr)
                printf("\n");
}

void BN::PrintDec(bool newstr)const
{
        BN bn=*this;
        int slen=ba.size() *bz*4+1;                // длина 10-ного числа не более чем в 2 раза больше 16-ричного
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

const vector<bt> BN::raw() const {
    return ba;
}

bool BN::is0() const
{
    if (ba.size() > 1 && ba.back() == 0)
        throw logic_error("wtf, ba.back() is 0");
    if (ba.size() > 1 || ba[0])
        return false;
    return true;
}

bool BN::isEven() const
{
    if(ba[0] & 1)
        return false;
    return true;
}

const BN BN::bn0() {
    static BN bn(1, 0);
    return bn;
}

const BN BN::bn1() {
    static BN bn(1, 0);
    bn.ba[0] = 1;
    return bn;
}
