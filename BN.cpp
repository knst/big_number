/*
 * File:   BN.cpp
 * Author: knst
 *
 * Created on 27 Май 2010 г., 3:16
 */

#include "BN.h"

#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stack>
#include <stdexcept>

using namespace std;

// anonymous namespace for private methods and constants.
namespace {

// Maximal size for cache-optimazed multiplication is:
// n < bt4_max * / (bt - 1) * bt_max^3
constexpr bt4 MaximalSizeForFastMul = numeric_limits<bt4>::max() / bmax / bmax / bmax * (bmax - 1) - 1;

// Minimal size of multiplicand, than karatsuba algorithm is most effective.
constexpr size_t karatsubaMinimalSize = 50;

// Normalization of BN: pop leading null
inline void Norm(vector<bt>& ba) noexcept
{
    while (ba.size() > 1 && ba.back() == 0)
        ba.pop_back();
}

const BN reductionBarrettPrecomputation(const BN& mod) {
    size_t rbc = mod.digitCount() * 2 + 1;
    vector<bt> ba(rbc);
    ba.back() = 1;
    return move(BN(ba, rbc) / mod);
}

constexpr size_t KarySize = 256;
constexpr bt KaryMask = 0xFF;
constexpr size_t KaryBits = 8;

} // namespace

BN::BN()
: ba(1)
{
}

BN::BN(uint64_t basecount, int)
: ba(basecount)
{
}

BN::BN(uint64_t x)
: ba((sizeof(uint64_t) + bz - 1) / bz)
{
    for(size_t i = 0; i < ba.size(); i++, x >>= bz8)
        ba[i] = static_cast<bt>(x);
    Norm(ba);
}

BN::BN(const BN& bn)
: ba(bn.ba)
{
}

BN::BN(BN&& bn) noexcept
: ba(move(bn.ba))
{
}

BN::BN(vector<bt>&& _ba) noexcept
: ba(move(_ba))
{
    Norm(ba);
}

BN::BN(const vector<bt>& _ba)
: ba(_ba)
{
    Norm(ba);
}


BN::BN(vector<bt>&& _ba, size_t rbc)
: ba(move(_ba))
{
    ba.resize(rbc);
}

BN::BN(const vector<bt>& _ba, size_t rbc)
: ba(_ba.begin(), _ba.begin() + rbc)
{
}

BN::BN(const string &str,const int &status)
{
    if(status == 1)
    {
        BN bn(1,0);
        for (auto i : str) {
            if (i < '0' || i > '9')
                continue;
            bn.mulbaseappr(10);
            bn += BN(i - '0');
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
    Norm(ba);
}

void BN::swap(BN& bn) noexcept {
    ba.swap(bn.ba);
}

BN & BN::operator = (const BN& bn)
{
    if (this == &bn)
        return *this;
    ba = bn.ba;

    return *this;
}

BN & BN::operator = (BN&& bn) noexcept
{
    if (this == &bn)
        return *this;
    ba = move(bn.ba);

    return *this;
}

const BN BN::operator + (const BN& bn)const {
    return move(BN(*this) += bn);
}

BN& BN::operator += (const BN& bn) {
    ba.resize(max(ba.size(), bn.ba.size()));

    bt2 sum = 0;
    for(size_t pos = 0; pos < bn.ba.size(); pos++) {
        sum = (sum >> bz8) + ba[pos] + bn.ba[pos];
        ba[pos] = sum;
    }

    for(size_t pos = bn.ba.size(); pos < ba.size(); pos++) {
        sum = (sum >> bz8) + ba[pos];
        ba[pos] = sum;
    }
    sum >>= bz8;
    if (sum)
        ba.emplace_back(sum);
    return *this;
}


BN & BN::operator ++()
{
    for (auto& i : ba) {
        ++i;
        if (i != 0)
            return *this;
    }
    ba.emplace_back(1);
    return *this;
}

const BN BN::operator - (const BN& bn) const
{
    return move(BN(*this) -= bn);
}

BN& BN::operator -= (const BN& bn)
{
    if (ba.size() < bn.ba.size())
        ba.resize(bn.ba.size());

    bool flag = 0;
    size_t pos = 0;

    for (; pos < bn.ba.size() && pos < ba.size(); ++pos) {
        bt2s res = static_cast<bt2s>(ba[pos]) - bn.ba[pos] - flag;
        ba[pos] = static_cast<bt>(res);
        flag = (res < 0);
    }

    for (; flag && pos < ba.size(); ++pos) {
        if (ba[pos])
            flag = false;
        --ba[pos];
    }

    Norm(ba);
    return *this;
}

BN & BN::operator --() noexcept
{
    for (auto& i : ba) {
        if (i) {
            --i;
            return *this;
        }
        --i;
    }
    Norm(ba);

    return *this;
}

BN BN::mulbt(size_t t) const
{
    if(t == 0)
        return *this;

    BN res(ba.size() + t, 0);
    for(size_t i = 0; i < ba.size(); ++i)
        res.ba[i + t] = ba[i];
    Norm(res.ba);
    return res;
}

BN BN::divbt(size_t t) const
{
    if(t >= ba.size())
        return BN::bn0();

    return move(vector<bt>(ba.begin() + t, ba.end()));
}

BN BN::modbt(size_t t) const
{
    if(t >= ba.size())
        return *this;

    return move(vector<bt>(ba.begin(), ba.begin() + t));
}

const BN BN::mulbase(const bt &multiplier)const
{
    return move(BN(*this).mulbaseappr(multiplier));
}

BN& BN::mulbaseappr(const bt &multiplier)
{
    if (!multiplier) {
        ba.resize(1);
        ba.front() = 0;
        return *this;
    }

    bt2 curr = 0;
    for(auto& i : ba) {
        i = curr += static_cast<bt2>(i) * multiplier;
        curr >>= bz8;
    }
    if (curr)
        ba.emplace_back(curr);
    return *this;
}

const BN BN::operator * (const BN& bn)const {
    return karatsubaMultiplication(bn);
}

const BN BN::classicMultiplication(const BN& bn) const {
    // classical O(n*n) multiplication.
    // Tested: b * a is faster than a * b
    const BN& b = ba.size() > bn.ba.size() ? *this : bn;
    const BN& a = ba.size() > bn.ba.size() ? bn : *this;

    if (a.ba.size() == 1)
        return b.mulbase(a.ba.front());


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
    Norm(result.ba);
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
    Norm(result.ba);
    return move(result);
}

inline vector<bt> karatsubaSum(const vector<bt>& u, size_t start, size_t n, size_t m) {
    vector<bt> result;
    result.reserve(m + 1);

    bt2 sum = 0;
    for (size_t pos = 0; pos < n; ++pos) {
        sum += u[start + pos] + u[start + pos + n];
        result.emplace_back(sum);
        sum >>= bz8;
    }
    if (n != m) {
        sum += u[start + n + n];
        result.emplace_back(sum);
        sum >>= bz8;
    }
    result.emplace_back(sum);
    return move(result);
}


vector<bt> karatsubaRecursive(
    const vector<bt>& U,
    const vector<bt>& V,
    size_t start,
    size_t count
) {
    const size_t len = count;
    const size_t n = len / 2;
    const size_t m = len - n;
    if (n < karatsubaMinimalSize) {
        vector<bt> result;
        result.reserve(len + len);

        bt4 t = 0;
        for(size_t s = 0; s < len + len - 1; s++) {

            size_t end_index = min(len - 1, s) + start;
            size_t start_index = s >= len ? s - len + 1 : 0;
            for(size_t i = start_index + start, j = s - start_index + start; i <= end_index; i++, --j)
                t += static_cast<bt2>(U[i]) * V[j];


            result.emplace_back(t);
            t = t >> bz8;
        }
        result.emplace_back(t);
        return move(result);
    }

    const vector<bt>& u01 = karatsubaSum(U, start, n, m);
    const vector<bt>& v01 = karatsubaSum(V, start, n, m);

    vector<bt> A = move(karatsubaRecursive(U, V, start + n, m));
    vector<bt> B = move(karatsubaRecursive(U, V, start, n));
    const vector<bt> C = move(karatsubaRecursive(u01, v01, 0, m + 1));

    vector<bt> result = B;
    result.resize(len + len);

    for (size_t i = 0; i < A.size(); ++i)
        result[i + n + n] = A[i];

    const size_t abcSize = m + m + 2;
    A.resize(abcSize);
    B.resize(abcSize);

    bt2s sum = 0;
    for (size_t i = 0; i < abcSize; ++i) {
        sum += result[i + n];
        sum += C[i];
        sum -= A[i];
        sum -= B[i];
        result[i + n] = sum;
        sum >>= bz8;
    }
    for (size_t i = n + abcSize; i < result.size(); ++i) {
        sum += result[i];
        result[i] = sum;
        sum >>= bz8;
    }
    return move(result);
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

const BN BN::divbase(const bt& diviser) const
{
    return move(BN(*this).divbaseappr(diviser));
}

BN& BN::divbaseappr(const bt &diviser)
{
    if(diviser == 0)
            throw "Div by 0";
    bt2 curr = 0;
    for (size_t i = ba.size(); i; --i) {
        curr = (curr % diviser << bz8) + ba[i - 1];
        ba[i - 1] = curr / diviser;
    }
    Norm(ba);
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

    BN result;
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

    BN delimoe(*this);
    BN delitel(bn);

    bt d = bsize  / (bn.ba.back() + 1);
    if (d != 1) {
        delimoe.mulbaseappr(d);
        delitel.mulbaseappr(d);
    }

    size_t n = delitel.ba.size();
    size_t m = delimoe.ba.size() - n + 1;

    delimoe.ba.resize(delimoe.ba.size() + 2);
    delitel.ba.resize(delitel.ba.size() + 1);

    div.ba.resize(m + 1);

    vector<bt> temp(n + 1);
    for (size_t j = m; j <= m; --j) {
        bt2 q = (delimoe.ba[j + n] * bsize + delimoe.ba[j + n - 1]) / delitel.ba[n-1];
        bt2 r = (delimoe.ba[j + n] * bsize + delimoe.ba[j + n - 1]) % delitel.ba[n-1];

        if (q == bsize || q * delitel.ba[n-2] > bsize * r + delimoe.ba[j + n - 2]) {
            --q;
            r += delitel.ba[n-1];
            if (r < bsize && q * delitel.ba[n-2] > bsize * r + delimoe.ba[j + n - 2])
                --q;
        }

        if (!q) {
            div.ba[j] = 0;
            continue;
        }

        bt4s x = 0;
        for (size_t i = 0; i < n; ++i) {
            x += delimoe.ba[j + i];
            x -= q * delitel.ba[i];
            delimoe.ba[j + i] = x;
            x >>= bz8;
        }
        x += delimoe.ba[j + n];
        delimoe.ba[j + n] = x;

        // If `x' is negative, than `q' is too large.
        // Decrement `q' and update `delimoe'.
        if (x < 0) {
            --q;
            x = 0;
            for (size_t i = 0; i < n; ++i) {
                x += delimoe.ba[j + i];
                x += delitel.ba[i];
                delimoe.ba[j + i] = x;
                x >>= bz8;
            }
            x += delimoe.ba[j + n];
            delimoe.ba[j + n] = x;
        }

        div.ba[j] = q;
    }
    Norm(div.ba);
    Norm(delimoe.ba);

    if (d != 1)
        delimoe.divbaseappr(d);
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

BN& BN::operator >>= (size_t shift) {
    if (shift == 0)
        return *this;
    size_t baseshift = shift / bz8;
    size_t realshift = shift - baseshift * bz8;

    if (realshift == 0)
        return *this = divbt(baseshift);

    if (baseshift >= ba.size()) {
        ba.resize(0);
        ba[0] = 0;
        return *this;
    }

    for (size_t i = 0; i < ba.size() - baseshift - 1; ++i) {
        ba[i] =
            (ba[i + baseshift] >> realshift) |
            (ba[i + baseshift + 1] << (bz8 - realshift));
    }
    ba[ba.size() - baseshift - 1] = ba.back() >> realshift;
    ba.resize(ba.size() - baseshift);

    Norm(ba);
    return *this;
}

const BN BN::operator >> (size_t shift) const {
    return move(BN(*this) >>= shift);
}

BN& BN::operator <<= (size_t shift) {
    if(shift == 0)
        return *this;

    size_t baseshift = shift / bz8;
    size_t realshift = shift - baseshift * bz8;

    if (realshift == 0)
        return *this = mulbt(baseshift);

    ba.resize(ba.size() + baseshift + 1, 0);
    ba.back() = ba.back() >> (bz8 - realshift);
    for (size_t i = ba.size() - 1; i; --i) {
        ba[i + baseshift] =
            (ba[i - 1] >> (bz8 - realshift)) |
            (ba[i] << realshift);
    }
    ba[baseshift] = ba[0] << realshift;
    Norm(ba);
    return *this;
}

const BN BN::operator << (size_t shift) const {
    return move(BN(*this) <<= shift);
}

bool BN::operator < (const BN& bn) const noexcept
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

bool BN::operator <= (const BN&bn) const noexcept
{
    return !(*this > bn);
}

bool BN::operator > (const BN& bn) const noexcept
{
    return bn < *this;
}

bool BN::operator >= (const BN& bn) const noexcept
{
    return !(*this < bn);
}

bool BN::operator == (const BN& bn) const noexcept
{
    return ba == bn.ba;
}

bool BN::operator !=(const BN& bn) const noexcept
{
    return !(*this == bn);
}

bt BN::operator [](size_t index) const noexcept
{
    return ba[index];
}

size_t BN::digitCount() const noexcept
{
    return ba.size();
}

size_t BN::bitCount() const noexcept
{
    size_t x = 0;
    bt value = ba.back();
    while (value) {
        ++x;
        value >>= 1;
    }
    return (ba.size() - 1) * bz8 + x;
}

BN BN::reductionBarrett(const BN& mod, const BN& mu) const {
    // m = m[k-1]...m[1]m[0], rbc = k
    size_t k = mod.ba.size();
    if(k * 2 < ba.size()) {
        return (*this) % mod;
    }

    BN q1 = divbt(k-1);
    BN q2 = q1*mu;
    BN q3 = q2.divbt(k+1);
    BN r1 = modbt(k+1);
    BN r2 = (q3 * mod).modbt(k+1);
    r1 -= r2;
    while (r1 >= mod)
        r1 -= mod;
    return r1;
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
    return expRightToLeft(power, mod);
}

BN BN::PowModBarrett(const BN& power, const BN& mod) const {
    if(power.is0())
        return BN::bn1();


    BN mu = reductionBarrettPrecomputation(mod);
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
            res = (res*t).reductionBarrett(mod, mu);

        if (i + 1 != len)
            t = t.Qrt().reductionBarrett(mod, mu);
        mask <<= 1;
    }
    return res;
}


BN BN::expRightToLeft(const BN& exponent, const BN& mod)const {
    if(exponent.is0())
        return BN::bn1();

    BN result(BN::bn1());
    BN S = *this % mod;

    size_t len = exponent.bitCount();
    bt mask = 1;
    const bt *curr = &*exponent.ba.begin();
    for(size_t i = 0; i < len; i++) {
        if(!mask) {
            mask = 1;
            ++curr;
        }
        if(*curr & mask)
            result = result * S % mod;

        if (i + 1 != len)
            S = S.Qrt() % mod;
        mask <<= 1;
    }
    return result;
}

vector <BN> BN::expLeftToRightK_aryPrecomputation(const BN& mod) const {
    BN g = *this % mod;
    vector <BN> garr(KarySize);
    garr[0] = BN::bn1();
    for(size_t i = 1; i < KarySize; i++) {
        garr[i] = garr[i-1] * g % mod;
    }
    return garr;
}

BN BN::expLeftToRightK_ary(const BN& exponent, const BN& mod, const vector<BN>& g) const {
    if(exponent.is0())
        return BN::bn1();

    BN A(BN::bn1());
    for(int i = exponent.ba.size() - 1; i >= 0; i--) {
        bt value = exponent.ba[i];
        for (size_t b = bz - 1; b < bz; --b) {
            for(size_t k = 0; k < KaryBits; k++)
                A = A.Qrt() % mod;
            A = A * g[(value >> KaryBits * b) & KaryMask] % mod;
        }
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
        A = A.Qrt() % mod;
        curr <<= 1;
        curr |= exponent.bitI(i);
    }
    return A * g[curr] % mod;
}

vector <BN> BN::expLeftToRightK_aryModifPrecomputation(BN mod) const {
    BN g = *this % mod;

    vector<BN> garr(KarySize);

    garr[0] = BN::bn1();
    garr[1] = g;
    garr[2] = g.Qrt() % mod;
    for(size_t i = 1; i < KarySize / 2; i++)
        garr[2 * i + 1] = garr[2 * i - 1] * garr[2] % mod;
    return garr;
}

BN BN::expLeftToRightK_aryMod(BN exponent, BN mod, vector <BN> g) const {
    if(exponent.is0())
        return BN::bn1();

    BN A(BN::bn1());
    for(int i = exponent.ba.size() - 1; i >= 0; i--) {
        for (size_t b = bz - 1; b < bz; --b) {
            bt ei = (exponent.ba[i] >> KaryBits * b) & KaryMask;

            int hi = 0;
            if(ei != 0) {
                while(! (ei & 1)) {
                    ei >>= 1;
                    hi++;
                }
            }

            for(size_t k = 0; k + hi < KaryBits; k++)
                A = A.Qrt() % mod;
            A = A * g[ei] % mod;
            for(int k = 0; k < hi; k++)
                A = A.Qrt() % mod;
        }
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
    vector <BN> garr (KarySize);
    BN mu = reductionBarrettPrecomputation(mod);
    BN g = this -> reductionBarrett(mod, mu);
    garr[0] = BN::bn1();
    garr[1] = g;
    garr[2] = g.Qrt().reductionBarrett(mod,mu);
    for(bt2 i = 1; i < KarySize/ 2; i++)
        garr[2 * i + 1] = (garr[2 * i - 1] * garr[2]).reductionBarrett(mod,mu);
    return garr;

}



BN BN::expBest_Slide(BN exponent, BN mod, vector <BN> g) const {
    BN A(BN::bn1());
    BN mu = reductionBarrettPrecomputation(mod);
    int i = exponent.bitCount() - 1;
    int k = KaryBits;
    while (i >= 0) {
        if(exponent.bitI(i) == 0) {
            A = A.Qrt().reductionBarrett(mod, mu);
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
            A = A.Qrt().reductionBarrett(mod, mu);
        A = (A * (g[gx])).reductionBarrett(mod, mu);
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
        x += *this / x;
        x >>= 1;
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
            bt2 m = static_cast<bt2>(ba[start_index]) * ba[end_index];
            t += m;
            t += m;
            ++start_index;
            --end_index;
        }
        if (start_index == end_index)
            t += static_cast<bt2>(ba[start_index]) * ba[end_index];

        result.ba[s] = t;
        t = t >> bz8;
    }

    result.ba[n + n - 1] = t;
    Norm(result.ba);
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
            cuv = static_cast<bt4>(res.ba[i + j]) +
                ((static_cast<bt4>(ba[i]) * ba[j]) << 1) +
                (cuv >> bz8);
            res.ba[i + j] = cuv;
        }
        cuv = res.ba[i + ba.size()] + (cuv >> bz8);
        res.ba[i + ba.size()] = cuv;
        res.ba[i + ba.size() + 1] += (cuv >> bz8);
    }
    Norm(res.ba);
    return res;
}

int BN::countzeroright() const noexcept
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

bool BN::bitI(size_t index) const noexcept {
    if(index >= bz8 * ba.size())
        return false;

    bt mask = 1;
    mask <<= (index % bz8);

    if (ba[index / bz8] & mask)
        return true;
    return false;
}

uint64_t BN::get64() const noexcept {
    uint64_t result = 0;
    for(size_t i = min(ba.size(), sizeof(uint64_t) / bz) - 1; i < sizeof(uint64_t); i--)
        result = (result << bz8) | ba[i];
    return result;
}

string to_hexstring(const BN& bn) {
    string result;

    const auto& raw = bn.raw();
    for (auto i = raw.rbegin(); i != raw.rend(); ++i) {
        stringstream stream;
        stream << hex << setfill('0') << setw(bz * 2) << static_cast<uint32_t>(*i);
        string group = stream.str();
        result = result + group;
    }
    return result;
}

string to_string(BN bn) {
    stack<char> chars;
    do {
        chars.push(bn.modbase(10).get64() + '0');
        bn.divbaseappr(10);
    } while (!bn.is0());

    string s;
    s.reserve(chars.size() + 1);
    while (!chars.empty()) {
        s = s + chars.top();
        chars.pop();
    }
    return s;
}

const vector<bt> BN::raw() const noexcept {
    return ba;
}

bool BN::is0() const noexcept
{
    if (ba.size() > 1 || ba[0])
        return false;
    return true;
}

bool BN::isEven() const noexcept
{
    if(ba[0] & 1)
        return false;
    return true;
}

const BN BN::bn0() noexcept {
    static BN bn(0);
    return bn;
}

const BN BN::bn1() noexcept {
    static BN bn(1);
    return bn;
}

BN BN::makeRandom(size_t byteCount) {
    vector<bt> result(byteCount / bz);
    if (result.empty())
        result.emplace_back(rand() & bmax);
    else
        for (auto& i : result)
            i = rand() & bmax;
    return move(BN(move(result)));
}
