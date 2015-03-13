/*
 * File:   BNsign.cpp
 * Author: knst
 *
 * Created on 27 Май 2010 г., 3:18
 */

#include "BNsign.h"

#include "BN.h"

using namespace std;

BNsign::BNsign()
: sign(false)
{
}

BNsign::BNsign(const BNsign& bn)
: value(bn.value)
, sign(bn.sign)
{
}

BNsign::BNsign(BNsign&& bn) noexcept
: value(move(bn.value))
, sign(bn.sign)
{
}

BNsign::BNsign(BN&& bn, const bool bnsign) noexcept
: value(move(bn))
, sign(bnsign)
{
}

BNsign::BNsign(const BN& bn, const bool bnsign)
: value(bn)
, sign(bnsign)
{
}

BNsign& BNsign::operator = (const BNsign& bn)
{
    if(this == &bn)
      return *this;

    value = bn.value;
    sign = bn.sign;

    return *this;
}

BNsign& BNsign::operator = (BNsign&& bn) noexcept
{
    if (this == &bn)
        return *this;
    value = move(bn.value);
    sign = bn.sign;

    return *this;
}

const BNsign BNsign::operator + (const BNsign& bn) const
{
    return move(BNsign(*this) += bn);
}

BNsign& BNsign::operator += (const BNsign& bn)
{
    bool firstIsGreater = (this->value >= bn.value);

    if(sign == bn.sign)
        value += bn.value;
    else {
        if (firstIsGreater)
            value -= bn.value;
        else
            value = bn.value - value;
    }

    sign = ((bn.sign && !firstIsGreater) || (sign && firstIsGreater));
    return *this;
}

const BNsign BNsign::operator - (const BNsign& bn) const
{
    return move(BNsign(*this) -= bn);
}

BNsign& BNsign::operator -= (const BNsign& bn)
{
    bool firstIsGreater = (this->value >= bn.value);
    bool newSign = !(bn.sign || firstIsGreater) || (sign && firstIsGreater);

    if (sign != bn.sign) {
        value += bn.value;
        sign = newSign;
        return *this;
    }

    sign = newSign;
    if (firstIsGreater) {
        value -= bn.value;
        return *this;
    }

    value = bn.value - value;
    return *this;
}

const BNsign BNsign::operator * (const BNsign& bn) const
{
    BNsign res;
    res.value = this->value * bn.value;
    res.sign = this->sign ^ bn.sign;
    return res;
}

//const BNsign BNsign::operator / (const BNsign&) const;
//bool BNsign::operator == (const BNsign&) const;
//bool BNsign::operator != (const BNsign&) const;
/*bool BNsign::operator >= (const BNsign& bn) const
{
   bool firstIsGreater = (this->value >= bn.value);
   return
      (sign && firstIsGreater) ||
         (!bn.sign && !firstIsGreater);
}*/


std::string signToChar(bool sign) {
    return sign ? "-" : "+";
}

std::string to_string(const BNsign& bn) {
    return signToChar(bn.sign) + to_string(bn.value);
}
std::string to_hexstring(const BNsign& bn) {
    return signToChar(bn.sign) + to_hexstring(bn.value);
}
