/*
 * File:   BNsign.cpp
 * Author: knst
 *
 * Created on 27 Май 2010 г., 3:18
 */

#include "BNsign.h"

#include "BN.h"

#include <cstdio>

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

BNsign::BNsign(BNsign&& bn)
: value(move(bn.value))
, sign(bn.sign)
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

BNsign& BNsign::operator = (BNsign&& bn)
{
    if (this == &bn)
        return *this;
    value = move(bn.value);
    sign = bn.sign;

    return *this;
}

const BNsign BNsign::operator + (const BNsign& bn) const
{
    bool firstIsGreater = (this->value >= bn.value);

    BNsign res;
    if(sign == bn.sign)
        res.value = this->value + bn.value;
    else {
        if (firstIsGreater)
            res.value = this->value - bn.value;
        else
            res.value = bn.value - this->value;
    }

    res.sign = ((bn.sign && !firstIsGreater) || (sign && firstIsGreater));
    return res;
}

const BNsign BNsign::operator - (const BNsign& bn) const
{
    bool firstIsGreater = (this->value >= bn.value);

    BNsign res;
    if (sign == bn.sign) {
        if (firstIsGreater)
            res.value = this->value - bn.value;
        else
            res.value = bn.value - this->value;
    } else
       res.value = this->value + bn.value;

    res.sign = !(bn.sign || firstIsGreater) || (sign && firstIsGreater);
    return res;
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


//debug
void BNsign::PrintSign() const
{
    if(sign)
        putchar('-');
    else
        putchar('+');
    value.PrintDec();
    return;
}
