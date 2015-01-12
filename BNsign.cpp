/* 
 * File:   BNsign.cpp
 * Author: knst
 * 
 * Created on 27 Май 2010 г., 3:18
 */

#include "BNsign.h"

#include "BN.h"
//#include <iostream>
#include <cstdio>
//
//using namespace std;

BNsign::~BNsign()
{
    return;
}

BNsign::BNsign()
{
	sign=false;
	return;
}

BNsign::BNsign(const BN &bn,const bool bnsign)
{
        value=bn;
	sign=bnsign;
	return;
}

BNsign::BNsign(const BNsign &bn)
{
	value=bn.value;
	sign=bn.sign;
	return;
}

BNsign & BNsign:: operator = (const BNsign& bn)
{
	if(this==&bn)
		return *this;

	value = bn.value;
	sign = bn.sign;
	return *this;
}

BNsign BNsign::operator + (const BNsign& bn)const
{
	BNsign res;
	if(sign == bn.sign)
		res.value = this->value + bn.value;
	else
		res.value = this->value - bn.value;
	bool first_is_greater = (this->value >= bn.value);
	res.sign = ((bn.sign && !first_is_greater) || (sign && first_is_greater));
	return res;
}

BNsign BNsign::operator - (const BNsign& bn)const
{
	BNsign res;
	if(sign == bn.sign)
		res.value = this->value - bn.value;
	else
		res.value = this->value + bn.value;
	bool first_is_greater = (this->value >= bn.value);
	res.sign = ((!bn.sign && !first_is_greater) || (sign && first_is_greater));
	return res;
}

BNsign BNsign::operator * (const BNsign& bn)const
{
	BNsign res;
	res.value = this->value * bn.value;
	res.sign = (this->sign ^ bn.sign);
	return res;
}

//BNsign BNsign::operator / (const BNsign&)const;
//bool BNsign::operator == (const BNsign&)const;
//bool BNsign::operator != (const BNsign&)const;
/*bool BNsign::operator >= (const BNsign& bn)const
{
	bool first_is_greater = (this->value >= bn.value);
	return
		(sign && first_is_greater) ||
			(!bn.sign && !first_is_greater);
}*/


//debug
void BNsign::PrintSign()const
{
    if(sign)
        putchar('-');
    else
        putchar('+');
    value.PrintDec();
    return;
}
