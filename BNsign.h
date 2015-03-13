/*
 * File:   BNsign.h
 * Author: knst
 *
 * Created on 27 Май 2010 г., 3:18
 */

#ifndef _BNSIGN_H
#define _BNSIGN_H

#include "BN.h"

class BNsign
{
public:
    BN value;

    // false: positive
    // true: negative
    bool sign;

public:
    BNsign();
    BNsign(const BNsign& bn);
    BNsign(BNsign&& bn) noexcept;
    BNsign(const BN& bn, bool negative = false);
    BNsign(BN&& bn, bool negative = false) noexcept;
    BNsign& operator = (const BNsign&);
    BNsign& operator = (BNsign&&) noexcept;

    const BNsign operator + (const BNsign&) const;
    BNsign& operator += (const BNsign& );

    const BNsign operator - (const BNsign&) const;
    BNsign& operator -= (const BNsign& );

    const BNsign operator * (const BNsign&) const;
    //const BNsign operator / (const BNsign&) const;
    //bool operator == (const BNsign&) const;
    //bool operator != (const BNsign&) const;
    //bool operator >= (const BNsign&) const;
};

std::string to_string(const BNsign& bn);
std::string to_hexstring(const BNsign& bn);

#endif /* _BNSIGN_H */

