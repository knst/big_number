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
    BNsign(const BN& bn, bool negative = false);
    BNsign(const BNsign &bn);
    BNsign & operator = (const BNsign&);
    const BNsign operator + (const BNsign&) const;
    const BNsign operator - (const BNsign&) const;
    const BNsign operator * (const BNsign&) const;
    //const BNsign operator / (const BNsign&) const;
    //bool operator == (const BNsign&) const;
    //bool operator != (const BNsign&) const;
    //bool operator >= (const BNsign&) const;

    void PrintSign() const;
};

#endif /* _BNSIGN_H */

