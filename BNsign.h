/*
 * File:   BNsign.h
 * Author: knst
 *
 * Created on 27 Май 2010 г., 3:18
 */

#ifndef _BNSIGN_H
#define    _BNSIGN_H

#include "BN.h"
class BNsign
{
    public:
        BN value;
        bool sign;
        ~BNsign();
        BNsign();
        BNsign(const BN &bn,const bool bnsign=false); //false==positive
        BNsign(const BNsign &bn);
        BNsign & operator = (const BNsign&);
        BNsign   operator + (const BNsign&)const;
        BNsign   operator - (const BNsign&)const;
        BNsign   operator * (const BNsign&)const;
        //BNsign operator / (const BNsign&)const;
        //bool operator == (const BNsign&)const;
        //bool operator != (const BNsign&)const;
        //bool operator >= (const BNsign&)const;

        //debug:
        void PrintSign()const;

};

#endif    /* _BNSIGN_H */

