/* 
 * File:   main.cpp
 * Author: knst
 *
 * Created on 27 Май 2010 г., 3:09
 */

#include <stdlib.h>
#include <iostream>
#include "BN.h"
#include "BNsign.h"
#include "testing.h"
#include "additional.h"
#include <cstdio>
#include <ctime>
using namespace std;
/*
 * 
 */

int main(int argc, char** argv) {
    BN x("4FEC11EDDFF42F604E79DA7E");
    BN y("5BE8E4A3D7E2908062D70280");
    BN m("DE8ED80881DC935AAE92C1DC");
    x.PowMod(y,m).PrintHex();
    BN z("9489D0CFF12992DB462CACF4");
    return 0;
    ull t = clock();
    try {
        testingExp();
        testingMul();
        testingBN();
        testing();
        resulttest(false);
    }
    catch (char const * str) {
        cout<<str;
        return -1;
    }
    cout << (clock() - t) << endl;
    return (EXIT_SUCCESS);
}

