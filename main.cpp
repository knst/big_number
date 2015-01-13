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

int main(int argc, char** argv) {
    uint64_t t = clock();
    try {
        testingExp();
        testingMul();
        testingBN();
        testing();
        if (argc > 1)
            resulttest();
    } catch (char const * str) {
        cout<<str;
        return -1;
    }
    cout.precision(3);
    cout << "Seconds: " << fixed << static_cast<double>(clock() - t) / CLOCKS_PER_SEC << endl;
    return 0;
}

