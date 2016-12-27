#include "BN.h"
#include "testing.h"

#include <ctime>
#include <iostream>

using namespace std;

int main(int argc, char**) {
    uint64_t t = clock();

    try {
        testingBN();
        testing();
        testingMul();
        testingExp();
        if (argc > 1)
            resulttest();
    } catch (const std::exception& exc) {
        cerr << exc.what() << endl;
        return -1;
    }
    cout.precision(3);
    cout << "Seconds: " << fixed << static_cast<double>(clock() - t) / CLOCKS_PER_SEC << endl;
    return 0;
}

