#include "testing.h"


#include <stdlib.h>
#include <iostream>
#include <vector>
#include <cstdio>

#include "testing.h"
#include "BN.h"
#include "additional.h"


int testingMul_ij(int max1,int max2, int i,int j) {
    BN bn1(rand()%max1+1,-1);
    BN bn2(rand()%max2+1,-1);
    BN mul = bn1*bn2;
    BN f_m = bn1.fast_mul(bn2);
    BN c_m = bn1.karatsuba(bn2);
    BN c_o = bn1.karatsuba_old(bn2);
    if(mul != f_m || mul != c_m || mul != c_o) {
        printf("bn1:\t");       bn1.PrintDec();
        printf("bn2:\t");       bn2.PrintDec();
        printf("bn1*bn2:\t");   mul.PrintDec();
        printf("bn1**bn2:\t");  f_m.PrintDec();
        printf("bn1c*bn2:\t");  c_m.PrintDec();
        printf("bn1cobn2:\t");  c_o.PrintDec();
        return 1;
    }
    return 0;
}

int testingExp_ij(int max1,int max2, int i,int j) {
    BN g(rand()%max1 + 1, -1);
    BN exp(rand()%max2 + 1, -1);
    BN mod(rand()%max2 + 1, -1);
    int vsize = rand() % (rand()%2 ? max1 : max2) + 1;
    if(vsize >  (int) sizeof(unsigned int) * 8 - 2)
        return 0;
    if(mod.is0())
        return 0;

    int k_slide = 8;          // TODO - rand!
    int K = 3;
    vector <BN> precomp = g.expLeftToRightK_aryPrecomputation(mod);
    vector <BN> precompVar = g.expLeftToRightK_aryVarPrecomputation(mod, K);
    vector <BN> precompModif = g.expLeftToRightK_aryModifPrecomputation(mod);
    vector <BN> precompSlide = g.expSlidingWindowPrecomputation(mod, k_slide);
    vector <BN> precompSlideU = g.expBest_SlidePrecomp(mod);
    vector <BN> precompFixWind = expFixedBaseWindowPrecomputation(g,exp.basecount()-1,mod);

    BN res1 = g.expLeftToRight(exp,mod);
    BN res2 = g.expRightToLeft(exp,mod);
    BN res3 = g.expLeftToRightK_ary(exp,mod,precomp);
    BN res3_1 = g.expLeftToRightK_aryVar(exp, mod, precompVar, K);
    BN res4 = g.expLeftToRightK_aryMod(exp,mod,precompModif);
    BN res5 = g.expSlidingWindow(exp, mod, precompSlide, k_slide);
    BN res6;
    if(gcdBinary(mod, (BN)bsize) == (BN) 1)
        res6 = g.expMontgomery(exp,mod);
    else {
        //cerr<<"!";
        res6 = res1;
    }

    BN res7 = expFixedBaseWindow(precompFixWind, exp, mod);
    BN res8;
    if(exp.basecount() < 3)
        res8 = res7;
    else
        res8 = expFixedBaseEuclidean(precompFixWind, exp, mod);

    int k = 5;
    BN res10 = expkaryStringReplacement(g,
            karyStringReplacementRepresentation(exp, k),
            mod, k);

    BN res12 = g.expBest_Slide(exp, mod, precompSlideU);

    if(
            res1 != res2 ||
            res1 != res3 ||
            res1 != res3_1 ||
            res1 != res4 ||
            res1 != res5 ||
            res1 != res6 ||
            res1 != res7 ||
//            res1 != res8 ||
            res1 != res10||
            res1 != res12) {
        printf("g:\t");         g.PrintDec();
        printf("exp:\t");       exp.PrintDec();
        printf("mod:\t");       mod.PrintDec();
        printf("bn1->bn2:\t");  res1.PrintDec();
        printf("bn1<-bn2:\t");  res2.PrintDec();
        printf("bn1=>bn2:\t");  res3.PrintDec();
        printf("bn1==bn2:\t");  res3_1.PrintDec();
        printf("bn1!>bn2:\t");  res4.PrintDec();
        printf("bn1>>bn2:\t");  res5.PrintDec();
        printf("bn1MMbn2:\t");  res6.PrintDec();
        printf("bn1FWbn2:\t");  res7.PrintDec();
        printf("bn1FEbn2:\t");  res8.PrintDec();
        printf("bn1SRbn2:\t");  res10.PrintDec();
        printf("bn1##bn2:\t");  res12.PrintDec();
        printf("precompFixWind:\n");
        for(vector <BN> :: iterator iter = precompFixWind.begin(); iter != precompFixWind.end(); iter++)
            iter->PrintDec();
        puts("");
        return 1;
    }
    return 0;
    if ( j % 5 == 0) {
        vector <BN> vg(vsize);
        vector <BN> ve(vsize);
        for(int i = 0; i < vsize; i++) {
            vg[i] = BN(rand() % max1 + 1, -1);
            ve[i] = BN(rand() % max2 + 1, -1);
        }
        vector <BN> G = expSimutaneousMulPrecomputation(vg,ve,mod);
        BN res1 = expSimutaneousMul(G,ve,mod);
        BN res2 = (BN) 1;
        for(int i = 0; i < vsize; i++) {
            res2 = res2 * vg[i].PowMod(ve[i], mod) % mod;
        }
        if(res1 != res2) {

            printf("res1:\t");  res1.PrintDec();
            printf("res2:\t");  res2.PrintDec();
            return 1;
        }

    }

    return 0;
}

int testingMul() {
    for(int i=0;i<=3;i++) {
        int max1=0;
        int max2=0;
        switch(i) {
            case 0:
                max1=3;
                max2=3;
                break;
            case 1:
                max1=100;
                max2=3;
                break;
            case 2:
                max1=3;
                max2=100;
                break;
            case 3:
                max1=100;
                max2=100;
                break;
        }
        for(int j=0;j<100000;j++) {
            int res=testingMul_ij(max1,max2,i,j);
            if(res != 0) {
                printf("Error int Mul-test: i=%d j=%d\ttest=%d\n",i,j,res);
                return res;
            }
        }
    }
    cout<<"Multiple test: OK" << endl;
    return 0;
}

int testingExp() {
    for(int i=0;i<=3;i++) {
        int max1=0;
        int max2=0;
        switch(i) {
            case 0:
                max1=3;
                max2=3;
                break;
            case 1:
                max1=100;
                max2=3;
                break;
            case 2:
                max1=3;
                max2=100;
                break;
            case 3:
                max1=100;
                max2=100;
                break;
        }
        int j_max = 20;
        for(int j=0;j< j_max;j++) {
            int res=testingExp_ij(max1,max2,i,j);
            if(j % 5 == 0)
                cerr << j*100/j_max << "% ...";
            if(res != 0) {
                printf("Error int Exp-test: i=%d j=%d\ttest=%d\n",i,j,res);
                return res;
            }
        }
        cerr<<"100 %" <<endl;
    }
    cerr << "Exp test: OK" << endl;
    return 0;
}

int testing_ij(int max1,int max2, int i,int j)
{
        BN bn_0 = 0;
        BN bn_1 = 1;
        BN bn1(rand()%max1+1,-1);
        BN bn2(rand()%max2+1,-1);
        BN mod(rand()%max2+1,-1);

        if(bn1!=bn_0&&bn2!=bn_0&&mod!=bn_0)
        {
                if((bn1+bn2-bn1)*bn1!=bn1*bn2)
                {
                        printf("bn1:\t");        bn1.PrintDec();
                        printf("bn2:\t");        bn2.PrintDec();
                        printf("bn1+bn2:\t");        (bn1+bn2).PrintDec();
                        printf("bn1*bn2:\t");        (bn1*bn2).PrintDec();
                        printf("f1:\t");        (bn1+bn2-bn1).PrintDec();
                        printf("f2:\t");        (bn1*bn2/bn1).PrintDec();
                        BN bn1bn2=bn1*bn2;
                        printf("bn1bn2:\t");        bn1bn2.PrintDec();
                        printf("f22:\t");        (bn1bn2/bn1).PrintDec();
                        return 1;
                }
                if((ull)(bn1*bn2/bn1/bn2)!=1)
                        return 2;
                if(bn1/bn2*bn2-bn1!=bn1%bn2+bn_0)
                        return 3;
                if((ull)(bn1*bn2)!=(ull)bn1*(ull)bn2)
                        return 4;
                if(bn1*bn1!=bn1.Qrt())
                        return 5;
                ull pow=(ull)bn1;
                BN powbn(pow);
                if(j % 10 == 0 && mod!=bn_0&&bn2.PowMod(pow,mod)!=bn2.PowMod(powbn,mod))
                {
                        printf("bn2:");bn2.PrintDec();
                        printf("pow:");printf("%lld\n",pow);
                        printf("pBN:");powbn.PrintDec();
                        printf("mod:");mod.PrintDec();
                        printf("Pll:");bn2.PowMod(pow,mod).PrintDec();
                        printf("Pbn:");bn2.PowMod(powbn,mod).PrintDec();
                        return 6;
        }
                BN gcd=gcdEuclidean(bn1,mod);
                if(gcd!=gcdBinary(bn1,mod))
                {
                        gcdEuclidean(bn1,mod).PrintDec();
                        gcdBinary(bn1,mod).PrintDec();
                        bn1.PrintDec();
                        mod.PrintDec();
                        return 7;
                }
                if(gcd!=gcdLehmer(bn1,mod))
                {
                        gcdEuclidean(bn1,mod).PrintDec();
                        gcdLehmer(bn1,mod).PrintDec();
                        bn1.PrintDec();
                        mod.PrintDec();
                        return 7;
                }
                BN gcdInverse=gcdInverseEuclidean(bn1,mod);
                BN gcdInverseBin=gcdInverseEuclideanBinary(bn1,mod);
                if(gcdInverse!=gcdInverseBin)
                {
                        printf("bn1:\t");        bn1.PrintDec();
                        printf("mod:\t");        mod.PrintDec();
                        printf("inverse:\t");        gcdInverse.PrintDec();
                        printf("inverseB:\t");        gcdInverseBin.PrintDec();
                        gcdExtendedEuclideanBinary(bn1,mod);//(.PrintDec();
                        return 8;
                }
                if(gcd!=bn_1||mod==bn_1)
                {
                        if(gcdInverse!=bn_0)
                        {
                                printf("bn1: ");        bn1.PrintDec();
                                printf("mod: ");        mod.PrintDec();
                                printf("inverse: ");        gcdInverseEuclidean(bn1,mod).PrintDec();
                                printf("nod: ");        gcdBinary(bn1,mod).PrintDec();
                                (gcdInverseEuclidean(bn1,mod)*bn1).PrintDec();
                                (gcdInverseEuclidean(bn1,mod)*bn1%mod).PrintDec();
                                return 9;
                        }
                }
        else
                {
                        if(gcdInverse*bn1%mod!=bn_1||gcdInverse>=mod)
                        {
                                printf("bn1: ");        bn1.PrintDec();
                                printf("mod: ");        mod.PrintDec();
                                printf("inverse: ");        gcdInverseEuclidean(bn1,mod).PrintDec();
                                printf("nod: ");        gcdBinary(bn1,mod).PrintDec();
                                (gcdInverseEuclidean(bn1,mod)*bn1).PrintDec();
                                (gcdInverseEuclidean(bn1,mod)*bn1%mod).PrintDec();
                                return 10;
                        }
                }
                if(i<2)
                {
                        BN sqrt = bn1.Sqrt();
                        BN sqrt1 = sqrt;
                        ++sqrt1;
                        BN odin(1);
                        if(sqrt*sqrt>bn1||sqrt1*sqrt1<=bn1)
                        {
                                sqrt.PrintDec();
                                bn1.PrintDec();
                                return 11;
                        }
        }
        }
        else
                if(mod!=bn_0)
                {
                        if(bn1*bn2!=bn_0||bn1+bn2!=max(bn1,bn2))
                                return 11;
                }
                else if(bn1*mod!=bn_0||bn1+mod!=max(bn1,mod))
                        return 12;
        return 0;
}

int testingBN()
{
        for(int i=0;i<=3;i++)
        {
                int max1=0;
                int max2=0;
                switch(i) {
                    case 0:
                        max1=3;
                        max2=3;
                        break;
                    case 1:
                        max1=100;
                        max2=3;
                        break;
                    case 2:
                        max1=3;
                        max2=100;
                        break;
                    case 3:
                        max1=100;
                        max2=100;
                        break;
                }
                for(int j=0;j<1000;j++)
                {
                        if(j%250==0)
                                cout<<"i="<<i<<"\tj="<<j<<'\n';
                        int res=testing_ij(max1,max2,i,j);
                        if(res!=0)
                        {
                                printf("Error: i=%d j=%d\ttest=%d\n",i,j,res);
                                return res;
                        }
                }
        }
        return 0;
}

void testing() {
        BN bn0(1,0);
        BN bn1(15000,-1);
        BN bn2(15000,-1);
        std::cout<<"Test 1:\t";
        cout.flush();
        if(bn1+bn2-bn1!=bn1*bn2/bn1)
                cout<<"FAIL\n";
        else
                cout<<"OK\n";
        cout<<"Test 2:\t";
        cout.flush();
        if((ull)(bn1*bn2/bn1/bn2)!=1)
                cout<<"FAIL\n";
        else
                cout<<"OK\n";
        cout<<"Test 3:\t";
        cout.flush();
        if(bn1/bn2*bn2-bn1!=bn1%bn2+bn0)
                cout<<"FAIL\n";
        else
                cout<<"OK\n";
        cout<<"Test 4:\t";
        cout.flush();
        if((ull)(bn1*bn2)!=(ull)bn1*(ull)bn2)
                cout<<"FAIL\n";
        else
                cout<<"OK\n";
        cout<<"Test sqrt:\t";
        cout.flush();
        BN sqrt;
        BN bnn(400,-1);
        sqrt=bnn.Sqrt();
        BN sqrt1=sqrt;
        ++sqrt1;
        if(bnn-sqrt*sqrt>bnn||(++sqrt)*(++sqrt)<=bnn)
                cout<<"FAIL\n";
        else
                cout<<"OK\n";
        cout<<"Test qrt:\t";
        cout.flush();
        if(bn1*bn1!=bn1.Qrt())
                cout<<"FAIL\n";
        else
                cout<<"OK\n";
        cout.flush();
        cout<<"Test gcd:\t";
        cout.flush();
        BN bnn1(2500,-1);
        BN bnn2(2500,-1);
        BN bngcd1=gcdEuclidean(bnn1,bnn2);
        cout<<"HaBePHo OK :)"<<endl;
        cout<<"Test binary gcd:\t";
        cout.flush();
        BN bngcd2=gcdBinary(bnn1,bnn2);
        if(bngcd1==bngcd2)
                cout<<"OK\n";
        else
                cout<<"FAIL\n";

        cout<<"Exiting...\n";
        cout.flush();

}

void multest(int base,int test)
{
        ull t;
        vector <BN> v1;
        vector <BN> v2;

        for(int i=0;i<test;i++) {
                v1.push_back( BN(base,-1) );
                v2.push_back( BN(base,-1) );
        }

        t = clock();
        for(vector <BN> :: iterator i = v1.begin(), j = v2.begin(); i != v1.end(); i++, j++)
                (*i) * (*j);
        float t1 = clock() - t;

        t = clock();
        for(vector <BN> :: iterator i = v1.begin(), j = v2.begin(); i != v1.end(); i++, j++)
                (*i).fast_mul(*j);
        float t2 = clock() - t;

//        t = clock();
//        for(vector <BN> :: iterator i = v1.begin(), j = v2.begin(); i != v1.end(); i++, j++)
//            (*i).karatsuba_old(*j);
//        float t3 = clock() - t;

        t = clock();
        for(vector <BN> :: iterator i = v1.begin(), j = v2.begin(); i != v1.end(); i++, j++)
            (*i).karatsuba(*j);
        float t3 = clock() - t;

        float diviser = test;
        t1 /= diviser;
        t2 /= diviser;
        t3 /= diviser;

        printf("%d\t%d\t%d\t%.2f\t\t%.2f\t\t%.2f\n",base,(int)(base*sizeof(bt)*8),(int)test, t1, t2, t3);
}

void modtest(int base,int test) {
    ull t;
    vector <BN> g;
    vector <BN> exp;
    vector <BN> mod;


    for(int i=0;i<test;i++) {
        g.push_back( BN(base,-1) );
        exp.push_back( BN(base,-1) );
        BN m = BN(base,-1);
        while (m.is0() || (gcdBinary(m, (BN)bsize) != (BN) 1))
            m = BN(base,-1);
        mod.push_back( m );
    }

    t = clock();
    for(vector <BN> :: iterator i = g.begin(), j = exp.begin(), k = mod.begin(); i != g.end(); i++, j++, k++)
        i->PowMod(*j, *k);
    float t1 = clock() - t;

    t = clock();
    for(vector <BN> :: iterator i = g.begin(), j = exp.begin(), k = mod.begin(); i != g.end(); i++, j++, k++)
        i->PowModBarrett(*j, *k);
    float t2 = clock() - t;

    t = clock();
    for(vector <BN> :: iterator i = g.begin(), j = exp.begin(), k = mod.begin(); i != g.end(); i++, j++, k++)
        i->expMontgomery(*j, *k);
    float t3 = clock() - t;

    float diviser = 1000.0 * test;
    t1 /= diviser;
    t2 /= diviser;
    t3 /= diviser;

    printf("%d\t%d\t%d\t%.2f\t\t%.2f\t\t%.2f\n",base,(int)(base*sizeof(bt)*8),test, t1, t2, t3);
}

void unitest(int base,int test) {
    ull t;

    BN g(base, -1);
    vector <BN> exp;
    BN mod(base, -1);


    for(int i=0;i<test;i++) {
        exp.push_back( BN(base,-1) );
    }

    t = clock();
    for(vector <BN> :: iterator e = exp.begin(); e != exp.end(); e++)
        g.expLeftToRight(*e, mod);
    float t1 = clock() - t;

    t = clock();
    for(vector <BN> :: iterator e = exp.begin(); e != exp.end(); e++)
        g.expRightToLeft(*e, mod);
    float t2 = clock() - t;

    t = clock();
    vector <BN> precomp = g.expLeftToRightK_aryPrecomputation(mod);
    for(vector <BN> :: iterator e = exp.begin(); e != exp.end(); e++)
        g.expLeftToRightK_ary(*e, mod, precomp);
    float t3 = clock() - t;

    t = clock();
    vector <BN> precompModif = g.expLeftToRightK_aryModifPrecomputation(mod);
    for(vector <BN> :: iterator e = exp.begin(); e != exp.end(); e++)
        g.expLeftToRightK_aryMod(*e, mod, precompModif);
    float t4 = clock() - t;

    int k = 8;
    t = clock();
    vector <BN> precompSlide = g.expSlidingWindowPrecomputation(mod, k);
    for(vector <BN> :: iterator e = exp.begin(); e != exp.end(); e++)
        g.expSlidingWindow(*e, mod, precompSlide, k);
    float t5 = clock() - t;

    float diviser = 1000.0 * test;
    t1 /= diviser;
    t2 /= diviser;
    t3 /= diviser;
    t4 /= diviser;
    t5 /= diviser;

    printf("%d\t%d\t%d\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\n",base,(int)(base*sizeof(bt)*8),test, t1, t2, t3, t4, t5);
}

void karytest(int K, int base,int test) {
    ull t;

    BN g(base, -1);
    vector <BN> exp;
    BN mod(base, -1);

    for(int i=0;i<test;i++) {
        exp.push_back( BN(base,-1) );
    }

    t = clock();
    if(1 || K == 8) {
        vector <BN> precomp = g.expLeftToRightK_aryPrecomputation(mod);
        for(vector <BN> :: iterator e = exp.begin(); e != exp.end(); e++)
            g.expLeftToRightK_ary(*e, mod, precomp);
    }
    float t1 = clock() - t;

    t = clock();
    vector <BN> precompVar = g.expLeftToRightK_aryVarPrecomputation(mod, K);
    for(vector <BN> :: iterator e = exp.begin(); e != exp.end(); e++)
        g.expLeftToRightK_aryVar(*e, mod, precompVar, K);
    float t2 = clock() - t;

    float diviser = 1000.0 * test;
    t1 /= diviser;
    t2 /= diviser;

    printf("%d\t%d\t%d\t%.2f\t\t%.2f\n",base,(int)(base*sizeof(bt)*8),test, t1, t2);
}

void slidetest(int base,int test) {
    ull t;

    BN g(base, -1);
    vector <BN> exp;
    BN mod(base, -1);


    for(int i=0;i<test;i++) {
        exp.push_back( BN(base,-1) );
    }

    int k;
    vector <BN> precompSlide;

    k = 4;
    t = clock();
    precompSlide = g.expSlidingWindowPrecomputation(mod, k);
    for(vector <BN> :: iterator e = exp.begin(); e != exp.end(); e++)
        g.expSlidingWindow(*e, mod, precompSlide, k);
    float t1 = clock() - t;

    k = 8;
    t = clock();
    precompSlide = g.expSlidingWindowPrecomputation(mod, k);
    for(vector <BN> :: iterator e = exp.begin(); e != exp.end(); e++)
        g.expSlidingWindow(*e, mod, precompSlide, k);
    float t2 = clock() - t;

    k = 9;
    t = clock();
    precompSlide = g.expSlidingWindowPrecomputation(mod, k);
    for(vector <BN> :: iterator e = exp.begin(); e != exp.end(); e++)
        g.expSlidingWindow(*e, mod, precompSlide, k);
    float t3 = clock() - t;

    k = 10;
    t = clock();
    precompSlide = g.expSlidingWindowPrecomputation(mod, k);
    for(vector <BN> :: iterator e = exp.begin(); e != exp.end(); e++)
        g.expSlidingWindow(*e, mod, precompSlide, k);
    float t4 = clock() - t;

    k = 11;
    t = clock();
    precompSlide = g.expSlidingWindowPrecomputation(mod, k);
    for(vector <BN> :: iterator e = exp.begin(); e != exp.end(); e++)
        g.expSlidingWindow(*e, mod, precompSlide, k);
    float t5 = clock() - t;

    float diviser = 1000.0;
    t1 /= diviser;
    t2 /= diviser;
    t3 /= diviser;
    t4 /= diviser;
    t5 /= diviser;

    printf("%d\t%d\t%d\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\n",base,(int)(base*sizeof(bt)*8),test, t1, t2, t3, t4, t5);
}

void simtest(int base,int test) {
    ull t;

    vector <BN> vg(test);
    vector <BN> ve(test);
    BN mod(base, -1);

    for(int i = 0; i < test; i++) {
        vg[i] = BN(base, -1);
        ve[i] = BN(base, -1);
    }

    t = clock();
    BN res2 = (BN) 1;
    for(int i = 0; i < test; i++) {
        res2 = res2 * vg[i].PowMod(ve[i], mod) % mod;
    }
    float t1 = clock() - t;

    t = clock();
    vector <BN> G = expSimutaneousMulPrecomputation(vg,ve,mod);
    expSimutaneousMul(G,ve,mod);
    float t2 = clock() - t;

    float diviser = 1000.0 * test;
    t1 /= diviser;
    t2 /= diviser;

    printf("%d\t%d\t%d\t%.2f\t\t%.2f\n",base,(int)(base*sizeof(bt)*8), test, t1, t2);
}

void fixtest(int base,int test) {
    ull t;

    BN g(base, -1);
    vector <BN> exp;
    BN mod(base, -1);

    for(int i=0;i<test;i++) {
        BN e = BN(base, -1);
        while(e.bitcount() != (int)(base*bz8))
            e = BN(base, -1);
        exp.push_back( e );
    }

    t = clock();
    BN res = 1;
    for(vector <BN> :: iterator e = exp.begin(); e != exp.end(); e++)
        res = res * g.expLeftToRight(*e, mod);
    float t1 = clock() - t;

    t = clock();
    vector <BN> precompFixWind = expFixedBaseWindowPrecomputation(g, base - 1,mod);
    float  t_2_3 = clock() - t;

    t = clock();
    for(vector <BN> :: iterator e = exp.begin(); e != exp.end(); e++)
        expFixedBaseWindow(precompFixWind, *e, mod);
    float t2 = clock() - t + t_2_3;

    t = clock();
    for(vector <BN> :: iterator e = exp.begin(); e != exp.end(); e++)
        expFixedBaseEuclidean(precompFixWind, *e, mod);
    float t3 = clock() - t + t_2_3;

    float diviser = 1000.0 * test;
    t1 /= diviser;
    t2 /= diviser;
    t3 /= diviser;

    printf("%d\t%d\t%d\t%.2f\t\t%.2f\t\t%.2f\n",base,(int)(base*sizeof(bt)*8), test, t1, t2, t3);
}


void exptest(int base,int test) {
    ull t;

    vector <BN> g;
    vector <BN> exp;
    vector <BN> mod;

    for(int i=0;i<test;i++) {
        exp.push_back(BN(base, -1));
        BN m(base, -1);
        BN G(base, -1);
        while(gcdBinary(G, m) != (BN)1) {
            m = BN(base, -1);
            G = BN(base, -1);
        }
        g.push_back( G );
        mod.push_back( m );
    }

    t = clock();
    for(vector <BN> :: iterator i = g.begin(), j = exp.begin(), k = mod.begin(); i != g.end(); i++, j++, k++)
        i->PowMod(*j, *k);
    float t1 = clock() - t;

    t = clock();
    for(vector <BN> :: iterator i = g.begin(), j = exp.begin(), k = mod.begin(); i != g.end(); i++, j++, k++)
        expkaryStringReplacement(*i, karyStringReplacementRepresentation(*j, 5),*k, 5);
    float t2 = clock() - t;

    float diviser = 1000.0 * test;
    t1 /= diviser;
    t2 /= diviser;

    printf("%d\t%d\t%d\t%.2f\t\t%.2f\n",base,(int)(base*sizeof(bt)*8),test,t1, t2);
}



void restest(int base,int test) {
    ull t;

    BN g(base, -1);
    vector <BN> exp;
    BN mod(base, -1);

    for(int i=0;i<test;i++) {
        exp.push_back( BN(base,-1) );
    }

    t = clock();
    for(vector <BN> :: iterator e = exp.begin(); e != exp.end(); e++)
        g.expLeftToRight(*e, mod);
    float t1 = clock() - t;

    t = clock();
    vector <BN> precomp_expBest_Slide = g.expBest_SlidePrecomp(mod);
    for(vector <BN> :: iterator e = exp.begin(); e != exp.end(); e++)
        g.expBest_Slide(*e, mod, precomp_expBest_Slide);
    float t2 = clock() - t;

//    t = clock();
//    for(vector <BN> :: iterator e = exp.begin(); e != exp.end(); e++)
//        expSignDigitRightToLeft(g, *e, mod);
//    float t3 = clock() - t;
//
    float diviser = 1000 * test;
    t1 /= diviser;
    t2 /= diviser;

    printf("%d\t%d\t%d\t%.2f\t\t%.2f\n",base,(int)(base*sizeof(bt)*8),test, t1, t2);
}


void resulttest() {
    cout<<"Test \"multiplication\":"<<endl;
    cout<<"Base\tBit\tTests\tClassic (µs)\tComma's (µs)\tCaracuba (µs)"<<endl;

    multest(16, 500000);
    multest(32, 100000);
    multest(64, 50000);
    multest(128, 20000);
    multest(256, 2000);
    multest(512, 1000);
    multest(1024, 200);
    multest(2048, 50);
    multest(4096, 20);

    cout<<"Test \"reduction\":"<<endl;
    cout<<"Base\tBit\tTests\tClassic (ms)\tBarrett (ms)\tMontgomery (ms)"<<endl;

    modtest(16, 50);
    modtest(32, 25);
    modtest(64, 5);
    modtest(128, 5);

    cout<<"Test \"Universal Methods\":"<<endl;
    cout<<"Base\tBit\tTests\tLeft-to-Right (ms)\tRight-to-Left (ms)\tk-ary (ms) \t k-ary [mod] (ms)\tsliding (ms)"<<endl;

    unitest(16, 50);
    unitest(32, 25);
    unitest(64, 5);
    unitest(128, 5);

    cout<<"Test \"k-ary method\":"<<endl;
    cout<<"Base\tBit\tTests\tk\tk-ary (ms) \t k-ary [Var] (ms)\t"<<endl;

    karytest(6, 16, 50);
    karytest(6, 32, 25);
    karytest(6, 64, 5);
    karytest(6, 128, 5);

    karytest(8, 16, 50);
    karytest(8, 32, 25);
    karytest(8, 64, 5);
    karytest(8, 128, 5);

    karytest(10, 16, 50);
    karytest(10, 32, 25);
    karytest(10, 64, 5);
    karytest(10, 128, 5);

    cout<<"Test \"k in Sliding Window\":"<<endl;
    cout<<"Base\tBit\tTests\tLeft-to-Right (ms)\tRight-to-Left (ms)\tk-ary (ms) \t k-ary [mod] (ms)\tsliding (ms)"<<endl;

    slidetest(16, 50);
    slidetest(32, 25);
    slidetest(64, 5);
    slidetest(128, 5);

    cout<<"Test \"Simultaneuos multiple exponentiation\":"<<endl;
    cout<<"Base\tBit\tTests\tClassic (ms)\tSimultaneous...(ms)"<<endl;

    simtest(16,10);
    simtest(32,10);
    simtest(64,10);
    simtest(128,10);

    simtest(16,5);
    simtest(32,5);
    simtest(64,5);
    simtest(128,5);

    simtest(16,2);
    simtest(32,2);
    simtest(64,2);
    simtest(128,2);

    cout<<"Test \"Method with Fixed-base\":"<<endl;
    cout<<"Base\tBit\tTests\tClassic (ms)\tWindow (ms)\tEuclidean (ms)"<<endl;

    fixtest(16, 50);
    fixtest(32, 50);
    fixtest(64, 50);


    fixtest(16, 10);
    fixtest(32, 10);
    fixtest(64, 10);
    fixtest(128, 10);

    fixtest(16, 5);
    fixtest(32, 5);
    fixtest(64, 5);
    fixtest(128, 5);

    fixtest(16, 2);
    fixtest(32, 2);
    fixtest(64, 2);
    fixtest(128, 2);

    cout<<"Test \"Method with exponent-recording\":"<<endl;
    cout<<"Base\tBit\tTests\tClassic (ms)\tSR(ms)"<<endl;

    exptest(16, 50);
    exptest(32, 25);
    exptest(64, 5);
    exptest(128, 5);

    cout<<"Test \"Super-Method\":"<<endl;
    cout<<"Base\tBit\tTests\tUniversal (ms)\tSlide mod(ms)\t SD mod"<<endl;

    restest(16, 50);
    restest(32, 25);
    restest(64, 5);
    restest(128, 5);
}

