#include <iostream>
#include <vector>
#include <algorithm>
#include "additional.h"
#include "BN.h"
#include "BNsign.h"

#include <cstdio>
bt2s abs(bt2s x) {
    return x<0?-x:x;
}

BN gcdEuclidean(BN a,BN b)
{
    while(!b.is0()) {
        BN temp=a%b;
        a=b;
        b=temp;
    }
    return a;
}

BN gcdInverseEuclidean(BN a, BN mod)
{
	BN start_mod=mod;
        a=a%mod;
        BN bn_1 = 1;
        BN bn_0(1,0);

        if(a.is0())
		return bn_0;
	BNsign x0;
	BNsign x1=bn_0;
	BNsign x2=bn_1;
	while(!(mod%a).is0())
	{
		BN q=mod/a;
		x0=x1;
		x1=x2;

		x2 = x0 - x1 * (BNsign)q;
		BN new_mod=a;
		a=mod%a;
		mod=new_mod;
	}
	if(a!=bn_1)
		return bn_0;
	if(x2.value.is0())
		x2.sign=false;
	if(!x2.sign)
		return x2.value;
	else
		return start_mod-x2.value;
}

BN gcdBinary(BN a,BN b)
{
	if(a.is0())
		return b;
	if(b.is0())
		return a;
	//BN g(1);

	int acount=a.countzeroright();
	int bcount=b.countzeroright();

	//g=g<<min(acount,bcount);
	a=a>>acount;
	b=b>>bcount;

	if(a<b)
	{
		BN t=a;
		a=b;
		b=t;
	}
	while(!a.is0())
	{
		/*if(a>=b)
			a=(a-b)>>(a-b).countzeroright();
		else
			b=(b-a)>>(b-a).countzeroright();*/
		if(a>=b)
		{
			BN ab=a.sub(b);
			a=ab>>ab.countzeroright();
		}
		else
		{
			BN ba=b.sub(a);
			b=ba>>ba.countzeroright();
		}
	}
	return b<<min(acount,bcount);
}

BN gcdExtendedEuclideanBinary(BN xx, BN yy)
{
	BN bn_1(1);
	BN bn_0(1,0);
	//if(yy.is0()||gcdEuclidean(xx,yy)!=bn_1)
		//return bn_0;
	BN g=bn_1;
	int xcount=xx.countzeroright();
	int ycount=yy.countzeroright();
	g=g<<min(xcount,ycount);
	xx=xx>>min(xcount,ycount);
	yy=yy>>min(xcount,ycount);

	BN u=xx;
	BN v=yy;
	BNsign x=xx;
	BNsign y=yy;

	BNsign a=bn_1;
	BNsign b=bn_0;
	BNsign c=bn_0;
	BNsign d=bn_1;

	do
	{
		/*
		u.PrintDec();
		v.PrintDec();
		a.PrintSign();
		b.PrintSign();
		c.PrintSign();
		d.PrintSign();
		printf("----next----\n");
		*/
		while(u.isEven())
		{
			u=u>>1;
			if(a.value.isEven() && b.value.isEven())
			{
				a.value=a.value>>1;
				b.value=b.value>>1;
			}
			else
			{
				a=(a+y);
				a.value=a.value>>1;
				b=(b-x);
				b.value=b.value>>1;
			}
		}
		while(v.isEven())
		{
			v=v>>1;
			if(c.value.isEven()&&d.value.isEven())
			{
				c.value=c.value>>1;
				d.value=d.value>>1;
			}
			else
			{
				c=(c+y);
				c.value=c.value>>1;
				d=(d-x);
				d.value=d.value>>1;
			}
		}
		if(u>=v)
		{
			u=u-v;
			a=a-c;
			b=b-d;
		}
		else
		{
			v=v-u;
			c=c-a;
			d=d-b;
		}
	}
	while(!u.is0());
	/*
	u.PrintDec();
	v.PrintDec();
	a.PrintSign();
	b.PrintSign();
	c.PrintSign();
	d.PrintSign();
	printf("----end----\n");
	*/
	//u.is0() == true;
	BNsign A=c;
	BNsign B=d;
	A.PrintSign();
	B.PrintSign();
	(g*v).PrintDec();
	return g*v;
	// A*x + B*y == v
}

BN gcdInverseEuclideanBinary(BN xx, BN mod)
{
	BN bn_1(1);
	BN bn_0(1,0);
	if(mod.is0())
		return bn_0;
//	if(mod.is0()||gcdEuclidean(xx,mod)!=bn_1)
//		return bn_0;

	BN g=bn_1;
	int xcount=xx.countzeroright();
	int ycount=mod.countzeroright();
	if(min(xcount,ycount))
		return bn_0;
	//min(xcount,ycount) = 0!
//	g=g<<min(xcount,ycount);
//	xx=xx>>min(xcount,ycount);
//	mod=mod>>min(xcount,ycount);

	BN u=xx;
	BN v=mod;
	BNsign x=xx;
	BNsign y=mod;

	BNsign a=bn_1;
	BNsign b=bn_0;
	BNsign c=bn_0;
	BNsign d=bn_1;
	if(mod.isEven())
		do
		{
			int ucount=u.countzeroright();
			u=u>>ucount;
			while(ucount)
			{
				int abcount=min(a.value.countzeroright(),b.value.countzeroright());
				if(abcount>ucount)
				{
					a.value=a.value>>ucount;
					b.value=b.value>>ucount;
					ucount=0;
				}
				else if(abcount)
				{
					a.value=a.value>>abcount;
					b.value=b.value>>abcount;
					ucount-=abcount;
				}
				else
				{
					a=(a+y);
					a.value=a.value>>1;
					b=(b-x);
					b.value=b.value>>1;
					ucount--;
				}
			}
			int vcount=v.countzeroright();
			v=v>>vcount;
			while(vcount)
			{
				int cdcount=min(c.value.countzeroright(),d.value.countzeroright());
				if(cdcount>vcount)
				{
					c.value=c.value>>vcount;
					d.value=d.value>>vcount;
					vcount=0;
				}
				else if(cdcount)
				{
					c.value=c.value>>cdcount;
					d.value=d.value>>cdcount;
					vcount-=cdcount;
				}
				else
				{
					c=(c+y);
					c.value=c.value>>1;
					d=(d-x);
					d.value=d.value>>1;
					vcount--;
				}
			}
			if(u>=v)
			{
				u=u-v;
				a=a-c;
				b=b-d;
			}
			else
			{
				v=v-u;
				c=c-a;
				d=d-b;
			}
		}
		while(!u.is0());
	else
		do
		{
			int ucount=u.countzeroright();
			u=u>>ucount;
			while(ucount)
			{
				int acount = a.value.countzeroright();
				if(acount)
				{
					a.value=a.value>>min(acount,ucount);
					ucount-=min(acount,ucount);
				}
				else
				{
					a=(a+y);
					a.value=a.value>>1;
					ucount--;
				}
			}
			int vcount=v.countzeroright();
			v=v>>vcount;
			while(vcount)
			{
				int ccount = c.value.countzeroright();
				if(ccount)
				{
					c.value=c.value>>min(ccount,vcount);
					vcount-=min(ccount,vcount);
				}
				else
				{
					c=(c+y);
					c.value=c.value>>1;
					vcount--;
				}
			}
			if(u>=v)
			{
				u=u-v;
				a=a-c;
			}
			else
			{
				v=v-u;
				c=c-a;
			}
		}
		while(!u.is0());

	if(v!=bn_1)
		return bn_0;

	if(c.value.is0())
		c.sign=false;
	if(c.sign)
		return mod-(c.value)%mod;
	else
		return c.value%mod;
}

vector <BN> multi_inverse(const vector <BN> &x, const BN &mod)
{
	int count=x.size();
	vector <BN> a;
	vector <BN> x_inverse;

	a.reserve(count);
	x_inverse.reserve(count);

	a.push_back(x[0]);

	for(	vector<BN>::const_iterator xi=x.begin()+1;
		xi!=x.end();
		xi++)
		a.push_back(a.back() * (*xi) % mod);

	BN curr = gcdInverseEuclideanBinary(a[0],mod);
	for(int i=x.size()-1;i>0;i--)
	{
		x_inverse.push_back(a[i-1] * curr % mod);	//index = i
		curr = x[i] * curr % mod;
		//x_inverse.push_back( (a[i-1] * curr).reduction_special(mod));
		//curr = (x[i] * curr).reduction_special(mod);
	}
	x_inverse.push_back(curr);				//index = 0
	reverse(x_inverse.begin(),x_inverse.end());
	return x_inverse;
}

//#define DEBUG_LEHMER
BN gcdLehmer(BN x,BN y)
{
	if(x<y)
	{
		BN t=x;
		x=y;
		y=t;
	}

	while(y.basecount()>1)
	{
		bt2s xp = x[x.basecount()-1];
		bt2s yp = y[y.basecount()-1];
		bt2s A = 1;
		bt2s B = 0;
		bt2s C = 0;
		bt2s D = 1;
		if(x.basecount()==y.basecount())
			while( (yp+C)!=0 && (yp+D)!=0 )
			{
				bt2s q = (xp+A)/(yp+C);
				bt2s qp= (xp+B)/(yp+D);
				if(q!=qp)
					break;
				bt2s t;
				t = A - q*C;	A = C;	C = t;
				t = B - q*D;	B = D;	D = t;
				t = xp- q*yp;	xp= yp;	yp= t;
			}
		if(B==0)
		{
			//TODO: BN T;T=x%y; - медленнее на 10 % (!)
			BN T = x%y;
			x = y;
			y = T;
		}
		else
		{
			BN T = (B<0 ? x.mulbase(abs(A)) - y.mulbase(abs(B)) : y.mulbase(abs(B)) - x.mulbase(abs(A)));
			BN U = (D<0 ? x.mulbase(abs(C)) - y.mulbase(abs(D)) : y.mulbase(abs(D)) - x.mulbase(abs(C)));
			//TODO: вариант с if-ами медленнее на 8-10 % (!)
			x=T;
			y=U;
		}
	}
	return gcdBinary(x,y);
}

BN Garner(vector <BN> m, vector <BN> v)
{
	if(m.size()!=v.size())
		throw "Garner: size's M and V are different!";
	int t = m.size();
	BN *C = new BN [t];
	BN u;
	for(int i=1;i<t;i++)
	{
		C[i] = 1;
		for(int j=0; j<i; j++)
		{
			u = gcdInverseEuclideanBinary(m[j],m[i]);
			C[i] = u * C[i] % m[i];
		}
	}
	u = v[0];
	BN x = u;
	BN Pm=1;
	for(int i=1;i<t;i++)
	{
		BN xmi = x%m[i];
		if(xmi >= m[i])
			cout<<"erunda!!!";
		BN vix = (v[i] >= xmi ? v[i] - xmi : m[i] + v[i] - xmi);
		//u = (v[i] - x) * C[i] % m[i];
		u = vix * C[i] % m[i];
		Pm = Pm * m[i-1];
		x = x + u * Pm;
	}
	return x;
}

BN CTO(vector <BN> m, vector <BN> v)
{
	if(m.size()!=v.size())
		throw "CTO: size's M and V are different!";

	int t = m.size();
	BN M = m[0];
	for (int i=1;i<t;i++)
		M = M * m[i];
	BN x(1,0);
	for (int i=0;i<t;i++)
	{
		BN Mi = M/m[i];
		x = x + v[i] * (M/m[i]) * gcdInverseEuclideanBinary(M/m[i],m[i]) % M;
	}
	return x%M;
}




vector <BN> expSimutaneousMulPrecomputation(vector <BN> g, vector <BN> exp, BN mod) {
    if(g.size() > sizeof(int)*8 - 2)
        throw "k is very long for this architecture\n";

    if(g.size() != exp.size())
        throw "incorrect k in expSimutaneousMulPrecomputation\n";

    int k = (int) g.size();
    int kpow = 1 << k;
    vector <BN> G(kpow);
    for(int i = 0; i < k; i++)
        g[i] = g[i] % mod;
    for(int i = 0; i < kpow; i++) {
        G[i] =  (BN) 1;
        unsigned int mask = 1;
        for(int j = 0; j < k; j++) {
            if(i & mask) {
                G[i] = G[i] * g[j] % mod;
            }
            mask <<= 1;
        }
    }
    return G;
}

BN expSimutaneousMul(vector <BN> G, vector <BN> exp, BN mod) {
    //TODO: k != bz8 possible

    //if(exp.size() != bz8)
    //    throw "incorect k\n";
    int k = exp.size();

    int t = exp[0].bitcount();
    for(int i = 1; i < k; i++) {
        if(t < exp[i].bitcount())
            t = exp[i].bitcount();
    }

    BN A = (BN) 1;
    for(int i = 1; i <= t; i++) {
        A = A.Qrt() % mod;
        unsigned int I = 0;
        for(int j = k-1; j >= 0 ; j--)
            I = (I << 1) | exp[j].bitI(t-i);
        A = A * G[I] % mod;
    }
    return A;
}

vector <BN> expFixedBaseWindowPrecomputation(BN g, int t, BN mod) {
    if(t < 0)
        throw "expFixedBaseWindow: t < 0\n";
    vector <BN> res(t+1);
    res[0] = g % mod;
    for(int i = 1; i <= t; i++) {
        res[i] = res[i-1];
        for(int j = 0; j < (int)bz8; j++)
            res[i] = res[i].Qrt() % mod;
    }
    return res;
}

BN expFixedBaseWindow(const vector <BN> & g, BN exp, BN mod) {
    BN A = 1;
    BN B = 1;
    int h = bsize;
    int t = g.size();
    for(int j = h - 1; j > 0; j--) {

        for(int i = 0; i < t; i++) {
            if(exp[i] == j)
                B = B * g[i] % mod;
        }
        
        A = A * B % mod;
    }
    return A;
}

pair <int,int> expFixedBaseEuclideanMN(vector <bt> x) {
    //if x[0] >= x[1]
    pair <int,int> res(0,0);    // (M, N)

    bt M_value = x[0];
    bt N_value;
    int t = x.size();
    for(int i = 1; i < t; i++) {
        if(x[i] >= M_value) {
            N_value = M_value;
            res.second = res.first;
            M_value = x[i];
            res.first = i;
        }
        else {
            if(x[i] >= N_value) {
                N_value = x[i];
                res.second = i;
            }
        }
    }
    return res;
}

BN expFixedBaseEuclidean(const vector<BN> &g, BN exp, BN mod) {
    int t = g.size();
    if(t < 3)
        throw "expFixedBaseEuclidean: t < 2\n";
    vector <bt> x(g.size());
    vector <BN> G(g.size());
    for(int i = 0; i < t; i++) {
        G[i] = g[i];
        x[i] = exp[i];
    }
    pair <int,int> MN = expFixedBaseEuclideanMN(x);
    int M = MN.first;
    int N = MN.second;
    while(x[N] != 0) {
        bt q = x[M] / x[N];
        for(int i = 0; i < q; i++) {
            G[N] = G[N] * G[M] % mod;
        }
        x[M] = x[M] % x[N];
        MN = expFixedBaseEuclideanMN(x);
        M = MN.first;
        N = MN.second;

    }
    BN res = 1;
    for(int i = 0; i < x[M]; i++) {
        res = res * G[M] % mod;
    }
    return res;
}


SDbn :: SDbn(const SDbn& sd) {
    lenght = sd.lenght;
    d = new char [lenght];
    for(int i = 0; i < lenght; i++)
        d[i] = sd.d[i];
}


SDbn::SDbn(const BN & bn) {
    int bncount = 0;
    int sdcount = 0;
    int t = bn.bitcount();
    lenght = t + 1;
    d = new char[lenght];
    char c0 = 0;
    char c1 = 0;
    for(int i = 0; i < lenght; i++) {
        c1 = (bn.bitI(i) + bn.bitI(i+1) + c0) / 2;
        d[i] = bn.bitI(i) + c0 - 2 * c1;
        if(bn.bitI(i))
            bncount ++;
        if(d[i])
            sdcount ++;
        c0 = c1;
    }
}

SDbn SDbn::operator = (const SDbn & sd) {
    if(this == &sd)
        return *this;
    delete [] d;
    lenght = sd.lenght;
    d = new char [lenght];
    for(int i = 0; i < lenght; i++)
        d[i] = sd.d[i];
    return *this;
}

char SDbn :: operator [] (int index) {
    return d[index];
}



int SDbn::Lenght() {
    return lenght;
}

BN expSignDigitRightToLeft(BN g, BN exponent, BN mod) {
    if(exponent.is0())
        return (BN) 1;

    BN A = 1;
    BN S = g % mod;
    // TODO: вернуть проверку! )
    if(gcdBinary(S, mod) != (BN) 1) {
        return g.PowModBarrett(exponent, mod);
    }

    BN s = gcdInverseEuclideanBinary(S, mod);         // g ^ (-1)
    SDbn exp = exponent;

    BN mu = mod.reduction_barrett_precomputation();
    for(int i = 0; i < exp.Lenght(); i++) {
        if(exp[i] == 1) {
            A = (A * S).reduction_barrett(mod, mu);
        }
        else if(exp[i] == -1) {
            A = (A * s).reduction_barrett(mod, mu);
        }
        S = S.Qrt().reduction_barrett(mod, mu);
        s = s.Qrt().reduction_barrett(mod, mu);
    }
    return A;
}


vector <int> karyStringReplacementRepresentation(BN exp, int k) {
    if(k > 30)
        throw "karyStringReplacementRepresentation: k > 30\n";
    if(k < 2)
        throw "karyStringReplacementRepresentation: k < 2\n";


    int t = exp.bitcount();
    vector <int> exp_sr;
    for(int i = 0; i < t; i++)
        exp_sr.push_back(exp.bitI(i));

    for(int i = k; i >= 2; i--)                         // для всех k

        for(int j = t - 1; j >= i - 1; j--) {
            int count = 0;                              // пытаемся найти начало последовательности

            while(exp_sr[j - count] == 1 && count < i)
                count++;

            if(count == i) {                            // нашли!
                for(count = 0; count < i - 1; count++)  // обнулили биты
                    exp_sr[j - count] = 0;
                exp_sr[j - i + 1] = i;
            }
        }

    return exp_sr;
}


BN expkaryStringReplacement(const BN &g, const vector <int> & exp, const BN & mod, int k) {
    vector <BN> G;
    G.push_back((BN) 1);
    
    G.push_back(g % mod);
    for(int i = 2; i <= k; i++) {
        BN ggg = (G[i-1].Qrt() % mod) * g % mod;
        G.push_back(ggg);
    }

    BN A = 1;
    int t = exp.size();

    for(int i = t - 1; i >= 0; i--) {
        A = A.Qrt() % mod;
        if(exp[i]) {
            A = A * G[exp[i]] % mod;
        }
    }
    return A;
}



int Fbc::BinaryToInt (int t, int *e) {
    int value = e[t - 1];
    for (int i = t - 2; i >= 0; i--)
        value = (value << 1) + e[i];
    return value;
}


Fbc::Fbc(BN g, BN mod, int bitcount) {
    m = mod;
    int t = bitcount;
    t--;

//    do h = rand() % min(8,(t + 2)); while (h == 0);
    h = 3;
    a = (t + h) / h;    // t+1 делить на h с округлением вверх

    int h2 = 1 << h;

//    do v = rand() % min(8, (a + 1)); while (v == 0);
    v = min(3, a);
    b = (a + v - 1) / v;

    G = new BN * [v];
    for (int j = 0; j < v; j++)
        G[j] = new BN [h2];

    BN * g1 = new BN [h];
    for (int i = 0; i < h; i++) {
        BN ia = 1;
        ia = ia << i*a;
        g1[i] = g.PowModBarrett(ia, m);
    }

    for (int i = 1; i < h2; i++) {
        BN f = i;
        int s = f.bitcount();
        BN A = 1;
        for (int j = 0; j < s; j++)
            if(f.bitI(j))
                A = A * g1[j] % m;
        G[0][i] = A;
        for (int j = 1; j < v; j++) {
            BN jb = 1;
            jb = jb << j*b;
            G[j][i] = G[0][i].PowModBarrett(jb, m);
        }
    }
    delete [] g1;
}

Fbc::~Fbc() {
    for (int i = 0; i < v; i++)
        delete [] G[i];
    delete [] G;
}

BN Fbc::fixed_base_comb(BN exp) {
    int ** EA = new int * [h];
    for (int i = 0; i < h; i++)
        EA[i] = new int [a];

    int ea_index = 0;
    for (int j = 0; j < h; j++) {
        for (int k = 0; k < a; k++) {
            EA[j][k] = exp.bitI(ea_index++);
        }
    }
    int * Y = new int [h];
    BN A = 1;
    for (int k = b - 1; k >= 0; k--) {
        A = A.Qrt() % m;
        for (int j = v - 1; j >= 0; j--) {
            int c = j * b + k;
            for (int i = 0; i < h; i++)
                if(c < a)
                    Y[i] = EA[i][c];
                else
                    Y[i] = 0;
            int Ijk = BinaryToInt(h, Y);


            if (Ijk != 0)
                A = A * G[j][Ijk] % m;
        }
    }

    delete [] Y;
    for (int i = 0; i < h; i++)
        delete [] EA[i];
    delete [] EA;

    return A;
}


