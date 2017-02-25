/////////////////////////////////////////////////////////
//
// aks.cpp
//
// Final project for CO685.
// Presented to Dr. D. Jao.
// Max Bennett, December 20, 2013
//
// The following C++ program uses the AKS algorithm
// to compute whether or not a given positive
// integer is prime or not. 
//
// Some of the steps of this program could be optimized
// or simplified, or in some cases removed altogether.
// Instead, the steps are meant to follow the outline
// of AKS explicitly, as defined by Agrawal, Kayal and
// Saxena themselves, in their paper PRIMES in P.
// Furthermore, readability was opted for, in leu of
// brevity.
// 
/////////////////////////////////////////////////////////

#include <stdio.h>
#include <time.h>
#include <math.h>
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>

using namespace std;
using namespace NTL;

/////////////////////////////////////////////////////////
//
// The NTL software package is used throughout the
// algorithm.  NTL is distributed under the  GNU LGPL.
// Although NTL does have some primality testing
// capabilities, none are used here.  NTL does however
// include a library of prime numbers, which are used
// in 'side steps' along the way.  In particular, in 
// calculating the totient of an integer, and elsewhere.
//
// NTL v6.0.0 was installed on the author's machine.
//
// More information about NTL may be found here:
//
// http://www.shoup.net/ntl
//
/////////////////////////////////////////////////////////

int step0(ZZ n);
int step1(ZZ n);
long step2(ZZ n);
int step3(ZZ n, long r);
int step4(ZZ n, long r);
int step5(ZZ n, long r);
int step6();
ZZ  gcd(ZZ m, ZZ n );
ZZ  order(ZZ  r, ZZ a);
ZZ phi(ZZ n);

/////////////////////////////////////////////////////////
//
// Main function
//
// This carries out the six different steps of the
// algorithm.
//
/////////////////////////////////////////////////////////

int main(void)
{

    int k;  // k is the truth value
    ZZ n;   // n is the integer being tested
    long r; // r is for step 2

    printf("Please input a positive integer: ");

    cin >> n;
    
    k = step0(n);
    if(k){
        if(k==1){
            cout << "Step 0, composite.\n";
            return 0;
        }
        cout << "Step 0, prime.\n";
        return 1;
    }
    cout << "Step 0 was not sure.\n";
    
    k = step1(n);
    if(k){
        cout << "Step 1, composite.\n";
        return 0;
    }
    cout << "Step 1 was not sure.\n";

    r = step2(n);

    k = step3(n,r);

    if(k){
        cout << "Step 3, composite.\n";
        return 0;
    }
    cout << "Step 3 was not sure.\n";

    k = step4(n,r);

    if(k){
        cout << "Step 4, prime.\n";
        return 1;
    }
    cout << "Step 4 was not sure.\n";

    k = step5(n,r);
    
    if(k){
        cout << "Step 5, composite.\n";
        return 1;
    }
    cout << "Step 5 was not sure.\n";
   
    k = step6();
    
    if(k){
        cout << "Step 6, prime.\n";
        return 1;
    }
}

/////////////////////////////////////////////////////////
//
// Extra functions.
//
// These include the Euclidean algorithm, the order of
// an element, and Euler's totient function.
//
/////////////////////////////////////////////////////////

/*
 * Euclidean algorithm, should be self explanitory.
 */

ZZ gcd( ZZ m, long n ){
    ZZ k, z;
    z = 0;
    z = n+z; // easiest way to promote n to ZZ.
    if (z<m) {
        swap(z,m);
    }
    
    while (m % z != 0) {
        k = m % z;
        m = z;
        z = k;
    }
    return z;
}

/*
 * The order of a modulo r is the smallest integer
 * k such that a^k = 1 mod r.  It's denoted o_r(a).
 */

ZZ order(ZZ a, long r){

    ZZ k,o,z;
    o = 1;
    z = 1;
    k = a;
    z = r + z; // promoting r to type ZZ

    while(k!=1){
        k*=a; k%=z; o++;
    }

    return o;
}

/*
 * Phi is Euler's totient function. This is not optimized
 * at all, and does not work if the input has a prime 
 * factor larger than the largest prime in the NTL
 * PrimeSeq library, which isn't usually a problem.
 */


ZZ phi(long x){

    ZZ n, ph,p; ph = 1;
    n = 0; n = n + x; // promoting x to type ZZ

    int k;

    for(PrimeSeq s; n>1; ){
        p = s.next();
        k = 0;
        while(n%p == 0){
            k++; n/=p;
        }
        if(k>0)
        ph *= power(p,k-1)*(p-1);
    }
    return ph;
}

/////////////////////////////////////////////////////////
//
// Steps.
//
// These are the steps of the algorithm as outlined
// in the paper.  A few of them could have been
// put in the main function due to their simplicity,
// but they were kept as steps to maintain consistency
// with the paper instead.
//
/////////////////////////////////////////////////////////

/*  Step 0:
 *  It is not difficult to test the first primes
 *  less than 2000 to see if they are divisors of
 *  n.  When we do this, step1 becomes much faster
 *  as we only need to check if n = a^b for a>2000,
 *  since if n = a^b for some a, then there exists
 *  a prime p <= a that divides n.
 *
 *  This code snippet was taken from Victor Shoup's
 *  website (the author of NTL).  It can be found
 *  here:
 *
 *  http://www.shoup.net/ntl/doc/tour-ex1.html
 *
 *  This function is not a part of the AKS algorithm.
 *
 *  Returns 1 if composite, 2 if prime, 0 if it's
 *  not sure.
 *
 */

int step0(ZZ n){

    long p;

    PrimeSeq s;
    
    p = s.next();  
    while (p && p < 2000) {
        if ((n % p) == 0){
            cout << p << " divides " << n << '\n';
            return (1 + (n == p)); // 1: composite
        }                          // 2: prime
        p = s.next();              
    }
    return 0; // 0: continue
}

/* Step 1:
 * If n=a^b for a in N and b>1, output composite.
 * The idea here is to test different values of b.
 * We start with b=2, and take the b-th root of n.
 * We call this number A, which is in the class RR.
 * If, when we round A down, we get the same thing,
 * we can test if A^b == n.  (We do this second step
 * to deal with rounding errors.  Just because when we
 * round down we get 0, our precision may not be perfect
 * so we test the exponetiation for these few cases.)  
 *
 * If A^b=n, our test is over, else we increase b by one.
 *
 * We keep doing this until b reaches log_2000 (n).
 * Since we tested divisibility of all primes less
 * than 2000, we can stop here.
 *
 * Returns 1 if composite, 0 if we don't know.
 */

int step1(ZZ n){

    RR A, B, N, log_2000_n;
    long b = 2;

    log_2000_n = log(n)/log(2000);
    N = to_RR(n);

    while(b<log_2000_n){

        B = 1.0/b;    // pow needs both args RR
        A = pow(N,B); // A is b-th root of N
        if((A - to_RR(FloorToZZ(A)))==0.0){ // A an int?
            if(power(FloorToZZ(A),b) == n)  // A^b==n?
                return 1; // 1: composite
        }
        b++;
    }

    return 0; // 0: ??

}

/* Step 2:
 * This is a function that searches for the 
 * smallest r such that ord_r(n) > log_2(n)^2.
 */

long step2(ZZ n){

    RR log2n, l;
    ZZ k;
    long r;

    log2n = NumBits(n);
    l = power(log2n,2);
    k = FloorToZZ(l);
    
    for(r=2;;r++){
        if(gcd(n,r)==1){
            if(order(n,r) > k){
                cout << "Step 2: r = " << r << '\n';
                return r; // result
            }
        }
    }
}

/* Step 3:
 * This computes gcd(a,n) for values less than
 * r, where r is the value computed in step 3.
 * Notice that since we have already determined that
 * for any value k<2000, k does not divide n, we
 * do not need to check values of gcd(a,n) for values
 * a<2000.
 *
 * Returns 1 if composite, 0 if we don't know.
 */

int step3(ZZ n, long r){

    long a;

    for(a=2002;++a<r;){
        if ( (gcd(n,a)%n)>1 ) // if 1 < (a,n) < n
            return 1; // 1: composite
    }
    
    return 0;
}

/* Step 4:
 * This is an extremely simple function and
 * is only laid out here to be consistent with the 
 * paper.  If n<=r, we know that n is prime.
 */

int step4(ZZ n, long r){
    if(n <= r)
        return 1; // 1: prime
    return 0;     // 0: ??
}

/* Step 5:
 * This is the heart of the algorithm.  It tests
 * for different values of a, whether or not
 * (x+a)^n = x^n + a, over the polynomial ring
 * F = Z_n[x] / x^r - 1.  If for any value of a the
 * equality doesn't hold, we output composite.  
 * Otherwise we continue.
 *
 * Unfortunately the input of polynomials using NTL is
 * not straight forward. Comments should provide insight
 * as to what is being defined and why.
 *
 * Inspiration for step 5 is taken from George Poulose's
 * AKS implementation.  Nothing was copied directly but
 * certain methods were immitated.  In particular, the 
 * use of polynomials, and their input.
 *
 * His program may be found here:
 *
 * http://www.gpoulose.com/gc/AKS_cpp.txt
 *
 * Function returns 1 if composite, 0 if we don't know.
 * Although, for all intents and purposes, this is the 
 * last step, so seeing a 0 means we have a prime.  But,
 * in keeping with the outline of the algorithm, we 
 * pretend that we don't know what happens when we see
 * a 0.
 */

int step5(ZZ n, long r){

    ZZ l;
    long a;

    NTL::ZZ_p::init(n); // Initializes modulus to n
                        // so we are working over Z_n.

    ZZ_pX polymod(r,1); // polymod = 1*x^r
    polymod -=1 ;       // polymod = x^r - 1
                        // Here we define what will
                        // be the modulus.

    ZZ_pXModulus mod(polymod); // mod is now a modulus.

    ZZ_pX RHS(1, 1);            // RHS = x
    PowerMod(RHS, RHS, n, mod); // RHS = x^n in F
    
    l = FloorToZZ(sqrt(to_RR(phi(r)))*NumBits(n));
                        // l is limit of the for loop

    for(a=1; a<=l; a++){
        ZZ_pX LHS(1,1); // LHS = x
        LHS += a;       // LHS = x + a
        PowerMod(LHS,LHS,n,mod); // LHS = (x+a)^n in F
        LHS -= a;       // LHS = (x+a)^n - a;
        if(LHS != RHS)  //     (x+a)^n - a != x^n
                        // iff (x+a)^n != x^n + a
          return 1;     // 1: composite
    }

    return 0; // 0: ??
}

/* Step 6:
 * If we've got this far, we know that n is prime.
 * Again, this function is merely a placeholder
 * to be consistent with the paper.
 */

int step6(){
    return 1; // 1: prime
}
