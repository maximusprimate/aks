#include <stdio.h>
#include <math.h>

/* Function declarations */

unsigned long long int power( unsigned long long int  base, unsigned long long int exp);

unsigned long long int euclid( unsigned long long int m, unsigned long long int n );

unsigned long long int order( unsigned long long int  element, unsigned long long int modulus);

unsigned long long int step1( unsigned long long int  n);

unsigned long long int step2( unsigned long long int  n, double log2n);

/* Main */

int main(void){

    unsigned long long int n;

    printf("Input unsigned long long integer: ");
    scanf("%lld",&n);

    double log2n = log(n)/log(2);
    printf("%lf\n",log2n);

    if(step1(n))
        printf("Continue..\n");
    else{
        printf("Composite.\n");
        return 0;
    }
    
    if(step2(n,log2n))
        printf("CONTINUE!");
    else
        return 0;



return 0;
}

/*
 * unsigned long long integer exponentiation:
 */

unsigned long long int power( unsigned long long int  base, unsigned long long int exp){

    unsigned long long int ans=base,i=1;
    while(i<exp){
        ans*=base;
        i++;
    }
    return ans;
}

unsigned long long int euclid( unsigned long long int m, unsigned long long int n );

void reduce(int numerator, unsigned long long int denominator,
            unsigned long long int *reduced_numerator,
            unsigned long long int *reduced_denominator);

/*
 * Euclidean algorithm
 */


unsigned long long int euclid( unsigned long long int m, unsigned long long int n ){
    unsigned long long int k, gcd;
    
    if (n<m) {
        k = n;
        n = m;
        m = k;
    }
    
    while (m % n != 0) {
        k = m % n;
        m = n;
        n = k;
    }
    
    gcd = n;
    
    return gcd;
}

/*
 * Order
 */

unsigned long long int order( unsigned long long int  r, unsigned long long int a){
    
    unsigned long long int k = 1;
    
    while(power(a,k)%r != 1)
        k++;

    return k;
}


/*
 * Step 1 of the AKS algorithm
 */

unsigned long long int step1( unsigned long long int  n){

    unsigned long long int a=2, b;
    for(a=2; a*a<(n+1); a++){
        for(b=2; pow(a,b)<(n+1); b++){
            //printf("%lld^%lld=%lld.\n",a,b,power(a,b));
            if(n==power(a,b))
                return 0;
        }
    }
    return 1;
}

unsigned long long int step2( unsigned long long int  n, double log2n){
    
    unsigned long long int r = 3;

    for(;r++;){
        if(euclid(r,n)!=1){
            printf("Composite! %lld | %lld\n",euclid(r,n),n);
            return 0;
        }
        if(order(r,n)>(log2n*log2n)){
            printf("Continue..");
            return 1;
        }
    }
    return 0;
}
