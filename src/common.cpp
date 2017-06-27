#include <FL.h>
#include "common.h"

char *AllocCharArray(int n)
{
    char *v=NULL;
    v = (char *) calloc(n,sizeof(char));
    if (v == NULL)
        Error(MSG1,"AllocCharArray");
    return(v);
}

uchar *AllocUCharArray(int n)
{
    uchar *v=NULL;
    v = (uchar *) calloc(n,sizeof(uchar));
    if (v == NULL)
        Error(MSG1,"AllocUCharArray");
    return(v);
}

ushort *AllocUShortArray(int n)
{
    ushort *v=NULL;
    v = (ushort *) calloc(n,sizeof(ushort));
    if (v == NULL)
        Error(MSG1,"AllocUShortArray");
    return(v);
}

int *AllocIntArray(int n)
{
    int *v=NULL;
    v = (int *) calloc(n,sizeof(int));
    if (v == NULL)
        Error(MSG1,"AllocIntArray");
    return(v);
}

float *AllocFloatArray(int n)
{
    float *v=NULL;
    v = (float *) calloc(n,sizeof(float));
    if (v == NULL)
        Error(MSG1,(char*)"AllocFloatArray");
    return(v);
}

double *AllocDoubleArray(int n)
{
    double *v=NULL;
    v = (double *) calloc(n,sizeof(double));
    if (v == NULL)
        Error(MSG1,(char*)"AllocDoubleArray");
    return(v);
}

void Error(char *msg,char *func){ /* It prints error message and exits
                                    the program. */
    fprintf(stderr,"Error:%s in %s\n",msg,func);
    exit(-1);
}

void Warning(char *msg,char *func){ /* It prints warning message and
                                       leaves the routine. */
    fprintf(stdout,"Warning:%s in %s\n",msg,func);

}

void Change(int *a, int *b){ /* It changes content between a and b */
    int c;
    c  = *a;
    *a = *b;
    *b = c;
}

int NCFgets(char *s, int m, FILE *f) {
    while(fgets(s,m,f)!=NULL)
        if (s[0]!='#') return 1;
    return 0;
}


/**
 * Gera um número inteiro aleatório no intervalo [low,high].
http://www.ime.usp.br/~pf/algoritmos/aulas/random.html
 **/

int RandomInteger (int low, int high){
    int k;
    double d;

    d = (double) rand () / ((double) RAND_MAX + 1);
    k = d * (high - low + 1);
    return low + k;
}

double randomNormalized(){
    double d;
    d = (double) rand () / ((double) RAND_MAX + 1);
    return d;
}

double randomNumber(double low,double high){
    float d;
    d = (float) rand () / ((float) RAND_MAX + 1);
    d = d * (high - low + 1);
    d += low;
    return d;
}



double generateGaussianNoise(double mean, double standardDeviation)
{
    const double epsilon = -DBL_MAX;
    const double two_pi = 2.0*M_PI;
    double sigma = standardDeviation;
    static double z0, z1;
    static bool generate;
    generate = !generate;

    if (!generate)
        return z1 * sigma + mean;

    double u1, u2;
    do
    {
        u1 = rand() * (1.0 / RAND_MAX);
        u2 = rand() * (1.0 / RAND_MAX);
    }
    while ( u1 <= epsilon );

    z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
    z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
    return z0 * sigma + mean;
}

char *copyString(const char *str) {
    char *copy = NULL;

    if (str != NULL) {
        copy = (char*)calloc(strlen(str)+1,sizeof(char));
        strcpy(copy, str);
    }

    return copy;
}
