#include "FL.h"
#include <time.h>
#include <iostream>
#include <vector>

void destroyFloatVector(void* data){
    float* aux = *( (float**) data);
    free(aux);
    aux = NULL;
}

int main(int argc, char **argv) {

    size_t n = 10000000;
    clock_t start, end;
    double time_normal;
    double time_vector;
    double time_vector_slow;
    const char* separator= "----------------------------------------------------------\n";
    printf("%s",separator);
    printf("workload is: %lu floats\n",n);
    printf("%s",separator);
    printf("Test: allocating data\n");
    start = clock();
    float* ref = (float*)calloc(n,sizeof(float));
    //float* refIncreased;
    end = clock();
    time_normal = ((double) (end - start)) / CLOCKS_PER_SEC;

    start = clock();
    GVector* vector = createNullVector(n,sizeof(float));
    end = clock();
    time_vector = ((double) (end - start)) / CLOCKS_PER_SEC;
    GVector* vector_slow = createNullVector(n,sizeof(float));

    printf("results:\nnormal float pointer: %f s\nvector: %f\n",time_normal,time_vector);
    printf("%s",separator);
    printf("Test: modefying data\n");
    start = clock();
    for (size_t i = 0; i < n; ++i) {
        float aux  = i;
        ref[i] = aux;
    }
    end = clock();
    time_normal = ((double) (end - start)) / CLOCKS_PER_SEC;

    start = clock();
    for (size_t i = 0; i < n; ++i) {
        float aux  = i;
        VECTOR_GET_ELEMENT_AS(float,vector,i) = aux;
    }
    end = clock();
    time_vector = ((double) (end - start)) / CLOCKS_PER_SEC;

    start = clock();
    for (size_t i = 0; i < n; ++i) {
        float aux  = i;
        *((float*)getElementInVectorAt(vector_slow,i)) = aux;
    }
    end = clock();
    time_vector_slow = ((double) (end - start)) / CLOCKS_PER_SEC;

    printf("results:\nnormal float pointer: %f s\nvector_macro: %f s\nvector_function: %f\n",time_normal,time_vector,time_vector_slow);
    printf("%s",separator);
    printf("Test: resize vector: increasing size\n");
    start = clock();
    ref = (float*) realloc (ref, 2*n * sizeof(float));
    end = clock();
    time_normal = ((double) (end - start)) / CLOCKS_PER_SEC;

    start = clock();
    resizeVector(vector,2*n);
    end = clock();
    time_vector = ((double) (end - start)) / CLOCKS_PER_SEC;

    printf("results:\nnormal float pointer: %f s\nvector_function: %f s\n",time_normal,time_vector);
    printf("%s",separator);
    printf("Test: resize vector: decreasing size\n");

    VECTOR_GET_ELEMENT_AS(float,vector,n+(n/2)) = 0;
    start = clock();
    ref = (float*) realloc (ref, n * sizeof(float));
    end = clock();
    time_normal = ((double) (end - start)) / CLOCKS_PER_SEC;


    start = clock();
    resizeVector(vector,n);
    end = clock();
    time_vector = ((double) (end - start)) / CLOCKS_PER_SEC;

    start = clock();
    shrinkToFit(vector);
    end = clock();
    time_vector_slow = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("results:\nnormal float pointer: %f s\nvector_function (resize): %f s\nvector_function (resize+shrink): %f s\n",
           time_normal,time_vector,time_vector+time_vector_slow);

    return 0;
}




