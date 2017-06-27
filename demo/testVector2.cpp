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
    size_t n_aux;
    clock_t start, end;
    double percentage = 0.05;
    double percentageInsertMiddle = 0.01;
    double vector_macro_time;
    double vector_function_time;
    n_aux = n*percentageInsertMiddle;

    const char* separator= "----------------------------------------------------------\n";
    GVector* vectorMacro = createVector(n,sizeof(float));
    GVector* vectorFunction = createVector(n,sizeof(float));

    printf("%s",separator);
    printf("workload is: %lu floats\n",n);
    printf("%s",separator);
    printf("Test: pushing back elements (no resize)\n");
    start = clock();
    for (size_t i = 0; i < n; ++i) {
        float aux = i;
        pushBackElementInVector(vectorFunction,&aux);
    }
    end = clock();
    vector_function_time = ((double) (end - start)) / CLOCKS_PER_SEC;

    start = clock();
    for (size_t i = 0; i < n; ++i) {
        float aux = i;
        VECTOR_PUSH_BACK(float,vectorMacro,aux+2);
    }
    end = clock();
    vector_macro_time = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("results:\nvector_function: %f s\nvector_macro: %f\n s",vector_function_time,vector_macro_time);
    printf("%s",separator);
    printf("Test: pushing back elements (with resize, starting from: %.2f%% of %lu)\n",percentage*100,n);
    destroyVector(&vectorFunction);
    destroyVector(&vectorMacro);
    vectorFunction = createVector(n*percentage,sizeof(float));
    vectorMacro = createVector(n*percentage,sizeof(float));
    start = clock();
    for (size_t i = 0; i < n; ++i) {
        float aux = i;
        pushBackElementInVector(vectorFunction,&aux);
    }
    end = clock();
    vector_function_time = ((double) (end - start)) / CLOCKS_PER_SEC;

    start = clock();
    for (size_t i = 0; i < n; ++i) {
        float aux = i;
        if(vectorMacro->size >= vectorMacro->capacity){
            setVectorCapacity(vectorMacro,vectorMacro->size*vectorMacro->growthFactor);
        }
        VECTOR_PUSH_BACK(float,vectorMacro,aux);
    }
    end = clock();
    vector_macro_time = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("results:\nvector_function: %f s\nvector_macro: %f s\n",vector_function_time,vector_macro_time);
    printf("%s",separator);


    destroyVector(&vectorFunction);
    destroyVector(&vectorMacro);
    vectorFunction = createVector(n_aux,sizeof(float));
    vectorMacro = createVector(n_aux,sizeof(float));
    printf("Test: insert element (one by one) at middle (no resize). Using  %.2f%% of %lu\n",percentageInsertMiddle*100,n);

    start = clock();
    for (size_t i = 0; i < n_aux; ++i) {
        float aux = i;
        insertElementInVectorAt(vectorFunction,&aux,vectorFunction->size/2);
    }
    end = clock();
    vector_function_time = ((double) (end - start)) / CLOCKS_PER_SEC;
    start = clock();
    for (size_t i = 0; i < n_aux; ++i) {
        float aux = i;
        VECTOR_INSERT_AT(float, vectorMacro,aux,vectorMacro->size/2);
    }
    end = clock();
    vector_macro_time = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("results:\nvector_function: %f s\nvector_macro: %f s\n",vector_function_time,vector_macro_time);
    printf("%s",separator);

    destroyVector(&vectorFunction);
    destroyVector(&vectorMacro);
    vectorFunction = createNullVector(n,sizeof(float));
    for (size_t i = 0; i < n; ++i) {
        float aux = i;
        VECTOR_GET_ELEMENT_AS(float,vectorFunction,i) = aux;
    }
    printf("Test: erase %lu elements (one by one) at middle from a vector of %lu elements (no resize).\n",
           n_aux/100,n);

    start = clock();
    for (size_t i = 0; i < n_aux/100; ++i) {
        removeElementInVectorAt(vectorFunction,vectorFunction->size/2);
    }
    end = clock();
    vector_function_time = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("results:\nvector_function: %f s\n",vector_function_time);
    destroyVector(&vectorFunction);
    printf("%s",separator);
    vectorFunction = createNullVector(n,sizeof(float));
    printf("Test: erase %lu elements (at same time) at middle from a vector of %lu elements (no resize).\n",
           n_aux/100,n);
    size_t middleIndex = vectorFunction->size/2;
    size_t low = middleIndex-(n_aux/(100*2));
    size_t high = middleIndex+(n_aux/(100*2));
    start = clock();
    removeElementsInVector(vectorFunction,low,high);
    end = clock();
    vector_function_time = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("results:\nvector_function: %f s\n",vector_function_time);
    printf("%s",separator);

    return 0;
}




