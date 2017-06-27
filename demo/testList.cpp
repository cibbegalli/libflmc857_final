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

    size_t n = 4000000;

    double percentageInsertMiddle = 0.05;
    clock_t start, end;
    double time_function = 0;
    double time_macro = 0;
    LinkedList *listFunction = createLinkedList(sizeof(float));
    LinkedList *listMacro = createLinkedList(sizeof(float));
    size_t n_aux = n*percentageInsertMiddle;
    const char* separator= "----------------------------------------------------------\n";
    printf("%s",separator);
    printf("workload is: %lu floats\n",n);
    printf("%s",separator);
    printf("Test: push back elements\n");

    start = clock();
    for (size_t i = 0; i < n; ++i) {
        float aux = i;
        pushBackElementInList(listFunction,&aux);
    }
    end = clock();
    time_function = ((double) (end - start)) / CLOCKS_PER_SEC;

    start = clock();
    for (size_t i = 0; i < n; ++i) {
        float aux = i;
        LIST_PUSH_BACK(float,listMacro,aux);
    }
    end = clock();
    time_macro = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("results:\nlist_function: %f s\nlist_macro: %f s\n",time_function,time_macro);
    printf("%s",separator);

    destroyLinkedList(&listFunction);
    destroyLinkedList(&listMacro);
    printf("Test: insert at middle. Using %.2f%% of %lu\n",percentageInsertMiddle*100,n);
    listFunction = createLinkedList(sizeof(float));
    start = clock();
    for (size_t i = 0; i < n_aux; ++i) {
        float aux = i;
        insertElementInListAt(listFunction, &aux, listFunction->length / 2);
    }
    end = clock();
    time_function = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("results:\nlist_function: %f s\n",time_function);
    printf("%s",separator);
    destroyLinkedList(&listFunction);


    printf("Test: deleting %lu at middle from a list of %lu elements. \n",n_aux,n);
    listFunction = createLinkedList(sizeof(float));
    for (size_t i = 0; i < n; ++i) {
        float aux = i;
        pushBackElementInList(listFunction,&aux);
    }

    start = clock();
    for (size_t i = 0; i < n_aux; ++i) {
        removeElementInListAt(listFunction, listFunction->length / 2);
    }
    end = clock();
    time_function = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("results:\nlist_function: %f s\n",time_function);
    printf("%s",separator);


    return 0;
}




