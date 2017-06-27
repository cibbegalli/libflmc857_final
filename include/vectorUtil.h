//
// Created by deangeli on 5/21/17.
//

#ifndef LIBFL_VECTORUTIL_H
#define LIBFL_VECTORUTIL_H

#include "vector.h"

inline double * myInsertionSort(double* vector, size_t n,size_t * indecesOrdered){
    double* auxVector = (double*)calloc(n,sizeof(double));
    double aux;
    size_t auxIndex;
    for (size_t i = 0; i < n; ++i) {
        auxVector[i] = vector[i];
        indecesOrdered[i] = i;
        for (size_t j = i; j > 0; --j) {

            if(auxVector[j-1] > auxVector[j]){
                aux = auxVector[j];
                auxVector[j] = auxVector[j-1];
                auxVector[j-1] =  aux;

                auxIndex = indecesOrdered[j];
                indecesOrdered[j] = indecesOrdered[j-1];
                indecesOrdered[j-1] = auxIndex;
            }else{
                break;
            }

        }
    }
    return auxVector;
}

inline void myInsertionSortInplace(double* vector, size_t n,size_t * indecesOrdered){
    double aux;
    size_t auxIndex;
    for (size_t i = 0; i < n; ++i) {
        indecesOrdered[i] = i;
        for (size_t j = i; j > 0; --j) {

            if(vector[j-1] > vector[j]){
                aux = vector[j];
                vector[j] = vector[j-1];
                vector[j-1] =  aux;

                auxIndex = indecesOrdered[j];
                indecesOrdered[j] = indecesOrdered[j-1];
                indecesOrdered[j-1] = auxIndex;
            }else{
                break;
            }

        }
    }
}

inline void myInsertionSortInplace(int* vector, size_t n){
    int aux;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i; j > 0; --j) {

            if(vector[j-1] > vector[j]){
                aux = vector[j];
                vector[j] = vector[j-1];
                vector[j-1] =  aux;
            }else{
                break;
            }

        }
    }
}



inline GVector* findUniquesInIntegerVector(GVector* vector,bool ordered){
    GVector* uniques = createNullVector(vector->size,sizeof(int));

    VECTOR_GET_ELEMENT_AS(int,uniques,0) = VECTOR_GET_ELEMENT_AS(int,vector,0);
    size_t sizeUniquesVector = 1;
    size_t countEquals;
    for (size_t vectorIndex = 1; vectorIndex < vector->size; ++vectorIndex) {
        countEquals = 0;
        for (size_t uniquesVectorIndex = 0; uniquesVectorIndex < sizeUniquesVector; ++uniquesVectorIndex) {
            if(VECTOR_GET_ELEMENT_AS(int,vector,vectorIndex) ==
               VECTOR_GET_ELEMENT_AS(int,uniques,uniquesVectorIndex)){
                countEquals++;
            }
        }
        if(countEquals == 0){
            VECTOR_GET_ELEMENT_AS(int,uniques,sizeUniquesVector) = VECTOR_GET_ELEMENT_AS(int,vector,vectorIndex);
            sizeUniquesVector++;
        }
    }
    uniques->size = sizeUniquesVector;
    shrinkToFit(uniques);
    if(ordered){
        myInsertionSortInplace((int*)uniques->data,uniques->size);
    }
    return uniques;
}

#endif //LIBFL_VECTORUTIL_H
