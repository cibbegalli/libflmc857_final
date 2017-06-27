#include "FL.h"
#include <time.h>


int main(int argc, char **argv) {

    Matrix* matrix1 = createMatrix(4,4,sizeof(float));
    Matrix* matrix2 = createMatrix(10,4,sizeof(float));
    for (size_t i = 0; i < matrix1->numberElements ; ++i) {
        MATRIX_GET_ELEMENT_BI_AS(float,matrix1,i) = i;
    }

    for (size_t i = 0; i < matrix2->numberRows ; ++i) {
        for (size_t j = 0; j < matrix2->numberColumns; ++j) {
            MATRIX_GET_ELEMENT_PO_AS(float,matrix2,i,j) = (i+j)+10;
        }
    }


    for (size_t i = 0; i < matrix1->numberRows ; ++i) {
        for (size_t j = 0; j < matrix1->numberColumns; ++j) {
            printf("%f ",MATRIX_GET_ELEMENT_PO_AS(float,matrix1,i,j));
        }
        printf("\n");
    }

    printf("\n");
    for (size_t i = 0; i < matrix2->numberRows ; ++i) {
        for (size_t j = 0; j < matrix2->numberColumns; ++j) {
            printf("%f ",MATRIX_GET_ELEMENT_PO_AS(float,matrix2,i,j));
        }
        printf("\n");
    }
    printf("\n");

    Matrix* matrix3 = stackVerticallyMatrices(matrix1,matrix2);

    for (size_t i = 0; i < matrix3->numberRows ; ++i) {
        for (size_t j = 0; j < matrix3->numberColumns; ++j) {
            printf("%f ",MATRIX_GET_ELEMENT_PO_AS(float,matrix3,i,j));
        }
        printf("\n");
    }
    printf("\n");

    stackVerticallyMatrices(&matrix1,matrix2);

    for (size_t i = 0; i < matrix1->numberRows ; ++i) {
        for (size_t j = 0; j < matrix1->numberColumns; ++j) {
            printf("%f ",MATRIX_GET_ELEMENT_PO_AS(float,matrix1,i,j));
        }
        printf("\n");
    }

    printf("\n");

    copyMatrixRow(matrix1,matrix2,1,4);

    for (size_t i = 0; i < matrix1->numberRows ; ++i) {
        for (size_t j = 0; j < matrix1->numberColumns; ++j) {
            printf("%f ",MATRIX_GET_ELEMENT_PO_AS(float,matrix1,i,j));
        }
        printf("\n");
    }

    printf("\n");

    return 0;
}




