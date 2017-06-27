//
// Created by deangeli on 5/15/17.
//

#ifndef _MATRIX_H
#define _MATRIX_H

#include "common.h"
#include "vector.h"

typedef struct _matrix{
    GVector* matrixData;
    size_t numberRows;
    size_t numberColumns;
    size_t numberElements;
}Matrix;

Matrix* createMatrix(size_t nrows, size_t ncols, size_t dataSize_bytes);
Matrix* createMatrix(size_t nrows, size_t ncols); //default double
Matrix* createMatrix(GVector* vector);
Matrix* stackVerticallyMatrices(Matrix* matrix1,Matrix* matrix2);
Matrix* stackVerticallyMatrices(Matrix* matrix1, GVector* vector);
void reshapeMatrix(Matrix* matrix, size_t nrows, size_t ncols);
void stackVerticallyMatrices(Matrix** matrix1,Matrix* matrix2);
void destroyMatrix(Matrix **pMatrix);
Matrix* copyMatrix(Matrix* matrix);

inline void copyMatrixRow(Matrix* matrix1,Matrix* matrix2, size_t rowIndex1,size_t rowIndex2){
    size_t nbytes =  matrix1->matrixData->elementSize*matrix1->numberColumns;
    size_t memoryShift1 = matrix1->numberColumns*rowIndex1*matrix1->matrixData->elementSize;
    size_t memoryShift2 = matrix2->numberColumns*rowIndex2*matrix2->matrixData->elementSize;
    memcpy( ((unsigned char*)matrix1->matrixData->data)+memoryShift1,
           ((unsigned char*)matrix2->matrixData->data)+memoryShift2,
           nbytes);
}

inline void setRowValueGivenVector(Matrix* matrix, GVector* vector, size_t rowIndex){
    size_t shift =  matrix->matrixData->elementSize*(matrix->numberColumns*rowIndex);
    size_t nbytes = vector->elementSize*vector->size;
    memcpy((unsigned char*)matrix->matrixData->data+shift,(unsigned char*)vector->data, nbytes);
}

//inline void getSubmetrix(Matrix* matrix,size_t startRow,size_t endRow, size_t startCol, size_t endRow){
//
//}

inline void copyMatrixRow(Matrix* matrix1,Matrix* matrix2, size_t rowIndex1){
    size_t rowIndex2 = rowIndex1;
    copyMatrixRow(matrix1,matrix2, rowIndex1,rowIndex2);
}





#define MATRIX_GET_ELEMENT_PO_AS(type, matrix, rowIndex,colIndex) VECTOR_GET_ELEMENT_AS(type,matrix->matrixData, MATRIX_COMPUTEINDEX(matrix,rowIndex,colIndex))

#define MATRIX_GET_ELEMENT_BI_AS(type, matrix, index) VECTOR_GET_ELEMENT_AS(type,matrix->matrixData, index)

#define MATRIX_COMPUTEINDEX(matrix,rowIndex,colIndex) ((rowIndex*matrix->numberColumns)+colIndex)

#define MATRIX_ADD_ROWS_AS(type,matrix1,matrix2,index1,index2) \
    { \
        type* aux1 = ( (type*) matrix1->matrixData->data) + (index1*matrix1->numberColumns); \
        type* aux2 = ( (type*) matrix2->matrixData->data) + (index2*matrix2->numberColumns); \
        for(size_t index=0; index < matrix1->numberColumns; index++) { \
            aux1[index] += aux2[index]; \
        } \
    } \

#define MATRIX_PRINT_AS(type,symbol,matrix) \
    for(size_t i = 0; i< matrix->numberRows; i++){\
        for(size_t j = 0; j< matrix->numberColumns; j++){\
            printf(symbol, MATRIX_GET_ELEMENT_PO_AS(type,matrix,i,j) ); \
        } \
        printf("\n"); \
    }
#endif //_MATRIX_H
