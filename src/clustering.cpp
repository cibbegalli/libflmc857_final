//
// Created by deangeli on 5/20/17.
//

#include "clustering.h"

Matrix* kmeansClustering(Matrix* featureMatrix, size_t nClusters,
                         size_t iterMax,
                         double tolerance, int seed,
                         DistanceFunction distanceFunction,
                         ArgumentList* argumentListDistanceFunction){
    Matrix* centroides_old = createMatrix(nClusters,featureMatrix->numberColumns,sizeof(float));
    if(nClusters >= featureMatrix->numberRows){
        float* minimuns = (float*)calloc(featureMatrix->numberColumns,sizeof(float));
        float* maximums = (float*)calloc(featureMatrix->numberColumns,sizeof(float));
        for (size_t col = 0; col < featureMatrix->numberColumns; ++col) {
            minimuns[col] = MATRIX_GET_ELEMENT_PO_AS(float,featureMatrix,0,col);
            maximums[col] = MATRIX_GET_ELEMENT_PO_AS(float,featureMatrix,0,col);
            for (size_t row = 0; row < featureMatrix->numberRows; ++row) {
                if(minimuns[col] > MATRIX_GET_ELEMENT_PO_AS(float,featureMatrix,row,col) ){
                    minimuns[col] =  MATRIX_GET_ELEMENT_PO_AS(float,featureMatrix,row,col);
                }
                if(maximums[col] < MATRIX_GET_ELEMENT_PO_AS(float,featureMatrix,row,col)){
                    maximums[col] = MATRIX_GET_ELEMENT_PO_AS(float,featureMatrix,row,col);
                }
            }
        }

        for (size_t row = 0; row < featureMatrix->numberRows; ++row) {
            for (size_t col = 0; col < featureMatrix->numberColumns; ++col) {
                MATRIX_GET_ELEMENT_PO_AS(float,centroides_old,row,col) =
                        MATRIX_GET_ELEMENT_PO_AS(float,featureMatrix,row,col);
            }
        }

        int diffRow = nClusters - featureMatrix->numberRows;
        for (int row = 0; row < diffRow; ++row) {
            for (size_t col = 0; col < featureMatrix->numberColumns; ++col) {
                MATRIX_GET_ELEMENT_PO_AS(float,centroides_old,row,col) = randomNumber(minimuns[col],maximums[col]);
            }
        }

        free(minimuns);
        free(maximums);
        return centroides_old;
    }

    size_t k = 0;
    bool *isUsed = (bool*)calloc(featureMatrix->numberRows,sizeof(*isUsed));
    while (k < nClusters) {
        int randomIndex = RandomInteger(0,featureMatrix->numberRows-1);
        if(isUsed[randomIndex] == false){
            copyMatrixRow(centroides_old,featureMatrix,k,randomIndex);
            isUsed[randomIndex] = true;
            k++;
        }
    }
    free(isUsed);
    unsigned int printStatusValue = iterMax*0.02;
    Matrix* centroides_new = createMatrix(centroides_old->numberRows, centroides_old->numberColumns,sizeof(float));
    int* count = (int*)calloc(nClusters,sizeof(int));
    for (size_t j = 0; j < nClusters; ++j) {
        count[j] = 0;
    }
    size_t iter = 0;
    double totalDiff = DBL_MAX;
    float currentDiff = 0;
    while( (iter < iterMax) && (totalDiff > tolerance) ){
        currentDiff = 0;
        totalDiff = 0;
        for (size_t featureVectorIndex = 0; featureVectorIndex < featureMatrix->numberRows; ++featureVectorIndex) {
            size_t nearestClusterIndex = findNearestRow(centroides_old,
                                                        featureMatrix,featureVectorIndex,
                                                        distanceFunction,argumentListDistanceFunction);
            MATRIX_ADD_ROWS_AS(float,centroides_new,featureMatrix,nearestClusterIndex,featureVectorIndex);
            count[nearestClusterIndex]++;
        }
        for (size_t i = 0; i < centroides_new->numberRows; ++i) {
            for (size_t j = 0; j < centroides_new->numberColumns; ++j) {
                MATRIX_GET_ELEMENT_PO_AS(float,centroides_new,i,j) =
                        (MATRIX_GET_ELEMENT_PO_AS(float,centroides_new,i,j) +
                         MATRIX_GET_ELEMENT_PO_AS(float,centroides_old,i,j))/(count[i]+1);
                currentDiff = MATRIX_GET_ELEMENT_PO_AS(float,centroides_new,i,j) -
                              MATRIX_GET_ELEMENT_PO_AS(float,centroides_old,i,j);
                currentDiff = fabsf(currentDiff);
                totalDiff += currentDiff;
                MATRIX_GET_ELEMENT_PO_AS(float,centroides_old,i,j) = MATRIX_GET_ELEMENT_PO_AS(float,centroides_new,i,j);
                MATRIX_GET_ELEMENT_PO_AS(float,centroides_new,i,j) = 0.0f;
            }
            count[i] = 0;
            totalDiff = totalDiff/centroides_new->numberRows;
        }
        if(iter % printStatusValue == 0 || iter == 0){
            printf("[kmeansClustering] At iteration %lu, the total displacement is %f\n",iter,totalDiff);
        }

        iter++;
    }
    printf("[kmeansClustering] At iteration %lu, the total displacement is %f\n",(iter-1),totalDiff);
    destroyMatrix(&centroides_new);

    free(count);
    return centroides_old;
}
