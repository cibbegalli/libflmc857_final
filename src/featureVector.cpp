//
// Created by deangeli on 3/27/17.
//

#include "featureVector.h"

FeatureVector* createFeatureVector(int size){
    FeatureVector* featureVector = (FeatureVector*)calloc(1,sizeof(FeatureVector));
    featureVector->size = size;
    featureVector->features = (float*)calloc((size_t)size, sizeof(float));
    return featureVector;
}
void destroyFeatureVector(FeatureVector** featureVector){
    free((*featureVector)->features);
    free((*featureVector));
    (*featureVector) = NULL;
}

FeatureVector* createRandomNormalizedFeatureVector(int size){
    FeatureVector* featureVector = (FeatureVector*)calloc(1,sizeof(FeatureVector));
    featureVector->size = size;
    featureVector->features = (float*)calloc((size_t)size, sizeof(float));
    for (int i = 0; i < size; ++i) {
        featureVector->features[i] = randomNormalized();
    }
    return featureVector;
}

void printFeatureVector(FeatureVector* featureVector){
    if(featureVector == NULL){
        printf("FeatureVector pointer is NULL\n");
        return;
    }
    for (int i = 0; i < featureVector->size; ++i) {
        printf("%f ",featureVector->features[i]);
    }
    printf("\n");
}

FeatureVector* mergeFeatureVectors(FeatureVector* vector1,FeatureVector* vector2){
    FeatureVector *mergedVector = NULL;
    if(vector1 == NULL || vector2 == NULL){
        printf("vector1 or/and vector2 are NULL\n");
        return mergedVector;
    }
    int mergedSize = vector1->size + vector2->size;
    mergedVector = createFeatureVector(mergedSize);
    for (int i = 0; i < vector1->size; ++i) {
        mergedVector->features[i] = vector1->features[i];
    }

    for (int i = 0; i < vector2->size; ++i) {
        mergedVector->features[i+vector1->size] = vector2->features[i];
    }
    return mergedVector;
}

FeatureVector* createFeatureVector(float* vec,int size){
    FeatureVector* featureVector = (FeatureVector*)calloc(1,sizeof(FeatureVector));
    featureVector->size = size;
    featureVector->features = (float*)calloc((size_t)size,sizeof(float));
    for (int i = 0; i < size; ++i) {
        featureVector->features[i] = vec[i];
    }
    return featureVector;
}


void writeFeatureVectors(FeatureVector** vectors, int nVectors, char *filename){
    FILE *fp = fopen(filename,"w");
    for (int i = 0; i < nVectors; ++i) {
        FeatureVector* vec = vectors[i];
        for (int j = 0; j < vec->size; ++j) {
            fprintf(fp,"%f",vec->features[j]);
            if(!(j == vec->size-1)){
                fprintf(fp," ");
            }
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
}


float vectorDifference(FeatureVector* vector1,FeatureVector* vector2){
    if(vector1->size != vector2->size){
        printf("vectors mismatch dimensions\n");
        return -1.0;
    }
    float difference = 0;
    FeatureVector* vector = createFeatureVector(vector1->size);
    float diff;
    for (int i = 0; i < vector1->size; ++i) {
        diff = (vector1->features[i]-vector2->features[i]);
        if(diff < 0){
            diff *= -1;
        }
        vector->features[i] = diff;
        //vector->features[i] = (vector1->features[i]-vector2->features[i])*(vector1->features[i]-vector2->features[i]);
        //vector->features[i] = sqrtf(vector->features[i]);
        difference += vector->features[i];
    }
    destroyFeatureVector(&vector);
    return difference;
}

FeatureVector* copyFeatureVector(FeatureVector* featureVector){
    return createFeatureVector(featureVector->features,featureVector->size);
}



double computeL1NormBetweenFeatureVectors(FeatureVector *featureVector1,FeatureVector *featureVector2){
    double totalRelativeDistance = 0;
    double featDistance = 0;
    for (int i = 0; i < featureVector1->size; ++i) {
        featDistance = (featureVector1->features[i]-featureVector2->features[i]);
        if(featDistance < 0){
            featDistance *= -1;
        }
        totalRelativeDistance += featDistance;
    }
    totalRelativeDistance /= featureVector1->size;
    return totalRelativeDistance;
}

void featureVectorAdditionInplace(FeatureVector *featureVector1,FeatureVector *featureVector2){
    for (int i = 0; i < featureVector2->size; ++i) {
        featureVector2->features[i] = featureVector1->features[i]+featureVector1->features[i];
    }
}

void setValueInFeatureVector(FeatureVector *featureVector,float value){
    for (int i = 0; i < featureVector->size; ++i) {
        featureVector->features[i] = value;
    }
}



FeatureMatrix* createFeatureMatrix(){
    FeatureMatrix* featureMatrix = NULL;
    featureMatrix = (FeatureMatrix*)calloc(1,sizeof(FeatureMatrix));
    featureMatrix->nFeaturesVectors = 0;
    featureMatrix->featureVector = NULL;
    return featureMatrix;
}

FeatureMatrix* createFeatureMatrix(int nFeaturesVectors){
    FeatureMatrix* featureMatrix = NULL;
    featureMatrix = (FeatureMatrix*)calloc(1,sizeof(FeatureMatrix));
    featureMatrix->nFeaturesVectors = nFeaturesVectors;
    featureMatrix->featureVector = (FeatureVector**)calloc((size_t)nFeaturesVectors,sizeof(FeatureVector*));
    return featureMatrix;
}

FeatureMatrix* createFeatureMatrix(int nFeaturesVectors,int vectorSize){
    FeatureMatrix* featureMatrix = NULL;
    featureMatrix = (FeatureMatrix*)calloc(1,sizeof(FeatureMatrix));
    featureMatrix->nFeaturesVectors = nFeaturesVectors;
    featureMatrix->featureVector = (FeatureVector**)calloc((size_t)nFeaturesVectors,sizeof(FeatureVector*));
    for (int i = 0; i < nFeaturesVectors; ++i) {
        featureMatrix->featureVector[i] = createFeatureVector(vectorSize);
    }
    return featureMatrix;
}

void destroyFeatureMatrix(FeatureMatrix** featureMatrix){
    if(*featureMatrix == NULL){
        return ;
    }
    for (int i = 0; i < (*featureMatrix)->nFeaturesVectors; ++i) {
        destroyFeatureVector( &((*featureMatrix)->featureVector[i]) );
    }
    free((*featureMatrix)->featureVector);
    free((*featureMatrix));
    *featureMatrix = NULL;
}

void writeFeatureMatrix(FeatureMatrix* featureMatrix, char *filename){
    FILE *fp = fopen(filename,"w");
    for (int i = 0; i < featureMatrix->nFeaturesVectors; ++i) {
        FeatureVector* vec = featureMatrix->featureVector[i];
        for (int j = 0; j < vec->size; ++j) {
            fprintf(fp,"%f",vec->features[j]);
            if(!(j == vec->size-1)){
                fprintf(fp," ");
            }
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
}

void addNewLines(FeatureMatrix** featureMatrix, int numberNewLines){
    FeatureMatrix* aux = *featureMatrix;

//    int numberLines = aux->nFeaturesVectors+numberNewLines;
//    FeatureVector** copy = (FeatureVector**)calloc(numberLines,sizeof(FeatureVector*));
//    for (int i = 0; i < aux->nFeaturesVectors; ++i) {
//        copy[i] = aux->featureVector[i];
//    }
//    free(aux->featureVector);
//    aux->featureVector = copy;
//    aux->nFeaturesVectors = numberLines;

    int numberLines = aux->nFeaturesVectors+numberNewLines;
    FeatureVector** newRows = (FeatureVector**)realloc(aux->featureVector,(numberLines)*sizeof(FeatureVector*));
    aux->featureVector = newRows;
    //aux->featureVector = (FeatureVector**)realloc(aux->featureVector,(numberLines)*sizeof(FeatureVector*));
    aux->nFeaturesVectors = numberLines;
}

void printFeatureMatrix(FeatureMatrix* featureMatrix){
    if(featureMatrix == NULL){
        printf("FeatureMatrix pointer is NULL\n");
        return;
    }
    for (int i = 0; i < featureMatrix->nFeaturesVectors; ++i) {
        printFeatureVector(featureMatrix->featureVector[i]);
    }
    printf("\n");
}

FeatureMatrix* copyFeatureMatrix(FeatureMatrix* featureMatrix){
    FeatureMatrix *copy = createFeatureMatrix(featureMatrix->nFeaturesVectors);
    for (int i = 0; i < featureMatrix->nFeaturesVectors; ++i) {
        copy->featureVector[i] = copyFeatureVector(featureMatrix->featureVector[i]);
    }
    return copy;
}

int findNearestVectorUsingL1Norm(FeatureMatrix* featureMatrix, FeatureVector *featureVector){
    int nearestClusterIndex = -1;
    double distance;
    double minDistance = DBL_MAX;
    for (int i = 0; i < featureMatrix->nFeaturesVectors; ++i) {
        distance = computeL1NormBetweenFeatureVectors(featureVector,
                                                      featureMatrix->featureVector[i]);

        if(minDistance > distance){
            minDistance = distance;
            nearestClusterIndex = i;
        }

    }
    return nearestClusterIndex;
}

void writeFeatureMatrix(FeatureMatrix *featureMatrix, char *filename, bool sameSize){
    if(featureMatrix == NULL){
        printf("[wirteFeatureMatrix] feature matrix is null");
        return ;
    }
    if(featureMatrix->nFeaturesVectors <= 0){
        printf("[wirteFeatureMatrix] invalid number of features vectors");
        return ;
    }

    FILE *fp = fopen(filename,"w");
    if(sameSize){
        fprintf(fp,",1\n"); //uniform
        fprintf(fp,"%d %d\n",featureMatrix->nFeaturesVectors, featureMatrix->featureVector[0]->size);
        for (int i = 0; i < featureMatrix->nFeaturesVectors; ++i) {
            FeatureVector* vec = featureMatrix->featureVector[i];
            for (int j = 0; j < vec->size; ++j) {
                fprintf(fp,"%f",vec->features[j]);
                if(!(j == vec->size-1)){
                    fprintf(fp," ");
                }
            }
            fprintf(fp,"\n");
        }
    }else{
        fprintf(fp,",2\n"); //not uniform
        fprintf(fp,"%d\n",featureMatrix->nFeaturesVectors);
        for (int i = 0; i < featureMatrix->nFeaturesVectors; ++i) {
            FeatureVector* vec = featureMatrix->featureVector[i];
            fprintf(fp,"%d\n",vec->size);
            for (int j = 0; j < vec->size; ++j) {
                fprintf(fp,"%f",vec->features[j]);
                if(!(j == vec->size-1)){
                    fprintf(fp," ");
                }
            }
            fprintf(fp,"\n");
        }
    }
}




