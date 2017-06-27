/*
 *Created by Deangeli Gomes Neves
 *
 * This software may be freely redistributed under the terms
 * of the MIT license.
 *
 */

#ifndef _FEATUREVECTOR_H_
#define _FEATUREVECTOR_H_

#include "common.h"

typedef struct _featureVector {
    float* features;
    int size;
} FeatureVector;

/*
 * estrutura para armazenar varios features vectors
 * */
typedef struct _featureMatrix {
    FeatureVector **featureVector;
    int nFeaturesVectors;
} FeatureMatrix;


FeatureVector* createFeatureVector(int size);
FeatureVector* createFeatureVector(float* vec,int size);
FeatureVector* createRandomNormalizedFeatureVector(int size);

/*
 * escreve uma arquivo txt onde caada linha no arquivo e um feature vector
 * */
void writeFeatureVector(FeatureVector* vector, FILE *fp);
void writeFeatureVectors(FeatureVector** vectors, int nVectors, char *filename);
void destroyFeatureVector(FeatureVector** vector);
void printFeatureVector(FeatureVector* featureVector);
FeatureVector* mergeFeatureVectors(FeatureVector* vector1,FeatureVector* vector2);
float vectorDifference(FeatureVector* vector1,FeatureVector* vector2);
FeatureVector* copyFeatureVector(FeatureVector* featureVector);
//compute the discrepancy between an exact value(featureVector1) and some approximation(featureVector2) to it
double computeL1NormBetweenFeatureVectors(FeatureVector *featureVector1,FeatureVector *featureVector2);
//add featureVector1 to featureVector2, and store the result in featureVector2
void featureVectorAdditionInplace(FeatureVector *featureVector1,FeatureVector *featureVector2);
//All elements in featureVector will be equal to a given value
void setValueInFeatureVector(FeatureVector *featureVector,float value);
/*
 * escreve uma arquivo txt onde cada linha no arquivo e um feature vector
 * */
FeatureMatrix* createFeatureMatrix();
FeatureMatrix* createFeatureMatrix(int nFeaturesVectors);
FeatureMatrix* createFeatureMatrix(int nFeaturesVectors,int vectorSize);
FeatureMatrix* copyFeatureMatrix(FeatureMatrix* featureMatrix);
void addNewLines(FeatureMatrix** featureMatrix, int numberNewLines);
void printFeatureMatrix(FeatureMatrix* featureMatrix);
void wirteFeatureMatrix(FeatureMatrix* featureMatrix, char *filename);
void destroyFeatureMatrix(FeatureMatrix** featureMatrix);
int findNearestVectorUsingL1Norm(FeatureMatrix* featureMatrix, FeatureVector *featureVector);
void writeFeatureMatrix(FeatureMatrix* featureMatrix, char *filename);
void writeFeatureMatrix(FeatureMatrix *featureMatrix, char *filename, bool sameSize);


void sortAt(FeatureVector featureVector, int lastIndex);

#endif //LIBFL_FEATUREVECTOR_H


