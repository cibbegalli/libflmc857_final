//
// Created by deangeli on 5/16/17.
//

#ifndef _DISTANCEFUNCTIONS_H_
#define _DISTANCEFUNCTIONS_H_

#include "vector.h"
#include "argumentList.h"

typedef double (*DistanceFunction)(float* vector1, float* vector2,size_t numberElements,ArgumentList* argumentList);

inline double  computeL1Norm(float* vector1,float* vector2,size_t  numberElements){
    float diff = 0;
    double total = 0;
    for (size_t i = 0; i < numberElements; ++i) {
        diff = vector1[i] - vector2[i];
        diff = (diff < 0) ? -diff : diff;
        total += diff;
    }
    return  total;
}

inline double  computeNormalizedL1Norm(float* vector1,float* vector2, size_t  numberElements){
    float diff = 0;
    double total = 0;
    for (size_t i = 0; i < numberElements; ++i) {
        diff = vector1[i] - vector2[i];
        diff = (diff < 0) ? -diff : diff;
        total += diff;
    }
    total /= (numberElements+0.000001);
    return  total;
}

inline double computeL1Norm(float* vector1,float* vector2, size_t  numberElements,ArgumentList* argumentList){
    return computeL1Norm(vector1,vector2,numberElements);
}

inline double computeNormalizedL1Norm(float* vector1,float* vector2, size_t  numberElements,ArgumentList* argumentList){
    return computeNormalizedL1Norm(vector1,vector2,numberElements);
}

inline double  computeL2Norm(float* vector1,float* vector2,size_t  numberElements){
    float diff = 0;
    double total = 0;
    for (size_t i = 0; i < numberElements; ++i) {
        diff = vector1[i] - vector2[i];
        diff = (diff*diff);
        total += diff;
    }
    total = sqrt(total);
    return  total;
}

inline double computeNormalizedL2Norm(float* vector1,float* vector2,size_t  numberElements){
    float diff = 0;
    double total = 0;
    for (size_t i = 0; i < numberElements; ++i) {
        diff = vector1[i] - vector2[i];
        diff = (diff*diff);
        total += diff;
    }
    total /= (numberElements+0.000001);
    total = sqrt(total);
    return  total;
}

inline double computeL2Norm(float* vector1,float* vector2, size_t  numberElements,ArgumentList* argumentList){
    return computeL2Norm(vector1,vector2,numberElements);
}

inline double computeNormalizedL2Norm(float* vector1,float* vector2, size_t  numberElements,ArgumentList* argumentList){
    return computeNormalizedL2Norm(vector1,vector2,numberElements);
}



#endif //LIBFL_DISTANCEFUNCTIONS_H
