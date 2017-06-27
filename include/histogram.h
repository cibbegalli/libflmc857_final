/*
 *Created by Deangeli Gomes Neves
 *
 * This software may be freely redistributed under the terms
 * of the MIT license.
 *
 */

#ifndef _HISTOGRAM_H_
#define _HISTOGRAM_H_

#include "image.h"
#include "common.h"
#include "adjacencyRelation.h"
#include "featureVector.h"
#include "vector.h"

typedef struct _histogram {
    float *val;
    int    n;
    int binSize;
} Histogram;


void writeHistogram(Histogram *hist, char *filename);
void destroyHistogram(Histogram **hist);
Histogram* createHistogram(int n);
Histogram* computeHistogram(Image *image,float binSize,bool normalization);
Image *ProbabilityDensityFunction(Image *image, double standardDeviation);
FeatureVector* createFeatureVector(Histogram *histogram);
FeatureVector* computeHistogramForFeatureVector(Image *image,float binSize, bool normalization);

//FeatureVector* computeHistogramForFeatureVectorGivenNBins(Image *image,int numberBinsPerChannel, bool normalization);
GVector* computeHistogramForFeatureVectorGivenNBins(Image *image,int numberBinsPerChannel, bool normalization);

GVector* computeHistogramOfOrientedGradient(Image *img, Image *mag, Image *phase, int bx, int by, int theta);


#endif //_HISTOGRAM_H
