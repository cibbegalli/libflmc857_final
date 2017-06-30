/*
 *Created by Deangeli Gomes Neves
 *
 * This software may be freely redistributed under the terms
 * of the MIT license.
 *
 */

#ifndef _BAGOFVISUALWORDS_H_
#define _BAGOFVISUALWORDS_H_

#include "common.h"
#include "featureVector.h"
#include "file.h"
#include "image.h"
#include "argumentList.h"
#include "vector.h"
#include "matrix.h"
#include "matrixUtil.h"
#include "distanceFunctions.h"
#include "sampling.h"
#include "featureExtractor.h"
#include "classifiers.h"
#include "clustering.h"
#include "featureVector.h"



typedef struct _bagOfVisualWordsManager BagOfVisualWordsManager;



typedef GVector* (*ImageSamplerFunction)(Image* image, Image *isf, BagOfVisualWordsManager* bagOfVisualWordsManager);
typedef Matrix* (*FeatureExtractorFunction)(GVector* outputSampler, Image* image, Image* isf, BagOfVisualWordsManager* bagOfVisualWordsManager);
typedef Matrix* (*ClusteringFunction)(Matrix* outputFeatureExtractor_allSamples, BagOfVisualWordsManager* bagOfVisualWordsManager);
typedef GVector* (*MountHistogramFunction) (Matrix* outputFeatureExtractor_singleSample,BagOfVisualWordsManager* bagOfVisualWordsManager);

typedef struct _bagOfVisualWordsManager {
        GVector* pathsToImages_dictionery;
        GVector* pathsToImages_train;
        GVector* pathsToImages_test;
        GVector* pathsToImages_isf;
        Matrix* dictionery;

        Matrix* histogramsTraining;
        GVector* labelsTraining;
        Matrix* histogramsPredictSamples;
        GVector* labelsPredicted;

        bool storeTrainData;
        bool storePredictedData;

        void* classifier;

        FreeFunction freeFunction2SamplerOutput;
        FreeFunction freeFunctionClassifier;

        ArgumentList* argumentListOfSampler;
        ArgumentList* argumentListOfFeatureExtractor;
        ArgumentList* argumentListOfClustering;
        ArgumentList* argumentListOfDistanceFunction;
        ArgumentList* argumentListOfHistogramMounter;

        FeatureExtractorFunction featureExtractorFunction;
        ImageSamplerFunction imageSamplerFunction;
        DistanceFunction distanceFunction;
        ClusteringFunction clusteringFunction;
        MountHistogramFunction mountHistogramFunction;
        FitFunction fitFunction;
        PredictFunction predictFunction;
}BagOfVisualWordsManager;

BagOfVisualWordsManager* createBagOfVisualWordsManager();
void destroyBagOfVisualWordsManager(BagOfVisualWordsManager** pBagOfVisualWordsManager);

//
GVector* gridSamplingBow(Image* image, Image* isf, BagOfVisualWordsManager* bagOfVisualWordsManager);
GVector* samplingByISFBow(Image* image, Image* isf, BagOfVisualWordsManager* bagOfVisualWordsManager);
GVector* samplingByISFColorBow(Image* image, Image* isf, BagOfVisualWordsManager* bagOfVisualWordsManager);

Matrix* computeColorHistogramBow(GVector* vector, Image* image, Image* isf,BagOfVisualWordsManager* bagOfVisualWordsManager);
//
Matrix* computeHogBow(GVector* vector, Image* image, Image* isf, BagOfVisualWordsManager* bagOfVisualWordsManager);

Matrix* computeColorAndHogBow(GVector* vector, Image* image, Image* isf, BagOfVisualWordsManager* bagOfVisualWordsManager);

//
Matrix* kmeansClusteringBow(Matrix* featureMatrix, BagOfVisualWordsManager* bagOfVisualWordsManager);
//
//
GVector* computeCountHistogram_bow(Matrix* featureMatrix,BagOfVisualWordsManager* bagOfVisualWordsManager);
//



void computeDictionery(BagOfVisualWordsManager* bagOfVisualWordsManager);
void trainClassifier(BagOfVisualWordsManager* bagOfVisualWordsManager);
GVector* predictLabels(BagOfVisualWordsManager* bagOfVisualWordsManager);






#endif //LIBFL_BAGOFVISUALWORDS_H
