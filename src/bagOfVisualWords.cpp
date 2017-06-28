//
// Created by deangeli on 4/7/17.
//

#include "bagOfVisualWords.h"


BagOfVisualWordsManager* createBagOfVisualWordsManager(){
    BagOfVisualWordsManager* bowManager = (BagOfVisualWordsManager*)calloc(1,sizeof(BagOfVisualWordsManager));
    bowManager->pathsToImages_dictionery = NULL;
    bowManager->pathsToImages_train = NULL;
    bowManager->pathsToImages_test = NULL;
    bowManager->dictionery = NULL;
    bowManager->featureExtractorFunction = NULL;
    bowManager->imageSamplerFunction = NULL;
    bowManager->clusteringFunction = NULL;
    bowManager->argumentListOfFeatureExtractor =NULL;
    bowManager->argumentListOfSampler = NULL;
    bowManager->argumentListOfClustering = NULL;
    bowManager->argumentListOfDistanceFunction = NULL;
    bowManager->freeFunction2SamplerOutput = NULL;
    bowManager->classifier = NULL;
    bowManager->freeFunctionClassifier = NULL;
    bowManager->histogramsTraining = NULL;
    bowManager->labelsTraining = NULL;
    bowManager->storeTrainData = false;
    bowManager->storePredictedData = false;
    return bowManager;
}

void destroyBagOfVisualWordsManager(BagOfVisualWordsManager** pBagOfVisualWordsManager){
    BagOfVisualWordsManager* aux = *pBagOfVisualWordsManager;
    if(aux == NULL){
        return;
    }

    destroyVector(&(aux->pathsToImages_dictionery));
    destroyVector(&(aux->pathsToImages_train));
    destroyVector(&(aux->pathsToImages_test));
    destroyVector(&(aux->labelsTraining));
    destroyVector(&(aux->labelsPredicted));

    destroyArgumentList(&(aux->argumentListOfSampler));
    destroyArgumentList(&(aux->argumentListOfFeatureExtractor));
    destroyArgumentList(&(aux->argumentListOfClustering));
    destroyArgumentList(&(aux->argumentListOfDistanceFunction));
    destroyArgumentList(&(aux->argumentListOfHistogramMounter));

    destroyMatrix(&(aux->dictionery));
    destroyMatrix(&(aux->histogramsTraining));
    destroyMatrix(&(aux->histogramsPredictSamples));


    if(aux->freeFunctionClassifier){
        if(aux->classifier){
            aux->freeFunctionClassifier(aux->classifier);
        }
    }


    free(*pBagOfVisualWordsManager);
}

//void setDefaultParameters(BagOfVisualWordsManager* bowManager){
//    if(bowManager == NULL){
//        printf("[setDefaultParameters] bowManager is null");
//        return;
//    }
//
//    ArgumentList* gridSamplingArguments = createArgumentList();
//    ARGLIST_PUSH_BACK_AS(size_t,gridSamplingArguments,64); //patch size X
//    ARGLIST_PUSH_BACK_AS(size_t,gridSamplingArguments,64); //patch size Y
//    bowManager->argumentListOfSampler = gridSamplingArguments;
//    bowManager->imageSamplerFunction = gridSampling;
//
//    //Neste exemplo eu usar o histograma de cores
//    ArgumentList* colorFeatureExtractorArguments = createArgumentList();
//    ARGLIST_PUSH_BACK_AS(size_t,colorFeatureExtractorArguments,8); //nBins per channel
//    ARGLIST_PUSH_BACK_AS(size_t,colorFeatureExtractorArguments,8*8*8); //total number of channels
//    bowManager->featureExtractorFunction = computeColorHistogram;
//
//}

GVector* gridSamplingBow(Image* image, Image* isf, BagOfVisualWordsManager* bagOfVisualWordsManager){
    ArgumentList* argumentList = bagOfVisualWordsManager->argumentListOfSampler;

    if(argumentList->length == 0){
        printf("[gridSampling] invalid argument list");
    }else if(argumentList->length == 1){
        size_t patchSize = ARGLIST_GET_ELEMENT_AS(size_t ,argumentList,0);
        Image* img = convertRGBtoYCbCr(image);
        GVector* patchs=gridSampling(img,patchSize,patchSize);
        destroyImage(&img);
        return patchs;
    }else if(argumentList->length == 2){
        size_t patchSizeX = ARGLIST_GET_ELEMENT_AS(size_t,argumentList,0);
        size_t patchSizeY = ARGLIST_GET_ELEMENT_AS(size_t,argumentList,1);
        Image* img = convertRGBtoYCbCr(image);
        GVector* patchs= gridSampling(img,patchSizeX,patchSizeY);
        destroyImage(&img);
        return patchs;
    }
    return NULL;
}

GVector* samplingByISFBow(Image* image, Image* isf, BagOfVisualWordsManager* bagOfVisualWordsManager){
    ArgumentList* argumentList = bagOfVisualWordsManager->argumentListOfSampler;

    if(argumentList->length == 0){
        printf("[samplingByISF] invalid argument list");
    }else if(argumentList->length == 1){
        size_t patchSize = ARGLIST_GET_ELEMENT_AS(size_t ,argumentList,0);
        return gridSampling(image,patchSize,patchSize);
    }else if(argumentList->length == 2){
        size_t patchSizeX = ARGLIST_GET_ELEMENT_AS(size_t,argumentList,0);
        size_t patchSizeY = ARGLIST_GET_ELEMENT_AS(size_t,argumentList,1);
        return gridSampling(image,patchSizeX,patchSizeY);
    }
    return NULL;
}

Matrix* computeColorHistogramBow(GVector* vector, Image* image, BagOfVisualWordsManager* bagOfVisualWordsManager){
    ArgumentList* argumentList = bagOfVisualWordsManager->argumentListOfFeatureExtractor;
    if(argumentList->length < 2){
        printf("[computeColorHistogram] invalid argument list");
        return NULL;
    }
    if(vector->size == 0){
        printf("[computeColorHistogram] vector has 0 elements");
        return NULL;
    }
    size_t nbinsPerChannel = ARGLIST_GET_ELEMENT_AS(size_t,argumentList,0);
    size_t totalBins = ARGLIST_GET_ELEMENT_AS(size_t,argumentList,1);
    return computeColorHistogram(vector,nbinsPerChannel,totalBins);
}
//Matrix* computeColorHistogram(GVector* vector_images,size_t nbinsPerChannel,size_t totalNumberBins)
//Matrix* computeHog(GVector* vector_images, Image *mag, Image *phase, int blocks_x, int blocks_y, int theta) 

Matrix* computeHogBow(GVector* vector, Image* image, BagOfVisualWordsManager* bagOfVisualWordsManager){
    ArgumentList* argumentList = bagOfVisualWordsManager->argumentListOfFeatureExtractor;
    if(argumentList->length < 3){
        printf("[computeHog] invalid argument list");
        return NULL;
    }
    if(vector->size == 0){
        printf("[computeHog] vector has 0 elements");
        return NULL;
    }
    Image* mag;
    Image* phase;
    computeGradient(image, &mag, &phase);

    int blocks_x = ARGLIST_GET_ELEMENT_AS(size_t,argumentList,0);
    int blocks_y = ARGLIST_GET_ELEMENT_AS(size_t,argumentList,1);
    int theta = ARGLIST_GET_ELEMENT_AS(size_t,argumentList,2);
    Matrix* m = computeHog(vector, mag, phase, blocks_x, blocks_y, theta);
    destroyImage(&mag);
    destroyImage(&phase);

    return m;
}


void computeDictionery(BagOfVisualWordsManager* bagOfVisualWordsManager){
    Matrix* allFeatures = NULL;
    if(!bagOfVisualWordsManager->imageSamplerFunction){
        printf("[computeDictionery] Sampler function not defined\n");
    }


    if(!bagOfVisualWordsManager->featureExtractorFunction){
        printf("[computeDictionery] Feature extractor function not defined\n");
        return;
    }
    if(!bagOfVisualWordsManager->clusteringFunction){
        printf("[computeDictionery] Clustering function not defined\n");
        return;
    }

    printf("[computeDictionery] Generating visual words...\n");
    for (size_t i = 0; i < bagOfVisualWordsManager->pathsToImages_dictionery->size; ++i) {
        char* imagePath = VECTOR_GET_ELEMENT_AS(char*,bagOfVisualWordsManager->pathsToImages_dictionery,i);
        Image* image = readImage(imagePath);

        Image* imageISF;

        if(bagOfVisualWordsManager->pathsToImages_isf != NULL) {
            char* isfPath = VECTOR_GET_ELEMENT_AS(char*, bagOfVisualWordsManager->pathsToImages_isf, i);
            imageISF = readImage(isfPath);
        }
        

        if(image == NULL){
            printf("[computeDictionery] invalid image path: %s",imagePath);
            continue;
        }

        //if(imageISF == NULL) {
        //    printf("[computeDictionery] invalid image path for superpixel: %s\n", isfPath);
        //    continue;
        //}

        GVector* samplingResults = NULL;
        if(bagOfVisualWordsManager->imageSamplerFunction){
            samplingResults = bagOfVisualWordsManager->imageSamplerFunction(image, imageISF,
                                                                            bagOfVisualWordsManager);
            samplingResults->freeFunction = bagOfVisualWordsManager->freeFunction2SamplerOutput;
        }else{
            samplingResults = createNullVector(1,sizeof(Image*));
            VECTOR_GET_ELEMENT_AS(Image*,samplingResults,0) = image;
        }
        Matrix* featureMatrix = bagOfVisualWordsManager->featureExtractorFunction(samplingResults, image,
                                                                                  bagOfVisualWordsManager);
        Matrix* newData = stackVerticallyMatrices(allFeatures,featureMatrix);
        destroyMatrix(&allFeatures);
        allFeatures = newData;
        destroyImage(&image);
        destroyVector(&samplingResults);
        destroyMatrix(&featureMatrix);
    }

    printf("[computeDictionery] Finding Visual words...\n");
    bagOfVisualWordsManager->dictionery = bagOfVisualWordsManager->clusteringFunction(allFeatures,
                                                bagOfVisualWordsManager);
    destroyMatrix(&allFeatures);
    printf("[computeDictionery] Dictioney computed\n");
}

void trainClassifier(BagOfVisualWordsManager* bagOfVisualWordsManager){
    if(!bagOfVisualWordsManager->imageSamplerFunction){
        printf("[trainClassifier] Sampler function not defined\n");
    }


    if(!bagOfVisualWordsManager->featureExtractorFunction){
        printf("[trainClassifier] Feature extractor function not defined\n");
        return;
    }
    if(bagOfVisualWordsManager->dictionery == NULL){
        printf("[trainClassifier] Dictionery is empty\n");
        return;
    }
    if(bagOfVisualWordsManager->mountHistogramFunction == NULL){
        printf("[trainClassifier] Mounter histogram function not defined\n");
        return;
    }
    if(bagOfVisualWordsManager->fitFunction == NULL){
        printf("[trainClassifier] Fit function not defined\n");
        return;
    }


    Matrix *bowHistograms = createMatrix(bagOfVisualWordsManager->pathsToImages_train->size,
                                         bagOfVisualWordsManager->dictionery->numberRows,
                                         sizeof(float));
    GVector *imagesLabels = createNullVector(bagOfVisualWordsManager->pathsToImages_train->size,sizeof(int));
    //Matrix *bowHistograms = NULL;
    printf("[trainClassifier] Generating histograms and labels from images\n");
    for (size_t index = 0; index < bagOfVisualWordsManager->pathsToImages_train->size; ++index) {
        char* imagePath = VECTOR_GET_ELEMENT_AS(char*,bagOfVisualWordsManager->pathsToImages_train,index);
        Image* image = readImage(imagePath);


        Image* imageISF;
        if(bagOfVisualWordsManager->pathsToImages_isf != NULL) {
            char* isfPath = VECTOR_GET_ELEMENT_AS(char*, bagOfVisualWordsManager->pathsToImages_isf, index);  
            imageISF = readImage(isfPath);  
        }
        
        


        if(image == NULL){
            printf("[computeDictionery] Invalid image path: %s",imagePath);
            continue;
        }
        GVector* samplingResults = NULL;
        if(bagOfVisualWordsManager->imageSamplerFunction){
            samplingResults = bagOfVisualWordsManager->imageSamplerFunction(image, imageISF,
                                                                            bagOfVisualWordsManager);
            samplingResults->freeFunction = bagOfVisualWordsManager->freeFunction2SamplerOutput;
        }else{
            samplingResults = createNullVector(1,sizeof(Image*));
            VECTOR_GET_ELEMENT_AS(Image*,samplingResults,0) = image;
        }
        Matrix* featureMatrix = bagOfVisualWordsManager->featureExtractorFunction(samplingResults, image, bagOfVisualWordsManager);

        GVector* histogram = bagOfVisualWordsManager->mountHistogramFunction(featureMatrix,bagOfVisualWordsManager);
        setRowValueGivenVector(bowHistograms,histogram,index);
        VECTOR_GET_ELEMENT_AS(int,imagesLabels,index) = findTrueLabelInName(imagePath);

        destroyImage(&image);
        destroyMatrix(&featureMatrix);
        destroyVector(&histogram);
        destroyVector(&samplingResults);
    }

    printf("[trainClassifier] Histograms and labels generated\n");
    printf("[trainClassifier] Classifier  trainning...\n");
    bagOfVisualWordsManager->fitFunction(bowHistograms,imagesLabels,bagOfVisualWordsManager->classifier);
    printf("[trainClassifier] Classifier trained...\n");

    if(bagOfVisualWordsManager->storeTrainData){
        if(bagOfVisualWordsManager->histogramsTraining){
            destroyMatrix(&(bagOfVisualWordsManager->histogramsTraining));
        }
        bagOfVisualWordsManager->histogramsTraining = bowHistograms;

        if(bagOfVisualWordsManager->labelsTraining){
            destroyVector(&(bagOfVisualWordsManager->labelsTraining));
        }
        bagOfVisualWordsManager->labelsTraining = imagesLabels;
    }else{
        destroyMatrix(&bowHistograms);
        destroyVector(&imagesLabels);
    }
}

GVector* predictLabels(BagOfVisualWordsManager* bagOfVisualWordsManager){
    if(!bagOfVisualWordsManager->imageSamplerFunction){
        printf("[predictLabels] Sampler function not defined\n");
    }


    if(!bagOfVisualWordsManager->featureExtractorFunction){
        printf("[predictLabels] Feature extractor function not defined\n");
        return NULL;
    }
    if(bagOfVisualWordsManager->dictionery == NULL){
        printf("[predictLabels] Dictionery is empty\n");
        return NULL;
    }
    if(bagOfVisualWordsManager->mountHistogramFunction == NULL){
        printf("[predictLabels] Mounter histogram function not defined\n");
        return NULL;
    }
    if(bagOfVisualWordsManager->predictFunction == NULL){
        printf("[predictLabels] Predict function not defined\n");
        return NULL;
    }
    Matrix *bowHistograms = createMatrix(bagOfVisualWordsManager->pathsToImages_test->size,
                                         bagOfVisualWordsManager->dictionery->numberRows,
                                         sizeof(float));
    printf("[predictLabels] Generating histograms and labels from images\n");
    for (size_t index = 0; index < bagOfVisualWordsManager->pathsToImages_test->size; ++index) {
        char* imagePath = VECTOR_GET_ELEMENT_AS(char*,bagOfVisualWordsManager->pathsToImages_test,index);
        Image* image = readImage(imagePath);

        Image* imageISF;
        if(bagOfVisualWordsManager->pathsToImages_isf != NULL) {
            char* isfPath = VECTOR_GET_ELEMENT_AS(char*, bagOfVisualWordsManager->pathsToImages_isf, index);
            imageISF = readImage(isfPath);
        }
        
        if(image == NULL){
            printf("[predictLabels] Invalid image path: %s",imagePath);
            continue;
        }
        GVector* samplingResults = NULL;
        if(bagOfVisualWordsManager->imageSamplerFunction){
            samplingResults = bagOfVisualWordsManager->imageSamplerFunction(image, imageISF,
                                                                            bagOfVisualWordsManager);
            samplingResults->freeFunction = bagOfVisualWordsManager->freeFunction2SamplerOutput;
        }else{
            samplingResults = createNullVector(1,sizeof(Image*));
            VECTOR_GET_ELEMENT_AS(Image*,samplingResults,0) = image;
        }
        Matrix* featureMatrix = bagOfVisualWordsManager->featureExtractorFunction(samplingResults, image, bagOfVisualWordsManager);

        GVector* histogram = bagOfVisualWordsManager->mountHistogramFunction(featureMatrix,bagOfVisualWordsManager);
        setRowValueGivenVector(bowHistograms,histogram,index);

        destroyImage(&image);
        destroyMatrix(&featureMatrix);
        destroyVector(&histogram);
        destroyVector(&samplingResults);
    }
    printf("[predictLabels] Histograms and labels generated\n");
    printf("[predictLabels] Predicting labels...\n");


    GVector* labelsPredicted = bagOfVisualWordsManager->predictFunction(bowHistograms,bagOfVisualWordsManager->classifier);
    printf("[predictLabels] Labels predicted...\n");
    if(bagOfVisualWordsManager->storePredictedData){
        if(bagOfVisualWordsManager->histogramsPredictSamples){
            destroyMatrix(&(bagOfVisualWordsManager->histogramsPredictSamples));
        }
        bagOfVisualWordsManager->histogramsPredictSamples = bowHistograms;
        if(bagOfVisualWordsManager->labelsPredicted){
            destroyVector(&(bagOfVisualWordsManager->labelsPredicted));
        }
        bagOfVisualWordsManager->labelsPredicted = labelsPredicted;
    }else{
        destroyMatrix(&bowHistograms);
    }
    return labelsPredicted;
}

GVector* computeCountHistogram_bow(Matrix* featureMatrix, BagOfVisualWordsManager* bagOfVisualWordsManager){
    GVector* bowHistogram = createNullVector(bagOfVisualWordsManager->dictionery->numberRows,sizeof(float));
    for (size_t patchIndex = 0; patchIndex < featureMatrix->numberRows; ++patchIndex) {
        size_t nearestClusterIndex = findNearestRow(bagOfVisualWordsManager->dictionery,
                                                    featureMatrix,patchIndex,
                                                    bagOfVisualWordsManager->distanceFunction,
                                                    bagOfVisualWordsManager->argumentListOfDistanceFunction);
        VECTOR_GET_ELEMENT_AS(float,bowHistogram,nearestClusterIndex) += 1.0;
    }
    for (size_t i = 0; i < bowHistogram->size; ++i) {
        VECTOR_GET_ELEMENT_AS(float,bowHistogram,i) /= featureMatrix->numberRows;
    }
    return bowHistogram;
}


Matrix* kmeansClusteringBow(Matrix* featureMatrix, BagOfVisualWordsManager* bagOfVisualWordsManager){


    ArgumentList* argumentList = bagOfVisualWordsManager->argumentListOfClustering;
    if(argumentList->length <= 0){
        printf("[kmeansClustering] invalid argument list\n");
        return NULL;
    }
    if(featureMatrix->numberRows == 0 || featureMatrix->numberColumns == 0){
        printf("[kmeansClustering] invalid matrix dimension: R%lu C%lu\n",featureMatrix->numberRows,featureMatrix->numberColumns);
        return NULL;
    }

    if(argumentList->length == 1){
        size_t numberOfCluster = ARGLIST_GET_ELEMENT_AS(size_t,argumentList,0);
        return kmeansClustering(featureMatrix,numberOfCluster);
    }
    else if(argumentList->length == 2){
        size_t numberOfCluster = ARGLIST_GET_ELEMENT_AS(size_t,argumentList,0);
        size_t maximumNumberIterations = ARGLIST_GET_ELEMENT_AS(size_t,argumentList,1);
        return kmeansClustering(featureMatrix,numberOfCluster,
                                maximumNumberIterations);
    }else if(argumentList->length == 3){
        size_t numberOfCluster = ARGLIST_GET_ELEMENT_AS(size_t,argumentList,0);
        size_t maximumNumberIterations = ARGLIST_GET_ELEMENT_AS(size_t,argumentList,1);
        double tolerance = ARGLIST_GET_ELEMENT_AS(double,argumentList,2);

        return kmeansClustering(featureMatrix,numberOfCluster,
                                maximumNumberIterations,tolerance);
    }else if(argumentList->length == 4){
        size_t numberOfCluster = ARGLIST_GET_ELEMENT_AS(size_t,argumentList,0);
        size_t maximumNumberIterations = ARGLIST_GET_ELEMENT_AS(size_t,argumentList,1);
        double tolerance = ARGLIST_GET_ELEMENT_AS(double,argumentList,2);
        int seed = ARGLIST_GET_ELEMENT_AS(int,argumentList,3);
        return kmeansClustering(featureMatrix,numberOfCluster,
                                maximumNumberIterations,tolerance,seed);
    }else if (argumentList->length == 6){
        size_t numberOfCluster = ARGLIST_GET_ELEMENT_AS(size_t,argumentList,0);
        size_t maximumNumberIterations = ARGLIST_GET_ELEMENT_AS(size_t,argumentList,1);
        double tolerance = ARGLIST_GET_ELEMENT_AS(double,argumentList,2);
        int seed = ARGLIST_GET_ELEMENT_AS(int,argumentList,3);
        DistanceFunction distanceFunction = ARGLIST_GET_ELEMENT_AS(DistanceFunction ,argumentList,4);
        ArgumentList* distanceFunctionArgs = ARGLIST_GET_ELEMENT_AS(ArgumentList* ,argumentList,5);
        return kmeansClustering(featureMatrix,numberOfCluster,
                                maximumNumberIterations,tolerance,seed,
                                distanceFunction, distanceFunctionArgs);
    }else{
        printf("[kmeansClusteringBow] invalid arguments for kmeans\n");
    }
    return NULL;

}





