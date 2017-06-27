//
// Created by deangeli on 3/22/17.
//
#include "filtering.h"

float* circularCrossCorrelation1D(float *X1,float *X2,int N1,int N2,bool normalization){
    float* signal = NULL;
    float* filter = NULL;
    float* output = NULL;
    int signalSize;
    int filterSize;

    if(N1 > N2){
        signal = X1;
        filter = X2;
        signalSize = N1;
        filterSize = N2;
    }else{
        filter = X1;
        signal = X2;
        signalSize = N2;
        filterSize = N1;
    }
    output = (float*)calloc(signalSize,sizeof(float));
    int filterCenter = signalSize/2;
    double sum;
    double signalSquared;
    double filterSquared;
    for (int i = 0; i < signalSize; ++i) {
        sum = 0;
        signalSquared = 0;
        filterSquared = 0;
        for (int j = 0; j < filterSize; ++j) {
            int index = i-filterCenter+j;
            if(index < 0){
                index = signalSize+index;
            }
            if(index >= signalSize){
                index %= signalSize;
            }
            sum += signal[index]*filter[j];
            if (normalization){
                signalSquared += signal[index]*signal[index];
                filterSquared += filter[j]*filter[j];
            }
        }
        if (normalization){
            double denominator = sqrt(signalSquared)*sqrt(filterSquared) + 0.000001;
            sum /= denominator;
        }
        output[i] = (float)sum;

    }
    return output;
}

FeatureVector* circularCrossCorrelation1D(FeatureVector*X1,FeatureVector*X2,bool normalization){
    FeatureVector* signal = NULL;
    FeatureVector* filter = NULL;
    FeatureVector* output = NULL;

    if(X1->size > X2->size){
        signal = X1;
        filter = X2;
    }else{
        filter = X1;
        signal = X2;
    }
    output = createFeatureVector(signal->size);

    int filterCenter = filter->size/2;
    double sum;
    double signalSquared;
    double filterSquared;
    for (int i = 0; i < signal->size; ++i) {
        sum = 0;
        signalSquared = 0;
        filterSquared = 0;
        for (int j = 0; j < filter->size; ++j) {
            int index = i-filterCenter+j;
            if(index < 0){
                index = signal->size+index;
            }
            if(index >= signal->size){
                index %= signal->size;
            }
            sum += signal->features[index]*filter->features[j];
            if (normalization){
                signalSquared += signal->features[index]*signal->features[index];
                filterSquared += filter->features[j]*filter->features[j];
            }
        }
        if (normalization){
            double denominator = sqrt(signalSquared)*sqrt(filterSquared) + 0.000001;
            sum /= denominator;
        }
        output->features[i] = sum;
    }
    return output;
}

FeatureVector* circularCrossCorrelation1D(FeatureVector* featureVector, Kernel* kernel ,bool normalization){
    FeatureVector* output = NULL;
    output = createFeatureVector(featureVector->size);

    int filterCenter = kernel->adjacencyRelation->n/2;
    float sum;
    float signalSquared;
    float filterSquared;
    for (int i = 0; i < featureVector->size; ++i) {
        sum = 0;
        signalSquared = 0;
        filterSquared = 0;
        for (int j = 0; j < kernel->adjacencyRelation->n; ++j) {
            int index = i-filterCenter+j;
            if(index < 0){
                index *= -1;
                index = index%featureVector->size;
                index = featureVector->size-index;
            }
            if(index >= featureVector->size){
                index %= featureVector->size;
            }
            sum += featureVector->features[index]*kernel->weight[j];
            if (normalization){
                signalSquared += featureVector->features[index]*featureVector->features[index];
                filterSquared += kernel->weight[j]*kernel->weight[j];
            }
        }
        if (normalization){
            float denominator = (sqrtf(signalSquared)*sqrtf(filterSquared)) + 0.000001;
            sum /= denominator;
        }
        output->features[i] = (float)sum;
    }
    return output;
}



FeatureVector* Convolution1D(FeatureVector*X1,FeatureVector*X2){
    int sizeOutputArray = 0;
    int windowSize = 0;
    float *signal = NULL;
    float *filter = NULL;
    FeatureVector *outputSignal = NULL;

    if(X1 == NULL || X2 == NULL){
        printf("X1 or/and X2 are NULL\n");
    }
    if (X1->size > X2->size){
        windowSize = X2->size;
        signal = X1->features;
        filter = X2->features;
    }else{
        windowSize = X1->size;
        signal = X2->features;
        filter = X1->features;
    }
    sizeOutputArray = X1->size + X2->size - 1;
    outputSignal = createFeatureVector(sizeOutputArray);
    float sum;
    int index;
    int windowPosition;
    int windowLastPositionIndex = windowSize-1;
    for (int i = 0; i < sizeOutputArray; ++i) {
        windowPosition = i-windowLastPositionIndex;
        sum = 0;
        for (int j = windowLastPositionIndex; j >= 0 ; j--) {
            index = j+windowPosition;
            if(index < 0 || index > X1->size){
                continue;
            }
            sum += signal[index]*filter[windowLastPositionIndex-j];
        }
        outputSignal->features[i] = sum;
    }
    return outputSignal;
}


FeatureVector* matchingWithCorrelation(FeatureVector*X1,FeatureVector*X2, int type){
    FeatureVector* signal = NULL;
    FeatureVector* filter = NULL;
    FeatureVector* output = NULL;

    if(X1->size > X2->size){
        signal = X1;
        filter = X2;
    }else{
        filter = X1;
        signal = X2;
    }
    output = createFeatureVector(signal->size);

    int filterCenter = filter->size/2;
    double sum;
    for (int i = 0; i < signal->size; ++i) {
        sum = 0;
        for (int j = 0; j < filter->size; ++j) {
            int index = i-filterCenter+j;
            if(index < 0){
                index = signal->size+index;
            }
            if(index >= signal->size){
                index %= signal->size;
            }
            if(type == 1){
                sum += fabs(signal->features[index]-filter->features[j]);
            }else if(type == 2){
                sum += (signal->features[index]-filter->features[j])*(signal->features[index]-filter->features[j]);
            }
        }
        output->features[i] = sum;
    }
    return output;
}

//filtragem
Image* convolveImageWithKernel(Image* image,Kernel* kernel,bool MultiBandFiltering, bool saturateValue){
    Image* outputImage = NULL;
    if(MultiBandFiltering){
        outputImage = convolveMultiBandImageWithKernel(image,kernel,saturateValue);
    }else{
        outputImage = convolveMonoBandImageWithKernel(image,kernel,saturateValue);
    }
    return outputImage;
}

Image* convolveMultiBandImageWithKernel(Image* image,Kernel* kernel,bool saturateValue){
    Image* outputImage = createImage(image->nx,image->ny,image->nchannels);
    outputImage->scalingFactor = image->scalingFactor;
#pragma omp parallel for
    for (int y = 0; y < image->ny; ++y) {
        for (int x = 0; x < image->nx; ++x) {
            for (int k = 0; k < kernel->adjacencyRelation->n; ++k) {
                int coordenateX = x + kernel->adjacencyRelation->dx[k];
                int coordenateY = y + kernel->adjacencyRelation->dy[k];
                if(isValidPixelCoordinate(image,coordenateX,coordenateY)){
                    for (int c = 0; c < image->nchannels; ++c) {
                        imageValCh(outputImage,x,y,c) += imageValCh(image,coordenateX,coordenateY,c)*kernel->weight[k];
                        if(saturateValue){
                            if(imageValCh(outputImage,x,y,c)>outputImage->scalingFactor){
                                imageValCh(outputImage,x,y,c) = outputImage->scalingFactor;
                            }
                            if(imageValCh(outputImage,x,y,c)<0){
                                imageValCh(outputImage,x,y,c) = 0;
                            }
                        }
                    }
                }
            }
        }
    }
    return outputImage;
}

Image* convolveMonoBandImageWithKernel(Image* image,Kernel* kernel,bool saturateValue){
    Image* outputImage = createImage(image->nx,image->ny,image->nchannels);
    outputImage->scalingFactor = image->scalingFactor;
#pragma omp parallel for
    for (int y = 0; y < image->ny; ++y) {
        for (int x = 0; x < image->nx; ++x) {
            for (int k = 0; k < kernel->adjacencyRelation->n; ++k) {
                int coordenateX = x + kernel->adjacencyRelation->dx[k];
                int coordenateY = y + kernel->adjacencyRelation->dy[k];
                if(isValidPixelCoordinate(image,coordenateX,coordenateY)){
                    imageVal(outputImage,x,y) += imageVal(image,coordenateX,coordenateY)*kernel->weight[k];
                    if(saturateValue){
                        if(imageVal(outputImage,x,y)>outputImage->scalingFactor){
                            imageVal(outputImage,x,y) = outputImage->scalingFactor;
                        }
                        if(imageVal(outputImage,x,y)<0){
                            imageVal(outputImage,x,y) = 0;
                        }
                    }

                }
            }
        }
    }
    return outputImage;
}


int comp (const void * elem1, const void * elem2)
{
    const float f = *((const float*)elem1);
    const float s = *((const float*)elem2);
    if (f > s) return  1;
    if (f < s) return -1;
    return 0;
}

Image* convolveMultiBandImageWithMedianFilter(Image* image,AdjacencyRelation* adjacencyRelation){
    Image* outputImage = createImage(image->nx,image->ny,image->nchannels);
    outputImage->scalingFactor= image->scalingFactor;
    FeatureMatrix *featureMatrix = createFeatureMatrix(image->nchannels,adjacencyRelation->n);
    bool isEvenSize = (adjacencyRelation->n%2 == 0);
    float median=0;
    int vecHalfSize = adjacencyRelation->n/2;
#pragma omp parallel for
    for (int y = 0; y < image->ny; ++y) {
        for (int x = 0; x < image->nx; ++x) {
            for (int k = 0; k < adjacencyRelation->n; ++k) {
                int coordenateX = x + adjacencyRelation->dx[k];
                int coordenateY = y + adjacencyRelation->dy[k];
                if(isValidPixelCoordinate(image,coordenateX,coordenateY)){
                    for (int c = 0; c < image->nchannels; ++c) {
                        featureMatrix->featureVector[c]->features[k] = imageValCh(image,coordenateX,coordenateY,c);
                    }
                }else{
                    for (int c = 0; c < image->nchannels; ++c) {
                        featureMatrix->featureVector[c]->features[k] = 0;
                    }
                }
            }

            for (int c = 0; c < image->nchannels; ++c) {
                qsort(featureMatrix->featureVector[c]->features,featureMatrix->featureVector[c]->size, sizeof(float ), comp);
                median = (isEvenSize)
                         ? (featureMatrix->featureVector[c]->features[vecHalfSize] + featureMatrix->featureVector[c]->features[vecHalfSize-1])*0.5f
                         : featureMatrix->featureVector[c]->features[vecHalfSize]
                        ;
                imageValCh(outputImage,x,y,c) = median;
            }

        }
    }

    destroyFeatureMatrix(&featureMatrix);
    return outputImage;

}

Image* convolveMonoBandImageWithMedianFilter(Image* image,AdjacencyRelation* adjacencyRelation){
    Image* outputImage = createImage(image->nx,image->ny,image->nchannels);
    outputImage->scalingFactor= image->scalingFactor;
    FeatureVector *featureVector = createFeatureVector(adjacencyRelation->n);
    bool isEvenSize = (adjacencyRelation->n%2 == 0);
    float median=0;
    int vecHalfSize = adjacencyRelation->n/2;
#pragma omp parallel for
    for (int y = 0; y < image->ny; ++y) {
        for (int x = 0; x < image->nx; ++x) {
            for (int k = 0; k < adjacencyRelation->n; ++k) {
                int coordenateX = x + adjacencyRelation->dx[k];
                int coordenateY = y + adjacencyRelation->dy[k];
                if(isValidPixelCoordinate(image,coordenateX,coordenateY)){
                    featureVector->features[k] = imageVal(image,coordenateX,coordenateY);
                }else{
                    featureVector->features[k] = 0;
                }
            }

            qsort(featureVector->features,featureVector->size, sizeof(float ), comp);
            if(isEvenSize){
                median = (featureVector->features[vecHalfSize] + featureVector->features[vecHalfSize-1])*0.5f;
            }else{
                median = featureVector->features[vecHalfSize];
            }

            imageVal(outputImage,x,y) = median;
        }
    }
    destroyFeatureVector(&featureVector);
    return outputImage;
}


Image* convolveImageWithMedianFilter(Image* image,AdjacencyRelation* adjacencyRelation,bool filterMultiBand){

    Image* outputImage = NULL;

    if(filterMultiBand){
        outputImage = convolveMultiBandImageWithMedianFilter(image,adjacencyRelation);
    }else{
        outputImage =convolveMonoBandImageWithMedianFilter(image,adjacencyRelation);
    }

    return outputImage;
}




