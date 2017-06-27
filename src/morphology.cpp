//
// Created by deangeli on 3/15/17.
//
#include "morphology.h"



Image *dilate(Image *image, AdjacencyRelation *AdjRel){
    Image* dilatedImage = createImage(image->nx,image->ny,1);
    dilatedImage->scalingFactor = image->scalingFactor;

#pragma omp parallel for
    for (int p=0; p < image->numberPixels; p++) {
        dilatedImage->channel[0][p] = image->channel[0][p];
        int pixelX = (p%image->nx);
        int pixelY = (p/image->ny);
        for (int i=0; i < AdjRel->n; i++) {
            int adjacentX = pixelX + AdjRel->dx[i];
            int adjacentY = pixelY + AdjRel->dy[i];

            if(isValidPixelCoordinate(image,adjacentX,adjacentY)){
                int adjacentIndex = (adjacentY*image->nx)+adjacentX;
                if(image->channel[0][adjacentIndex] > dilatedImage->channel[0][p]){
                    dilatedImage->channel[0][p] = image->channel[0][adjacentIndex];
                }
            }
        }
    }


    return(dilatedImage);
}



Image *erode(Image *image, AdjacencyRelation *AdjRel){
    Image* erodedImage = createImage(image->nx,image->ny,1);
    erodedImage->scalingFactor = image->scalingFactor;
#pragma omp parallel for
    for (int p=0; p < image->numberPixels; p++) {
        //iftVoxel u = iftGetVoxelCoord(img,p);
        erodedImage->channel[0][p] = image->channel[0][p];
        int pixelX = (p%image->nx);
        int pixelY = (p/image->nx);
        for (int i=0; i < AdjRel->n; i++) {
            int adjacentX = pixelX + AdjRel->dx[i];
            int adjacentY = pixelY + AdjRel->dy[i];

            if(isValidPixelCoordinate(image,adjacentX,adjacentY)){
                int adjacentIndex = (adjacentY*image->nx)+adjacentX;
                if(image->channel[0][adjacentIndex] < erodedImage->channel[0][p]){
                    erodedImage->channel[0][p] = image->channel[0][adjacentIndex];
                }
            }
        }
    }
    return erodedImage;
}


Image *open(Image *image, AdjacencyRelation *AdjRel){
    Image* outputImage = NULL;
    Image* erodedImage = NULL;
    erodedImage = erode(image,AdjRel);
    outputImage = dilate(erodedImage,AdjRel);
    destroyImage(&erodedImage);
    return outputImage;
}


Image *close(Image *image, AdjacencyRelation *AdjRel){
    Image* outputImage = NULL;
    Image* dilatedImage = NULL;
    dilatedImage = dilate(image,AdjRel);
    outputImage = erode(dilatedImage,AdjRel);
    destroyImage(&dilatedImage);
    return outputImage;
}




Image *topHat(Image *image,AdjacencyRelation *AdjRel){
    Image* outputImage = NULL;
    Image* openedImage = open(image,AdjRel);
    outputImage = imageSubtraction(image,openedImage,false);
    destroyImage(&openedImage);
    return outputImage;
}


Image *bottomHat(Image *image,AdjacencyRelation *AdjRel){
    Image* outputImage = NULL;
    Image* closedImage = close(image,AdjRel);
    outputImage = imageSubtraction(closedImage,image,false);
    destroyImage(&closedImage);
    return outputImage;
}


Image *morphologicGradient(Image *image,AdjacencyRelation *AdjRel){
    Image* outputImage = NULL;
    Image* dilatedImage = dilate(image,AdjRel);
    Image* erodedImage = erode(image,AdjRel);
    outputImage = imageSubtraction(dilatedImage,erodedImage,false);
    destroyImage(&dilatedImage);
    destroyImage(&erodedImage);
    return outputImage;
}

Image* transformAdjacencyRelation2Image(AdjacencyRelation *adjRel,int nx,int ny,int centerX,int centerY){
    Image* image = createImage(nx,ny,1);
    image->scalingFactor = 255;
    int coordinateX = 0;
    int coordinateY = 0;
    int index;
    for (int i = 0; i < adjRel->n; ++i) {
        coordinateX = adjRel->dx[i] + centerX;
        coordinateY = adjRel->dy[i] + centerY;
        index = (coordinateY*nx) +  coordinateX;
        if(isValidPixelCoordinate(image,coordinateX,coordinateY)){
            image->channel[0][index] = image->scalingFactor;
        }
    }
    return image;
}

AdjacencyRelation* transformImage2AdjacencyRelation(Image *image, float thresholding,int centerX,int centerY){

    int counter = 0;
    for (int p = 0; p < image->numberPixels; ++p) {
        if(image->channel[0][p] > thresholding){
            counter++;
        }
    }
    AdjacencyRelation *adjacencyRelation = createAdjacencyRelation(counter);
    int k = 0;

    for (int p = 0; p < image->numberPixels; ++p) {

        if(image->channel[0][p] > thresholding) {
            int dx = (p % image->nx);
            int dy = (p / image->nx);
            dx -= centerX;
            dy -= centerY;


            if (abs(dx) > adjacencyRelation->maxDx) {
                adjacencyRelation->maxDx = abs(dx);
            }
            if (abs(dy) > adjacencyRelation->maxDy) {
                adjacencyRelation->maxDy = abs(dy);
            }
            int k_copy;

            k_copy = k;
            k++;

            adjacencyRelation->dx[k_copy] = dx;
            adjacencyRelation->dy[k_copy] = dy;
        }

    }
    return adjacencyRelation;
}

AdjacencyRelation* dilate(AdjacencyRelation *adjacencyRelation1, AdjacencyRelation *adjacencyRelation2){
    int nx = ((adjacencyRelation1->maxDx*2 + 1) + (adjacencyRelation2->maxDx*2 + 1)) - 1;
    int ny = ((adjacencyRelation1->maxDy*2 + 1) + (adjacencyRelation2->maxDy*2 + 1)) -1;
    int centerX = nx/2;
    int centerY = ny/2;
    Image* imageStructElement = transformAdjacencyRelation2Image(adjacencyRelation2,nx,ny,centerX,centerY);
    Image* imageDilated = dilate(imageStructElement,adjacencyRelation1);
    AdjacencyRelation* adjacencyRelationDilated = transformImage2AdjacencyRelation(imageDilated, imageDilated->scalingFactor/2,centerX,centerY);
    destroyImage(&imageStructElement);
    destroyImage(&imageDilated);
    return adjacencyRelationDilated;
}

FeatureVector *getMorphologicalPdf(Image *image, AdjacencyRelation* adjacencyRelation, int k_times){

    FeatureVector *featureVector = NULL;
    if(k_times < 0){
        return featureVector;
    }
    featureVector = createFeatureVector(k_times+1);
    featureVector->features[0] = 0;
    float v0 = sumUpAllPixelsValues(image,true);
    v0 += 0.00000001;//somo isso para evitar divisao por zero
    float vk = 0;
    float feature = 0;
    Image* openedImage = NULL;
    AdjacencyRelation* currentStructorElement = copyAdjcencyRelation(adjacencyRelation);
    AdjacencyRelation* auxStructorElement = NULL;

    for (int i = 1; i <= k_times; ++i) {
        openedImage = open(image,currentStructorElement);
        vk = sumUpAllPixelsValues(openedImage,true);
        feature = 1 - (vk/v0);
        featureVector->features[i] = feature;
        auxStructorElement = dilate(adjacencyRelation,currentStructorElement);
        destroyAdjacencyRelation(&currentStructorElement);
        currentStructorElement = auxStructorElement;
        destroyImage(&openedImage);
    }

    destroyAdjacencyRelation(&currentStructorElement);
    //destroyAdjacencyRelation(&auxStructorElement);
    //FeatureVector *vector = createFeatureVector(featureVector,k);
    return featureVector;
}

