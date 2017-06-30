//
// Created by deangeli on 5/20/17.
//

#include "featureExtractor.h"
#include "../demo/parameters.cpp"

Matrix* computeColorHistogram(GVector* vector_images,size_t nbinsPerChannel,size_t totalNumberBins){
    Matrix* matrix = createMatrix(vector_images->size,totalNumberBins,sizeof(float));
    int k = 0;
    for (size_t i = 0; i < vector_images->size; ++i) {
        Image* patch = VECTOR_GET_ELEMENT_AS(Image*,vector_images,i);
        GVector* featureVector = computeHistogramForFeatureVectorGivenNBins(patch,nbinsPerChannel,true);
        for (size_t j = 0; j < matrix->numberColumns; ++j) {
            MATRIX_GET_ELEMENT_BI_AS(float,matrix,k) = VECTOR_GET_ELEMENT_AS(float,featureVector,j);
            k++;
        }
        destroyVector(&featureVector);
    }
    return matrix;
}

//GVector* computeHistogramOfOrientedGradient(Image *img, Image *mag, Image *phase, int bx, int by, int theta)

Matrix* computeHog(GVector* vector_images, int blocks_x, int blocks_y, int theta) {
    
    //parameters for hog
    if(360 % theta != 0) {
        printf(" [computeHog] parameter theta not divide 360\n");
        return NULL;
    }
    
    Image* img = VECTOR_GET_ELEMENT_AS(Image*,vector_images,0);
    int cx = img->nx/blocks_x;
    int cy = img->ny/blocks_y;

    if(img->nx % cx != 0 || img->ny % cy != 0) {
        printf(" [computeHog] ceil size not divide dimensions of image\n");
        return NULL;
    }

    int nBins = 360/theta;
    int totalNumberBins = nBins * blocks_x * blocks_y;
    Matrix* matrix = createMatrix(vector_images->size,totalNumberBins,sizeof(float));
    int k = 0;
    for (size_t i = 0; i < vector_images->size; ++i) {
        Image* patch_color = VECTOR_GET_ELEMENT_AS(Image*,vector_images,i);
        Image* patch = convertRGBtoYCbCr(patch_color);

        Image* patch_mag;
        Image* patch_phase;
        computeGradient(patch_color, &patch_mag, &patch_phase);

        GVector* featureVector = computeHistogramOfOrientedGradient(patch, patch_mag, patch_phase, blocks_x, blocks_y, theta);
        for (size_t j = 0; j < matrix->numberColumns; ++j) {
            MATRIX_GET_ELEMENT_BI_AS(float,matrix,k) = VECTOR_GET_ELEMENT_AS(float,featureVector,j);
            k++;
        }
        destroyVector(&featureVector);
        destroyImage(&patch_mag);
        destroyImage(&patch_phase);
        destroyImage(&patch);
    }

    return matrix;
}

Matrix* computeColorHistogramAndHog(GVector* vector_images, int blocks_x, int blocks_y, 
                                    int theta, size_t nbinsPerChannel,size_t totalNumberBins) {
    Matrix* m1 = computeHog(vector_images, blocks_x, blocks_y, theta);
    Matrix* m2 = computeColorHistogram(vector_images, nbinsPerChannel, totalNumberBins);

    int totalNumberBinsHOG = (360/theta) * blocks_x * blocks_y;
    int numberColumns = totalNumberBins + totalNumberBinsHOG;

    Matrix* matrix = createMatrix(vector_images->size,numberColumns,sizeof(float));
    int k = 0, k1=0, k2=0;
    for (size_t i = 0; i < vector_images->size; ++i) {
        for(int j=0; j<totalNumberBinsHOG; j++) {
            MATRIX_GET_ELEMENT_BI_AS(float,matrix,k) = MATRIX_GET_ELEMENT_BI_AS(float,m1,k1);
            k++;
            k1++;
        }
        for(int j=totalNumberBinsHOG; j<numberColumns; j++) {
            MATRIX_GET_ELEMENT_BI_AS(float,matrix,k) = MATRIX_GET_ELEMENT_BI_AS(float,m2,k2);
            k++;
            k2++;
        }
    }
    destroyMatrix(&m1);
    destroyMatrix(&m2);

    return matrix;
}

void computeGradient(Image* image, Image** p_mag, Image** p_phase) {

    AdjacencyRelation *adjRel = createCircularAdjacency(RADIUS_GRADIENT);

    Kernel *Kx = createKernel(adjRel);
    Kernel *Ky = createKernel(adjRel);

    double c = 2*pow((RADIUS_GRADIENT/3),2.0);
    double epsilon = 0.0000001;

    for(int i=0; i<adjRel->n; i++) {

        double dist = sqrt((double)adjRel->dx[i]*adjRel->dx[i] + adjRel->dy[i]*adjRel->dy[i])+epsilon;
    
        int dist_x = adjRel->dx[i];
        double Wx = exp( (-1) * dist * dist/c) * (dist_x/dist);

        int dist_y = adjRel->dy[i];
        double Wy = exp( (-1) * dist * dist/c) * (dist_y/dist);
        
        Kx->weight[i] = (float)Wx;
        Ky->weight[i] = (float)Wy;
    }

    Image* img = convertRGBtoYCbCr(image);
    Image* img_mag = createImage(img->nx, img->ny, 1);
    Image* img_phase = createImage(img->nx, img->ny, 1);

    for(int i=0; i<img->nx; i++) { // Percorrendo pixels da imagem
        for(int j=0; j<img->ny; j++) {

            float val = imageVal(img, i, j); //intensidade do pixel atual 
            int index = (j*img->nx) + i; //índice do pixel atual

            float Gx = 0.0, Gy = 0.0;
            
            for(int k=0; k<adjRel->n; k++) { //Percorrendo relação de adjacência

                int x = i+adjRel->dx[k]; //coordenada x do pixel adjacente
                int y = j+adjRel->dy[k]; //coordenada y do pixel adjacente

                if(isValidPixelCoordinate(img, x, y)) {

                    float val_adj = imageVal(img, x, y);
                    Gx += (val_adj - val) * Kx->weight[k];
                    Gy += (val_adj - val) * Ky->weight[k];
                }               
            }
            
            double G = sqrt((double)Gx*Gx + Gy*Gy+epsilon); //obtendo magnitude 
            img_mag->channel[0][index] = G;

            Gx = Gx/G; Gy = Gy/G;

            double angle;

            if(Gy >= 0.0){
                angle = (180/PI)*acos(Gx);
            } else {
                angle = 360 - (180/PI)*acos(Gx);
            }
            img_phase->channel[0][index] = angle;

        }

    }
    destroyKernel(&Kx);
    destroyKernel(&Ky);
    destroyAdjacencyRelation(&adjRel);
    destroyImage(&img);
    (*p_mag) = img_mag;
    (*p_phase) = img_phase;
}