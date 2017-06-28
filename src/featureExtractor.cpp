//
// Created by deangeli on 5/20/17.
//

#include "featureExtractor.h"

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

Matrix* computeHog(GVector* vector_images, Image *mag, Image *phase, int blocks_x, int blocks_y, int theta) {
    // sampling mag
    GVector* vector_mag = gridSampling(mag, 64, 64); 

    // sampling phase
    GVector* vector_phase = gridSampling(phase, 64, 64); 
    
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
        Image* patch = VECTOR_GET_ELEMENT_AS(Image*,vector_images,i);
        Image* patch_mag = VECTOR_GET_ELEMENT_AS(Image*,vector_mag,i);
        Image* patch_phase = VECTOR_GET_ELEMENT_AS(Image*,vector_phase,i);

        GVector* featureVector = computeHistogramOfOrientedGradient(patch, patch_mag, patch_phase, blocks_x, blocks_y, theta);
        for (size_t j = 0; j < matrix->numberColumns; ++j) {
            MATRIX_GET_ELEMENT_BI_AS(float,matrix,k) = VECTOR_GET_ELEMENT_AS(float,featureVector,j);
            k++;
        }
        destroyVector(&featureVector);
        destroyImage(&patch_mag);
        destroyImage(&patch_phase);
    }

    destroyVector(&vector_mag);
    destroyVector(&vector_phase);

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
    destroyImage(&img);
    (*p_mag) = img_mag;
    (*p_phase) = img_phase;
}
