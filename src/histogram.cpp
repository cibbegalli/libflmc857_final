#include "histogram.h"

void writeHistogram(Histogram *hist, char *filename) {
    FILE *fp = fopen(filename, "w");
    int i;

    for (i = 0; i < hist->n; i++) {
        fprintf(fp, "%d %f\n", i, hist->val[i]);
    }

    fclose(fp);
}

void destroyHistogram(Histogram **hist)
{
    if (*hist != NULL) {
        free((*hist)->val);
        free(*hist);
    }
    *hist = NULL;
}

FeatureVector* createFeatureVector(Histogram *histogram){
    FeatureVector* featureVector = (FeatureVector*)calloc(1,sizeof(FeatureVector));
    featureVector->size = histogram->n;
    featureVector->features = (float *)calloc(histogram->n,sizeof(float));
    for (int i = 0; i < histogram->n; ++i) {
        featureVector->features[i] = histogram->val[i];
    }
    return featureVector;
}

Histogram* computeHistogram(Image *image,float binSize, bool normalization){

    int numberBinsPerChannel = ceil((image->scalingFactor+1)/binSize) ;
    int totalNumberBins = pow(numberBinsPerChannel,image->nchannels);
    Histogram *histogram = createHistogram(totalNumberBins);
    histogram->binSize = binSize;
    for (int yp=0; yp < image->ny; yp++){
        for (int xp=0; xp < image->nx; xp++) {
            int index = (yp*image->nx) + xp;
            int binIndex = 0;
            int w = 1;
            for (int cp=0; cp < image->nchannels; cp++) {
                int binIndexOnChannel = image->channel[cp][index]/binSize;
                binIndex += (binIndexOnChannel)*w;
                w *= numberBinsPerChannel;
            }
            histogram->val[binIndex] += 1.0f;
        }
    }
    if(normalization){
        int numberPixels = image->nx*image->ny;
        for (int i = 0; i < histogram->n; ++i) {
            histogram->val[i] /= numberPixels;
        }
    }
    return histogram;
}

FeatureVector* computeHistogramForFeatureVector(Image *image,float binSize, bool normalization){

    int numberBinsPerChannel = ceil((image->scalingFactor+1)/binSize) ;
    int totalNumberBins = pow(numberBinsPerChannel,image->nchannels);
    FeatureVector *histogram = createFeatureVector(totalNumberBins);
    for (int yp=0; yp < image->ny; yp++){
        for (int xp=0; xp < image->nx; xp++) {
            int index = (yp*image->nx) + xp;
            int binIndex = 0;
            int w = 1;
            for (int cp=0; cp < image->nchannels; cp++) {
                int binIndexOnChannel = image->channel[cp][index]/binSize;
                binIndex += (binIndexOnChannel)*w;
                w *= numberBinsPerChannel;
            }
            histogram->features[binIndex] += 1.0f;
        }
    }

    if(normalization){
        int numberPixels = image->nx*image->ny;
        for (int i = 0; i < histogram->size; ++i) {
            histogram->features[i] /= numberPixels;
        }
    }
    return histogram;
}

GVector* computeHistogramForFeatureVectorGivenNBins(Image *image,int numberBinsPerChannel, bool normalization) {

    float binSize = ( (image->scalingFactor+1) / numberBinsPerChannel);

    int totalNumberBins = pow(numberBinsPerChannel, image->nchannels);
    GVector *histVec = createNullVector(totalNumberBins,sizeof(float));

    for (int yp = 0; yp < image->ny; yp++) {
        for (int xp = 0; xp < image->nx; xp++) {
            int index = (yp * image->nx) + xp;
            int binIndex = 0;
            int w = 1;
            for (int cp = 0; cp < image->nchannels; cp++) {
                int binIndexOnChannel = (image->channel[cp][index]) / (binSize+1);
                binIndex += (binIndexOnChannel) * w;
                w *= numberBinsPerChannel;
            }
            VECTOR_GET_ELEMENT_AS(float,histVec,binIndex) += 1.0f;
        }
    }
    if (normalization) {
        int numberPixels = image->nx * image->ny;
        for (size_t i = 0; i < histVec->size; ++i) {
            VECTOR_GET_ELEMENT_AS(float,histVec,i) /= numberPixels;
        }
    }

    return histVec;
}


GVector* computeHistogramOfOrientedGradient(Image *img, Image *mag, Image *phase, int bx, int by, int theta) {
    
    int nBins = 360/theta;

    int nCeils = bx*by;
    int cx = img->nx/bx;
    int cy = img->ny/by;

    FeatureMatrix *histograms = createFeatureMatrix(nCeils, nBins);
    for(int i=0; i<nCeils; i++){
        setValueInFeatureVector(histograms->featureVector[i], 0.0);
    }

    FeatureVector *sum = createFeatureVector(nCeils);
    setValueInFeatureVector(sum, 0.0);


    for(int i=0; i<img->nx; i++) {
        for(int j=0; j<img->ny; j++) {

            // Valores magnitude e fase
            float w = imageVal(mag, i, j);
            float z = imageVal(phase, i, j);

            int b1 = ((int)(z/theta)) % nBins;
            int b2 = ((int)(z/theta + 1)) % nBins;
            float z1 = (b1+1)*theta - theta/2.0;
            float z2 = b2*theta + theta/2.0;

            int column = i/cx, line = j/cy;
            int ceil = column + line*(img->nx/cx);

            float w1 = w * (z-z1)/theta;
            float w2 = w * (z2-z)/theta;

            if(w1 < 0) w1 *= -1.0;
            if(w2 < 0) w2 *= -1.0;
     
            histograms->featureVector[ceil]->features[b1]+=w1;
            histograms->featureVector[ceil]->features[b2]+=w2;

            sum->features[ceil] += (w1*w1+w2*w2);
        }
    }

    GVector *histVec = createNullVector(nCeils*nBins,sizeof(float));
    float epsilon = 0.0000001;
    int k=0;
    for(int i=0; i<nCeils; i++) {
        float norm = sqrt(sum->features[i]) + epsilon;
        for(int j=0; j<nBins; j++) {
            float val = histograms->featureVector[i]->features[j]/norm;
            VECTOR_GET_ELEMENT_AS(float,histVec,k) = val;
            k++;
        }
    }
    destroyFeatureVector(&sum);
    destroyFeatureMatrix(&histograms);

    return histVec;
}

Image *ProbabilityDensityFunction(Image *image, double standardDeviation)
{
    Image *pdf = createImage(image->nx, image->ny,1);
    pdf->scalingFactor = 255; //8-bits
    float  *val = (float *)calloc(image->nx*image->ny,sizeof(float));
    double  K = 2.0*standardDeviation*standardDeviation, maxdist = 3*standardDeviation*3*standardDeviation;
    float   Imax = image->scalingFactor, maxval, minval;
    maxval = INT_MIN; minval = INT_MAX;
#pragma omp parallel for
    for (int p=0; p < image->numberPixels; p++) {
        float *vec_p = (float*)calloc(image->nchannels,sizeof(float));
        float *vec_q = (float*)calloc(image->nchannels,sizeof(float));
        for (int i = 0; i < image->nchannels; ++i) {
            vec_p[i] = image->channel[i][p]/Imax;
        }

        for (int q=0; q < image->numberPixels; q++) {
            if (p != q) {
                double dist = 0;
                for (int i = 0; i < image->nchannels; ++i) {
                    vec_q[i] = image->channel[i][q]/Imax;
                    dist += (vec_p[i] - vec_q[i])*(vec_p[i] - vec_q[i]);
                }
                dist = sqrt(dist);
                if (dist <= maxdist){
                    val[p] += exp(-dist/K);
                }
            }
        }
        free(vec_p);
        free(vec_q);

        if (val[p] > maxval)
            maxval = val[p];
        if (val[p] < minval)
            minval = val[p];
    }

    if (minval != maxval)
        for (int p=0; p < image->numberPixels; p++) {
            pdf->channel[0][p] = (int)(255*(val[p]-minval)/(maxval-minval));
        }
    free(val);

    return(pdf);

}

Histogram* createHistogram(int n){
    Histogram *histogram = (Histogram *) calloc(1,sizeof(Histogram));
    histogram->n = n;
    histogram->val = (float *)calloc(histogram->n,sizeof(float));
    histogram->binSize = 1;
    return histogram;
}


/*
GVector* computeHOG(Image *image) {

    //Validação tamanho da janela/imagem/patch
    if(image->nx % HOG_N2 != 0 || image>ny % HOG_M2 != 0) {
        printf("Imagem não alinhada número de células\n");
        return NULL;
    }

    return hog(image);
}
*/