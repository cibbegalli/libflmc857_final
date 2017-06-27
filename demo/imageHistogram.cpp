#include "FL.h"

int main(int argc, char **argv) {
    char             ext[10],*pos;
    FeatureVector       *hist=NULL;

    if (argc != 4) {
        printf("histogram <image.[pgm,ppm]> <nbins> <histogram.txt>\n");
        exit(-1);
    }

    pos = strrchr(argv[1],'.') + 1;
    sscanf(pos,"%s",ext);

    if (strcmp(ext,"pgm")==0){
        Image *grayImage = readImage(argv[1]);
        hist = computeHistogramForFeatureVectorGivenNBins(grayImage,atoi(argv[2]),true);
        destroyImage(&grayImage);
    } else { /* ppm */
        Image *colorImage =  readImage(argv[1]);
        hist = computeHistogramForFeatureVectorGivenNBins(colorImage,atoi(argv[2]),true);
        destroyImage(&colorImage);
    }
    Image* image = readImage(argv[1]);
    Histogram *histogram = computeHistogram(image,43,true);
    writeHistogram(histogram,"histogram.txt");
    destroyHistogram(&histogram);


    writeHistogram(hist,argv[3]);
    destroyHistogram(&hist);

    return 0;
}

