#include "FL.h"

int main(int argc, char **argv) {

    //gray image
    Image *grayImage = readImage("../data/lena.pgm");

    int index=0;
    for (int i = 0; i < grayImage->ny; ++i) {
        for (int j = 0; j < grayImage->nx; ++j) {
            index = (i*grayImage->nx) + j;
            grayImage->channel[0][index] = 255 - grayImage->channel[0][index];
        }
    }
    writeImage(grayImage,"outputGray.pgm");

    //color image
    Image *colorImage = readImage("../data/lenaColor.ppm");
    for (int i = 0; i < colorImage->ny; ++i) {
        for (int j = 0; j < colorImage->nx; ++j) {
            //int index = (i*colorImage->nx) + j;

            imageValCh(colorImage,i,j,0) = 255 - imageValCh(colorImage,i,j,0);
            imageValCh(colorImage,i,j,1) = 255 - imageValCh(colorImage,i,j,1);
            imageValCh(colorImage,i,j,2) = 255 - imageValCh(colorImage,i,j,2);
//
//          colorImage->channel[0][index] = 255 - colorImage->channel[0][index];
//            colorImage->channel[1][index] = 255 - colorImage->channel[1][index];
//            colorImage->channel[2][index] = 255 - colorImage->channel[2][index];

        }
    }
    writeImage(colorImage,"outputColor.ppm");

    destroyImage(&grayImage);
    destroyImage(&colorImage);
    return 0;
}
