#include "FL.h"



int main(int argc, char **argv) {

    Image* image = readImage("../data/airplane.png");
    for (int y = 0; y < image->ny; ++y) {
        for (int x = 0; x < image->nx; ++x) {
            printf("Pixel at position [x:%d - y:%d] has values: ",x,y);
            for (int c = 0; c < image->nchannels; ++c) {
                printf("%d ", (int)imageValCh(image,x,y,c));
                imageValCh(image,x,y,c) = imageValCh(image,x,y,c)*(c/4.0f);
            }
            printf("\n");
        }
    }
    printf("nx:%d ny:%d chanelDepth:%d nChannels:%d colorSpaceId:%d \n",
           image->nx,image->ny,image->channelDepth,image->nchannels,image->colorSpace);
    writeImage(image,"airplaneModified.png");

    return 0;
}
