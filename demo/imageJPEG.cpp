#include "FL.h"



int main(int argc, char **argv) {

    Image* image = readImageJPEG("../data/flower.jpg");
    writeImage(image,"outputRGB.jpg");
    image->colorSpace = YCbCr;
    writeImage(image,"outputYCbCr.jpg");
    destroyImage(&image);
    return 0;
}
