#include "FL.h"



int main(int argc, char **argv) {

    Image*image = readImageTIFF("../data/MARBLES.tif");
    writeImageTIFF(image,"output.tif");
    return 0;
}
