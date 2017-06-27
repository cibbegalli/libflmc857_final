#include "FL.h"

int main(int argc, char **argv) {

    /*
    * To clear all images in demo folder, just type "make cleanImages"
    * */

    int patchSizeX = 64;
    int patchSizeY = 64;
    Image* image = readImage("../data/lenaColor.ppm");
    Image* subImage = NULL;
    char number[15];
    char filename[80];

    int i = 0;
    for (int y = 0; y < image->ny; y += patchSizeY) {
        for (int x = 0; x < image->nx; x += patchSizeX) {
            subImage = extractSubImage(image,x,y,patchSizeX,patchSizeY,true);
            sprintf(number,"%d",i);
            memset(filename,0,sizeof(filename));
            strcat(filename,"patch");
            strcat(filename,number);
            strcat(filename,".ppm");
            memset(number,0,sizeof(number));
            writeImage(subImage,filename);
            destroyImage(&subImage);
            i++;
        }
    }
    destroyImage(&image);
    return 0;
}


