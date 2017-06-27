#include "FL.h"

int main(int argc, char **argv) {


    /*
     * To clear all images in demo folder, just type "make cleanImages"
     * */

    DirectoryManager* directoryManager = loadFilesFromDirBySuffix("../data/object6/", "ppm");
    Image* images = packImagesFromDirectory(directoryManager);
    //printf("%d %d")
    char number[15];
    char filename[80];
    for (int i = 0; i < images->nz; ++i) {
        Image* imageSlice = getSlice(images,i);
        sprintf(number,"%d",i);
        memset(filename,0,sizeof(filename));
        strcat(filename,"slice");
        strcat(filename,number);
        strcat(filename,".ppm");
        memset(number,0,sizeof(number));
        writeImage(imageSlice,filename);
    }

    destroyDirectoryManager(&directoryManager);
    destroyImage(&images);
//    destroyImage(&imageSlice);
//    destroyImage(&images);
    return 0;
}


