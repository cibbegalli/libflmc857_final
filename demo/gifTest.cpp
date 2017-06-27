#include "FL.h"

int main(int argc, char **argv) {


    char number[15];
    char filename[80];
    Image* image = NULL;
    Image* imageAlpha = NULL;
    uint8_t* imageVec = NULL;
    uint32_t delay = 10;
    GifWriter* gifWriter = (GifWriter*)calloc(1,sizeof(GifWriter));


    for (int i = 0; i < 360; i = i + 5) {
        sprintf(number,"%d",i);
        memset(filename,0,sizeof(filename));
        strcat(filename,"../data/object6/obj6__");
        strcat(filename,number);
        strcat(filename,".ppm");
        printf("%s\n",filename);
        image = readImage(filename);
        imageAlpha = createAlphaChannel(image,255);
        if(i == 0){
            GifBegin(gifWriter,"test.gif",image->nx,image->ny,delay);
        }
        imageVec = convertImage2IntergerArray8bits(imageAlpha);
        GifWriteFrame(gifWriter,imageVec,image->nx,image->ny,delay);
        memset(number,0,sizeof(number));
        destroyImage(&image);
        destroyImage(&imageAlpha);
        free(imageVec);
    }
    GifEnd(gifWriter);
    return 0;
}
