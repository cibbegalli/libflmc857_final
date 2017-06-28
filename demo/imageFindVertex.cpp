#include "FL.h"

#define DIFF_SUPER_PIXELS 3
#define SPARSE_HIST 5


int main(int argc, char **argv) {

    /*
    * To clear all images in demo folder, just type "make cleanImages"
    * */

    int patchSizeX = 64;
    int patchSizeY = 64;
    Image* image = readImagePGM("urso.pgm");
    Image* subImage = NULL;
    Image* realImage = readImage("image1.png");
    Histogram* histo;
    char number[15];
    char filename[80];

    int k = 0;
    int count=0;
    int max;
    int patch = 0;
    float aux;
    for (int y = 0; y < image->ny; y += patchSizeY) {
        for (int x = 0; x < image->nx; x += patchSizeX) {
            
            count = 0;
            aux = (image->channel[0][k]);
            //printf("image: %d value: %f\n", k, aux);
            for (int j = y; j < y + patchSizeY; j++) {
                for (int i = x; i < x + patchSizeX; i++) {
                    if ((image->channel[0][k]) != aux)
                        count ++;
                }
                k++;
            }
            
            //imprime patch
            if (count >= DIFF_SUPER_PIXELS){
                subImage = extractSubImage(realImage, x, y, patchSizeX, patchSizeY,true);
                histo = computeHistogram(subImage,128, true);
                max = 0;
                
                for(int i = 0; i<histo->n; i++){
                    if (histo->val[i] != 0)
                        max++;
                }
                
                if (max > SPARSE_HIST){ // quantidade de valores no histograma -> diferentes cores na imagem
                    sprintf(number,"%d",patch);
                    memset(filename,0,sizeof(filename));
                    strcat(filename,"patch");
                    strcat(filename,number);
                    strcat(filename,".ppm");
                    memset(number,0,sizeof(number));
                    writeImage(subImage,filename);
                    destroyImage(&subImage);
                    patch++;
                }
                destroyHistogram(&histo);
            }
            
        }
    }
    destroyImage(&image);
    return 0;
}


