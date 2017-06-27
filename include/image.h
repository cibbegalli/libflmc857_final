/*
 *Created by Deangeli Gomes Neves
 *
 * This software may be freely redistributed under the terms
 * of the MIT license.
 *
 */


#ifndef _IMAGE_H_
#define _IMAGE_H_

#include "common.h"
#include "file.h"
#include "png.h"
#include "jpeglib.h"
//#include  "tiffio.h"




enum ColorSpace {
    UNKNOWN = 0,
    GRAYSCALE = 1,
    GRAYSCALE_ALPHA = 2,
    RGB = 3,
    RGBA = 4,
    YCbCr = 5,
    YCbCrK = 6,
    HSV = 7,
    HSVA = 8,
    CMYK = 9
};

enum DataType {
    FLOAT,
    DOUBLE
};


typedef struct _image {
    int   nx,ny,nz;    /* dimensoes da imagem */
    float dx,dy, dz;    /* tamanho do pixel em unidades de comprimento */
    unsigned  char channelDepth;
    int numberPixels;
    float **channel; /*might be RGB, or Ycbcr*/
    int scalingFactor;
    char  unid[10]; /* unidade de comprimento */
    int nchannels;
    ColorSpace colorSpace;
    DataType dataTypeId;
} Image;


//macro para facil acesso aos pixels da imagem. a macro computa o indice para acessar a posicao (x,y) da imagem
#define imageVal(image, x, y) image->channel[0][(y*image->nx) + x]
#define imageValCh(image, x, y, c) image->channel[c][(y*image->nx) + x]
#define imageVolume(image, x, y,z) image->channel[0][(z*image->ny*image->nx) + (y*image->nx) + x]
#define imageVolumeCh(image, x, y, z,c) image->channel[c][(z*image->ny*image->nx) + (y*image->nx) + x]


Image* createImage(int nx, int ny,int nchannels);
Image* createImage(int nx, int ny,int nz, int nchannels);
Image* createImage(int nx, int ny);
void findAppropriateColorSpace(Image* image);
void destroyImage(Image**image );
void destroyImageVoidPointer(void* data);
Image *readImagePGM(char *filename);
Image *readImagePPM(char *filename);
Image *readImageTIFF(char *filename);
Image *readImagePNG(char *filename);
Image *readImageJPEG(char *filename);
Image* readImage(char *filename);
void writeImageP2(Image *image,char *filename);
void writeImageP3(Image *image, char *filename);
void writeImageP5(Image *image, char *filename);
void writeImageP6(Image *image, char *filename);
void writeImagePGM(Image* image, char *filename,  char* magicNumber);
void writeImagePPM(Image* image, char *filename,  char* magicNumber);
void writeImagePNG(Image *image,char *filename);
void writeImageJPEG(Image *image,char *filename);
void writeImageTIFF(Image *image,char *filename);
void writeImage(Image* image, char *filename);
bool isValidPixelCoordinate(Image *image,int pixelCoordinateX,int pixelCoordinateY);
Image *imageSubtraction(Image *image1, Image *image2, bool saturation);
bool isImagesSameDomain(Image *image1,Image *image2);
Image *convertRGBtoYCbCr(Image *rgbImage);
float sumUpAllPixelsValues(Image *image, bool normalize);
Image* extractImageChannel(Image* image, int channel);
void printImage(Image* image);
Image* createAlphaChannel(Image* image,float alpha);
Image* convertGrayImage2RGBImage(Image* image_ppm);
Image* mergeImages(Image* image1, Image* image2);
Image* getSlice(Image* image, int sliceIndex);
void putSlice(Image* volume, Image* image, int sliceIndex);
Image* packImagesFromDirectory(DirectoryManager* directoryManager);
Image* extractSubImage(Image*image, int xCoord,int yCoord, int xsize, int ysize, bool zeroPadding);



uint8_t* convertImage2IntergerArray8bits(Image* imageRGBA);
Image* copyImage(Image* image, bool copyChannels);
void addUniformNoise(Image* image, float uniformValue, double probability);
void addSaltAndPepperNoise(Image* image, double probability);
void addAdditiveGaussianNoise(Image* image, double mean, double variance);

#endif
