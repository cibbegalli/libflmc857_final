//
// Created by deangeli on 5/20/17.
//

#include "sampling.h"

GVector* gridSampling(Image* image, size_t patchSizeX,size_t patchSizeY) {
    size_t nPatchs_X = image->nx/patchSizeX;
    size_t nPatchs_Y = image->ny/patchSizeY;
    size_t nPatchs = nPatchs_X*nPatchs_Y;
    GVector* vector_images = createNullVector(nPatchs,sizeof(Image*));
    int k = 0;
    for (size_t y = 0; y <= (size_t)image->ny-patchSizeY; y +=patchSizeY) {
        for (size_t x = 0; x <= (size_t)image->nx-patchSizeX; x += patchSizeX) {
            VECTOR_GET_ELEMENT_AS(Image*,vector_images,k) = extractSubImage(image,x,y,patchSizeX,patchSizeY,true);
            k++;
        }
    }
    return vector_images;
}

GVector* samplingByISF(Image* image, float threshold, char *path) {

}

