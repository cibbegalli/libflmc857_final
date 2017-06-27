//
// Created by deangeli on 5/20/17.
//

#ifndef LIBFL_FEATUREEXTRACTOR_H
#define LIBFL_FEATUREEXTRACTOR_H

#include "vector.h"
#include "matrix.h"
#include "image.h"
#include "histogram.h"
#include "sampling.h"
#include "kernel.h"

#define RADIUS_GRADIENT 3.0

void computeGradient(Image* img, Image** p_mag, Image** p_phase);

Matrix* computeColorHistogram(GVector* vector_images,size_t nbinsPerChannel,size_t totalNumberBins);
Matrix* computeHog(GVector* vector_images, Image *mag, Image *phase, int blocks_x, int blocks_y, int theta);



#endif //LIBFL_FEATUREEXTRACTOR_H
