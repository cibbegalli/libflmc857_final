//
// Created by deangeli on 5/20/17.
//

#ifndef _SAMPLING_H
#define _SAMPLING_H

#include "vector.h"
#include "image.h"
#include "featureVector.h"
#include "adjacencyRelation.h"
#include "featureExtractor.h"

GVector* gridSampling(Image* image, size_t patchSizeX,size_t patchSizeY);
GVector* samplingBySuperpixelAndGradient(Image* image, Image* imageISF, float threshold, int patchSizeX, int patchSizeY, int width_adj_rel, int height_adj_rel);
GVector* samplingBySuperPixel(Image* realImage, Image* superPixel, size_t patchSizeX,size_t patchSizeY, size_t DIFF_SUPER_PIXELS, size_t binSize, size_t SPARSE_HIST);

#endif //LIBFL_SAMPLING_H
