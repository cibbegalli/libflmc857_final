//
// Created by deangeli on 5/20/17.
//

#ifndef _CLUSTERING_H
#define _CLUSTERING_H

#include "matrix.h"
#include "matrixUtil.h"

Matrix* kmeansClustering(Matrix* featureMatrix, size_t nClusters,
                         size_t iterMax = 1000,
                         double tolerance = 1, int seed = 0,
                         DistanceFunction distanceFunction = computeNormalizedL1Norm,
                         ArgumentList* argumentListDistanceFunction = NULL);



#endif //LIBFL_CLUSTERING_H
