/*
 *Created by Deangeli Gomes Neves
 *
 * This software may be freely redistributed under the terms
 * of the MIT license.
 *
 */

#ifndef _ADJACENCYRELATION_H_
#define _ADJACENCYRELATION_H_

#include "common.h"

typedef struct _adjrel {
    /* displacements to achieve the n adjacent voxels. */
    int *dx, *dy;
    /* number of adjacent pixels/voxels. */
    int n;
    /* useful for operations between two AdjacencyRelations. */
    int maxDx;
    int maxDy;

} AdjacencyRelation;

AdjacencyRelation * createAdjacencyRelation(int n);
void  destroyAdjacencyRelation(AdjacencyRelation **adjRel);
int countNumberAdjacentsInCircle(float radius);
AdjacencyRelation *createCircularAdjacency(float radius);
AdjacencyRelation *createRectangularAdjacency(int height,int width);
AdjacencyRelation* copyAdjcencyRelation(AdjacencyRelation* adjacencyRelation);
void copyAdjcencyRelationInPlace(AdjacencyRelation* adjacencyRelationSource, AdjacencyRelation** adjacencyRelationTarget);
void getAdjacentsInCircle(float radius,AdjacencyRelation* adjRel);
void resizeCircularAdjacency(AdjacencyRelation **adjRel, float radius);
AdjacencyRelation *createLosangeAdjacency(int k);
void resizeLosangeAdjacency(AdjacencyRelation **adjRel,int k);
void printAdjacencyRelation(AdjacencyRelation* adjacencyRelation);

#endif //LIBFL_ADJACENCYRELATION_H
