//
// Created by deangeli on 3/14/17.
//
#include "adjacencyRelation.h"


AdjacencyRelation * createAdjacencyRelation(int n)
{
    AdjacencyRelation *AdjRel =(AdjacencyRelation *)calloc(1,sizeof(AdjacencyRelation));

    AdjRel->dx = (int *)calloc(n, sizeof(int));
    AdjRel->dy = (int *)calloc(n, sizeof(int));
    AdjRel->n  = n;

    return(AdjRel);
}



void  destroyAdjacencyRelation(AdjacencyRelation **AdjRel) {
    AdjacencyRelation *aux = *AdjRel;

    if (aux != NULL){
        if (aux->dx != NULL) free(aux->dx);
        if (aux->dy != NULL) free(aux->dy);
        free(aux);
        *AdjRel = NULL;
    }
}

AdjacencyRelation *createCircularAdjacency(float radius)
{

    AdjacencyRelation *AdjRel=NULL;
    int numberAdjacents = countNumberAdjacentsInCircle(radius);
    AdjRel = createAdjacencyRelation(numberAdjacents);
    //int indexCerter = AdjRel->n/2;
    getAdjacentsInCircle(radius, AdjRel);


    return AdjRel;
}

void resizeCircularAdjacency(AdjacencyRelation **adjRel,float radius){
    AdjacencyRelation *AdjRelOut = NULL;
    destroyAdjacencyRelation(adjRel);
    AdjRelOut = createCircularAdjacency(radius);
    *adjRel = AdjRelOut;
}

AdjacencyRelation *createLosangeAdjacency(int k)
{
    AdjacencyRelation *adjRel = NULL;
    int numberAdjacents = 0;
    if(k == 0){
        adjRel = createCircularAdjacency(0.0f);
        return adjRel;
    }

    for (int dy = -k; dy <= k; ++dy) {
        for (int dx = -k; dx <= k; ++dx) {
            if(abs(dx) + abs(dy) <= k){
                numberAdjacents++;
            }
        }
    }

    adjRel = createAdjacencyRelation(numberAdjacents);
    int index = 0;
    for (int dy = -k; dy <= k; ++dy) {
        for (int dx = -k; dx <= k; ++dx) {
            if(abs(dx) + abs(dy) <= k){
                adjRel->dx[index] = dx;
                adjRel->dy[index] = dy;
                index++;

                if(abs(dx) > adjRel->maxDx){
                    adjRel->maxDx = abs(dx);
                }
                if(abs(dy) > adjRel->maxDy){
                    adjRel->maxDy = abs(dy);
                }
            }
        }
    }

    return adjRel;
}

void resizeLosangeAdjacency(AdjacencyRelation **adjRel,int k){
    AdjacencyRelation *AdjRelOut = NULL;
    destroyAdjacencyRelation(adjRel);
    AdjRelOut = createLosangeAdjacency(k);
    *adjRel = AdjRelOut;
}

int countNumberAdjacentsInCircle(float radius){
    int radius_int = radius;
    int counter = 0;
    float radius_squared = (radius*radius)+0.1;
    for(int dy= -radius_int;dy <= radius_int; dy++){
        for(int dx= -radius_int;dx <= radius_int; dx++){
            int distance = (dx*dx) + (dy*dy);
            if( distance <= radius_squared){
                counter++;
            }
        }
    }
    return counter;
}

void getAdjacentsInCircle(float radius,AdjacencyRelation* AdjRel){
    if(AdjRel == NULL){
        return;
    }
    int radius_int = radius;

    float radius_squared = (radius*radius)+0.1;
    int i=0;
    for(int dy= -radius_int;dy <= radius_int; dy++){
        for(int dx= -radius_int;dx <= radius_int; dx++){
            int distance = (dx*dx) + (dy*dy);
            if( distance <= radius_squared){
                AdjRel->dx[i] = dx;
                AdjRel->dy[i] = dy;
                i++;
                if(abs(dx) > AdjRel->maxDx){
                    AdjRel->maxDx = abs(dx);
                }
                if(abs(dy) > AdjRel->maxDy){
                    AdjRel->maxDy = abs(dy);
                }
            }
        }
    }
}

AdjacencyRelation* copyAdjcencyRelation(AdjacencyRelation* adjacencyRelation){
    AdjacencyRelation* adjacencyRelationCopy = createAdjacencyRelation(adjacencyRelation->n);
    adjacencyRelationCopy->maxDx = adjacencyRelation->maxDx;
    adjacencyRelationCopy->maxDy = adjacencyRelation->maxDy;
    for (int i = 0; i < adjacencyRelation->n; ++i) {
        adjacencyRelationCopy->dx[i] = adjacencyRelation->dx[i];
        adjacencyRelationCopy->dy[i] = adjacencyRelation->dy[i];
    }
    return adjacencyRelationCopy;
}

void copyAdjcencyRelationInPlace(AdjacencyRelation* adjacencyRelationSource, AdjacencyRelation** adjacencyRelationTarget){
    if(*adjacencyRelationTarget != NULL){
        destroyAdjacencyRelation(adjacencyRelationTarget);
    }

    *adjacencyRelationTarget = copyAdjcencyRelation(adjacencyRelationSource);
}


AdjacencyRelation *createRectangularAdjacency(int height,int width){
    AdjacencyRelation *adjacencyRelation=NULL;
    int numberAdjacents = height*width;
    adjacencyRelation = createAdjacencyRelation(numberAdjacents);
    int centerY = height/2;
    int centerX = width/2;
    int k =0;
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            adjacencyRelation->dx[k] = x - centerX;
            adjacencyRelation->dy[k] = y - centerY;
            k++;
        }
    }
    return adjacencyRelation;
}

void printAdjacencyRelation(AdjacencyRelation* adjacencyRelation){
    for (int i = 0; i < adjacencyRelation->n; ++i) {
        printf("dx: %d  dy: %d\n",adjacencyRelation->dx[i],adjacencyRelation->dy[i]);
    }
}



