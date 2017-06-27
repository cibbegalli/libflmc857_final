/*
 *Created by Deangeli Gomes Neves
 *
 * This software may be freely redistributed under the terms
 * of the MIT license.
 *
 */


#ifndef _KERNEL_H_
#define _KERNEL_H_


#include "common.h"
#include "adjacencyRelation.h"
#include "featureVector.h"



typedef struct _kernel {
  AdjacencyRelation *adjacencyRelation;
  float *weight;
} Kernel;

//aloca memoria para a estrutura do kernel e cria uma relacao de adjacencia e inicializa os pesos
Kernel *createKernel(AdjacencyRelation *adjRel);
Kernel *createKernel(FeatureVector* weights);

//destroi uma estrutura de kernel
void destroyKernel(Kernel **kernel);

//define os pesos do kernel
void setKernelWeights(Kernel *kernel,float* weights, int n);
void setKernelWeights(Kernel *kernel,FeatureVector* weights);

//rotationa o kernel
void rotateKernel2D(Kernel *kernel, float angle_rad);
Kernel *copyKernel(Kernel* kernel);
void printKernel(Kernel* kernel);


//filtros
Kernel* createMeanKernel(AdjacencyRelation* adjacencyRelation);
Kernel* createGaussianKernel(AdjacencyRelation* adjacencyRelation,double variance);
Kernel* createGaussianKernel(AdjacencyRelation* adjacencyRelation,double varianceX,double varianceY);
Kernel* createSobelHorizontalGradientKernel(int height,int width);
Kernel* createSobelHorizontalGradientKernel(AdjacencyRelation* adjacencyRelation);
Kernel* createSobelVerticalGradientKernel(int height,int width);
Kernel* createSobelVerticalGradientKernel(AdjacencyRelation* adjacencyRelation);







#endif
