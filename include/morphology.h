/*
 *Created by Deangeli Gomes Neves
 *
 * This software may be freely redistributed under the terms
 * of the MIT license.
 *
 */

#include "image.h"
#include "adjacencyRelation.h"


#ifndef _MORPHOLOGY_H_
#define _MORPHOLOGY_H_

#include "adjacencyRelation.h"
#include "featureVector.h"


Image *dilate(Image *image, AdjacencyRelation *AdjRel);
Image *erode(Image *image, AdjacencyRelation *AdjRel);
Image *close(Image *image, AdjacencyRelation *AdjRel);
Image *open(Image *image, AdjacencyRelation *AdjRel);
Image *close(Image *image, AdjacencyRelation *AdjRel);
Image *topHat(Image *image,AdjacencyRelation *AdjRel);
Image *bottomHat(Image *image,AdjacencyRelation *AdjRel);
Image *morphologicGradient(Image *image,AdjacencyRelation *AdjRel);

/*
 * nx - numero de pixels na horizontal da imagem
 * ny - numero de pixels na vertical da imagem
 * centerX - coordenada X (na imagem) onde o centro da relacao de adjacencia ira ser colocado
 * centerY - coordenada Y (na imagem) onde o centro da relacao de adjacencia ira ser colocado
 * */
Image* transformAdjacencyRelation2Image(AdjacencyRelation *adjRel,int nx,int ny,int centerX,int centerY);

/*
 * thresholding - pixels com valor abaixo de thresholding nao
 * ny - numero de pixels na vertical da imagem
 * centerX - coordenada X (na imagem) onde o centro da relacao de adjacencia ira ser colocado
 * centerY - coordenada Y (na imagem) onde o centro da relacao de adjacencia ira ser colocado
 * */
AdjacencyRelation* transformImage2AdjacencyRelation(Image *image, float thresholding,int centerX,int centerY);
/*usa o elemento estruturante de adjacencyRelation1 para dilatar o elemento estruturante em adjacencyRelation2
 * */
AdjacencyRelation* dilate(AdjacencyRelation *adjacencyRelation1, AdjacencyRelation *adjacencyRelation2);


FeatureVector *getMorphologicalPdf(Image *image, AdjacencyRelation* adjacencyRelation, int k_times);

#endif //LIBFL_MORPHOLOGY_H
