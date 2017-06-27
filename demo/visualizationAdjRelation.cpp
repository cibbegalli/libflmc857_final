#include "FL.h"

int main(int argc, char **argv) {
    int  k = 30;
    AdjacencyRelation *adjacencyRelationLosango = createLosangeAdjacency(k);
    int paddingFactor = 2;
    int losangoImageSize = (((k*2) +1)*paddingFactor)-1;
    int centerX = losangoImageSize/2;
    int centerY = losangoImageSize/2;
    Image* imageLosango = transformAdjacencyRelation2Image(adjacencyRelationLosango,losangoImageSize,losangoImageSize,centerX,centerY);
    writeImage(imageLosango,"adjRelLosango.pgm");
    destroyAdjacencyRelation(&adjacencyRelationLosango);
    destroyImage(&imageLosango);

    float radius = 30;
    AdjacencyRelation *adjacencyRelationCircle = createCircularAdjacency(radius);
    int circleImageSize = (((radius*2) +1)*paddingFactor)-1;
    centerX = circleImageSize/2;
    centerY = circleImageSize/2;
    Image* imageCircle = transformAdjacencyRelation2Image(adjacencyRelationCircle,circleImageSize,circleImageSize,centerX,centerY);
    writeImage(imageCircle,"adjRelCircle.pgm");
    destroyAdjacencyRelation(&adjacencyRelationCircle);
    destroyImage(&imageCircle);

    return 0;
}
