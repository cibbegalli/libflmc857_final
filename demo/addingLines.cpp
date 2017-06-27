#include "FL.h"

int main(int argc, char **argv) {



    FeatureMatrix* featureMatrix = createFeatureMatrix();
    for (int i = 0; i < 5; ++i) {
        FeatureVector *featureVector = createRandomNormalizedFeatureVector(i+1);
        addNewLines(&featureMatrix,1);
        featureMatrix->featureVector[i] = featureVector;
        printFeatureMatrix(featureMatrix);
    }

    destroyFeatureMatrix(&featureMatrix);

    return 0;
}


