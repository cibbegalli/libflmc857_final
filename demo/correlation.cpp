#include "FL.h"

int main(int argc, char **argv) {

    float v1[16] = {3,2,4,1,3,8,4,0,3,8,0,7,7,7,1,2};
    float v2[3] = {3,7,5};
//    float v1[16] = {10,50,60,10,20,40,30};
//    float v2[3] = {1./3.,1./3.,1./3.};

    FeatureVector* featureVector1 = createFeatureVector(v1,16);
    FeatureVector* featureVector2 = createFeatureVector(v2,3);
    Kernel *kernel = createKernel(featureVector2);

    FeatureVector* output1 = circularCrossCorrelation1D(featureVector1,kernel,true);
    FeatureVector* output2 = circularCrossCorrelation1D(featureVector1,featureVector2,true);
    FeatureVector* output3 = matchingWithCorrelation(featureVector1,featureVector2,1);
    printFeatureVector(output1);
    printFeatureVector(output2);
    printFeatureVector(output3);
    destroyFeatureVector(&featureVector1);
    destroyFeatureVector(&featureVector2);
    destroyFeatureVector(&output1);
    destroyFeatureVector(&output2);
    destroyFeatureVector(&output3);
    destroyKernel(&kernel);
    return 0;

}
