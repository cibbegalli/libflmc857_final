#include "FL.h"

int main(int argc, char **argv) {



//    if (argc!=2){
//        printf("Operation types:\n");
//        printf("0: Mean Filter\n");
//        printf("1: Median Filter\n");
//        printf("2: Gaussian Filter\n");
//        printf("3: Sobel Filter\n");
//        printf("Usage: <operation type> \n");
//        //adj. radius tipically  1  to sqrt(8)
//        exit(-1);
//    }

    //gray image
    Image *image = readImage("../data/lena.pgm");

    //addUniformNoise(image,40,0.03);
    //addSaltAndPepperNoise(image, 0.03);
    //addAdditiveGaussianNoise(image,0,0.01);
    AdjacencyRelation* adjacencyRelation = createRectangularAdjacency(3,3);
    Kernel* kernel = NULL;
    Image *outputImage = NULL;
    int type = 3;
    switch (type){
        case 0:
        {
            kernel = createMeanKernel(adjacencyRelation);
            outputImage = convolveImageWithKernel(image,kernel,false,false);
            writeImage(outputImage,"outputMeanFilter.pgm");
        }
            break;
        case 1:
        {
            kernel = createGaussianKernel(adjacencyRelation,2);
            outputImage = convolveImageWithKernel(image,kernel,false,false);
            writeImage(outputImage,"outputGaussianFilter.pgm");
        }
            break;
        case 2:
        {
            outputImage = convolveImageWithMedianFilter(image,adjacencyRelation,false);
        }

            break;

        case 3:
        {

            printf("%f\n", M_PI);
            Kernel *kernelX = createSobelHorizontalGradientKernel(adjacencyRelation);
            Kernel *kernelY = createSobelVerticalGradientKernel(adjacencyRelation);
            Image *sobelX = convolveImageWithKernel(image, kernelX, false, false);
            Image *sobelY = convolveImageWithKernel(image, kernelY, false, false);
            Image* gradientMagnitude = createImage(sobelX->nx,sobelX->ny,sobelX->nchannels);
            Image* gradientphase = createImage(sobelX->nx,sobelX->ny,sobelX->nchannels);
            //compute the gradient magnitude and  phase
            for (int k = 0; k < sobelX->numberPixels; ++k) {
                double gradientMagnitudeValue = (sobelX->channel[0][k]*sobelX->channel[0][k]) +
                        (sobelY->channel[0][k]*sobelY->channel[0][k]);

                gradientMagnitudeValue = (float)sqrt(gradientMagnitudeValue);
                double angle = atan2(sobelY->channel[0][k],sobelX->channel[0][k]+0.000001);
                if(angle < 0){
                    angle = 2.0*M_PI - angle;
                }
                angle /= 2.0*M_PI;
                angle *= sobelX->scalingFactor;

                gradientMagnitude->channel[0][k] = gradientMagnitudeValue;
                gradientphase->channel[0][k] = angle;
            }
            writeImage(gradientMagnitude, "gradMagnitude.pgm");
            writeImage(gradientphase, "gradPhase.pgm");
            //outputImage = copyImage(gradientMagnitude);
            destroyImage(&gradientMagnitude);
            destroyImage(&gradientphase);
            destroyKernel(&kernelX);
            destroyKernel(&kernelY);
            destroyImage(&sobelX);
            destroyImage(&sobelY);
            
        }
            break;

        default:
            printf("Invalid Operation type:\n");
            exit(-1);
    }
//    diffImage = imageSubtraction(image,outputImage,true);
    //writeImageP2(outputImage,"output.pgm");
    writeImageP2(image,"input.pgm");
//    writeImageP2(diffImage,"diff.pgm");
    destroyImage(&image);
    //destroyImage(&outputImage);
    destroyAdjacencyRelation(&adjacencyRelation);
    destroyKernel(&kernel);
//    destroyImage(&diffImage);

    return 0;
}
