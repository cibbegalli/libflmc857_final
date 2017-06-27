#include "kernel.h"

Kernel *createKernel(AdjacencyRelation *adjRel){
    Kernel* kernel = NULL;
    size_t size = adjRel->n;
    kernel  = (Kernel*) calloc(1,sizeof(Kernel));
    kernel->adjacencyRelation = copyAdjcencyRelation(adjRel);
    kernel->weight = (float*)calloc( size , sizeof(float) );
    return kernel;
}

Kernel *createKernel(FeatureVector* weights){
    Kernel* kernel = NULL;
    size_t size = (size_t)weights->size;
    kernel  = (Kernel*) calloc(1,sizeof(Kernel));
    kernel->weight = (float*)calloc(size,sizeof(float));
    AdjacencyRelation* adjacencyRelation = createAdjacencyRelation(weights->size);
    int centerX = weights->size/2;
    adjacencyRelation->maxDx = centerX;
    adjacencyRelation->maxDy = 0;
    for (int i = 0; i < weights->size; ++i) {
        adjacencyRelation->dx[i] = i - centerX;
        adjacencyRelation->dy[i] = 0;
        kernel->weight[i] = weights->features[i];
    }
    kernel->adjacencyRelation = adjacencyRelation;
    return kernel;
}


void destroyKernel(Kernel **kernel){
    if((*kernel) == NULL){
        return;
    }
    destroyAdjacencyRelation(&((*kernel)->adjacencyRelation));
    free((*kernel)->weight);
    free(*kernel);
    *kernel = NULL;
}

void setKernelWeights(Kernel* kernel, float* weights, int n){
    if(n != kernel->adjacencyRelation->n){
        printf("dimensions mismatch\n");
    }
    for (int i = 0; i < kernel->adjacencyRelation->n; ++i) {
        kernel->weight[i] = weights[i];
    }
}

void setKernelWeights(Kernel *kernel,FeatureVector* weights){
    if(weights->size != kernel->adjacencyRelation->n){
        printf("dimensions mismatch\n");
    }
    for (int i = 0; i < kernel->adjacencyRelation->n; ++i) {
        kernel->weight[i] = weights->features[i];
    }
}

void rotateKernel2D(Kernel *kernel, float angle_rad){
    for (int i = 0; i < kernel->adjacencyRelation->n; ++i) {
        kernel->adjacencyRelation->dx[i] = round((cos(angle_rad)*kernel->adjacencyRelation->dx[i]) +
                (-sin(angle_rad)* kernel->adjacencyRelation->dy[i]));

        kernel->adjacencyRelation->dy[i] = round((sin(angle_rad)*kernel->adjacencyRelation->dx[i]) +
                                           (cos(angle_rad)* kernel->adjacencyRelation->dy[i]));

        if(abs(kernel->adjacencyRelation->dx[i]) >  kernel->adjacencyRelation->maxDx){
            kernel->adjacencyRelation->maxDx = abs(kernel->adjacencyRelation->dx[i]);
        }
        if(abs(kernel->adjacencyRelation->dy[i]) >  kernel->adjacencyRelation->maxDy){
            kernel->adjacencyRelation->maxDy = abs(kernel->adjacencyRelation->dy[i]);
        }
    }
}

void rotateKernel2D_degree(Kernel *kernel, float angle_degree){
    float angle_rad = (angle_degree * M_PI)/180.0;
    rotateKernel2D(kernel,angle_rad);
}

Kernel *copyKernel(Kernel* kernel){
    Kernel* kernelCopy = createKernel(kernel->adjacencyRelation);
    setKernelWeights(kernelCopy,kernel->weight,kernel->adjacencyRelation->n);
    return kernelCopy;
}

//filtros
Kernel* createMeanKernel(AdjacencyRelation* adjacencyRelation){
    Kernel *kernel = createKernel(adjacencyRelation);
    float factor = 1.0f/(adjacencyRelation->n);
    for (int i = 0; i < adjacencyRelation->n; ++i) {
        kernel->weight[i] = factor;
    }
    return kernel;
}

Kernel* createGaussianKernel(AdjacencyRelation* adjacencyRelation,double variance){
    Kernel *kernel = createKernel(adjacencyRelation);
    double arg;
    double amplitude = 1.0/(2.0*M_PI*variance);
    float normalizingFactor = 0;
    for (int i = 0; i < adjacencyRelation->n; ++i) {
        //numerator
        arg = (adjacencyRelation->dx[i]*adjacencyRelation->dx[i]) + (adjacencyRelation->dy[i]*adjacencyRelation->dy[i]);
        arg /= (2*variance);
        kernel->weight[i] = (float)(exp(-arg)*amplitude);
        normalizingFactor += kernel->weight[i];
    }

    //normalizing Kernel
    for (int i = 0; i < adjacencyRelation->n; ++i) {
        kernel->weight[i] /= normalizingFactor;
    }

    return kernel;
}

Kernel* createGaussianKernel(AdjacencyRelation* adjacencyRelation,double varianceX,double varianceY){
    Kernel *kernel = createKernel(adjacencyRelation);
    double arg;
    double crossVariance = sqrt(varianceX*varianceY);
    double amplitude = 1.0/(2.0*M_PI*crossVariance);
    float normalizingFactor = 0;
    for (int i = 0; i < adjacencyRelation->n; ++i) {

        arg = (adjacencyRelation->dx[i]*adjacencyRelation->dx[i])/(2*varianceX);
        arg += (adjacencyRelation->dy[i]*adjacencyRelation->dy[i])/(2*varianceY);
        kernel->weight[i] = (float)(exp(-arg)*amplitude);
        normalizingFactor += kernel->weight[i];
    }

    //normalizing Kernel
    for (int i = 0; i < adjacencyRelation->n; ++i) {
        kernel->weight[i] /= normalizingFactor;
    }

    return kernel;
}

Kernel* createSobelHorizontalGradientKernel(int height,int width){
    AdjacencyRelation *adjacencyRelation = createRectangularAdjacency(height,width);
    Kernel* kernel = createKernel(adjacencyRelation);
    int increment = 1;
    int startValue = width/2;
    int currentValue;
    for (int dx = 0; dx < width; ++dx) {
        currentValue = startValue;

        if(dx == width/2){
            startValue -= 1;
            continue;
        }
        if(dx >= width/2){
            increment = -1;
        }else{
            increment = 1;
        }
        for (int dy = 0; dy < height; ++dy) {
            int index = (dy*width) + dx;
            kernel->weight[index] = currentValue;
            if(dy == height/2){
                increment *= -1;
            }
            currentValue += increment;
        }
        startValue -= 1;
    }

    return kernel;
}

Kernel* createSobelHorizontalGradientKernel(AdjacencyRelation* adjacencyRelation){
    //assuming rectangular symmetric kernel
    int height = (adjacencyRelation->maxDy*2)+1;
    int width =  (adjacencyRelation->maxDx*2)+1;
    return createSobelHorizontalGradientKernel(height,width);
}

Kernel* createSobelVerticalGradientKernel(int height,int width){
    AdjacencyRelation *adjacencyRelation = createRectangularAdjacency(height,width);
    Kernel* kernel = createKernel(adjacencyRelation);
    int increment = 1;
    int startValue = height/2;
    int currentValue;
    for (int dy = 0; dy < height; ++dy) {
        currentValue = startValue;

        if(dy == height/2){
            startValue -= 1;
            continue;
        }
        if(dy >= height/2){
            increment = -1;
        }else{
            increment = 1;
        }
        for (int dx = 0; dx < height; ++dx) {
            int index = (dy*width) + dx;
            kernel->weight[index] = currentValue;
            if(dx == width/2){
                increment *= -1;
            }
            currentValue += increment;
        }
        startValue -= 1;
    }

    return kernel;
}

Kernel* createSobelVerticalGradientKernel(AdjacencyRelation* adjacencyRelation){
    //assuming rectangular symmetric kernel
    int height = (adjacencyRelation->maxDy*2)+1;
    int width =  (adjacencyRelation->maxDx*2)+1;
    return createSobelVerticalGradientKernel(height,width);
}

void printKernel(Kernel* kernel){
    if(kernel == NULL){
        return;
    }
    for (int i = 0; i < kernel->adjacencyRelation->n; ++i) {
        printf("dx: %d  dy: %d w: %f\n",kernel->adjacencyRelation->dx[i],kernel->adjacencyRelation->dy[i],kernel->weight[i]);
    }
    printf("\n");
}



