#include "sampling.h"

GVector* gridSampling(Image* image, size_t patchSizeX,size_t patchSizeY) {
    size_t nPatchs_X = image->nx/patchSizeX;
    size_t nPatchs_Y = image->ny/patchSizeY;
    size_t nPatchs = nPatchs_X*nPatchs_Y;
    GVector* vector_images = createNullVector(nPatchs,sizeof(Image*));
    int k = 0;
    for (size_t y = 0; y <= (size_t)image->ny-patchSizeY; y +=patchSizeY) {
        for (size_t x = 0; x <= (size_t)image->nx-patchSizeX; x += patchSizeX) {
            VECTOR_GET_ELEMENT_AS(Image*,vector_images,k) = extractSubImage(image,x,y,patchSizeX,patchSizeY,true);
            k++;
        }
    }
    return vector_images;
}

int numberLabels(Image* img, int x0, int y0, AdjacencyRelation* adjRel, int scalingFactor) {
    FeatureVector* labels = createFeatureVector(scalingFactor);
    setValueInFeatureVector(labels, 0.0);

    for(int i=0; i<adjRel->n; i++) {
        int x = x0+adjRel->dx[i];
        int y = y0+adjRel->dy[i];
        if(isValidPixelCoordinate(img, x, y)) {
            float label = imageVal(img, x, y);
            labels->features[(int)(label-1.0)] = 1.0; // imagem superpixel rotulada de 1 a maxLabel
        } else {
            destroyFeatureVector(&labels);
            return 0;
        }
    }

    int count = 0;
    for(int k=0; k<scalingFactor; k++) {
        if(labels->features[k] == 1.0)
            count++;
    }
    destroyFeatureVector(&labels);
    return count;
}

GVector* samplingBySuperpixelAndGradient(Image* image, Image* imageISF, float threshold, int patchSizeX, int patchSizeY, int width_adj_rel, int height_adj_rel) {
    
    Image *mag, *phase;
    computeGradient(image , &mag, &phase); 
    //writeImage(magOrig, "mag4_orig.pgm");

    //normalizar imagem magnitude
    /*float max = magOrig->scalingFactor;
    float min = 0.0f;
    float L = 255.0;
    Image* mag = createImage(magOrig->nx, magOrig->ny, magOrig->nchannels);
    for(int i=0; i<image->nx; i++) {
        for(int j=0; j<image->ny; j++) {
            mag->channel[0][(j*image->nx) + i] = L*imageVal(magOrig, i, j)/max;
        }
    }
    //writeImage(mag, "mag4_norm.pgm");
    destroyImage(&magOrig);*/

    int scalingFactor = imageISF->scalingFactor;
    /*for(int i=0; i<image->nx; i++) {
        for(int j=0; j<image->ny; j++) {
            if(((int)imageVal(imageISF, i, j)) > scalingFactor) 
                scalingFactor = imageVal(imageISF, i, j);
        }
    }*/

    AdjacencyRelation *adjRel = createRectangularAdjacency(height_adj_rel, width_adj_rel);
    int nPatchs = 0;
    for(int i=0; i<image->nx; i+= width_adj_rel/2) {
        for(int j=0; j<image->ny; j+=height_adj_rel/2) {
            int nLabels = numberLabels(imageISF, i, j, adjRel, scalingFactor);
            if(nLabels >= 3) {

                if(imageVal(mag, i, j) > threshold) {
                    nPatchs++;
                    i += width_adj_rel;
                    j += height_adj_rel;
                    //i += patchSizeX/2 - width_adj_rel;
                    //j += patchSizeY/2 - height_adj_rel;
                }

            }
        }
    }
    if(nPatchs == 0) {
        GVector* vector_images = createNullVector(1,sizeof(Image*));
        int xc = image->nx/2, yc = image->ny/2; //coordenada do centro
        int x = xc - patchSizeX/2, y = yc - patchSizeY/2;

        VECTOR_GET_ELEMENT_AS(Image*,vector_images,0) = extractSubImage(image,x,y,patchSizeX,patchSizeY,true);
        
        destroyImage(&mag);
        destroyImage(&phase);
        destroyAdjacencyRelation(&adjRel);

        return vector_images;
    }

    GVector* vector_images = createNullVector(nPatchs,sizeof(Image*));
    int k=0;
    for(int i=0; i<image->nx; i+= width_adj_rel/2) {
        for(int j=0; j<image->ny; j+=height_adj_rel/2) {
            int nLabels = numberLabels(imageISF, i, j, adjRel, scalingFactor);
            if(nLabels >= 3) {

                if(imageVal(mag, i, j) > threshold) {

                    int x = i-patchSizeX/2;
                    int y = j-patchSizeY/2;
                    if(x < 0) x = 0;
                    if(y < 0) y = 0;
                    if(x + patchSizeX > image->nx) {
                        x = image->nx - patchSizeX;
                    }
                    if(y + patchSizeY > image->ny) {
                        y = image->ny - patchSizeY;
                    }

                    
                    
                    Image *sub = extractSubImage(image,x,y,patchSizeX,patchSizeY,true);


                    VECTOR_GET_ELEMENT_AS(Image*,vector_images,k) = sub;
                    k++;
                    
                    char name[20];
                    sprintf(name, "%d_sub.ppm", i + j*image->nx);
                    writeImage(sub, name);
                    destroyImage(&sub);
                    
                    i += width_adj_rel;
                    j += height_adj_rel;
                    //i += patchSizeX/2 - width_adj_rel;
                    //j += patchSizeY/2 - height_adj_rel;

                }
            }
        }
    }


    destroyImage(&mag);
    destroyImage(&phase);
    destroyAdjacencyRelation(&adjRel);

    return vector_images;
}

GVector* samplingBySuperPixel(Image* realImage, Image* superPixel, size_t patchSizeX,size_t patchSizeY, size_t DIFF_SUPER_PIXELS, size_t binSize, size_t SPARSE_HIST) {
    
    //int patchSizeX = 20; // mudar de acordo com base utilizada
    //int patchSizeY = 20;
    //Image* image = readImagePGM("out2.pgm");
    Image* subImage = NULL;
    //Image* realImage = readImage("obj2__140.png");
    Histogram* histo;
    char number[15];
    char filename[80];
    
    int k = 0;
    int count=0;
    int max;
    int patch = 0;
    float aux;
    //Cria vetor de imagens para 1 imagem
    GVector* vector_images = NULL; //createNullVector(1,sizeof(Image*));
        
    for (int y = 0; y < superPixel->ny; y += (int)(patchSizeY)) {
        for (int x = 0; x < superPixel->nx; x += (int)(patchSizeX)) {
            
            //conta quantos super pixels com valor diferente
            count = 0;
            aux = (superPixel->channel[0][k]);
            for (int j = y; j < y + (int)(patchSizeY); j++) {
                for (int i = x; i < x + (int)(patchSizeX); i++) {
                    if ((superPixel->channel[0][k]) != aux)
                        count ++;
                    k++;
                }
                
            }
            
            //imprime patch
            if (count >= DIFF_SUPER_PIXELS){
                
                subImage = extractSubImage(realImage, x, y, patchSizeX, patchSizeY,true);
                histo = computeHistogram(subImage,binSize, true);
                max = 0;
                for(int i = 0; i<histo->n; i++){
                    if (histo->val[i] != 0)
                        max++;
                }
                
                if (max > SPARSE_HIST){ // quantidade de valores no histograma -> diferentes cores na imagem
                    patch++;
                }
                destroyImage(&subImage);
                destroyHistogram(&histo);
                
            }
            
        }
    }
    k=0;
    vector_images = createNullVector(patch,sizeof(Image*));
    patch = 0;
    for (int y = 0; y < superPixel->ny; y += (int)(patchSizeY)) {
        for (int x = 0; x < superPixel->nx; x += (int)(patchSizeX)) {
            
            //conta quantos super pixels com valor diferente
            count = 0;
            aux = (superPixel->channel[0][k]);
            for (int j = y; j < y + (int)(patchSizeY); j++) {
                for (int i = x; i < x + (int)(patchSizeX); i++) {
                    if ((superPixel->channel[0][k]) != aux)
                        count ++;
                    k++;
                }
                
            }
            
            //imprime patch
            if (count >= DIFF_SUPER_PIXELS){
                
                subImage = extractSubImage(realImage, x, y, patchSizeX, patchSizeY,true);
                histo = computeHistogram(subImage,binSize, true);
                max = 0;
                for(int i = 0; i<histo->n; i++){
                    if (histo->val[i] != 0)
                        max++;
                }
                
                if (max > SPARSE_HIST){ // quantidade de valores no histograma -> diferentes cores na imagem
                    VECTOR_GET_ELEMENT_AS(Image*,vector_images,patch) = subImage;
                    //VECTOR_GET_ELEMENT_AS(Image*,vector_images,patch) = subImage;
                    // printa os patchs
                    /*sprintf(number,"%d",patch);
                    memset(filename,0,sizeof(filename));
                    strcat(filename,"patch/patch");
                    strcat(filename,number);
                    strcat(filename,".ppm");
                    memset(number,0,sizeof(number));
                    writeImage(subImage,filename);*/
                    //destroyImage(&subImage);
                    patch++;
                } else {
                    destroyImage(&subImage);
                }
                destroyHistogram(&histo);
                
            }
            
        }
    }
    //destroyImage(&superPixel);
    //destroyImage(&realImage);
    //Image *val = (Image)(vector_images->data[0]);
    
    
    return vector_images;
}
/*
GVector* samplingBySuperPixel(Image* realImage, Image* superPixel, size_t patchSizeX,size_t patchSizeY, size_t DIFF_SUPER_PIXELS, size_t binSize, size_t SPARSE_HIST) {
    
    //int patchSizeX = 20; // mudar de acordo com base utilizada
    //int patchSizeY = 20;
    //Image* image = readImagePGM("out2.pgm");
    Image* subImage = NULL;
    //Image* realImage = readImage("obj2__140.png");
    Histogram* histo;
    //char number[15];
    //char filename[80];
    
    int k = 0;
    int count=0;
    int max;
    int patch = 0;
    float aux;
    int index = 0;
    //Cria vetor de imagens para 1 imagem
    GVector* vector_images = createNullVector(1,sizeof(Image*));
    
    for (int y = 0; y < superPixel->ny; y += (int)patchSizeY) {
        for (int x = 0; x < superPixel->nx; x += (int)patchSizeX) {
            
            //conta quantos super pixels com valor diferente
            count = 0;
            aux = (superPixel->channel[0][k]);
            for (int j = y; j < y + (int)patchSizeY; j++) {
                for (int i = x; i < x + (int)patchSizeX; i++) {
                    if ((superPixel->channel[0][k]) != aux)
                        count ++;
                }
                k++;
            }
            //imprime patch
            if (count >= DIFF_SUPER_PIXELS){
                
                subImage = extractSubImage(realImage, x, y, patchSizeX, patchSizeY,true);
                histo = computeHistogram(subImage,binSize, true);
                max = 0;
                for(int i = 0; i<histo->n; i++){
                    if (histo->val[i] != 0)
                        max++;
                }
                
                if (max > SPARSE_HIST){ // quantidade de valores no histograma -> diferentes cores na imagem
                    setVectorCapacity(vector_images, patch+1); // amplia capacidade do vetor
                    
                    VECTOR_GET_ELEMENT_AS(Image*,vector_images,index) = subImage;
                    index++;
                     // printa os patchs
                    //sprintf(number,"%d",patch);
                    //memset(filename,0,sizeof(filename));
                    //strcat(filename,"patch/patch");
                    //strcat(filename,number);
                    //strcat(filename,".ppm");
                    //memset(number,0,sizeof(number));
                    //writeImage(subImage,filename);
                    //destroyImage(&subImage);
                    patch++;
                } else {
                    destroyImage(&subImage);    
                }
                destroyHistogram(&histo);
                
            }
            
        }
    }
    //destroyImage(&superPixel);
    //destroyImage(&realImage);
    
    return vector_images;
}*/