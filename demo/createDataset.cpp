#include "FL.h"
#include <iostream>
#include <string>


using namespace std;

//COIL

int main(){
	float binSize = 32;
	DirectoryManager* directoryManager = loadDirectory("/home/deangeli/Desktop/coil-100", 0);
	FeatureMatrix* featureMatrix = createFeatureMatrix(directoryManager->nfiles-2);
	char number[15];
	char number2[15];
    char filename[80];
    size_t k = 0;
    float progress;
    string str;
    string strNew;
    int k_times = 20;
    AdjacencyRelation *losangeAdjacency = createLosangeAdjacency(1);
    int label;
    int *countClass = (int*)calloc(1000,sizeof(int));
	for (size_t i = 0; i < directoryManager->nfiles; ++i)
	{
		//cout << directoryManager->files[i]->path << endl;
		// str = directoryManager->files[i]->path;
		int first = str.find("obj");
		int last = str.find("_");
		if(first == -1 || last == -1){
			continue;
		}
		strNew = str.substr(first+3,last-first-3);
		const char* strC = strNew.c_str();
		label = atoi(strC);
		// //cout << strNew << " " << label << endl;
		// Image* image = readImage(directoryManager->files[i]->path);
		// Image *imageYcbcr = convertRGBtoYCbCr(image);
  //       Image *imageY =  extractImageChannel(imageYcbcr,0);
		// FeatureVector* featureVector = getMorphologicalPdf(imageY,losangeAdjacency,k_times);
		// featureMatrix->featureVector[k] = featureVector;
		// float* aux = (float*)calloc(featureMatrix->featureVector[k]->size+1,sizeof(float));
		// for (int j = 0; j < featureMatrix->featureVector[k]->size; ++j)
		// {
		// 	aux[j] = featureMatrix->featureVector[k]->features[j];
		// }
		// aux[featureMatrix->featureVector[k]->size] = label;
		// free(featureMatrix->featureVector[k]->features);
		// featureMatrix->featureVector[k]->features = aux;
		// featureMatrix->featureVector[k]->size += 1;
		int currentIndex = countClass[label];
		sprintf(number,"%lu",currentIndex);
		sprintf(number2,"%lu",label);
        memset(filename,0,sizeof(filename));
        strcat(filename,number);
        strcat(filename,"_");
        strcat(filename,number2);
        strcat(filename,".png");
        memset(number,0,sizeof(number));
        memset(number2,0,sizeof(number2));
		writeImage(image,filename);
		//k++;
		//break;
		//progress = ((float)i)/directoryManager->nfiles;
		// if(k % 100 == 0){
		// 	//printf("%f %%\n", progress*100);
		// 	printf("%d \n",k);
		// }
		destroyImage(&image);
		// destroyImage(&imageYcbcr);
		// destroyImage(&imageY);
		//free((char*)strC);

	}
	printFeatureMatrix(featureMatrix);
	writeFeatureMatrix(featureMatrix,"datasetCoil_descriptorColor.txt");
	return 0;

}