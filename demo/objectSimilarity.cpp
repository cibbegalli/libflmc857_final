#include "FL.h"

int main(int argc, char **argv) {

    int k_times = 30;

    if (argc!=4){
        printf("Usage: <image1.[pgm/ppm]> <image2.[pgm/ppm]> <image3.[pgm/ppm]>\n");
        exit(-1);
    }
    FeatureVector **featureMatrix = (FeatureVector **)calloc(3,sizeof(FeatureVector *));
    AdjacencyRelation *losangeAdjacency = createLosangeAdjacency(1);

    /*
     * computa um feature vector para cada imagem
     * */
    for (int i = 0; i < 3; ++i) {
        Image *imageRGB = readImage(argv[i+1]);
        Image *imageYcbcr = convertRGBtoYCbCr(imageRGB);
        Image *imageY =  extractImageChannel(imageYcbcr,0);
        FeatureVector *vec2 = getMorphologicalPdf(imageY,losangeAdjacency,k_times);
        featureMatrix[i] = copyFeatureVector(vec2);
        destroyImage(&imageRGB);
        destroyImage(&imageYcbcr);
        destroyImage(&imageY);
        destroyFeatureVector(&vec2);
        printf("foi computado o feature vector da imagem%d\n",i);
    }


    /*
     * escreve uma arquivo txt onde caada linha no arquivo e um feature vector
     * */
    wirteFeatureVectors(featureMatrix,3, "featureVectors.txt");

    /*
     * computa a correlacao cruzada para entre os features vector da imagem0 e imagem1, e entre os
     * features vector da imagem0 e imagem2.
     * */
    //FeatureVector* vec1 = circularCrossCorrelation1D(featureMatrix[0],featureMatrix[1],true);
    //FeatureVector* vec2 = circularCrossCorrelation1D(featureMatrix[0],featureMatrix[2],true);
    FeatureVector* vec1 = matchingWithCorrelation(featureMatrix[0],featureMatrix[1],1);
    FeatureVector* vec2 = matchingWithCorrelation(featureMatrix[0],featureMatrix[2],1);

    /*
     * desaloca a memoria dos features vectors obtidos para cada imagem
     * */
    for (int i = 0; i < 3; ++i) {
        destroyFeatureVector(&featureMatrix[i]);
    }

    /*
     * Escreve os resultados obtidos do matching obtido por correlacao. Para plotar os dados obtidos va na pasta "pythonUtils" e la
     * tera um script chamado "plotVec". Entao copia e cole o arquivo R_xy dentro da pastar pythonUtils e execute
     * o script. Para executar o script basta fazer
     *
     * python plotVec.py R_xy4.txt
     *
     * Para poder executar esse script em python voce deve ter instalado os seguintes pacotes
     * *NumPy
     * *matplotlib
     *
     *
     * Se voce conseguir plotar o grafico, e esperado que o matching tenha menores valores para
     * objetos com estruturas fisicas similares.   
     * */
    FeatureVector **correlations = (FeatureVector **)calloc(3,sizeof(FeatureVector *));
    correlations[0] = vec1;
    correlations[1] = vec2;
    wirteFeatureVectors(correlations,2, "R_xy.txt");

    destroyFeatureVector(&vec1);
    destroyFeatureVector(&vec2);
    free(correlations);
    return 0;
}
