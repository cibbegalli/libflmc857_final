#include "FL.h"

#include "parameters.cpp"

int main(int argc, char **argv) {
    //size_t numberOfVisualWords = 500;
    size_t numberOfVisualWords = NUMBER_OF_VISUAL_WORDS;


    //Caminhos onde esta o arquivo txt gerado pelo o script python "selec_samples2.py"
    /*char const* const fileName_createDict = "/home/cibelle/Documentos/Databases/COIL100/train_paths.txt";
    char const* const fileName_createTrain = "/home/cibelle/Documentos/Databases/COIL100/train_paths.txt";
    char const* const fileName_createTest = "/home/cibelle/Documentos/Databases/COIL100/test_paths.txt";*/
    char const* const fileName_createDict = FILENAME_DICT;
    char const* const fileName_createTrain = FILENAME_TRAIN;
    char const* const fileName_createTest = FILENAME_TEST;
    char const* const fileName_superpixels = FILENAME_SUPERPIXEL;
    
    //cada posicao do vetor tem uma string para o caminho de uma imagem
    GVector* vectorSamplesUsed2CreateDict =  splitsLinesInTextFile(fileName_createDict);
    GVector* vectorSamplesUsed2TrainClassifier =  splitsLinesInTextFile(fileName_createTrain);
    GVector* vectorSamplesUsed2TestClassifier =  splitsLinesInTextFile(fileName_createTest);
    GVector* vectorSamplesUsed2Superpixel =  splitsLinesInTextFile(fileName_superpixels);

    //apenas checkando se o vetor vazio. Caso o vetor esteja vazio, talvez seu caminho ate o arquivo
    //txt nao esteja correto
    if(vectorSamplesUsed2CreateDict->size == 0){
        printf("no path found");
        return -1;
    }

    if(vectorSamplesUsed2TrainClassifier->size == 0){
        printf("no path found");
        return -1;
    }

    if(vectorSamplesUsed2TestClassifier->size == 0){
        printf("no path found");
        return -1;
    }

    printf("Escolha descritor:\n");
    printf("\t1: Histograma de cor\n");
    printf("\t2: HOG\n");
    printf("\t3: Histograma de cor + HOG\n");
    int option;
    do {
        scanf("%d", &option);
    } while(option < 1 && option > 3);

    //pipeline para a construncao do dicionario. Para mais detalhes olhe a imagem que esta
    //em data/bowArquiteturaImplementada.png

    //a estrutura bow manager e encarregada de fazer o processo do bow
    BagOfVisualWordsManager* bowManager = createBagOfVisualWordsManager();

    ////////////////////////////////////////////////////////////////////////
    //Passando os vetores que contem os caminhos das imagens para...
    bowManager->pathsToImages_dictionery = vectorSamplesUsed2CreateDict;//criar o dicionario
    bowManager->pathsToImages_train = vectorSamplesUsed2TrainClassifier;//treinar o classificador
    bowManager->pathsToImages_test = vectorSamplesUsed2TestClassifier;//testar o classificador
    bowManager->pathsToImages_isf = vectorSamplesUsed2Superpixel;
    //////////////////////////////////////////////////////////////////////

    printf("Escolha método grid:\n");
    printf("\t1: Grid simples\n");
    printf("\t2: Grid superpixel com gradiente\n");
    printf("\t3: Grid superpixel com histograma de cor\n");
    int optionGrid;
    do {
        scanf("%d", &optionGrid);
    } while(optionGrid < 1 && optionGrid > 3);

    if(optionGrid == 1) {

        //////////////////////////////////////////////////////////////////////
        //metodo de sampling que vai ser usado para criar os patchs. Se vc passar NULL aqui o estrutura
        // do bow vai criar um vetor de tamanho 1 onde o unico elemento desse vetor vai ser a imagem.
        bowManager->imageSamplerFunction = gridSamplingBow;//ponteiro da funcao para o sampling

        //Nesta demo o metodo de sampling  usado é o grid. Entao eu vou criar um argument list
        //para colocar os parametros do metodo de grinding que eu fiz.
        //Note que o cabecalho geral para a funcao de sammpling e
        //GVector* minhaFuncaoDeSampling(Image* image, BagOfVisualWordsManager* bagOfVisualWordsManager);
        ArgumentList* gridSamplingArguments = createArgumentList();
        ARGLIST_PUSH_BACK_AS(size_t,gridSamplingArguments,PATCH_SIZE_X); //patch size X
        ARGLIST_PUSH_BACK_AS(size_t,gridSamplingArguments,PATCH_SIZE_Y); //patch size Y
        bowManager->argumentListOfSampler = gridSamplingArguments;//passando a lista de argumentos para o bow manager
        //////////////////////////////////////////////////////////
    } else if(optionGrid == 2) {
        bowManager->imageSamplerFunction = samplingByISFBow;//ponteiro da funcao para o sampling
        ArgumentList* gridSamplingArguments = createArgumentList();
        ARGLIST_PUSH_BACK_AS(size_t,gridSamplingArguments,GRID_GRADIENT_THRESOULD); //patch size X
        ARGLIST_PUSH_BACK_AS(size_t,gridSamplingArguments,PATCH_SIZE_Y); //patch size Y
        ARGLIST_PUSH_BACK_AS(size_t,gridSamplingArguments,PATCH_SIZE_X); //patch size X
        ARGLIST_PUSH_BACK_AS(size_t,gridSamplingArguments,GRID_WIDTH_ADJ_REL); //patch size Y
        ARGLIST_PUSH_BACK_AS(size_t,gridSamplingArguments,GRID_HEIGHT_ADJ_REL); //patch size X
        bowManager->argumentListOfSampler = gridSamplingArguments;//passando a lista de argumentos para o bow manager
    } else if(optionGrid == 3) {
        bowManager->imageSamplerFunction = samplingByISFColorBow;//ponteiro da funcao para o sampling
        ArgumentList* gridSamplingArguments = createArgumentList();
        ARGLIST_PUSH_BACK_AS(size_t,gridSamplingArguments,PATCH_SIZE_Y); //patch size Y
        ARGLIST_PUSH_BACK_AS(size_t,gridSamplingArguments,PATCH_SIZE_X); //patch size X
        ARGLIST_PUSH_BACK_AS(size_t,gridSamplingArguments,DIFF_SUPER_PIXELS); //patch size Y
        ARGLIST_PUSH_BACK_AS(size_t,gridSamplingArguments,BIN_SIZE); //patch size X
        ARGLIST_PUSH_BACK_AS(size_t,gridSamplingArguments,SPARSE_HIST); //patch size X
        bowManager->argumentListOfSampler = gridSamplingArguments;//passando a lista de argumentos para o bow manager
    }

    /////////////////////////////////////////////////////////////////////////////
    //Essa função serve como um garbage collector para o metodo do sampling. Ao final de
    //cada iteracao, ela limpa da memoria os patchs gerados.
    //Se por acaso seu metodo de sampling nao gerar um vetor de imagens (GVector* de Image*),
    //voce pode passar NULL, porém fique consciente que vai ter um pouco de memory leak.
    //Ao final do programa seu sistema operional vai limpar toda a sujeira.
    bowManager->freeFunction2SamplerOutput = destroyImageVoidPointer;
    //bowManager->freeFunction2SamplerOutput = NULL;
    /////////////////////////////////////////////////


    ///////////////////////////////////////
    if(option == 1) { // histograma de cor
        
        bowManager->featureExtractorFunction = computeColorHistogramBow;//ponteiro da funcao para a extracao de features
        ArgumentList* colorFeatureExtractorArguments = createArgumentList();
        size_t nbins = COLOR_HISTOGRAM_NBINS;
        ARGLIST_PUSH_BACK_AS(size_t,colorFeatureExtractorArguments,nbins); //nBins per channel
        ARGLIST_PUSH_BACK_AS(size_t,colorFeatureExtractorArguments,nbins*nbins*nbins); //total number of channels
        bowManager->argumentListOfFeatureExtractor = colorFeatureExtractorArguments; //passando a lista de argumentos do feature extractor para o bow manager

    } else if(option == 2) { // hog
        
        bowManager->featureExtractorFunction = computeHogBow;
        ArgumentList* hogFeatureExtractorArguments = createArgumentList();
        ARGLIST_PUSH_BACK_AS(size_t,hogFeatureExtractorArguments,HOG_BLOCKS_X); //número de blocos na horizontal
        ARGLIST_PUSH_BACK_AS(size_t,hogFeatureExtractorArguments,HOG_BLOCKS_Y); //número de blocos na vertical
        ARGLIST_PUSH_BACK_AS(size_t,hogFeatureExtractorArguments,HOG_THETA); //theta 
        bowManager->argumentListOfFeatureExtractor = hogFeatureExtractorArguments; //passando a lista de argumentos do feature extractor para o bow manager

    } else { // histograma de cor + hog
        bowManager->featureExtractorFunction = computeColorAndHogBow;
        ArgumentList* colorAndhogFeatureExtractorArguments = createArgumentList();
        ARGLIST_PUSH_BACK_AS(size_t,colorAndhogFeatureExtractorArguments,HOG_BLOCKS_X); //número de blocos na horizontal
        ARGLIST_PUSH_BACK_AS(size_t,colorAndhogFeatureExtractorArguments,HOG_BLOCKS_Y); //número de blocos na vertical
        ARGLIST_PUSH_BACK_AS(size_t,colorAndhogFeatureExtractorArguments,HOG_THETA); //theta
        size_t nbins = COLOR_HISTOGRAM_NBINS;
        ARGLIST_PUSH_BACK_AS(size_t,colorAndhogFeatureExtractorArguments,nbins); //nBins per channel
        ARGLIST_PUSH_BACK_AS(size_t,colorAndhogFeatureExtractorArguments,nbins*nbins*nbins); //total number of channels
        bowManager->argumentListOfFeatureExtractor = colorAndhogFeatureExtractorArguments; //passando a lista de argumentos do feature extractor para o bow manager
    
    }
        
    ///////////////////////////////////////////////////////
    //Existem muitas maneiras de computar distancias entre pontos e vetores. A mais comum delas talvez
    //seja a distancia Euclidianda (norma l2). Neste exemplo eu vou usar a norma l1.
    //Quando vc implementar seu metodo de sampling, feature extraxtion, ou clustering, voce
    //pode usar essa funcao distancia.
    bowManager->distanceFunction = computeNormalizedL1Norm;
    bowManager->argumentListOfDistanceFunction = NULL;
    ////////////////////////////////////////////////////


    /////////////////////////////////////////////////////
    //Aqui precisamos definir qual funcao de clsutering vamos usar para encontrar as palavras
    //do dicionario. Eu optei de usar o kmeans clustering. Meu metodo do kmeans recebe 6 parametros,
    //desta forma eu preciso criar uma ArgumentList com 6 parametros.
    //Note que o cabecalho geral para a funcao de clustering e
    //typedef Matrix* minhaFuncaoDeClustering(Matrix* outputFeatureExtractor_allSamples, BagOfVisualWordsManager* bagOfVisualWordsManager);
    bowManager->clusteringFunction = kmeansClusteringBow;
    ArgumentList* clusteringMethodArguments = createArgumentList();
    ARGLIST_PUSH_BACK_AS(size_t,clusteringMethodArguments,numberOfVisualWords); //number of words
    ARGLIST_PUSH_BACK_AS(size_t,clusteringMethodArguments,100); //maximum number of iterations
    ARGLIST_PUSH_BACK_AS(double,clusteringMethodArguments,0.0001); //tolerance
    ARGLIST_PUSH_BACK_AS(int,clusteringMethodArguments,0); //seed
    ARGLIST_PUSH_BACK_AS(DistanceFunction,clusteringMethodArguments,computeNormalizedL1Norm); //distance function
    ARGLIST_PUSH_BACK_AS(ArgumentList*,clusteringMethodArguments,NULL); //
    bowManager->argumentListOfClustering = clusteringMethodArguments;
    ///////////////////////////////////////////////////////////////

    ////////////
    //computa o dicionario
    computeDictionery(bowManager);
    /////////////

    //////////////////////
    //define a funcao para montar o histograma
    bowManager->mountHistogramFunction = computeCountHistogram_bow;
    bowManager->argumentListOfHistogramMounter = NULL;
    ///////////////////////


    /////////////////////////////////////////////////
    //criar um classificador e define os parametros
    //do classficiador. Em seguida, o bow manager recebe
    //o ponteiro do classificador. Desta forma o classificador
    //podera ser usado internamente dentro do bow manager.

    //SVM Classifier
    SVM_Classifier* classifiersvm = createSVMClassifier();
    classifiersvm->param.kernel_type = RBF;
    classifiersvm->param.gamma = 3.5;
    bowManager->classifier = (void*)classifiersvm;
    bowManager->fitFunction = svm_Classifier_fit;
    bowManager->storeTrainData = false;
    bowManager->predictFunction = svm_Classifier_predict;
    bowManager->storePredictedData = false;
    bowManager->freeFunctionClassifier = destroySVMClassifierForVoidPointer;
    //////////////////////////////////////

    ///////
    //monta os histogramas, le os label e em seguida treina o classificador
    trainClassifier(bowManager);
    //////////

    /////////////////////////////////////////////////////
    //monta os histogramas e usa o classificador treinado para
    //classificar as amostras do conjunto de teste
    GVector* labelsPredicted = predictLabels(bowManager);
    //////////////////////////

    //////////////////////////
    //Le os true labels das imagens e checa com os labels predizidos pelo o classificador.
    //computa uma simples acuracia (numero de amostras rotuladas corretamente / numero de amostras do conjunto)
    
    

    GVector* trueLabels = createNullVector(bowManager->pathsToImages_test->size,sizeof(int));
 
    int hit = 0;
    printf("file | predicted true\t\tcorrect\n");
    char symbol;
    for (size_t index = 0; index < bowManager->pathsToImages_test->size; ++index) {
        symbol = 'X';
        char * path = VECTOR_GET_ELEMENT_AS(char*,bowManager->pathsToImages_test,index);
        VECTOR_GET_ELEMENT_AS(int,trueLabels,index) = findTrueLabelInName(path);
        int trueLabel = VECTOR_GET_ELEMENT_AS(int,trueLabels,index);
        int predictedLabel = VECTOR_GET_ELEMENT_AS(int,labelsPredicted,index);
        if(trueLabel == predictedLabel){
            hit++;
            symbol = 'O';

        } 
        printf("%s | %03d %03d\t\t%c\n",
               path,
               VECTOR_GET_ELEMENT_AS(int,labelsPredicted,index),
               VECTOR_GET_ELEMENT_AS(int,trueLabels,index),symbol
        );
    }

    printf("Acucaria:%f\n", ((double)hit)/bowManager->pathsToImages_test->size);
    /////////////////////////////////////
    destroyBagOfVisualWordsManager(&bowManager);
    destroyVector(&trueLabels);
    destroyVector(&labelsPredicted);
    
    return 0;

}


