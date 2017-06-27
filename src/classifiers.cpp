//
// Created by deangeli on 5/19/17.
//

#include "classifiers.h"

Knn_Classifier* createKnnClassifier(){
    Knn_Classifier* knnClassifier = (Knn_Classifier*)calloc(1,sizeof(Knn_Classifier));
    knnClassifier->dataMatrixTrain = NULL;
    knnClassifier->labelsTrain = NULL;
    knnClassifier->k = 0;
    knnClassifier->fitFunction = knn_Classifier_fit;
    knnClassifier->predictFunction = knn_Classifier_predict;
    knnClassifier->distanceFunction = computeNormalizedL1Norm;
    knnClassifier->distanceFunctionParameters = NULL;
    knnClassifier->nlabels = 0;
    return knnClassifier;
}
Knn_Classifier* createKnnClassifier(Matrix* dataMatrixTrain,GVector* labelsTrain){
    Knn_Classifier* knnClassifier = createKnnClassifier();
    knnClassifier->dataMatrixTrain = dataMatrixTrain;
    knnClassifier->labelsTrain = labelsTrain;
    knnClassifier->k = dataMatrixTrain->numberRows;
    knnClassifier->fitFunction = knn_Classifier_fit;
    knnClassifier->predictFunction = knn_Classifier_predict;
    knnClassifier->distanceFunction = computeNormalizedL1Norm;
    knnClassifier->distanceFunctionParameters = NULL;
    knnClassifier->nlabels = 0;
    return knnClassifier;
}

void destroyKnnClassifier(Knn_Classifier** pknnClassifier){
    Knn_Classifier* aux = *pknnClassifier;
    destroyMatrix(&(aux->dataMatrixTrain));
    destroyVector(&(aux->labelsTrain));
    destroyArgumentList(&(aux->distanceFunctionParameters));
    free(*pknnClassifier);
}

void destroyKnnClassifierForVoidPointer(void* pknnClassifier){
    Knn_Classifier* aux = ((Knn_Classifier*)pknnClassifier);
    destroyKnnClassifier(&aux);
    pknnClassifier = NULL;
}

void knnClassifierFit(Knn_Classifier* knnClassifier, Matrix* dataMatrixTrain,GVector* labelsTrain){
    knnClassifier->dataMatrixTrain =   copyMatrix(dataMatrixTrain);
    knnClassifier->labelsTrain = copyVector(labelsTrain);
}

void knn_Classifier_fit(Matrix* matrix, GVector* vector,void* knnClassifier){
    Knn_Classifier* classifier = (Knn_Classifier*)knnClassifier;
    knnClassifierFit(classifier,matrix,vector);
}

GVector* knnClassifierPredict(Matrix* matrix,Knn_Classifier* knnClassifier){
    if(knnClassifier->nlabels < 1){
        printf("[knnClassifeirPredict] nlabels is less than one\n");
    }
    if(knnClassifier->dataMatrixTrain == NULL){
        printf("[knnClassifeirPredict] classifier is not trained\n");
    }
    if(knnClassifier->labelsTrain == NULL){
        printf("[knnClassifeirPredict] classifier is not trained\n");
    }

    GVector* predictedLabels = createNullVector(matrix->numberRows,sizeof(int));
    size_t *indices = (size_t *)calloc(knnClassifier->dataMatrixTrain->numberRows,sizeof(size_t));
    size_t *count = (size_t *)calloc(knnClassifier->nlabels,sizeof(size_t));
    int maxLabelIndex = 0;
    size_t maxCount = 0;
    for (size_t testSampleIndex = 0; testSampleIndex < matrix->numberRows; ++testSampleIndex) {
        for (size_t i = 0; i < knnClassifier->nlabels; ++i) {
            count[i] = 0;
        }
        double* distances = computeAllDistancesBetweenRowAndMatrix(knnClassifier->dataMatrixTrain,
                                                                   matrix,
                                                                   testSampleIndex,
                                                                   knnClassifier->distanceFunction,
                                                                   knnClassifier->distanceFunctionParameters);
        myInsertionSortInplace(distances,knnClassifier->dataMatrixTrain->numberRows,indices);

        maxCount = 0;
        for (size_t i = 0; i < knnClassifier->k; ++i) {
            int label = VECTOR_GET_ELEMENT_AS(int,knnClassifier->labelsTrain,indices[i]);
            count[label-1]++;
            if(count[label-1] > maxCount){
                maxCount = count[label-1];
                maxLabelIndex = label;
            }
        }
        VECTOR_GET_ELEMENT_AS(int,predictedLabels,testSampleIndex) = maxLabelIndex;
        free(distances);
    }
    free(indices);
    free(count);
    return predictedLabels;
}

GVector* knn_Classifier_predict(Matrix* matrix,void* knnClassifier){
    Knn_Classifier* classifier = (Knn_Classifier*)knnClassifier;
    return knnClassifierPredict(matrix,classifier);
}

//kmeans-classifier
Kmeans_Classifier* createKmeansClassifier(){
    Kmeans_Classifier* kmeansClassifier = (Kmeans_Classifier*)calloc(1,sizeof(Kmeans_Classifier));
    kmeansClassifier->labelsTrain = NULL;
    kmeansClassifier->centroides = NULL;
    kmeansClassifier->distanceFunction = computeNormalizedL1Norm;
    kmeansClassifier->distanceFunctionParameters = NULL;
    kmeansClassifier->fitFunction = kmeans_Classifier_fit;
    kmeansClassifier->predictFunction = kmeans_Classifier_predict;
    return kmeansClassifier;
}

void destroyKmeansClassifier(Kmeans_Classifier** pkmeansClassifier){
    Kmeans_Classifier* aux = *pkmeansClassifier;
    destroyMatrix(&(aux->centroides));
    destroyVector(&(aux->labelsTrain));
    destroyArgumentList(&(aux->distanceFunctionParameters));
    free(*pkmeansClassifier);
}

void destroyKmeansClassifierForVoidPointer(void* pkmeansClassifier){
    Kmeans_Classifier* aux = ((Kmeans_Classifier*)pkmeansClassifier);
    destroyKmeansClassifier(&aux);
    pkmeansClassifier = NULL;
}

void kmeansClassifierFit(Kmeans_Classifier* kmeansClassifier, Matrix* dataMatrixTrain,GVector* labelsTrain){
    if(kmeansClassifier->nlabels < 1){
        printf("[kmeansClassifierFit] nlabels is less than one\n");
    }
    Matrix* matrix = createMatrix(kmeansClassifier->nlabels,dataMatrixTrain->numberColumns,sizeof(float));
    size_t *count = (size_t*)calloc(kmeansClassifier->nlabels,sizeof(size_t));

    for (size_t currentRow = 0; currentRow < dataMatrixTrain->numberRows; ++currentRow) {
        int label = VECTOR_GET_ELEMENT_AS(int,labelsTrain,currentRow);
        MATRIX_ADD_ROWS_AS(float,matrix,dataMatrixTrain,(label-1),currentRow);
        count[label-1]++;
    }
    for (size_t i = 0; i < matrix->numberRows; ++i) {
        for (size_t j = 0; j< matrix->numberColumns; ++j) {
            MATRIX_GET_ELEMENT_PO_AS(float,matrix,i,j) /= count[i];
        }
    }
    free(count);
    kmeansClassifier->centroides = matrix;
    kmeansClassifier->labelsTrain = copyVector(labelsTrain);
}

void kmeans_Classifier_fit(Matrix* matrix, GVector* vector,void* kmeansClassifier){
    Kmeans_Classifier* classifier = (Kmeans_Classifier*)kmeansClassifier;
    kmeansClassifierFit(classifier,matrix,vector);
}

GVector* kmeansClassifierPredict(Matrix* matrix,Kmeans_Classifier* kmeansClassifier){
    GVector* predictedLabels = createNullVector(matrix->numberRows,sizeof(int));
    if(kmeansClassifier->nlabels < 1){
        printf("[kmeansClassifierPredict] nlabels is less than one\n");
    }
    for (size_t currentRow = 0; currentRow < matrix->numberRows; ++currentRow) {
        size_t nearestRowIndex = findNearestRow(kmeansClassifier->centroides,matrix,currentRow,
                                                kmeansClassifier->distanceFunction,
                                                kmeansClassifier->distanceFunctionParameters);
        VECTOR_GET_ELEMENT_AS(int,predictedLabels,currentRow) = nearestRowIndex+1;
    }
    return predictedLabels;
}

GVector* kmeans_Classifier_predict(Matrix* matrix,void* kmeansClassifier){
    Kmeans_Classifier* classifier = (Kmeans_Classifier*)kmeansClassifier;
    return kmeansClassifierPredict(matrix,classifier);
}

//SVM
SVM_Classifier* createSVMClassifier(){
    SVM_Classifier* svm = (SVM_Classifier*)calloc(1,sizeof(SVM_Classifier));
    setSVMDefaultParameters(svm);
    return svm;
}

void setSVMDefaultParameters(SVM_Classifier* svmClassifier){
    svmClassifier->param.svm_type = C_SVC;
    svmClassifier->param.kernel_type = RBF;
    svmClassifier->param.degree = 3;
    svmClassifier->param.gamma = 0;	// 1/num_features
    svmClassifier->param.coef0 = 0;
    svmClassifier->param.nu = 0.5;
    svmClassifier->param.cache_size = 100;
    svmClassifier->param.C = 1;
    svmClassifier->param.eps = 1e-3;
    svmClassifier->param.p = 0.1;
    svmClassifier->param.shrinking = 1;
    svmClassifier->param.probability = 0;
    svmClassifier->param.nr_weight = 0;
    svmClassifier->param.weight_label = NULL;
    svmClassifier->param.weight = NULL;
    svmClassifier->crossValidation = false;
    svmClassifier->nfolds_crossValidation = 2;
    svmClassifier->useSparseDataRepresentation = false;
    svmClassifier->isClassification = true;
    svmClassifier->predict_probability = false;
}

int convertOriginalDataFormat2SvmFormatTrain(double* X,
                                             size_t nrows,size_t ncols, double* Y,
                                             struct svm_problem* prob,
                                             struct svm_node **x_data){
    size_t nElements = 0;
    size_t k;
    nElements = nrows*ncols;
    prob->l = nrows;
    prob->y = (double*)calloc(prob->l,sizeof(double));
    prob->x = (struct svm_node **)calloc(prob->l,sizeof(struct svm_node *));
    //every line must include -1 for this reason it is added prob.l
    (*x_data) = (struct svm_node*)calloc(nElements+prob->l,sizeof(struct svm_node));
    int j = 0;
    k = 0;
    int inst_max_index = 0;
    int max_index = 0;
    for (size_t row = 0; row < nrows; ++row) {
        prob->x[row] = &(*x_data)[j];
        prob->y[row] = Y[row];
        for (size_t col = 0; col < ncols; ++col) {
            (*x_data)[j].index = col+1;
            inst_max_index =  col+1;
            (*x_data)[j].value = (double) X[k];
            //printf("%d %f %f\n",k,X[k],(*x_data)[j].value);
            j++;
            k++;
        }
        if(inst_max_index > max_index){
            max_index = inst_max_index;
        }
        (*x_data)[j].index = -1;
        j++;
    }

    return max_index;
}




int convertOriginalDataFormat2SvmFormatTrain_usingSparseRepresentation(double* X,
                                                                       size_t nrows,size_t ncols, double* Y,
                                                                       struct svm_problem* prob,
                                                                       struct svm_node **x_data){
//    int max_index = 0;
//    size_t nElements = nrows*ncols;
//    size_t nzeros = 0;
//    size_t k;
//    for (k = 0; k < nrows*ncols; ++k) {
//        if( isAlmostZero(X[k]) ){
//            nzeros++;
//        }
//    }
//    prob->l = nrows;
//    nElements -= nzeros;
//    prob->y = (double*)calloc(prob->l,sizeof(double));
//    prob->x = (struct svm_node **)calloc(prob->l,sizeof(struct svm_node *));
//    //every line must include -1 for this reason it is added prob.l
//    (*x_data) = (struct svm_node*)calloc(nElements+prob->l,sizeof(struct svm_node));
//    int j = 0;
//    k = 0;
//    int inst_max_index = 0;
//
//    for (size_t row = 0; row < nrows; ++row) {
//        prob->x[row] = &((*x_data)[j]);
//        prob->y[row] = Y[row];
//        for (size_t col = 0; col < ncols; ++col) {
//            if( isAlmostZero(X[k]) ){
//                k++;
//                continue;
//            }
//            printf("%lu:%f\n",col+1,X[k]);
//            (*x_data)[j].index = col+1;
//            inst_max_index = col+1;
//            (*x_data)[j].value = X[k];
//            printf("%f %f %lu:%f\n",(*x_data)[j].index,(*x_data)[j].value,col+1,X[k]);
//            j++;
//            k++;
//        }
//        if(inst_max_index > max_index){
//            max_index = inst_max_index;
//        }
//        (*x_data)[j].index = -1;
//        j++;
//    }
//
//
//    return max_index;
    return 0;
}

//based on svm-train.c
void svmClassifierFit(SVM_Classifier* svmClassifier,double* X, size_t nrows,size_t ncols,
                      double* Y){

    /////////////// read in a problem (in svmlight format)
    int max_index = 0;
    int inst_max_index=0;
    int nElements;
    int k = 0;
    int j = 0;
    if(svmClassifier->x_space){
        free(svmClassifier->x_space);
    }

    if(svmClassifier->useSparseDataRepresentation){
        nElements = nrows*ncols;
        int nzeros=0;
        for (int element = 0; element < nElements; ++element) {
            if(isAlmostZero(X[element])){
                nzeros++;
            }
        }
        svmClassifier->prob.l = nrows;
        svmClassifier->prob.y = (double*)calloc(nrows,sizeof(double));
        svmClassifier->prob.x = (struct svm_node **)calloc(nrows,sizeof(struct svm_node *));
        svmClassifier->x_space = (struct svm_node*)calloc(nElements+nrows,sizeof(struct svm_node));
        k=0;
        j=0;
        for (size_t row = 0; row < nrows; ++row) {
            inst_max_index = -1;
            svmClassifier->prob.x[row] = &(svmClassifier->x_space[j]);
            svmClassifier->prob.y[row] = Y[row];
            for (size_t col = 0; col < ncols; ++col) {
                if(isAlmostZero(X[k])){
                    k++;
                    continue;
                }
                svmClassifier->x_space[j].index = col+1;
                svmClassifier->x_space[j].value = X[k];
                inst_max_index = col+1;
                k++;
                j++;
            }
            if(inst_max_index > max_index){
                max_index = inst_max_index;
            }
            svmClassifier->x_space[j].index = -1;
            svmClassifier->x_space[j].value = 0.0;
            j++;
        }

    }else{
        nElements = nrows*ncols;
        svmClassifier->prob.l = nrows;
        svmClassifier->prob.y = (double*)calloc(nrows,sizeof(double));
        svmClassifier->prob.x = (struct svm_node **)calloc(nrows,sizeof(struct svm_node *));
        svmClassifier->x_space = (struct svm_node*)calloc(nElements+nrows,sizeof(struct svm_node));
        k=0;
        j=0;
        for (size_t row = 0; row < nrows; ++row) {
            inst_max_index = -1;
            svmClassifier->prob.x[row] = &(svmClassifier->x_space[j]);
            svmClassifier->prob.y[row] = Y[row];
            for (size_t col = 0; col < ncols; ++col) {
                svmClassifier->x_space[j].index = col+1;
                svmClassifier->x_space[j].value = X[k];
                k++;
                j++;
            }
            svmClassifier->x_space[j].index = -1;
            svmClassifier->x_space[j].value = 0.0;
            j++;
        }

        max_index = ncols+1;
    }
    if(svmClassifier->param.gamma == 0 && max_index > 0){
        svmClassifier->param.gamma = 1.0/max_index;
    }

    if(svmClassifier->param.kernel_type == PRECOMPUTED){
        for(size_t i=0;i<nrows;i++)
        {
            if (svmClassifier->prob.x[i][0].index != 0)
            {
                fprintf(stderr,"Wrong input format: first column must be 0:sample_serial_number\n");
                exit(1);
            }
            if ((int)svmClassifier->prob.x[i][0].value <= 0 || (int)svmClassifier->prob.x[i][0].value > max_index)
            {
                fprintf(stderr,"Wrong input format: sample_serial_number out of range\n");
                exit(1);
            }
        }
    }
    /////////////// end read in a problem (in svmlight format)
    if(svmClassifier->crossValidation){
        do_cross_validation_SVM(svmClassifier);
    }


    svmClassifier->model = svm_train(&svmClassifier->prob,&svmClassifier->param);
    free(svmClassifier->prob.y);
    free(svmClassifier->prob.x);
    //free(svmClassifier->x_space);


}

double* predictLabelsLineByLine_SVM(SVM_Classifier* svmClassifier,
                                    double* X, size_t nrows,size_t ncols){
    double* prob_estimates = NULL;
    int svm_type=svm_get_svm_type(svmClassifier->model);
    int nr_class=svm_get_nr_class(svmClassifier->model);
    double *predictedLabels = (double*)calloc(nrows,sizeof(double));
    if(svmClassifier->predict_probability)
    {
        if (svm_type==NU_SVR || svm_type==EPSILON_SVR){
            printf("[svmClassifierPredict] Prob. model for test data: target value = predicted value + z,\nz: Laplace distribution e^(-|z|/sigma)/(2sigma),sigma=%g\n",svm_get_svr_probability(svmClassifier->model));
        }
        else
        {
            //int *labels=(int *) malloc(nr_class*sizeof(int));
            //svm_get_labels(svmClassifier->model,labels);
            prob_estimates = (double *) malloc(nr_class*sizeof(double));
            //free(labels);
        }
    }
    struct svm_node *x_dataRow;
    size_t k;
    double predict_label;
    //the last node should be -1 for this reason it is added a column
    x_dataRow = (struct svm_node*)calloc(ncols+1,sizeof(struct svm_node));
    int j = 0;
    k = 0;
    for (size_t row = 0; row < nrows; ++row) {
        j = 0;
        //insert row information to node
        for (size_t col = 0; col < ncols; ++col) {
            x_dataRow[j].index = col+1;
            x_dataRow[j].value = (double) X[k];
            j++;
            k++;
        }
        x_dataRow[j].index = -1;
        if (svmClassifier->predict_probability && (svm_type == C_SVC || svm_type == NU_SVC)){
            predict_label = svm_predict_probability(svmClassifier->model,x_dataRow,prob_estimates);
            predictedLabels[row] = predict_label;
        }else{
            predict_label = svm_predict(svmClassifier->model,x_dataRow);
            predictedLabels[row] = predict_label;
        }
    }
    if(svmClassifier->predict_probability){
        free(prob_estimates);
    }
    free(x_dataRow);
    return predictedLabels;
}

double* predictLabelsLineByLineUsingSparseRepresentation_SVM(SVM_Classifier* svmClassifier,
                                                             double* X, size_t nrows,size_t ncols){
    double* prob_estimates = NULL;
    int svm_type=svm_get_svm_type(svmClassifier->model);
    int nr_class=svm_get_nr_class(svmClassifier->model);
    double *predictedLabels = (double*)calloc(nrows,sizeof(double));
    if(svmClassifier->predict_probability)
    {
        if (svm_type==NU_SVR || svm_type==EPSILON_SVR){
            printf("[svmClassifierPredict] Prob. model for test data: target value = predicted value + z,\nz: Laplace distribution e^(-|z|/sigma)/(2sigma),sigma=%g\n",svm_get_svr_probability(svmClassifier->model));
        }
        else
        {
            //int *labels=(int *) malloc(nr_class*sizeof(int));
            //svm_get_labels(svmClassifier->model,labels);
            prob_estimates = (double *) malloc(nr_class*sizeof(double));
            //free(labels);
        }
    }
    struct svm_node *x_dataRow;
    size_t k;
    double predict_label;
    //the last node should be -1 for this reason it is added a column
    x_dataRow = (struct svm_node*)calloc(ncols+1,sizeof(struct svm_node));
    int j = 0;
    k = 0;
    for (size_t row = 0; row < nrows; ++row) {
        j = 0;
        //insert row information to node
        for (size_t col = 0; col < ncols; ++col) {
            if(isAlmostZero(X[k])){
                k++;
                continue;
            }
            x_dataRow[j].index = col+1;
            x_dataRow[j].value = (double) X[k];
            j++;
            k++;
        }
        x_dataRow[j].index = -1;

        if (svmClassifier->predict_probability && (svm_type == C_SVC || svm_type == NU_SVC)){
            predict_label = svm_predict_probability(svmClassifier->model,x_dataRow,prob_estimates);
            predictedLabels[row] = predict_label;
        }else{
            predict_label = svm_predict(svmClassifier->model,x_dataRow);
            predictedLabels[row] = predict_label;
        }
    }
    if(svmClassifier->predict_probability){
        free(prob_estimates);
    }
    free(x_dataRow);
    return predictedLabels;

}

//based on svm-predict.c
double* svmClassifierPredict(SVM_Classifier* svmClassifier,double* X, size_t nrows,size_t ncols){


    if(svmClassifier->predict_probability)
    {
        if(svm_check_probability_model(svmClassifier->model)==0)
        {
            printf("[svmClassifierPredict] Model does not support probabiliy estimates\n");
            return NULL;
        }
    }
    else
    {
        if(svm_check_probability_model(svmClassifier->model)!=0) {
            printf("[svmClassifierPredict] Model supports probability estimates, but disabled in prediction.\n");
        }
    }
    double* predictedData = NULL;
    if(svmClassifier->useSparseDataRepresentation){
        predictedData = predictLabelsLineByLineUsingSparseRepresentation_SVM(svmClassifier,X,nrows,ncols);
    }else{
        predictedData = predictLabelsLineByLine_SVM(svmClassifier,X,nrows,ncols);
    }
    return predictedData;
}

void destroySVMClassifier(SVM_Classifier** svmClassifier){
    if( (*svmClassifier)->model ){
        svm_free_and_destroy_model(&(*svmClassifier)->model);
    }
    if( (*svmClassifier)->param.weight  ||  (*svmClassifier)->param.weight_label){
        svm_destroy_param(&(*svmClassifier)->param);
    }
    if( (*svmClassifier)->x_space){
        free((*svmClassifier)->x_space);
    }
    free((*svmClassifier));
    (*svmClassifier) = NULL;
}

void destroySVMClassifierForVoidPointer(void* psvmClassifier){
    SVM_Classifier* aux = ((SVM_Classifier*)psvmClassifier);
    destroySVMClassifier(&aux);
    psvmClassifier = NULL;
}

void svm_Classifier_fit(Matrix* matrix, GVector* vector,void* svmClassifier){
    SVM_Classifier* classifier = (SVM_Classifier*)svmClassifier;
    double* data = (double*)calloc(matrix->numberElements,sizeof(double));
    double* labels = (double*)calloc(matrix->numberRows,sizeof(double));
    int k = 0;
    for (size_t row = 0; row < matrix->numberRows; ++row) {
        labels[row] = VECTOR_GET_ELEMENT_AS(int,vector,row);
        for (size_t col = 0; col < matrix->numberColumns; ++col) {
            data[k] = MATRIX_GET_ELEMENT_PO_AS(float,matrix,row,col);
            k++;
        }
    }

    svmClassifierFit(classifier,data,matrix->numberRows,matrix->numberColumns,labels);
    free(data);
    free(labels);
}

GVector* svm_Classifier_predict(Matrix* matrix,void* svmClassifier){
    SVM_Classifier* classifier = (SVM_Classifier*)svmClassifier;
    double* data = (double*)calloc(matrix->numberElements,sizeof(double));
    for (size_t k = 0; k < matrix->numberElements; ++k) {
        data[k] = (double)MATRIX_GET_ELEMENT_BI_AS(float,matrix,k);
    }

    double* predictedLabels = svmClassifierPredict(classifier,data,matrix->numberRows,matrix->numberColumns);
    GVector* predictedLabels_int = NULL;
    if(predictedLabels){
        predictedLabels_int = createNullVector(matrix->numberRows,sizeof(int));
        for (size_t row = 0; row < matrix->numberRows; ++row) {
            VECTOR_GET_ELEMENT_AS(int,predictedLabels_int,row) = (predictedLabels[row]+0.1);
        }

    }
    free(data);
    free(predictedLabels);
    return predictedLabels_int;
}

void do_cross_validation_SVM(SVM_Classifier* svmClassifier)
{
    int i;
    int total_correct = 0;
    double total_error = 0;
    double sumv = 0, sumy = 0, sumvv = 0, sumyy = 0, sumvy = 0;
    double *target = Malloc(double,svmClassifier->prob.l);

    svm_cross_validation(&svmClassifier->prob,&svmClassifier->param,
                         svmClassifier->nfolds_crossValidation,target);
    if(svmClassifier->param.svm_type == EPSILON_SVR ||
            svmClassifier->param.svm_type == NU_SVR)
    {
        for(i=0;i<svmClassifier->prob.l;i++)
        {
            double y = svmClassifier->prob.y[i];
            double v = target[i];
            total_error += (v-y)*(v-y);
            sumv += v;
            sumy += y;
            sumvv += v*v;
            sumyy += y*y;
            sumvy += v*y;
        }
        printf("Cross Validation Mean squared error = %g\n",total_error/svmClassifier->prob.l);
        printf("Cross Validation Squared correlation coefficient = %g\n",
               ((svmClassifier->prob.l*sumvy-sumv*sumy)*(svmClassifier->prob.l*sumvy-sumv*sumy))/
               ((svmClassifier->prob.l*sumvv-sumv*sumv)*(svmClassifier->prob.l*sumyy-sumy*sumy))
        );
    }
    else
    {
        for(i=0;i<svmClassifier->prob.l;i++)
            if(target[i] == svmClassifier->prob.y[i])
                ++total_correct;
        printf("Cross Validation Accuracy = %g%%\n",100.0*total_correct/svmClassifier->prob.l);
    }
    free(target);
}




