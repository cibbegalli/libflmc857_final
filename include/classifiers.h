//
// Created by deangeli on 5/19/17.
//

#ifndef _CLASSIFIERS_H
#define _CLASSIFIERS_H

#include "matrix.h"
#include "matrixUtil.h"
#include "vector.h"
#include "vectorUtil.h"
#include "distanceFunctions.h"
#include "argumentList.h"
#include "svm.h"

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

typedef void (*FitFunction) (Matrix* matrix,GVector* vector,void* classifierObj);
typedef GVector* (*PredictFunction) (Matrix* matrix,void* classifierObj);


typedef struct _knn_Classifier{
    Matrix* dataMatrixTrain;
    GVector* labelsTrain;
    FitFunction fitFunction;
    PredictFunction predictFunction;
    DistanceFunction distanceFunction;
    ArgumentList* distanceFunctionParameters;
    size_t k;
    size_t nlabels;
}Knn_Classifier;


typedef struct _kmeans_Classifier{
    Matrix* centroides;
    GVector* labelsTrain;
    FitFunction fitFunction;
    PredictFunction predictFunction;
    DistanceFunction distanceFunction;
    ArgumentList* distanceFunctionParameters;
    size_t nlabels;
}Kmeans_Classifier;


//"param.svm_type : set type of SVM (default 0)\n"
//"	0 -- C-SVC		(multi-class classification)\n"
//"	1 -- nu-SVC		(multi-class classification)\n"
//"	2 -- one-class SVM\n"
//"	3 -- epsilon-SVR	(regression)\n"
//
//"param.kernel_type : set type of kernel function (default 2)\n"
//"	0 -- linear: u'*v\n"
//"	1 -- polynomial: (gamma*u'*v + coef0)^degree\n"
//"	2 -- radial basis function: exp(-gamma*|u-v|^2)\n"
//"	3 -- sigmoid: tanh(gamma*u'*v + coef0)\n"
//" param.degree : set degree in kernel function (default 3)\n"
//" param.gamma : set gamma in kernel function (default 1/num_features)\n"
//" param.coef0 : set coef0 in kernel function (default 0)\n"
//" param.cost : set the parameter C of C-SVC, epsilon-SVR, and nu-SVR (default 1)\n"
//" param.nu : set the parameter nu of nu-SVC, one-class SVM, and nu-SVR (default 0.5)\n"
//" param.epsilon : set the epsilon in loss function of epsilon-SVR (default 0.1)\n"
//" param.cachesize : set cache memory size in MB (default 100)\n"
//" param.epsilon : set tolerance of termination criterion (default 0.001)\n"
//" param.shrinking : whether to use the shrinking heuristics, 0 or 1 (default 1)\n"
//" param.probability_estimates : whether to train a SVC or SVR model for probability estimates, 0 or 1 (default 0)\n"
//" param.weight : set the parameter C of class i to weight*C, for C-SVC (default 1)\n"
//go to ./external/libsvm/include/svm.h for for more details
typedef struct _svm_Classifier{
    struct svm_parameter param;		// set by parse_command_line
    struct svm_model *model;
    struct svm_problem prob;		// set by read_problem
    struct svm_node *x_space;
    bool crossValidation;
    int nfolds_crossValidation;
    bool useSparseDataRepresentation;
    bool isClassification;
    bool predict_probability;
}SVM_Classifier;



//knn-classifier
Knn_Classifier* createKnnClassifier();
Knn_Classifier* createKnnClassifier(Matrix* dataMatrixTrain,GVector* labelsTrain,
                                    DistanceFunction distanceFunction = computeNormalizedL1Norm,
                                    ArgumentList* argumentListDistanceFunction = NULL);
void destroyKnnClassifier(Knn_Classifier** pknnClassifier);
void destroyKnnClassifierForVoidPointer(void* pknnClassifier);
void knnClassifierFit(Knn_Classifier* knnClassifier, Matrix* dataMatrixTrain,GVector* labelsTrain);
void knn_Classifier_fit(Matrix* matrix, GVector* vector,void* knnClassifier);
GVector* knnClassifierPredict(Matrix* matrix,Knn_Classifier* knnClassifier);
GVector* knn_Classifier_predict(Matrix* matrix,void* knnClassifier);

//kmeans-classifier
Kmeans_Classifier* createKmeansClassifier();
void kmeansClassifierFit(Kmeans_Classifier* kmeansClassifier, Matrix* dataMatrixTrain,GVector* labelsTrain);
void kmeans_Classifier_fit(Matrix* matrix, GVector* vector,void* kmeansClassifier);
GVector* kmeansClassifierPredict(Matrix* matrix,Kmeans_Classifier* kmeansClassifier);
GVector* kmeans_Classifier_predict(Matrix* matrix,void* kmeansClassifier);
void destroyKmeansClassifier(Kmeans_Classifier** pkmeansClassifier);
void destroyKmeansClassifierForVoidPointer(void* pkmeansClassifier);


//SVM-Classifier
SVM_Classifier* createSVMClassifier();
void destroySVMClassifier(SVM_Classifier** svmClassifier);
void destroySVMClassifierForVoidPointer(void* psvmClassifier);
void setSVMDefaultParameters(SVM_Classifier* svmClassifier);
void svmClassifierFit(SVM_Classifier* svmClassifier,double* X, size_t nrows,size_t ncols, double* Y);

int convertOriginalDataFormat2SvmFormatTrain_usingSparseRepresentation(double* X, size_t nrows,size_t ncols, double* Y,
                                                                   struct svm_problem* prob,
                                                                   struct svm_node **x_space);

int convertOriginalDataFormat2SvmFormatTrain(double* X, size_t nrows,size_t ncols, double* Y,
                                         struct svm_problem* prob,
                                         struct svm_node **x_space);

double* convertOriginalDataFormat2SvmFormatPredict(SVM_Classifier* svmClassifier,
                                                   double* X, size_t nrows,size_t ncols);

double* convertOriginalDataFormat2SvmFormatPredict_usingSparseRepresentation(SVM_Classifier* svmClassifier,
                                                   double* X, size_t nrows,size_t ncols,
                                                   double** prob_estimates);


double* svmClassifierPredict(SVM_Classifier* svmClassifier,double* X, size_t nrows,size_t ncols);
double* predictLabelsLineByLine_SVM(SVM_Classifier* svmClassifier,
                                    double* X, size_t nrows,size_t ncols);
double* predictLabelsLineByLineUsingSparseRepresentation_SVM(SVM_Classifier* svmClassifier,
                                                             double* X, size_t nrows,size_t ncols);

void svm_Classifier_fit(Matrix* matrix, GVector* vector,void* svmClassifier);
GVector* svm_Classifier_predict(Matrix* matrix,void* svmClassifier);
void do_cross_validation_SVM(SVM_Classifier* svmClassifier);
#endif //LIBFL_CLASSIFIERS_H
