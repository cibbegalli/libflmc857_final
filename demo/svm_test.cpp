#include "FL.h"

int main(int argc, char **argv) {

    SVM_Classifier* svm = createSVMClassifier();
    int nrows = 3;
    int ncols = 4;

    double *data = (double*)calloc(nrows*ncols,sizeof(double));
    double *labels = (double*)calloc(nrows,sizeof(double));

    int k = 0;
    for (int i = 0; i < nrows; ++i) {
        labels[i] = (i%3) + 1;
        printf("%f",labels[i]);
        for (int j = 0; j < ncols; ++j) {
            data[k] = (double)((i*ncols)+j)*10;
            printf(" %d:%f",j+1,data[k]);
            k++;
        }
        printf("\n");
    }
    svmClassifierFit(svm,data,nrows,ncols,labels);
    double* labelsPredicted = svmClassifierPredict(svm,data,nrows,ncols);
    for (int row = 0; row < nrows; ++row) {
        printf("%f\n",labelsPredicted[row]);
    }

    destroySVMClassifier(&svm);
    free(labelsPredicted);
    free(data);
    free(labels);
    return 0;
}


