#include "FL.h"

bool compareIntegers(void *value1,void *value2)
{
    int n1 = (*(int*)value1);
    int n2 = (*(int*)value2);
    return (n1 == n2);
}

void destroyImage(void *image){
    Image** aux = ((Image**)image);
    destroyImage(aux);
}


int main(int argc, char const *argv[]) {

    //generating numbers
    LinkedList *list = createLinkedList(sizeof(int),true,NULL);
    for(int i = 1; i <= 10; i++) {
        appendElementInList(list, &i);
    }
    //iterate over the list nodes (foward)
    for (unsigned  int i = 0; i < list->length*2; ++i) {
        int a = *((int*)getElementInList(list,i));
        printf("%d ",a);
    }
    printf("\n");
    //iterate over the list nodes (backwards)
    for ( int i = list->length*2 - 1; i >= 0; --i) {
        int a = *((int*)getElementInList(list,i));
        printf("%d ",a);
    }
    printf("\n");

    removeListHead(list);
    for (unsigned  int i = 0; i < list->length*2; ++i) {
        int a = *((int*)getElementInList(list,i));
        printf("%d ",a);
    }
    printf("\n");
//    resetIterator(list);
    removeListTail(list);
    for (unsigned  int i = 0; i < list->length*2; ++i) {
        int a = *((int*)getElementInList(list,i));
        printf("%d ",a);
    }
    printf("\n");
    int n1 = 30;
    int n2 = 60;
    insertElementInListAt(list, &n1, 3);
    insertElementInListAt(list, &n2, 6);
    resetIterator(list);
    for (unsigned  int i = 0; i < list->length*2; ++i) {
        int a = *((int*)getElementInList(list,i));
        printf("%d ",a);
    }
    printf("\n");
    removeElementInListAt(list, 6);
    removeElementInListAt(list, 3);
    for (unsigned  int i = 0; i < list->length*2; ++i) {
        int a = LIST_GET_ELEMENT_AS(int,list,i);
        printf("%d ",a);
    }
    printf("\n");
    int n = *((int*) getElementInList(list, 5));
    printf("%d\n",n);
    removeElementInListByReference(list, &n);
    for (unsigned  int i = 0; i < list->length*2; ++i) {
        int a = LIST_GET_ELEMENT_AS(int,list,i);
        printf("%d ",a);
    }
    printf("\n");
    list->comparatorFunction = compareIntegers;
    removeElementInListGivenValue(list, &n);
    for (unsigned  int i = 0; i < list->length*2; ++i) {
        int a = *((int*)getElementInList(list,i));
        printf("%d ",a);
    }
    printf("\n");




    Image* image1 = readImage("../data/lena.pgm");
    Image* image2 = readImage("../data/lena.pgm");
    Image* image3 = readImage("../data/lena.pgm");

    LinkedList *list_image = createLinkedList(sizeof(Image*),false,destroyImage);
    appendElementInList(list_image, &image1);
    appendElementInList(list_image, &image2);
    appendElementInList(list_image, &image3);


    Image* aux = *((Image**)getNextElement(list_image));
    writeImage(aux,"img1.pgm");

    aux = *((Image**)getNextElement(list_image));
    addUniformNoise(aux,30,0.05);
    writeImage(aux,"img2.pgm");

    aux = *((Image**)getNextElement(list_image));
    addSaltAndPepperNoise(aux,0.05);
    writeImage(aux,"img3.pgm");

    destroyLinkedList(&list);
    destroyLinkedList(&list_image);


    return 0;
}




