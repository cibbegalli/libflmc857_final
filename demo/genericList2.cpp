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
    LinkedList *list = createLinkedList(sizeof(float));
    for(int i = 1; i <= 25; i++) {
        float aux = (i*i)/100.0f;
        LIST_PUSH_BACK(float,list,aux);
        //appendElementInList(list, &i); //same as LIST_PUSH_BACK
    }

    ListIterator* iterator = getListIteratorBegin(list);

    while(hasNextInListIterator(iterator)){//only for non-circular lists
        printf("%f ", *((float*)getNextValueInListIterator(iterator)));
        //printf("%f\n",   LISTITERATOR_GET_NEXT_AS(float,iterator));
    }
    printf("\n");
    iterator = getListIteratorEnd(list);
    for (size_t j = 0; j < list->length; ++j) { //may works in both, circular and non-circluar lists
        printf("%f ", LISTITERATOR_GET_PREVIOUS_AS(float,iterator));
    }
    printf("\n");


    destroyLinkedList(&list);

    return 0;
}




