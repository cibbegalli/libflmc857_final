//
// Created by deangeli on 4/30/17.
//

#include "linkedlist.h"

LinkedList* createLinkedList(size_t dataSize,  bool isCircular, FreeFunction freeFunction){
    if(dataSize == 0){
        printf("[createList] data size is 0\n");
        return NULL;
    }
    LinkedList* list = (LinkedList*)calloc(1,sizeof(LinkedList));
    list->dataSize = dataSize;
    list->head = list->tail = NULL;
    list->freeFunction = freeFunction;
    list->length = 0;
    list->currentNode = NULL;
    list->currentNodeIndex = -1;
    list->isCircular = isCircular;
    list->comparatorFunction = NULL;
    return list;
}

void removeListHead(LinkedList *list){
    if(list->head == NULL){
        return;
    }
    LinkedListNode *currentNode = NULL;
    currentNode = list->head;
    list->head = currentNode->next;
    destroyNodeList(list,&currentNode);
    setCursorPosition(list,0,list->head);
    if(list->isCircular){
        list->tail->next = list->head;
        list->head->previous = list->tail;
    }
}

void removeListTail(LinkedList *list){
    if(list->head == NULL){
        return;
    }

    LinkedListNode *previousNode = list->tail->previous;
    LinkedListNode *currentTail = list->tail;
    destroyNodeList(list,&currentTail);
    list->tail = previousNode;
    setCursorPosition(list,list->length,list->tail);
    //setCursorPosition(list,-1,NULL);
    if(list->isCircular){
        list->tail->next = list->head;
        list->head->previous = list->tail;
    }

    return;
}



void insertElementInListAt(LinkedList *list, void *element, size_t index){
    if(index > list->length){
        printf("[insertElementInListAt] invalid position %lu. The list length is %lu (indexing start from 0)\n", index,list->length);
        return;
    }

    if (index == 0){
        prependElementInList(list, element);
        return;
    }
    if(index == list->length){
        appendElementInList(list, element);
        return;
    }

    LinkedListNode *node = (LinkedListNode*)malloc(sizeof(LinkedListNode));
    node->data = malloc(list->dataSize);
    //memcpy(node->data, element, list->dataSize);
    TRANSFER_DATA_COPY(node->data, element, list->dataSize);

    int distance2Head = index;
    int distance2Tail = list->length -index;
    int distance2Current = index - list->currentNodeIndex;
    int currentDirection = 0; //foward
    if(distance2Current <= 0){
        currentDirection = 1;//backward
        distance2Current = -distance2Current;
    }

    if(distance2Head <= distance2Tail && distance2Head <= distance2Current){//head 2 element
        LinkedListNode *currentNode = list->head;
        for (size_t i = 0; i < list->length; ++i) {
            if(i == index){
                node->previous = currentNode->previous;
                node->next = currentNode;
                currentNode->previous->next = node;
                currentNode->previous = node;
                list->length++;
                setCursorPosition(list,i,node);
                return;
            }else{
                currentNode = currentNode->next;
            }

        }
    }else if(distance2Tail <= distance2Current){//tail 2 element
        LinkedListNode *currentNode = list->tail;
        for (int i = list->length-1; i >= 0; --i) {
            if(i == (int)index){
                node->previous = currentNode->previous;
                node->next = currentNode;
                currentNode->previous->next = node;
                currentNode->previous = node;
                list->length++;
                setCursorPosition(list,i,node);
                return;
            }else{
                currentNode = currentNode->previous;
            }
        }
    }else{//current 2 element
        if(currentDirection){//element is back
            LinkedListNode *currentNode = list->currentNode;
            for (int i = list->currentNodeIndex; i >= 0; --i) {
                if(i == (int)index){
                    node->previous = currentNode->previous;
                    node->next = currentNode;
                    currentNode->previous->next = node;
                    currentNode->previous = node;
                    list->length++;
                    setCursorPosition(list,i,node);
                    return;
                }else{
                    currentNode = currentNode->previous;
                }
            }
        }else{//element is front
            LinkedListNode *currentNode = list->currentNode;
            for (size_t i = list->currentNodeIndex; i < list->length; ++i) {
                if(i == index){
                    node->previous = currentNode->previous;
                    node->next = currentNode;
                    currentNode->previous->next = node;
                    currentNode->previous = node;
                    list->length++;
                    setCursorPosition(list,i,node);
                    return;
                }else{
                    currentNode = currentNode->next;
                }
            }
        }

    }
}

void removeElementInListAt(LinkedList *list, size_t index){
    if(index >= list->length){
        printf("[insertElementInListAt] invalid position %lu. The list length is %lu (indexing start from 0)\n", index,list->length);
        return;
    }

    if (index == 0){
        removeListHead(list);
        return;
    }
    if(index == list->length-1){
        removeListTail(list);
        return;
    }

    int distance2Head = index;
    int distance2Tail = list->length -index;
    int distance2Current = index - list->currentNodeIndex;
    int currentDirection = 0; //foward
    if(distance2Current <= 0){
        currentDirection = 1;//backward
        distance2Current = -distance2Current;
    }
    if(distance2Head <= distance2Tail && distance2Head <= distance2Current){//head 2 element
        LinkedListNode *currentNode = list->head;
        for (size_t i = 0; i < list->length; ++i) {
            if(i == index){
                currentNode->previous->next = currentNode->next;
                currentNode->next->previous = currentNode->previous;
                setCursorPosition(list,i,currentNode->next);
                destroyNodeList(list,&currentNode);//already decreasing list size
                return;
            }else{
                currentNode = currentNode->next;
            }

        }
    }else if(distance2Tail <= distance2Current) {//tail 2 element
        LinkedListNode *currentNode = list->tail;
        for (int i = list->length-1; i >= 0; --i) {
            if(i == (int)index){
                currentNode->previous->next = currentNode->next;
                currentNode->next->previous = currentNode->previous;
                //setCursorPosition(list,-1,NULL);
                destroyNodeList(list,&currentNode);//already decreasing list size
                return;
            }else{
                currentNode = currentNode->previous;
            }
        }
    }else{//current 2 element
        if(currentDirection){//element is back
            LinkedListNode *currentNode = list->currentNode;
            for (int i = list->currentNodeIndex; i >= 0; --i) {
                if(i == (int)index){
                    currentNode->previous->next = currentNode->next;
                    currentNode->next->previous = currentNode->previous;
                    setCursorPosition(list,i,currentNode->next);
                    destroyNodeList(list,&currentNode);//already decreasing list size
                    return;
                }else{
                    currentNode = currentNode->previous;
                }
            }
        }else{//element is front
            LinkedListNode *currentNode = list->currentNode;
            for (size_t i = list->currentNodeIndex; i < list->length; ++i) {
                if(i == index){
                    currentNode->previous->next = currentNode->next;
                    currentNode->next->previous = currentNode->previous;
                    setCursorPosition(list,i,currentNode->next);
                    destroyNodeList(list,&currentNode);//already decreasing list size
                    return;
                }else{
                    currentNode = currentNode->next;
                }
            }
        }
    }

}

void removeElementInListByReference(LinkedList *list, void *element){
    if(list->head->data == element){
        removeListHead(list);
        return;
    }
    if(list->tail->data == element){
        removeListTail(list);
        return;
    }

    LinkedListNode *currentNode = list->head;
    for (size_t i = 0; i < list->length; ++i) {
        if(currentNode->data == element){
            currentNode->previous->next = currentNode->next;
            currentNode->next->previous = currentNode->previous;
            setCursorPosition(list,i,currentNode->next);
            destroyNodeList(list,&currentNode);//already decreasing list size
            return;
        }else{
            currentNode = currentNode->next;
        }
    }
    printf("[removeElement] could not find the pointer in the list\n");
}

void removeElementInListGivenValue(LinkedList *list, void *element){
    if(!list->comparatorFunction){
        printf("[removeElementInListGivenValue] compare function not defined\n");
    }

    if(list->head->data == element){
        removeListHead(list);
        return;
    }
    if(list->tail->data == element){
        removeListTail(list);
        return;
    }

    LinkedListNode *currentNode = list->head;
    for (size_t i = 0; i < list->length; ++i) {
        if( list->comparatorFunction(currentNode->data,element)){
            currentNode->previous->next = currentNode->next;
            currentNode->next->previous = currentNode->previous;
            setCursorPosition(list,i,currentNode->next);
            destroyNodeList(list,&currentNode);//already decreasing list size
            return;
        }else{
            currentNode = currentNode->next;
        }
    }
}

void removeElementsInListGivenValue(LinkedList *list, void *element){
    if(!list->comparatorFunction){
        printf("[removeElementInListGivenValue] compare function not defined\n");
    }

    if(list->head->data == element){
        removeListHead(list);
        return;
    }
    if(list->tail->data == element){
        removeListTail(list);
        return;
    }

    LinkedListNode *currentNode = list->head;
    for (size_t i = 0; i < list->length; ++i) {
        if( list->comparatorFunction(currentNode->data,element)){
            currentNode->previous->next = currentNode->next;
            currentNode->next->previous = currentNode->previous;
            setCursorPosition(list,i,currentNode->next);
            destroyNodeList(list,&currentNode);//already decreasing list size
        }else{
            currentNode = currentNode->next;
        }
    }
}

LinkedListNode* getLinkedListNode(LinkedList *list, size_t  index){
    if(index >= list->length){
        if(list->isCircular){
            index = index % list->length;
        }else{
            printf("[getElement] invalid position %lu. The list length is %lu (indexing start from 0)\n", index,list->length);
            return NULL;
        }
    }

    if (index == 0){
        return list->head;
    }
    if(index == list->length-1){
        return list->tail;
    }

    int distance2Head = index;
    int distance2Tail = list->length -index;
    int distance2Current = index - list->currentNodeIndex;
    int currentDirection = 0; //foward
    if(distance2Current <= 0){
        currentDirection = 1;//backward
        distance2Current = -distance2Current;
    }

    if(distance2Head <= distance2Tail && distance2Head <= distance2Current){//head 2 element
        LinkedListNode *currentNode = list->head;
        for (size_t i = 0; i < list->length; ++i) {
            if(i == index){
                setCursorPosition(list,i,currentNode);
                return currentNode;
            }else{
                currentNode = currentNode->next;
            }

        }
    }else if(distance2Tail <= distance2Current) {//tail 2 element
        LinkedListNode *currentNode = list->tail;
        for (int i = list->length-1; i >= 0; --i) {
            if(i == (int)index){
                setCursorPosition(list,i,currentNode);
                return currentNode;
            }else{
                currentNode = currentNode->previous;
            }
        }
    }else{//current 2 element
        if(currentDirection){//element is back
            LinkedListNode *currentNode = list->currentNode;
            for (int i = list->currentNodeIndex; i >= 0; --i) {
                if(i == (int)index){
                    setCursorPosition(list,i,currentNode);
                    return currentNode;
                }else{
                    currentNode = currentNode->previous;
                }
            }
        }else{//element is front
            LinkedListNode *currentNode = list->currentNode;
            for (size_t i = list->currentNodeIndex; i < list->length; ++i) {
                if(i == index){
                    setCursorPosition(list,i,currentNode);
                    return currentNode;
                }else{
                    currentNode = currentNode->next;
                }
            }
        }
    }

    //unlikely hit this bit
    return NULL;
}

void destroyNodeList(LinkedList* list,LinkedListNode **node){
    LinkedListNode *aux = *node;
    if(aux == NULL){
        return;
    }
    if(list->freeFunction) {
        list->freeFunction(aux->data);
        aux->data = NULL;
    }
    free(aux);
    aux = NULL;
    list->length--;
}

void clearLinkedList(LinkedList* list){
    LinkedListNode *currentNode = NULL;
    if(list->isCircular){
        list->tail->next = NULL;
        list->head->previous = NULL;
    }
    while(list->head != NULL) {
        currentNode = list->head;
        list->head = currentNode->next;
        destroyNodeList(list,&currentNode);
    }
    setCursorPosition(list,-1,NULL);
}

void destroyLinkedList(LinkedList **list)
{

    LinkedList* aux = *list;
    LinkedListNode *currentNode = NULL;
    if(aux == NULL){
        return;
    }

    aux->tail->next = NULL;
    aux->head->previous = NULL;

    while(aux->head != NULL) {
        currentNode = aux->head;
        aux->head = currentNode->next;
        destroyNodeList(aux,&currentNode);
    }
    free(aux);
    aux = NULL;
}

ListIterator* getListIteratorBegin(LinkedList *list){
    ListIterator* iterator = (ListIterator*)calloc(1,sizeof(ListIterator));
    iterator->node = list->head;
    return iterator;
}

ListIterator* getListIteratorEnd(LinkedList *list){
    ListIterator* iterator = (ListIterator*)calloc(1,sizeof(ListIterator));
    iterator->node = list->tail;
    return iterator;
}

size_t getIteratorIndexInList(LinkedList* list, ListIterator* iterator){
    LinkedListNode* node = list->head;
    for (size_t i = 0; i < list->length; ++i) {
        if(iterator->node == node){
            return i;
        }else{
            node = node->next;
        }
    }
    return -1;
}
