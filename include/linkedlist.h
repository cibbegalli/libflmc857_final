//
// Created by deangeli on 5/1/17.
//

#ifndef _LINKEDLIST_H_
#define _LINKEDLIST_H_

#include "common.h"
#include "dataTransference.h"

typedef void (*FreeFunction)(void *);
typedef bool (*ComparatorFunction)(void *,void *);


typedef struct _LinkedListNode {
    void *data;
    struct _LinkedListNode *next;
    struct _LinkedListNode *previous;
} LinkedListNode;

typedef struct _linkedList{
    size_t length;
    size_t dataSize;
    LinkedListNode *head;
    LinkedListNode *tail;
    LinkedListNode *currentNode;
    int currentNodeIndex;
    bool isCircular;
    FreeFunction freeFunction;
    ComparatorFunction comparatorFunction;
} LinkedList;

typedef struct _listIterator {
    LinkedListNode *node;
} ListIterator;

inline void  setCursorPosition(LinkedList* list, int index,_LinkedListNode* node){
    list->currentNodeIndex = index;
    list->currentNode = node;
}



//LinkedList* createLinkedList(size_t dataSize);
LinkedList* createLinkedList(size_t dataSize, bool isCircular = false, FreeFunction freeFunction = NULL);
inline void appendElementInList(LinkedList *list, void *element){
    LinkedListNode *node = (LinkedListNode*)calloc(1,sizeof(LinkedListNode));
    node->data = calloc(1,list->dataSize);
    node->next = NULL;

    TRANSFER_DATA_COPY(node->data, element, list->dataSize);
    //memcpy(node->data, element, list->dataSize);

    if(list->length == 0) {
        list->head = list->tail = node;
    } else {
        list->tail->next = node;
        node->previous = list->tail;
        list->tail = node;
    }
    if(list->isCircular){
        list->tail->next = list->head;
        list->head->previous = list->tail;
    }

    list->length++;
}
#define LIST_PUSH_BACK(type,list,element) \
    LinkedListNode *node = (LinkedListNode*)calloc(1,sizeof(LinkedListNode)); \
    node->data = calloc(1,list->dataSize); \
    node->next = NULL; \
    ((type*)node->data)[0] = element; \
    if(list->length == 0) {\
        list->head = list->tail = node;\
    }\
    else {\
        list->tail->next = node;\
        node->previous = list->tail;\
        list->tail = node;\
    }\
    if(list->isCircular){\
        list->tail->next = list->head; \
        list->head->previous = list->tail;\
    }\
    list->length++;

#define LIST_GET_ELEMENT_AS(type,list,index) (*((type*)getElementInList(list,index)))

inline void pushBackElementInList(LinkedList *list, void *element){
    appendElementInList(list,element);
}
inline void prependElementInList(LinkedList *list, void *element){
    LinkedListNode *node = (LinkedListNode*)malloc(sizeof(LinkedListNode));
    node->data = malloc(list->dataSize);
    TRANSFER_DATA_COPY(node->data, element, list->dataSize);

    if(list->length == 0) {
        list->head = node;
        list->tail = list->head;
    }else{
        node->next = list->head;
        list->head->previous = node;
        list->head = node;
    }
    if(list->isCircular){
        list->tail->next = list->head;
        list->head->previous = list->tail;
    }
    list->length++;
}

inline void pushFrontElementInList(LinkedList *list, void *element){
    prependElementInList(list,element);
}
void insertElementInListAt(LinkedList *list, void *element, size_t index);
void removeListHead(LinkedList *list);
void removeListTail(LinkedList *list);
void removeElementInListAt(LinkedList *list, size_t index);
void removeElementInListByReference(LinkedList *list, void *element);
void removeElementInListGivenValue(LinkedList *list, void *element);
void removeElementsInListGivenValue(LinkedList *list, void *element);
void resetIterator(LinkedList *list);
void* getNextElement(LinkedList *list);
void* getPreviousElement(LinkedList *list);
inline void* getElementInList(LinkedList *list, size_t index){
    if(index >= list->length){
        if(list->isCircular){
            index = index % list->length;
        }else{
            printf("[getElement] invalid position %lu. The list length is %lu (indexing start from 0)\n", index,list->length);
            return NULL;
        }
    }

    if (index == 0){
        return list->head->data;
    }
    if(index == list->length-1){
        return list->tail->data;
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
                return currentNode->data;
            }else{
                currentNode = currentNode->next;
            }

        }
    }else if(distance2Tail <= distance2Current) {//tail 2 element
        LinkedListNode *currentNode = list->tail;
        for (int i = list->length-1; i >= 0; --i) {
            if(i == (int)index){
                setCursorPosition(list,i,currentNode);
                return currentNode->data;
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
                    return currentNode->data;
                }else{
                    currentNode = currentNode->previous;
                }
            }
        }else{//element is front
            LinkedListNode *currentNode = list->currentNode;
            for (size_t i = list->currentNodeIndex; i < list->length; ++i) {
                if(i == index){
                    setCursorPosition(list,i,currentNode);
                    return currentNode->data;
                }else{
                    currentNode = currentNode->next;
                }
            }
        }
    }

    //unlikely hit this bit
    return NULL;
}



LinkedListNode* getLinkedListNode(LinkedList *list, size_t index);
void destroyLinkedList(LinkedList **list);
void destroyNodeList(LinkedList* list,LinkedListNode **node);
void clearLinkedList(LinkedList* list);

//ITERATOR
ListIterator* getListIteratorBegin(LinkedList *list);
ListIterator* getListIteratorEnd(LinkedList *list);

#define LISTITERATOR_GET_NEXT_AS(type,iterator) (*((type*)getNextValueInListIterator(iterator)))
#define LISTITERATOR_GET_PREVIOUS_AS(type,iterator) (*((type*)getPreviousValueInListIterator(iterator)))
#define LISTITERATOR_GET_CURRENT_AS(type,iterator) (*((type*)getValueInListIterator(iterator)))

inline void* getValueInListIterator(ListIterator* iterator){
    return iterator->node->data;
}

inline void* getNextValueInListIterator(ListIterator* iterator){
    void* data = iterator->node->data;
    iterator->node = iterator->node->next;
    return data;
}

inline void* getPreviousValueInListIterator(ListIterator* iterator){
    void* data = iterator->node->data;
    iterator->node = iterator->node->previous;
    return data;
}

inline void* hasNextInListIterator(ListIterator* iterator){
    return iterator->node;
}



size_t getIteratorIndexInList(LinkedList* list, ListIterator* iterator);


#endif //_LIST_H_
