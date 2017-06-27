#include "FL.h"
#include <time.h>
#include <iostream>
#include <vector>

void destroyFloatVector(void* data){
    float* aux = *( (float**) data);
    free(aux);
    aux = NULL;
}

int main(int argc, char **argv) {
    int number = 40;
    float anotherNumber = 1.5;
    Image* image = createImage(100,100);
    ArgumentListNode* arg0 = createArgumentNode(&number,sizeof(int),"inteiro",NULL);
    ArgumentListNode* arg1 = createArgumentNode(&anotherNumber,sizeof(float),"float_teste",NULL);
    ArgumentListNode* arg2 = createArgumentNode(&image,sizeof(Image*),"imagem_buga_buga",NULL);
    ArgumentList* arglist = createArgumentList();
    pushBackElementInArgumentList(arglist,arg0);
    pushBackElementInArgumentList(arglist,arg1);
    pushBackElementInArgumentList(arglist,arg2);

    int arg0Value = ARGLIST_GET_ELEMENT_AS(int,arglist,0);
    float arg1Value = ARGLIST_GET_ELEMENT_AS(float,arglist,1);
    Image* arg2Value = ARGLIST_GET_ELEMENT_AS(Image*,arglist,2);

    printf("arg0: %d\n",arg0Value);
    printf("arg1: %f\n",arg1Value);
    printf("arg2: %d %d\n",arg2Value->nx,arg2Value->ny);
    printf("\n");

    //another way to do the same
    ArgumentList* arglist2 = createArgumentList();
    pushBackElementInArgumentList(arglist2,&number, sizeof(int));
    pushBackElementInArgumentList(arglist2,&anotherNumber, sizeof(int));
    pushBackElementInArgumentList(arglist2,&image, sizeof(Image*));
    arg0Value = ARGLIST_GET_ELEMENT_AS(int,arglist,0);
    arg1Value = ARGLIST_GET_ELEMENT_AS(float,arglist,1);
    arg2Value = ARGLIST_GET_ELEMENT_AS(Image*,arglist,2);

    printf("arg0: %d\n",arg0Value);
    printf("arg1: %f\n",arg1Value);
    printf("arg2: %d %d\n",arg2Value->nx,arg2Value->ny);
    printf("\n");

    //another way to do the same
    ArgumentListNode* arg0_1 = createRawArgumentNode(sizeof(int));
    ArgumentListNode* arg1_1 = createRawArgumentNode(sizeof(float));
    ArgumentListNode* arg2_1 = createRawArgumentNode(sizeof(Image*));
    ArgumentList* arglist3 = createArgumentList();
    pushBackElementInArgumentList(arglist3,arg0_1);
    pushBackElementInArgumentList(arglist3,arg1_1);
    pushBackElementInArgumentList(arglist3,arg2_1);
    ARGLIST_GET_ELEMENT_AS(int,arglist3,0) = number;
    ARGLIST_GET_ELEMENT_AS(float,arglist3,1) = anotherNumber;
    ARGLIST_GET_ELEMENT_AS(Image*,arglist3,2) = image;

    printf("arg0: %d\n", ARGLIST_GET_ELEMENT_AS(int,arglist3,0) );
    printf("arg1: %f\n", ARGLIST_GET_ELEMENT_AS(float,arglist3,1));
    printf("arg2: %d %d\n", (ARGLIST_GET_ELEMENT_AS(Image*,arglist3,2))->nx,(ARGLIST_GET_ELEMENT_AS(Image*,arglist,2))->ny );
    printf("\n");

    ArgumentList* arglist4 = createArgumentList();
    ARGLIST_PUSH_BACK_AS(int,arglist4,number);
    ARGLIST_PUSH_BACK_AS(float,arglist4,anotherNumber);
    ARGLIST_PUSH_BACK_AS(Image*,arglist4,image);
    printf("arg0: %d\n", ARGLIST_GET_ELEMENT_AS(int,arglist4,0) );
    printf("arg1: %f\n", ARGLIST_GET_ELEMENT_AS(float,arglist4,1));
    printf("arg2: %d %d\n", (ARGLIST_GET_ELEMENT_AS(Image*,arglist4,2))->nx,(ARGLIST_GET_ELEMENT_AS(Image*,arglist,2))->ny );
    printf("\n");

    destroyArgumentList(&arglist);
    destroyArgumentList(&arglist2);
    destroyArgumentList(&arglist3);
    return 0;
}




