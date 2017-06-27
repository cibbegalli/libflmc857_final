#include "FL.h"

int main(int argc, char **argv) {

    AdjacencyRelation *adjRel = NULL;
    Image *grayImage = NULL;
    Image *outputImage = NULL;
    char *pos = NULL;
    char fileExtension[10];

    if (argc!=5){
        printf("Operation types:\n");
        printf("0: Dilation\n");
        printf("1: Erosion\n");
        printf("2: Opening\n");
        printf("3: Closing\n");
        printf("4: Top Hat \n");
        printf("5: Bottom Hat \n");
        printf("6: Morphologic Gradient\n");
        printf("Usage: <image.[pgm]> <operation type> <adj. radius> <output_filename>.[pgm]\n");
        //adj. radius tipically  1  to sqrt(8)
        exit(-1);
    }

    pos = strrchr(argv[1],'.') + 1;
    sscanf(pos,"%s",fileExtension);

    if (strcmp(fileExtension,"pgm")==0){
        grayImage = readImage(argv[1]);
        adjRel = createCircularAdjacency(atof(argv[3]));

        switch (atoi(argv[2])){
            case 0:
                //printf("%d\n",adjRel->n);
                outputImage = dilate(grayImage,adjRel);
                break;
            case 1:
                outputImage = erode(grayImage,adjRel);
                break;
            case 2:
                outputImage = open(grayImage,adjRel);
                break;
            case 3:
                outputImage = close(grayImage,adjRel);
                break;
            case 4:
                outputImage = topHat(grayImage,adjRel);
                break;
            case 5:
                outputImage = bottomHat(grayImage,adjRel);
                break;
            case 6:
                outputImage = morphologicGradient(grayImage,adjRel);
                break;
            default:
                printf("Invalid Operation type:\n");
                exit(-1);
        }
        writeImage(outputImage,argv[4]);
        destroyImage(&grayImage);
        destroyImage(&outputImage);
        destroyAdjacencyRelation(&adjRel);



    } else { /* ppm */
        printf("Invalide image format\n");
    }

    //free(pos);
    return 0;
}

