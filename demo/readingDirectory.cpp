#include "FL.h"

int main(int argc, char **argv) {


//    if (argc != 2){
//        printf("Usage: <directory path [such as ../data/]> \n");
//        exit(-1);
//    }
    //list all files in directory
    DirectoryManager* directoryManager = loadDirectory("../data/", 1);
    printf("reading:\n");
    for (int i = 0; i < (int)directoryManager->nfiles; ++i) {
        printf("%s\n",directoryManager->files[i]->path);
    }

    printf("\n\n");
    //list all files in directory with extension .ppm
    DirectoryManager* directoryManager2 = loadFilesFromDirBySuffix(argv[1], "ppm");
    printf("reading:\n");
    for (int i = 0; i < (int)directoryManager2->nfiles; ++i) {
        printf("%s\n",directoryManager2->files[i]->path);
    }

    destroyDirectoryManager(&directoryManager);
    destroyDirectoryManager(&directoryManager2);

    return 0;
}
