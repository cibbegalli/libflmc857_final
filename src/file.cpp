#include "file.h"

DirectoryManager *loadDirectory(const char *dir_pathname, size_t hier_levels) {
    DirectoryManager *dir = NULL;


    if (pathnameExists(dir_pathname)) {
        // it is really a directory and it exists
        if (directoryExists(dir_pathname)) {
            dir = (DirectoryManager*) calloc(1, sizeof(DirectoryManager));
            dir->path = (char*)calloc(strlen(dir_pathname) + 2, sizeof(char)); // one more char to put the separation '/'
            strcpy(dir->path, dir_pathname);

            // puts the '/' at the end of the pathname
            if (dir->path[strlen(dir->path) - 1] != SEP_C[0])
                strcat(dir->path, SEP_C);

            listDirectory(dir, hier_levels);
        }
            // it is a File instead of a Directory
        else {
            printf("Pathname \"%s\" is a File", dir_pathname);
        }
    }
    else {
        printf("Pathname \"%s\" does not exist!", dir_pathname);
    }

    return dir;
}

int findTrueLabelInName(char *pathname){
    char *posLastSlash = NULL;
    char *posLastUnderScore = NULL;
    //char *posLastPoint = NULL;
    char labelNumber[10];
    int k = 0;
    size_t nChars;
    posLastSlash = strrchr(pathname,'/') +1;
    posLastUnderScore = strrchr(pathname,'_');

    if(posLastSlash == NULL && posLastUnderScore == NULL){
        return -1;
    }if(posLastSlash == NULL){
        nChars = posLastUnderScore-pathname;
        for (size_t i = 0; i < nChars; ++i) {
            labelNumber[i] = pathname[i];
            k++;
        }
        labelNumber[k] = '\0';
        int label = atoi(labelNumber);
        return label;
    }

    nChars = posLastUnderScore - posLastSlash;
    for (size_t i = 0; i < nChars; ++i) {
        labelNumber[i] = posLastSlash[i];
        k++;
    }
    labelNumber[k] = '\0';
    int label = atoi(labelNumber);
    return label;
}

void findTrueLabelInCurrentDirectory(DirectoryManager *directoryManager){

    for (size_t j = 0; j < directoryManager->nfiles; ++j) {
        directoryManager->files[j]->label = findTrueLabelInName(directoryManager->files[j]->path);
        printf("%d\n",directoryManager->files[j]->label);
    }
}

bool pathnameExists(const char *pathname) {
    struct stat buffer;
    return (stat (pathname, &buffer) == 0);
}

bool directoryExists(const char *pathname) {
    struct stat st;

    if (stat(pathname, &st) == 0) {
        if (S_ISDIR(st.st_mode))
            return true; //it's a directory
    }
    return false ;
}

void listDirectory(DirectoryManager *root_dir, size_t hier_levels) {
    if (root_dir == NULL) {
        printf("Directory is NULL");
        return ;
    }
    else {
        if (hier_levels == 0)
            hier_levels = INT_MAX; // trick to set the hier_levels as the possible maximum

        root_dir->nfiles   = 0;
        root_dir->nsubdirs = 0;
        root_dir->files    = NULL;
        root_dir->subdirs  = NULL;

        size_t curr_level = 1;
        listDirectoryRecursively(root_dir, hier_levels, curr_level);
    }
}

char *joinPathnames(const char *pathname1, const char *pathname2){
    if (pathname1 == NULL) {
        printf("First Pathname is NULL");
        return NULL;
    }
    if (pathname2 == NULL){
        printf("Second Pathname is NULL");
        return NULL;
    }


    char *joined_pathname = (char*)calloc(strlen(pathname1)+strlen(pathname2)+2,sizeof(char)); // one char for '/' (separation char) and other for '\0'

    if (pathname2[0] == SEP_C[0])
        pathname2++; // skip the first char

    if (endsWith(pathname1, SEP_C))
        sprintf(joined_pathname, "%s%s", pathname1, pathname2);
    else
        sprintf(joined_pathname, "%s%s%s", pathname1, SEP_C, pathname2);

    return joined_pathname;
}

bool endsWith(const char *str, const char *suffix) {
    if (str == NULL) {
        printf("String is NULL\n");
        return false;
    }
    if (suffix == NULL){
        printf("Suffix is NULL\n");
        return false;
    }

    size_t len_suffix = strlen(suffix);
    size_t len_str    = strlen(str);

    if (len_suffix <= len_str) {
        size_t shift = len_str - len_suffix;
        return (strncmp(str+shift, suffix, len_suffix) == 0);
    }
    else
        return false;
}


void countFilesInDirectory(const char *dir_pathname, size_t *nfiles, size_t *nsubdirs) {
    //http://pubs.opengroup.org/onlinepubs/007908799/xsh/dirent.h.html
    //http://www.delorie.com/gnu/docs/glibc/libc_270.html
    DIR *system_dir;
    struct dirent *entry;
    char *pathname = NULL;
    *nfiles = 0;
    *nsubdirs = 0;

    system_dir = opendir(dir_pathname);
    if (system_dir) {
        while ((entry = readdir(system_dir)) != NULL)
            // it excludes the system_dir . and ..
            if ((strcmp(entry->d_name, ".") != 0) && (strcmp(entry->d_name, "..") != 0)) {
                pathname = joinPathnames(dir_pathname, entry->d_name);

                if (directoryExists(pathname)) {
                    (*nsubdirs)++;
                }
                else {
                    (*nfiles)++;
                }

                free(pathname);
                pathname = NULL;
            }
        closedir(system_dir);
    }
    else {
        printf("Error opening directory path: \"%s\"\n", dir_pathname);
    }

}

void listDirectoryRecursively(DirectoryManager *dir, size_t hier_levels, size_t curr_level) {
    DIR *system_dir;
    struct dirent *entry;
    char *pathname = NULL;

    dir->files   = NULL;
    dir->subdirs = NULL;

    if (curr_level <= hier_levels) {
        system_dir = opendir(dir->path);
        if (system_dir) {
            countFilesInDirectory(dir->path, &dir->nfiles, &dir->nsubdirs);
            if (dir->nfiles != 0)
                dir->files   = (FileManager**) calloc(dir->nfiles, sizeof(FileManager*));
            if (dir->nsubdirs != 0)
                dir->subdirs = (DirectoryManager**) calloc(dir->nsubdirs, sizeof(DirectoryManager*));

            size_t i = 0, j = 0;
            while ((entry = readdir(system_dir)) != NULL) {
                // it excludes the dir . and ..
                if ((strcmp(entry->d_name, ".") != 0) && (strcmp(entry->d_name, "..") != 0)) {
                    pathname = joinPathnames(dir->path, entry->d_name);

                    if (directoryExists(pathname)) { // it is a directory
                        DirectoryManager *subdir = (DirectoryManager*) calloc(1, sizeof(DirectoryManager));
                        subdir->path = pathname;

                        subdir->nfiles   = 0;
                        subdir->nsubdirs = 0;
                        subdir->files    = NULL;
                        subdir->subdirs  = NULL;

                        listDirectoryRecursively(subdir, hier_levels, curr_level+1);
                        dir->subdirs[j++] = subdir;
                    }
                    else { // it is a File
                        FileManager *f = (FileManager*) calloc(1, sizeof(FileManager));
                        f->path = pathname;
                        dir->files[i++] = f;
                        f->suffix = NULL;
                    }
                }
            }
            closedir(system_dir);

            /* sorts the pathnames using qsort functions */
            qsort(dir->files, dir->nfiles, sizeof(FileManager*), compareFiles);
            qsort(dir->subdirs, dir->nsubdirs, sizeof(DirectoryManager*), compareDirectories);
        }
        else {
            printf("Error opening directory path: \"%s\"\n", dir->path);
        }
    }
}

int compareFiles(const void *a, const void *b) {
    FileManager **f1 = (FileManager**) a;
    FileManager **f2 = (FileManager**) b;

    return strcmp((*f1)->path, (*f2)->path);
}

int compareDirectories(const void *a, const void *b) {
    FileManager **dir1 = (FileManager**) a;
    FileManager **dir2 = (FileManager**) b;

    return strcmp((*dir1)->path, (*dir2)->path);
}

DirectoryManager *loadFilesFromDirBySuffix(const char *dir_pathname, const char *suffix) {

    DirectoryManager *dir = NULL;

    if (pathnameExists(dir_pathname)) {
        // it is really a directory and it exists
        if (directoryExists(dir_pathname)) {
            dir = (DirectoryManager*) calloc(1, sizeof(DirectoryManager));
            dir->path = (char*)calloc(strlen(dir_pathname) + 2,sizeof(char)); // one more char to put the separation '/'
            strcpy(dir->path, dir_pathname);

            // puts the '/' at the end of the pathname
            if (dir->path[strlen(dir->path) - 1] != SEP_C[0])
                strcat(dir->path, SEP_C);

            listFilesFromDirectory(dir, suffix);
        }
        else { // it is really a directory and it exists
            printf("\"%s\" is a file!\n", dir->path);
        }
    }
    else {
        printf("Pathname \"%s\" does not exist!\n", dir_pathname);
    }

    return dir;
}

void listFilesFromDirectory(DirectoryManager *dir, const char *suffix) {
    DIR *system_dir;
    struct dirent *entry;
    char *pathname = NULL;

    if (dir != NULL) {
        dir->files    = NULL;
        dir->subdirs  = NULL;
        dir->nfiles   = 0;
        dir->nsubdirs = 0;

        system_dir = opendir(dir->path);
        if (system_dir) {
            dir->nfiles = countOnlyFilesInDirectory(dir->path, suffix);

            if (dir->nfiles != 0) {
                dir->files = (FileManager**) calloc(dir->nfiles, sizeof(FileManager*));

                size_t i = 0;
                while ((entry = readdir(system_dir)) != NULL) {
                    pathname = joinPathnames(dir->path, entry->d_name);

                    if (!directoryExists(pathname) && (endsWith(pathname, suffix))) {
                        FileManager *f = createFileManager(pathname);
                        dir->files[i++] = f;
                    }
                }
                closedir(system_dir);

                /* sorts the pathnames using qsort functions */
                sortDirectory(dir);
            }
        }
        else {
            printf("Error opening directory path: \"%s\"\n", dir->path);
        }
    }
}

int countOnlyFilesInDirectory(const char *dir_pathname, const char *extension) {
    DIR *dir;
    struct dirent *entry;

    size_t file_count = 0;

    dir = opendir(dir_pathname);
    if (dir) {
        while ((entry = readdir(dir)) != NULL)
            if (!directoryExists(entry->d_name) && (endsWith(entry->d_name, extension)))
                file_count++;
        closedir(dir);
    }
    else {
        printf("Error opening directory path: \"%s\"\n", dir_pathname);
    }

    return file_count;
}

FileManager *createFileManager(const char *pathname) {
    FileManager *f = NULL;

    // it is a Directory instead of a File
    if (directoryExists(pathname)) {
        printf("Pathname \"%s\" is a directory\n", pathname);
        return NULL;
    }
    else {
        f = (FileManager*) calloc(1, sizeof(FileManager));
        f->path = copyString(pathname);
        f->suffix = NULL;
        //updateDatasetFileInfo(f);
    }

    return f;
}

void sortDirectory(DirectoryManager *files) {
    qsort(files->files, files->nfiles, sizeof(FileManager*), compareFiles);
}

void destroyFileManager(FileManager **fileManager) {
    if (fileManager != NULL) {
        FileManager *f_aux = *fileManager;

        if (f_aux != NULL) {
            if (f_aux->path != NULL) {
                free(f_aux->path);
                f_aux->path = NULL;
            }
            if(f_aux->suffix != NULL) {
                free(f_aux->suffix);
                f_aux->suffix = NULL;
            }
            free(f_aux);
            *fileManager = NULL;
        }
    }
}

void destroyDirectoryManager(DirectoryManager **directoryManager) {
    if (*directoryManager != NULL) {
        DirectoryManager *dir_aux = *directoryManager;

        if (dir_aux != NULL) {
            if (dir_aux->path != NULL)
                free(dir_aux->path);

            // deallocates the files
            if (dir_aux->files != NULL) {
                for (size_t i = 0; i < dir_aux->nfiles; i++)
                    destroyFileManager(&dir_aux->files[i]);

                free(dir_aux->files);
                dir_aux->files = NULL;
            }

            if (dir_aux->subdirs != NULL) {
                // deallocates the subdirs
                for (size_t j = 0; j < dir_aux->nsubdirs; j++)
                    destroyDirectoryManager(&dir_aux->subdirs[j]);

                free(dir_aux->subdirs);
                dir_aux->subdirs = NULL;
            }
            free(dir_aux);
            *directoryManager = NULL;
        }
    }
}

void destroyFileSet(FileSet **farr) {
    if (farr != NULL) {
        FileSet *faux = *farr;

        if (faux != NULL) {
            if (faux->files != NULL) {
                for (size_t i = 0; i < faux->n; i++)
                    destroyFileManager(&(faux->files[i]));
                free(faux->files);
            }
            free(faux);
            *farr = NULL;
        }
    }
}

//void updateDatasetFileInfo(FileManager* file) {
//    if (isvalidDataSetFilenameConvention(file->path)) {
//        file->sample = imageSampleId(file->path);
//        file->label  = imageLabel(file->path);
//        file->suffix = imageSuffix(file->path);
//    }
//}

//bool isvalidDataSetFilenameConvention(const char *path) {
//    if (path == NULL){
//        printf("Pathname is NULL\n");
//        return false;
//    }
//
//    char *filename = getFilename(path, NULL);
//
//    bool valid_dataset_file = iftRegexMatch(filename, "^[0-9]+_[0-9]+.+$");
//
//    iftFree(filename);
//
//    return valid_dataset_file;
//}

//char *getFilename(const char *pathname, const char *suffix) {
//    if (pathname == NULL){
//        printf("Pathname is NULL\n");
//        return NULL;
//    }
//
//
//    char *base = iftSplitStringAt(pathname, SEP_C, -1);
//
//    if ((suffix != NULL) && (!iftCompareStrings(suffix, ""))) {
//        char *out_base = iftRemoveSuffix(base, suffix);
//        iftFree(base);
//        base = out_base;
//    }
//
//    return base;
//}

void destroyCharPointer(void* data){
    char* aux = *((char**) data);
    if(aux == NULL){
        return;
    }
    free(aux);
    data = NULL;
    aux = NULL;
}

GVector* splitsLinesInTextFile(const char* filename){
    size_t initalSize = 10;
    GVector* vector = createVector(initalSize,sizeof(char*));
    vector->freeFunction = destroyCharPointer;
    FILE* file = fopen(filename, "r");
    char line[256];
    int index;
    while (fgets(line, sizeof(line), file)) {
        char* path = (char*)calloc(sizeof(line)+1,sizeof(char));
        index = 0;
        while(line[index] != '\n'){
            path[index] = line[index];
            index++;
        }
        path[index] = '\0';
        VECTOR_PUSH_BACK(char*,vector,path);
    }
    shrinkToFit(vector);
    return vector;
}

