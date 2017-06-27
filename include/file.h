/*
 *Created by Deangeli Gomes Neves
 *
 * This software may be freely redistributed under the terms
 * of the MIT license.
 *
 */

#ifndef _FILE_H_
#define _FILE_H_

#include "common.h"
#include <sys/stat.h>
#include "vector.h"

#define SEP_C "/"

typedef struct _fileManager {
    /** File path */
    char *path;
    /**Label for machine learning dataset files*/
    int label;
    /** Sample id following the new convention for IFT datasets (Thiago Vallin Spina, Mar 09, 2016) */
    int sample;
    /** Suffix following the new convention for IFT datasets (Thiago Vallin Spina, Mar 09, 2016) */
    char *suffix;

    /**Status for machine learning dataset files
     * iftFile can be marked as TRAIN, TEST, OUTLIER as the iftSample*/
    char status;
} FileManager;


typedef struct _fileset {
    /** Number of Files */
    size_t n;
    /** Array of Files */
    FileManager **files;
} FileSet;


typedef struct _directoryManager {
    char *path;
    /** Number of files from the directory */
    size_t nfiles;
    /** Vector of iftFiles with the children */
    FileManager **files;
    /** Number of files from the directory */
    size_t nsubdirs;
    /** Vector of iftDir with the subdirs */
    struct _directoryManager **subdirs;
} DirectoryManager;

DirectoryManager *loadDirectory(const char *dir_pathname, size_t hier_levels);
bool pathnameExists(const char *pathname);
bool directoryExists(const char *pathname);
void listDirectory(DirectoryManager *root_dir, size_t hier_levels);
void listDirectoryRecursively(DirectoryManager *dir, size_t hier_levels, size_t curr_level);
void countFilesInDirectory(const char *dir_pathname, size_t *nfiles, size_t *nsubdirs);
bool endsWith(const char *str, const char *suffix);
void countFilesInDirectory(const char *dir_pathname, size_t *nfiles, size_t *nsubdirs);
int compareFiles(const void *a, const void *b);
int compareDirectories(const void *a, const void *b);
int countOnlyFilesInDirectory(const char *dir_pathname, const char *extension);
void updateDatasetFileInfo(FileManager* file);
bool isvalidDataSetFilenameConvention(const char *path);
char *getFilename(const char *pathname, const char *suffix);
FileManager *createFileManager(const char *pathname);
void sortDirectory(DirectoryManager *files);
void listFilesFromDirectory(DirectoryManager *dir, const char *suffix);
char *joinPathnames(const char *pathname1, const char *pathname2);
DirectoryManager *loadFilesFromDirBySuffix(const char *dir_pathname, const char *suffix);

void destroyDirectoryManager(DirectoryManager **directoryManager);
void destroyFileManager(FileManager **fileManager);
void destroyFileSet(FileSet **farr);

void findTrueLabelInCurrentDirectory(DirectoryManager *directoryManager);
int findTrueLabelInName(char *pathname);


GVector* splitsLinesInTextFile(const char* filename);

#endif //LIBFL_FILE_H
