# makefile for zlib using gcc (generic, static library)
#
# This code is released under the zlib license.
# For conditions of distribution and use, see the disclaimer
# and license in zlib.h

LIB=./lib
INCLUDESDIR = ./include
SRC=./source
OBJ=./obj
BIN=./bin

INCLUDE = -I$(INCLUDESDIR)
#LOC = -DASMV
#LOC = -DZLIB_DEBUG -g
FLAGS = $(LOC) -W -Wall -O3


zlib: $(LIB)/zlib.a
	echo "zlib.a built..."

$(LIB)/zlib.a: \
$(OBJ)/adler32.o \
$(OBJ)/crc32.o \
$(OBJ)/deflate.o \
$(OBJ)/gzclose.o \
$(OBJ)/gzlib.o \
$(OBJ)/gzread.o \
$(OBJ)/gzwrite.o \
$(OBJ)/infback.o \
$(OBJ)/inffast.o \
$(OBJ)/inflate.o \
$(OBJ)/inftrees.o \
$(OBJ)/trees.o \
$(OBJ)/uncompr.o \
$(OBJ)/zutil.o \

	ar csr $(LIB)/zlib.a \
$(OBJ)/*.o

$(OBJ)/%.o: $(SRC)/%.c
	$(CC) $(FLAGS) -c $< $(INCLUDE) -o $@



clean:
	rm $(LIB)/zlib*.a; rm $(OBJ)/*.o;