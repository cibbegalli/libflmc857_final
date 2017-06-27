# This code is released under the libpng license.
# For conditions of distribution and use, see the disclaimer
# and license in png.h

# Location of the zlib library and include files
ZLIBINC = ../zlib/include
ZLIBLIB = ../zlib/lib

LIB=./lib
INCLUDESDIR = ./include
SRC=./source
OBJ=./obj
BIN=./bin

INCLUDE = -I$(INCLUDESDIR)

WARNMORE = -Wwrite-strings -Wpointer-arith -Wshadow \
	-Wmissing-declarations -Wtraditional -Wcast-align \
	-Wstrict-prototypes -Wmissing-prototypes # -Wconversion
CPPFLAGS = -I$(ZLIBINC) # -DPNG_DEBUG=5
CFLAGS = -W -Wall -O2 # $(WARNMORE) -g
FLAGS = $(CPPFLAGS) $(CFLAGS)
LIBS = -lz -lm

# Pre-built configuration
# See scripts/pnglibconf.mak for more options
#PNGLIBCONF_H_PREBUILT = pnglibconf.h.prebuilt


libpng: zlib $(LIB)/libpng.a
	echo "libpng.a built..."

$(LIB)/libpng.a: \
$(OBJ)/png.o \
$(OBJ)/pngerror.o \
$(OBJ)/pngget.o \
$(OBJ)/pngmem.o \
$(OBJ)/pngpread.o \
$(OBJ)/pngread.o \
$(OBJ)/pngrio.o \
$(OBJ)/pngrtran.o \
$(OBJ)/pngrutil.o \
$(OBJ)/pngset.o \
$(OBJ)/pngtrans.o \
$(OBJ)/pngwio.o \
$(OBJ)/pngwrite.o \
$(OBJ)/pngwtran.o \
$(OBJ)/pngwutil.o \

	ar csr $(LIB)/libpng.a \
$(OBJ)/*.o

$(OBJ)/%.o: $(SRC)/%.c
	$(CC) $(FLAGS) -c $< $(INCLUDE) -o $@

zlib:
	cd ../zlib; $(MAKE) ; cd -\

clean:
	rm $(LIB)/lib*.a; rm $(OBJ)/*.o;