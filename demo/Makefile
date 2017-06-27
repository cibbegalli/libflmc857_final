COMPILER = $(CXX)
#COMPILER = clang++

#FLAGS= -g -O0 -Wall -Wno-unused-result -std=gnu99 -pedantic
FLAGS=   -Wall -Wno-unused-result  -pedantic -Wno-write-strings -fsanitize=address -g -fopenmp
#-std=gnu99

FLLIBPNGINCLUDES = -I../externals/libpng/include
FLLIBJPEGINCLUDES = -I../externals/libjpeg/include
FLLIBTIFFINCLUDES = -I../externals/libtiff/include
FLLIBSVMINCLUDES = -I../externals/libsvm/include

FL = -L../lib -lFL
LINKER =  $(FL)

FLFLAGS  = -I../include
FLAGSIM =  $(FLLIBPNGINCLUDES) $(FLLIBJPEGINCLUDES) $(FLLIBSVMINCLUDES) -I../include $(FLAGS)


libFL:
	cd ../; $(MAKE) ; cd demo\

$@.cpp: libFL $@.cpp
	make
	$(COMPILER) $(FLAGSIM) $@.cpp -o $@.bin  $(LINKER); \

#Removing everything without an extension
#clean:
#	find . -type f -executable -exec rm '{}' \;

clean:
	rm *.bin ;


cleanImages:
	rm *.ppm *.pgm *.gif *.png *.jpg *.jpeg

cleanText:
	rm *.txt
