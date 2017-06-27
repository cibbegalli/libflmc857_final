
LIB=./lib
INCLUDESDIR = ./include
SRC=./source
OBJ=./obj
BIN=./bin

INCLUDE = -I$(INCLUDESDIR)
CFLAGS = -Wall -Wconversion -O3 -fPIC
FLAGS =  $(CFLAGS)

libsvm: $(LIB)/libsvm.a
	echo "libsvm.a built..."

$(LIB)/libsvm.a: \
$(OBJ)/svm.o \

	ar csr $(LIB)/libsvm.a \
$(OBJ)/*.o

$(OBJ)/%.o: $(SRC)/%.cpp
	$(CC) $(FLAGS) -c $< $(INCLUDE) -o $@

clean:
	rm $(LIB)/lib*.a; rm $(OBJ)/*.o;