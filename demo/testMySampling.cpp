#include "FL.h"

int main(int argc, char **argv) {
	Image *img = readImage("../data/000001_00000805.ppm");
	Image *imgISF = readImage("../data/000001_000008.pgm");


	GVector* r = samplingBySuperpixelAndGradient(img, imgISF, 80.0, 40, 40, 4, 4);

	destroyImage(&img);
	destroyImage(&imgISF);
	
	destroyVector(&r);
	return 0;
}