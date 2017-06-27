//#include "FL.h"
#include "image.h"
#include "common.h"
#include "adjacencyRelation.h"
#include "kernel.h"
#include "featureVector.h"
#include "vector.h"


#define HOG_r_gradient 3.0

// tamanho janela de detecção (em pixels)
#define HOG_N1 12
#define HOG_M1 16

// tamanho célula (em pixels)
#define HOG_N2 4
#define HOG_M2 4

// tamanho bloco (em células)
#define HOG_N3 2
#define HOG_M3 2

#define HOG_epsilon 0.000001

#define HOG_THETA 45.0
#define HOG_B1 22.5
#define HOG_B2 67.5
#define HOG_B3 112.5
#define HOG_B4 157.5
#define HOG_B5 202.5
#define HOG_B6 247.5
#define HOG_B7 292.5
#define HOG_B8 337.5

#ifndef _HOG_H_
#define _HOG_H_


int get_bin1(float z, float mag) {
	if(mag == 0.0) return 0;
	if(z >= B1 && z < B2) return 1;
	if(z >= B2 && z < B3) return 2;
	if(z >= B3 && z < B4) return 3;
	if(z >= B4 && z < B5) return 4;
	if(z >= B5 && z < B6) return 5;
	if(z >= B6 && z < B7) return 6;
	if(z >= B7 && z < B8) return 7;
	return 8;
}

int get_bin2(float z, float mag) {
	if(mag == 0.0) return 0;
	if(z >= B1 && z < B2) return 2;
	if(z >= B2 && z < B3) return 3;
	if(z >= B3 && z < B4) return 4;
	if(z >= B4 && z < B5) return 5;
	if(z >= B5 && z < B6) return 6;
	if(z >= B6 && z < B7) return 7;
	if(z >= B7 && z < B8) return 8;
	return 1;
}

void celulas_adjacentes(int i, int j, int *celula1, int *celula2, int *celula3, int *celula4) {

	int celula =  i/HOG_N2 + (j/HOG_N2)*(HOG_N1/HOG_N2);
	int x_centro = ((celula % (HOG_N1/HOG_N2)) * HOG_N2 + HOG_N2/2);
	int y_centro = ((celula /(HOG_N1/HOG_N2)) * HOG_M2 + HOG_M2/2);

	if(i < x_centro) {
		if(j < y_centro) {

			(*celula4) = celula;
			(*celula3) = (*celula4) - 1;
			(*celula2) = (*celula4) - (HOG_N1/HOG_N2);
			(*celula1) = (*celula2) - 1;

			if(((*celula2) % (HOG_N1/HOG_N2)) == 0)
				(*celula1) = -1;

			if(((*celula4) % (HOG_N1/HOG_N2)) == 0)
				(*celula3) = -1;

		} else {

			(*celula2) = celula;
			(*celula1) = (*celula2) - 1;                                         
			(*celula4) = (*celula2) + (HOG_N1/HOG_N2);
			(*celula3) = (*celula4) - 1;

			if(((*celula2) % (HOG_N1/HOG_N2)) == 0)
				(*celula1) = -1;

			if(((*celula4) % (HOG_N1/HOG_N2)) == 0)
				(*celula3) = -1;
		}

	} else {
		if(j < y_centro) {

			(*celula3) = celula;
			(*celula4) = (*celula3) + 1;
			(*celula1) = (*celula3) - (HOG_N1/HOG_N2);
			(*celula2) = (*celula1) + 1;

			if(((*celula2) % (HOG_N1/HOG_N2)) == 0)
				(*celula2) = -1;
			if(((*celula4) % (HOG_N1/HOG_N2)) == 0)
				(*celula4) = -1;
		} else {

			(*celula1) = celula;
			(*celula2) = (*celula1) + 1;
			(*celula3) = (*celula1) + (HOG_N1/HOG_N2);
			(*celula4) = (*celula3) + 1;

			if(((*celula2) % (HOG_N1/HOG_N2)) == 0)
				(*celula2) = -1;                

			if(((*celula4) % (HOG_N1/HOG_N2)) == 0)
				(*celula4) = -1;

		}
	}

	if((*celula1) < 0) (*celula1) = -1;
	if((*celula2) < 0) (*celula2) = -1;
	if((*celula3) < 0) (*celula3) = -1;
	if((*celula4) < 0) (*celula4) = -1;

	int numCelulas = (HOG_N1*HOG_M1)/(HOG_N2*HOG_M2);

	if((*celula1) >= numCelulas) (*celula1) = -1;
	if((*celula2) >= numCelulas) (*celula2) = -1;
	if((*celula3) >= numCelulas) (*celula3) = -1;
	if((*celula4) >= numCelulas) (*celula4) = -1;

	if((*celula1) == -1 && (*celula2) == -1 && (*celula4) == -1) { //canto superior direito
		celulas_adjacentes(i-HOG_N2/2, j+HOG_M2/2, celula1, celula2, celula3, celula4);
	}

	else if((*celula1) == -1 && (*celula2) == -1 && (*celula3) == -1) { //canto superior esquerdo
		celulas_adjacentes(i+HOG_N2/2, j+HOG_M2/2, celula1, celula2, celula3, celula4);
	}

	else if((*celula1) == -1 && (*celula3) == -1 && (*celula4) == -1) { //canto inferior esquerdo
		celulas_adjacentes(i+HOG_N2/2, j-HOG_M2/2, celula1, celula2, celula3, celula4);
	}

	else if((*celula2) == -1 && (*celula3) == -1 && (*celula4) == -1) { //canto inferior direito
		celulas_adjacentes(i-HOG_N2/2, j-HOG_M2/2, celula1, celula2, celula3, celula4);
	}

	else if((*celula1) == -1 && (*celula2) == -1) { //linha superior
		celulas_adjacentes(i, j+HOG_M2/2, celula1, celula2, celula3, celula4);
	}

	else if((*celula3) == -1 && (*celula4) == -1) { //linha inferior
		celulas_adjacentes(i, j-HOG_M2/2, celula1, celula2, celula3, celula4);
	}

	else if((*celula1) == -1 && (*celula3) == -1) { //coluna esquerda
		celulas_adjacentes(i+HOG_N2/2, j, celula1, celula2, celula3, celula4);
	}

	else if((*celula2) == -1 && (*celula4) == -1) { //coluna direita
		celulas_adjacentes(i-HOG_N2/2, j, celula1, celula2, celula3, celula4);
	}
}

int coordenada_centro_eixox( int celulaA, int celulaB) {
	if(celulaA < 0) {
		if(celulaB < 0)
			return 0;
		return (celulaB % (HOG_N1/HOG_N2)) * HOG_N2 + HOG_N2/2;
	} 
	return (celulaA % (HOG_N1/HOG_N2)) * HOG_N2 + HOG_N2/2;	
}

int coordenada_centro_eixoy(int celulaA, int celulaB) {
	if(celulaA < 0) {
		if(celulaB < 0)
			return 0;
		return ((celulaB /(HOG_N1/HOG_N2)) * HOG_M2 + HOG_M2/2);
	}
	return ((celulaA /(HOG_N1/HOG_N2)) * HOG_M2 + HOG_M2/2);
}


void gradient(Image* img, Image **g_mag, Image **g_orient);

GVector* hog(Image* img);

#endif //LIBFL_MORPHOLOGY_H
