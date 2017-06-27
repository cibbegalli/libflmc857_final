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

#include "FL.h"

#define HOG_THETA 45.0
#define HOG_B1 22.5
#define HOG_B2 67.5
#define HOG_B3 112.5
#define HOG_B4 157.5
#define HOG_B5 202.5
#define HOG_B6 247.5
#define HOG_B7 292.5
#define HOG_B8 337.5

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


int main(int argc, char **argv) {

	//////////////////////////////////////////////////////////////////////////////
	// Calculando pesos do kernel para calcular imagem magnitude e orientação
	AdjacencyRelation *adjRel = createCircularAdjacency(HOG_r_gradient);

	Kernel *Kx = createKernel(adjRel);
	Kernel *Ky = createKernel(adjRel);

	double c = 2*pow((HOG_r_gradient/3),2.0);

	if(DEBUG) {
		printf("%f\n", c);
		printf("|adjRel|=%d\n", adjRel->n);
	}
	
	for(int i=0; i<adjRel->n; i++) {

		double dist = sqrt(adjRel->dx[i]*adjRel->dx[i] + adjRel->dy[i]*adjRel->dy[i])+epsilon;

		
		int dist_x = adjRel->dx[i];
		double Wx = exp( (-1) * dist * dist/c) * (dist_x/dist);

		int dist_y = adjRel->dy[i];
		double Wy = exp( (-1) * dist * dist/c) * (dist_y/dist);
		

		Kx->weight[i] = Wx;
		Ky->weight[i] = Wy;

		if(DEBUG) {
			printf("%d %d\n", adjRel->dx[i], adjRel->dy[i]);
			printf("dist=%f\n", dist);
			printf("Wx=%f Wy=%f\n", Kx->weight[i], Ky->weight[i]);
		}
	
	}
	//////////////////////////////////////////////////////////////////////////////

	Image *img = readImage("../data/lena.pgm");

	Image* img_mag = createImage(img->nx, img->ny, 1);
	Image* img_orient = createImage(img->nx, img->ny, 1);

	for(int i=0; i<img->nx; i++) { // Percorrendo pixels da imagem
		for(int j=0; j<img->ny; j++) {

			float val = imageVal(img, i, j); //intensidade do pixel atual 
			int index = (j*img->nx) + i; //índice do pixel atual

			float Gx = 0.0, Gy = 0.0;
			
			for(int k=0; k<adjRel->n; k++) { //Percorrendo relação de adjacência

				int x = i+adjRel->dx[k]; //coordenada x do pixel adjacente
				int y = j+adjRel->dy[k]; //coordenada y do pixel adjacente

	
				if(isValidPixelCoordinate(img, x, y)) {

					float val_adj = imageVal(img, x, y);
					Gx += (val_adj - val) * Kx->weight[k];
					Gy += (val_adj - val) * Ky->weight[k];
				}				
			}
			
			double G = sqrt((double)Gx*Gx + Gy*Gy); //obtendo magnitude	
			img_mag->channel[0][index] = G;

			Gx = Gx/G; Gy = Gy/G;

			double angle;

			if(Gy >= 0.0){
				angle = (180/PI)*acos(Gx);
			} else {
				angle = 360 - (180/PI)*acos(Gx);
			}
			img_orient->channel[0][index] = angle;

		}

	}

	//////////////////////////////////////////////////////////////////////////////
	int numCelulas1 = (HOG_N1*HOG_M1)/(HOG_N2*HOG_M2);

	//Validação tamanho da janela/imagem/patch
	if(img->nx % HOG_N2 != 0 || img->ny % HOG_M2 != 0) {
		printf("Imagem não alinhada número de células\n");
		//return 0;
	}
	
	int numCelulas = (HOG_N1*HOG_M1)/(HOG_N2*HOG_M2);
	int numBins = 9;

	FeatureMatrix *histogramas = createFeatureMatrix(numCelulas, numBins);
	for(int i=0; i<numCelulas; i++){
		setValueInFeatureVector(histogramas->featureVector[i], 0.0);
	}
	
	for(int i=0; i<HOG_N1; i++){
		for(int j=0; j<HOG_M1; j++) {
	//for(int i=0; i<img->nx; i++) {//for(int i=0; i<HOG_N1; i++){
	//	for(int j=0; j<img->ny; j++) {//for(int j=0; j<HOG_M1; j++)

			// Valores magnitude e orientação
			float w = imageVal(img_mag, i, j);
			float z = imageVal(img_orient, i, j);

			int b1 = get_bin1(z, w);
			int b2 = get_bin2(z, w);

			float z1 = b1*THETA - THETA/2.0;
			float z2 = b2*THETA - THETA/2.0;

			int celula1, celula2, celula3, celula4;
			celulas_adjacentes(i,j, &celula1, &celula2, &celula3, &celula4);


			int x_q1q3 = coordenada_centro_eixox(celula1, celula3);
			int x_q2q4 = coordenada_centro_eixox(celula2, celula4);
			int y_q1q2 = coordenada_centro_eixoy(celula1, celula2);
			int y_q3q4 = coordenada_centro_eixoy(celula3, celula4);

			if(DEBUG) {
				printf("[%d %d] %d %d %d %d\n", i, j, celula1, celula2, celula3, celula4);
				printf("(%d %d) (%d %d) (%d %d) (%d %d)\n", x_q1q3, y_q1q2, x_q2q4, y_q1q2, x_q1q3, y_q3q4, x_q2q4, y_q3q4);
			}
			
			int x = i, y = j;

			double w1 = w*(x_q2q4-x)/HOG_N2;
			double w2 = w*(x-x_q1q3)/HOG_N2;


			double w3 = w1*(y-y_q1q2)/HOG_M2;
			double w4 = w1*(y_q3q4-y)/HOG_M2;
			double w5 = w2*(y-y_q1q2)/HOG_M2;
			double w6 = w2*(y_q3q4-y)/HOG_M2;

			int dif_z = z2 - z;
			if(dif_z < 0)
				dif_z += 360;
			double w7 = w3*(dif_z)/THETA;
			double w8 = w4*(dif_z)/THETA;
			double w9 = w6*(dif_z)/THETA;
			double w10 = w5*(dif_z)/THETA;

			dif_z = z - z1;
			if(dif_z < 0)
				dif_z += 360;
			double w11 = w3*(dif_z)/THETA;
			double w12 = w4*(dif_z)/THETA;
			double w13 = w6*(dif_z)/THETA;
			double w14 = w5*(dif_z)/THETA;


			histogramas->featureVector[celula1]->features[b1]+=w8;
			histogramas->featureVector[celula1]->features[b2]+=w12;
			
			//printf("%f %f\n", w8, w12);

			histogramas->featureVector[celula2]->features[b1]+=w9;
			histogramas->featureVector[celula2]->features[b2]+=w13;
			
			//printf("%f %f\n", w9, w13);

			histogramas->featureVector[celula3]->features[b1]+=w7;
			histogramas->featureVector[celula3]->features[b2]+=w11;
			
			//printf("%f %f\n", w7, w11);
			histogramas->featureVector[celula4]->features[b1]+=w10;
			histogramas->featureVector[celula4]->features[b2]+=w14;
			
			//printf("%f %f\n", w14, w10);			

		}
	}
	
    

	// Stride blocos e normalização
	int numBlocos = (HOG_N1/HOG_N2-HOG_N3+1)*(HOG_M1/HOG_M2-HOG_M3+1);
	
	FeatureMatrix *histogramasConcatNormalizados = createFeatureMatrix(numBlocos, HOG_N3*HOG_M3*numBins);
	for(int i=0; i<numBlocos; i++){
		setValueInFeatureVector(histogramas->featureVector[i], 0.0);
	}
	
	int numCelulasHorizontal = (HOG_N1/HOG_N2);
	int numBlocosHorizontal = (HOG_N1/HOG_N2-HOG_N3+1);

	
	FeatureVector *soma = createFeatureVector(numBlocos);
	setValueInFeatureVector(soma, 0.0);

		for(int k=0; k<numBlocos; k++) {
		int ind = 0;
		for(int i=0; i<HOG_N3; i++) {
			int celula = (k/numBlocosHorizontal)*numCelulasHorizontal + i*numCelulasHorizontal + k%(numCelulasHorizontal-HOG_N3+1);
			for(int j=0; j<HOG_M3; j++) {
				printf("Bloco %d, celula %d \n", k, celula+j);
				for(int b=0; b<numBins; b++) {
					float val = histogramas->featureVector[celula+j]->features[b];
					histogramasConcatNormalizados->featureVector[k]->features[ind] = val; 
					soma->features[k] += (val*val);
					ind++;
				}
			}
		}
	}
	

	GVector *histVec = createNullVector(numBlocos*HOG_N3*HOG_M3*numBins,sizeof(float));
	for(int k=0; k<numBlocos; k++) {
		float norm = sqrt(soma->features[k]) + epsilon;
		for(int i=0; i<HOG_N3*HOG_M3*numBins; i++) {
			float val = histogramasConcatNormalizados->featureVector[k]->features[i]/norm;
			VECTOR_GET_ELEMENT_AS(float,histVec,i) = val;
		}
	}


    //return histVec;

	writeImageP2(img_mag, "gradMagnitude.pgm");
    writeImageP2(img_orient, "gradPhase.pgm");

	destroyKernel(&Kx);
	destroyKernel(&Ky);
	destroyAdjacencyRelation(&adjRel);
	destroyImage(&img);
	destroyImage(&img_mag);
	destroyImage(&img_orient);
	destroyFeatureVector(&soma);

	destroyFeatureMatrix(&histogramas);
	destroyFeatureMatrix(&histogramasConcatNormalizados);
	destroyVector(&histVec);
	return 0;
}
/*
	*g_mag = img_mag;
	*g_orient = img_orient;
}*/