#include "hog.h"

void gradient(Image* img, Image **g_mag, Image **g_orient) {

	//////////////////////////////////////////////////////////////////////////////
	// Calculando pesos do kernel para calcular imagem magnitude e orientação
	AdjacencyRelation *adjRel = createCircularAdjacency(HOG_r_gradient);

	Kernel *Kx = createKernel(adjRel);
	Kernel *Ky = createKernel(adjRel);

	double c = 2*pow((HOG_r_gradient/3),2.0);

	
	for(int i=0; i<adjRel->n; i++) {

		double dist = sqrt(adjRel->dx[i]*adjRel->dx[i] + adjRel->dy[i]*adjRel->dy[i])+epsilon;

		
		int dist_x = adjRel->dx[i];
		double Wx = exp( (-1) * dist * dist/c) * (dist_x/dist);

		int dist_y = adjRel->dy[i];
		double Wy = exp( (-1) * dist * dist/c) * (dist_y/dist);
		

		Kx->weight[i] = Wx;
		Ky->weight[i] = Wy;
	
	}
	//////////////////////////////////////////////////////////////////////////////

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
	iftDestroyAdjRel(&adjRel);
	iftDestroyKernel(&Kx);
	iftDestroyKernel(&Ky);
	*g_mag = img_mag;
	*g_orient = img_orient;
}

GVector* hog(Image* img) {
	
	Image *img_mag;
	Image *img_orient;
	gradient(img, &img_mag, &img_orient);

	int numCelulas = (HOG_N1*HOG_M1)/(HOG_N2*HOG_M2);
	int numBins = 9;


	////////////////////////////////////////////////////////////////////////
	//Calculando histogramas das células
	FeatureMatrix *histogramas = createFeatureMatrix(numCelulas, numBins);
	for(int i=0; i<numCelulas; i++){
		setValueInFeatureVector(histogramas->featureVector[i], 0.0);
	}


	for(int i=0; i<img->nx; i++) {//for(int i=0; i<HOG_N1; i++){
		for(int j=0; j<img->ny; j++) {//for(int j=0; j<HOG_M1; j++)

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
			

			histogramas->featureVector[celula2]->features[b1]+=w9;
			histogramas->featureVector[celula2]->features[b2]+=w13;
			

			histogramas->featureVector[celula3]->features[b1]+=w7;
			histogramas->featureVector[celula3]->features[b2]+=w11;
			
			histogramas->featureVector[celula4]->features[b1]+=w10;
			histogramas->featureVector[celula4]->features[b2]+=w14;
			
		}
	}
	//////////////////////////////////////////////////////////////////////
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

	destroyImage(&img_mag);
	destroyImage(&img_orient);
	destroyFeatureVector(&soma);
	destroyFeatureMatrix(&histogramas);
	destroyFeatureMatrix(&histogramasConcatNormalizados);

	return histVec;
	
}

	
