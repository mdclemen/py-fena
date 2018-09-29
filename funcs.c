// 20180923 MDC, C function implementations to make stiffness and mass matrices for FEM example
// the speedup gained by using these functions instead of pure python is REAL

int isInArray(int val, int width, int* array) {
  int i;
  for (i = 0; i < width; i++) {
	if (array[i] == val) {
	  return 1;
	}
  }
  return 0;
}

void makeK(int nn, int mm, double* K, double* As, double* phis) {
  int i, j, m;
  for (i = 0; i < nn; i++) {
	for (j = 0; j < nn; j++) {
	  for (m = 0; m < mm; m++) {
		K[i*nn + j] += As[m]*(phis[2*m*nn+2*i]*phis[2*m*nn+2*j] +
							  phis[2*m*nn+2*i+1]*phis[2*m*nn+2*j+1]);
	  }
	}
  }
}

void makeM(int nn, int mm, int width, double* M, int* elements, double* As) {
  int i, j, m;
  double val;
  for (i = 0; i < nn; i++) {
	for (j = 0; j < nn; j++) {
	  val = 0.;
	  for (m = 0; m < mm; m++) {
		if (isInArray(i, width, &elements[m*width]) && isInArray(j, width, &elements[m*width])) {
		  val += As[m];
		}
	  }
	  if (i == j)
		M[i*nn + j] = (1./6.)*val;
	  else
		M[i*nn + j] = (1./12.)*val;
	}
  }
}
