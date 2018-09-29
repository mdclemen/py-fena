//funcs.cl
//20180926 MDC, OpenCL functions to create stiffness and mass matrices

int isInArray(int val, int width, __global const int* array) {
  for (int i = 0; i < width; i++) {
	if (array[i] == val) {
	  return 1;
	}
  }
  return 0;
}

__kernel void makeK(int nn, int mm, __global double* K, __global const double* As, __global const double* phis) {

  int tid = get_global_id(0);
  for (int j = 0; j < nn; j++) {
	for (int m = 0; m < mm; m++) {
	  K[tid*nn + j] += As[m]*(phis[2*m*nn+2*tid]*phis[2*m*nn+2*j] +
							  phis[2*m*nn+2*tid+1]*phis[2*m*nn+2*j+1]);
	}
  }
}

__kernel void makeM(int nn, int mm, int width, __global double* M, __global const int* elements, __global const double* As) {

  int tid = get_global_id(0);
  double val = 0.;
  for (int j = 0; j < nn; j++) {
	val = 0.;
	for (int m = 0; m < mm; m++) {
	  if (isInArray(tid, width, &elements[m*width]) && isInArray(j, width, &elements[m*width])) {
		val += As[m];
	  }
	}
	if (tid == j)
	  M[tid*nn + j] = (1./6.)*val;
	else
	  M[tid*nn + j] = (1./12.)*val;
  }
}
