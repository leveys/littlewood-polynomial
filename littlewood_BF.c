#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <stdlib.h>
#include <gsl/gsl_poly.h>

#define DEFAULT_RES 2000

// --- brute force algorithm ---

// recursive function to get the roots of all Littlewood polynomials
int roots_rec(double* poly, int degree, int max_degree, double complex* roots, int index) {

	/* for debugging

	printf("degree = %d\n", degree);

	for (int i = 0; i < max_degree+1; i++) {
		printf("a_%i = %.1f\n", i, poly[i]);
	}
	printf("index =  %i\n", index);
	printf("\n");
	*/

	// find roots of given polynomial using gsl

	double* z = (double*)malloc(degree*2*sizeof(double));

	gsl_poly_complex_workspace* w = gsl_poly_complex_workspace_alloc(degree+1);
	gsl_poly_complex_solve(poly, degree+1, w, z);
	gsl_poly_complex_workspace_free(w);

	for (int i=0; i<degree*2; i+=2) {

		double complex z_i = z[i] + z[i+1]*I;
		roots[index] = z_i;
		index++;
	}
	free(z);

	// check if max degree is reached
	if (degree == max_degree) {
		return index;
	}

	// find roots of next polynomials

	degree++;

	for (int c=-1; c < 2; c += 2) {

		poly[degree] = c;
		index = roots_rec(poly, degree, max_degree, roots, index);
	}

	return index;
}

// find all roots of Littlewood polynomials, using the recursive function
int all_roots(double complex* roots, int max_degree) {

    if (max_degree < 0) {
		return 0;
	}

    // initialize polynomial with first 2 coefficients
    double* poly = (double*)calloc(max_degree+1, sizeof(double));

	if (!poly) {
		return 0;
	}

	// set constant term to 1
	poly[0] = 1.0;

	int index = 0;

	// call the recursive function starting with a polynomial of degree 1
    for (int c1 = -1; c1 < 2; c1 += 2) {

		poly[1] = (double)c1;

		index = roots_rec(poly, 1, max_degree, roots, index);

    }

    return 1;
}

// get the expected number of roots for the given maximum degree using a formula
int num_roots(int max_degree) {
	return 2*(1 + (max_degree - 1) * pow(2, max_degree));
}


// --- create image from array of roots ---

// compare function to sort complex values according to their position in a text file
// the text file will be created top to bottom, left to right; this means the roots are sorted first by imaginary part, descending, and then by real part, ascending.
int cmpfunc(const void* a, const void* b) {

	int complex z_a = *(int complex*) a;
	int complex z_b = *(int complex*) b;

	int cmp_imag = cimag(z_b) - cimag(z_a);

	if (cmp_imag) {
		return cmp_imag;
	}

	return creal(z_a) - creal(z_b);
}

// create a PBM text file which can later be converted to a PNG file
void make_image(char* filename, double xmin, double xmax, double ymin, double ymax, int res, double complex* roots, int roots_size) {

	// create the file and set the correct header
	FILE* fp = fopen(filename, "w");
    fprintf(fp, "P1\n%i %i\n", res, res);

	// convert the complex roots to a format which can be easily written to the text file 

	int complex* roots_int = (int complex*)malloc(sizeof(int complex)*roots_size);

	for (int i = 0; i < roots_size; i++) {

		double z_real = creal(roots[i]);
		double z_imag = cimag(roots[i]);

		if (xmin <= z_real && z_real <= xmax && ymin <= z_imag && z_imag <= ymax) {
			int x = (int) ((z_real - xmin) / (xmax - xmin) * res);
			int y = (int) ((z_imag - ymin) / (ymax - ymin) * res);
			roots_int[i] = x + y*I;
		} else { 
			// if the root is not whithin the specified bounds, make a dummy value
			roots_int[i] = 0 - I;
		}

		/* for debugging
		printf("z = %.5f + %.5fi\n", z_real, z_imag);
		printf("z = %.5f + %.5fi\n", creal(roots_int[i]), cimag(roots_int[i]));
		printf("\n");
		*/
	}

	// sort the converted roots
	qsort(roots_int, roots_size, sizeof(int complex), cmpfunc);

	/* for debuggin
	printf("after sorting:\n\n");
	for (int i = 0; i < roots_size; i++) {
		
		printf("z = %.5f + %.5fi\n", creal(roots_int[i]), cimag(roots_int[i]));
	}
	*/
	
	// write the roots to the file

	int index_real = 0;
	int index_imag = res;

	for (int i = 0; i < roots_size; i++) {

		if (i > 0 && roots_int[i] == roots_int[i-1]) {
			continue;
		}

		while (index_imag > (int)cimag(roots_int[i]) && index_imag >= 0) {

			while (index_real < res) {
				fprintf(fp, "1 "); // 1 is a black pixel
				index_real++;
			}

			index_real = 0;
			index_imag--;
			fputc('\n', fp);
		}

		if (index_imag < 0) {
			break;
		}

		while (index_real < (int)creal(roots_int[i])-1) {
			fprintf(fp, "1 ");
			index_real++;
		}

		fprintf(fp, "0 "); // 0 is a white pixel
		index_real++;
	}

	while (index_imag > 0) {

		while (index_real < res) {

			fprintf(fp, "1 ");
			index_real++;
		}

		index_real = 0;
		index_imag--;
		fputc('\n', fp);
	}

	free(roots_int);

	fclose(fp);
}



int main(int argc, char* argv[]) {

	struct timeval start, stop;
	gettimeofday(&start, NULL);

	if (argc < 3) {
		printf("at least two arguments are required: max_degree and filename\n");
		printf("usage: %s <max_degree> <filename> [resolution (default=%i)]\n", argv[0], DEFAULT_RES);
		return 1;
	}

	if (argc > 4) {
		printf("at most three arguments are permitted\n");
		printf("usage: %s <max_degree> <filename> [resolution (default=%i)]\n", argv[0], DEFAULT_RES);
		return 1;
	}

	int max_degree = atoi(argv[1]);

	if (max_degree < 1) {
		printf("max_degree argument must be > 0\n");
		return 1;
	}

	int res = DEFAULT_RES;

	if (argc > 3) {
		res = atoi(argv[3]);
		if (res < 1) {
			printf("resolution argument should be > 0\n");
			return 1;
		}
	}


	printf("max_degree = %i\n", max_degree);
	printf("number of roots = %i\n", num_roots(max_degree));

	// find roots of littlewood functions
	
	int roots_size = num_roots(max_degree);
	double complex* roots = (double complex*)malloc(sizeof(double complex)*roots_size);

	if (!roots) {
		return 1;
	}

	all_roots(roots, max_degree);

	/* for debugging
	for (int i=0; i<roots_size; i++) {

		printf("z = %.10f + %.10fi\n", creal(roots[i]), cimag(roots[i]));
	}
	*/

	// create the image

	make_image(argv[2], -1.6, 1.6, -1.6, 1.6, res, roots, roots_size);

	free(roots);

	gettimeofday(&stop, NULL);
	double time = stop.tv_sec - start.tv_sec + (double)(stop.tv_usec - start.tv_usec) / 1000000.0;
	printf("total time = %.3f sec\n\n", time);

	return 0;
}