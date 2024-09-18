#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>

#define TOL 1e-3
#define DEFAULT_RES 2000
#define MIN(x, y) (((x) < (y)) ? (x) : (y))


// representation of Littlewood polynomials

// because we constantly have to extend the polynomials, arrays are not the best choice to represent them
// instead, i chose to represent them as bitvectors, because each coefficient in a Littlewood polynomial is one of 2 values: 1 or -1
// this can be easily represented as a bitvector:

//   x^3 - x^2 + x + 1
//   1     0     1   1
//  = 27

// the only problem here is that if the leading coefficient is -1, the first bit in the bitvector has to be 0, which would be the same as there not being anything
// we can solve this be adding a leading one to every bitvector:

//  1   -x^3 - x^2 + x + 1
//  1    0     0     1   1
//  = 19

// it is now easy to extend a given polynomial:

//       x^4 - x^3 - x^2 + x + 1
//  1    1     0     0     1   1
//  original bitvector + 2 * the highest power of 2 of original bitvector
//  19 + 2*16 = 51
// or:
//      -x^4 - x^3 - x^2 + x + 1
//  1    0     0     0     1   1
//  original bitvector + the highest power of 2 of original bitvector
//  19 + 16 = 35


// evaluation of a Littlewood polynomial

double complex eval_Littlewood(int poly, double complex z) {

    double complex res = 0.0;

    double complex z_pow = 1.0;

    while (poly > 1) {

        int coeff = poly%2 * 2 - 1;     // 1 -> 1 and 0 -> -1

        res += coeff * z_pow;

        poly = poly >> 1;

        z_pow *= z;
    }

    return res;
}


// check if a complex value z is a root of any Littlewood polynomial
// this is checked by recursively increasing the degree of the polynomials
int root_rec(int poly, double complex z, int max_degree, double tol) {

    double l_z_abs = cabs(eval_Littlewood(poly, z));

    //printf("poly = %i\n", poly);
    //printf("|l(z)| = %.16f\n", l_z_abs);
    //printf("\n");

    // check if poly(z) is close enough to 0
    if (l_z_abs < tol) {

        //printf("poly: %i is 0 for z = %.1f + %.1fI\n", poly, creal(z), cimag(z));
        
        // z is a solution of poly
        return 1;
    }

    // degree + 1 of poly = (int)log2(poly)
    // e.g. 101101 (= 45) -> -x^4 + x^3 + x^2 - x + 1, degree + 1 = 5
    // (int)log2(45) = 5

    int degplus1 = (int)log2(poly);

    if (degplus1 >= max_degree+1) {

        // max degree reached
        return 0;
    }

    if (l_z_abs > pow(cabs(z), degplus1)/(1-cabs(z))) {

        // z is not a root of any extension of poly
        return 0;
    }

    // extend poly to next degree

    int highestpowof2 = pow(2, degplus1);

    int polyplus = poly + 2*highestpowof2;
    int polymin = poly + highestpowof2;

    return root_rec(polyplus, z, max_degree, tol) || root_rec(polymin, z, max_degree, tol);
}

int is_root(double complex z, int max_degree, double tol) {

    //printf("z = %.1f + %.1fI", creal(z), cimag(z));

    // start recursive search with degree 1 
    // => starting polynomials: x + 1, -x + 1
    // =>                       111     101
    // (the polynomials with -1 as the constant term have the exact same roots as their counterparts)

    for (int poly=5; poly < 8; poly+=2) {
        if (root_rec(poly, z, max_degree, tol)) {
            return 1;
        }
    }

    return 0;
}

// get the expected size of the array to store the roots for the given maximum degree using a formula
int max_roots_size(int max_degree, int res) {
	return MIN(4*(1 + (max_degree - 1) * pow(2, max_degree)), res*res);
}

void loop_pixels(int max_degree, int res, double complex roots[]) {

    int roots_index = 0;

    for (int real=0; real < res; real++) {
        for (int imag=0; imag < res; imag++) {

            double complex z = (double)real/res + (double)imag/res * I;

            if (cabs(z) <= 1 && is_root(z, max_degree, TOL)) {

                roots[roots_index] = z;
                roots_index++;

                if (z != 0) {
                    roots[roots_index] = -z;
                    roots_index++;

                    if (imag != 0) {
                        roots[roots_index] = conj(z);
                        roots[roots_index+1] = -conj(z);
                        roots_index += 2;
                    }

                    if (1 - cabs(z) > TOL) {
                        double complex inv = 1/z;
                        roots[roots_index] = inv;
                        roots[roots_index+1] = -inv;
                        roots_index += 2;

                        if (imag != 0) {
                            roots[roots_index] = conj(inv);
                            roots[roots_index+1] = -conj(inv);
                            roots_index += 2;
                        }
                    }
                }
            }
        }
    }
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

	/* for debugging
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

	int roots_size;

	if (max_degree < 5) {
		roots_size = 300;
	} else {
		roots_size = max_roots_size(max_degree, res);
	}

	double complex* roots = (double complex*)malloc(sizeof(double complex)*roots_size);

	if (!roots) {
		return 1;
	}

	loop_pixels(max_degree, (int)res/3.2, roots);

	// create the image
	make_image(argv[2], -1.6, 1.6, -1.6, 1.6, res, roots, roots_size);

	free(roots);
    
	gettimeofday(&stop, NULL);
	double time = stop.tv_sec - start.tv_sec + (double)(stop.tv_usec - start.tv_usec) / 1000000.0;
	printf("total time = %.3f sec\n\n", time);

	return 0;
}