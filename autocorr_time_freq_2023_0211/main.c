#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define DebugMode 0
#define SWAP(a, b) tempr=(a);(a)=(b);(b)=tempr;
typedef double RealType;
typedef int intType;

void four1(RealType data[], unsigned long nn, int isign);
void autocorr_freq_domain(RealType *, unsigned long);

void autocorr_time_domain(int, int, RealType *, RealType *, int);

void get_rand_sig(RealType *sig, int n);

int main(void) {

	int i, RunTimes = 10000;
	double t1, t2, dt1;
	// IO
	srand(1);
	unsigned long n = 150;
	unsigned long nNextPow = 256;
	int nlags = n - 1;
	// input
	RealType *input = (RealType *)calloc(n, sizeof(RealType));
	get_rand_sig(input, n);
	// output
	RealType *y1 = (RealType *)calloc(nlags, sizeof(RealType));
	RealType *y2 = (RealType *)calloc(nlags, sizeof(RealType));

	// freq
	RealType *acf1 = (RealType *)calloc(4 * nNextPow, sizeof(RealType));
	t1 = clock();

	for (int c = 0; c < RunTimes; c++) {
		for (i = 0; i < n; i++)
			acf1[2 * i] = input[i];
		autocorr_freq_domain(acf1, 2 * nNextPow);
		for (i = 0; i < nlags; i++)
			y1[i] = acf1[2 * i] / acf1[0];
		memset(&acf1[0], 0, 4 * nNextPow * sizeof(RealType));
	}

	t2 = clock();
	dt1 = (t2 - t1) / CLOCKS_PER_SEC;
	printf("autocorr_freq_domain : t1 = %lf\tt2 = %lf\tdt1 = %lf\n", t1, t2, dt1);

	// time
	RealType *acf2 = (RealType *)calloc(nlags, sizeof(RealType));
	t1 = clock();

	for (int c = 0; c < RunTimes; c++) {
		autocorr_time_domain(nlags - 1, n, input, acf2, 1);
		for (i = 0; i < nlags; i++)
			y2[i] = acf2[i] / acf2[0];
		memset(&acf2[0], 0, nlags * sizeof(RealType));
	}

	t2 = clock();
	dt1 = (t2 - t1) / CLOCKS_PER_SEC;
	printf("autocorr_time_domain : t1 = %lf\tt2 = %lf\tdt1 = %lf\n", t1, t2, dt1);

#if(DebugMode)
	printf("input : \n\n");
	for (i = 0; i < n; i++) 
		printf("%d\t %f\n", i, input[i]);
	printf("\n\n");

	printf("autocorr in freq, autocorr in time :\n");
	for (i = 0; i < nlags; i++) 
		printf("%d\t %f\t %f\n", i, y1[i], y2[i]);

	printf("\nerror :\n");
	for (i = 0; i < nlags; i++) 
		printf("%d\t %f\n", i, y1[i] - y2[i]);
#endif
	free(acf2);
	free(acf1);
	free(y2);
	free(y1);
	free(input);
	return(0);
}

void get_rand_sig(RealType *sig, int n) {

	for (int i = 0; i < n; i++) {
		//input[i] = (double)rand() / (RAND_MAX + 1.0); // 0~1
		sig[i] = ((RealType)rand() * 2 / RAND_MAX - 1.0); // -1~1
	}
}

void autocorr_freq_domain(RealType *data, unsigned long n) {
	
	int i;
	four1(data - 1, n, 1);
	for (i = 0; i < 2 * n; i += 2) {
		data[i] = data[i] * data[i] + data[i + 1] * data[i + 1];
		data[i + 1] = 0;
		//printf("%f\n", data[i]);
	}
	//printf("\n\n");
	four1(data - 1, n, -1);
	for (int i = 0; i < 2 * n; i += 2 ) 
		data[i] /= n;
	//for (int i = 0; i < 2 * n; i += 2) printf("%f\n", data[i]);
}

// NUMERICAL RECIPES IN C, SECOND EDITION Page.507
void four1(RealType data[], unsigned long nn, int isign)
/* Replaces data[1..2*nn] by its discrete Fourier transform,if isign is input as 1;
or replaces data[1..2*nn] by nn times its inverse discrete Fourier transform, if isgin is inputs as -1.
data is a complex array of length nn or, equivalently, a real array of length 2*nn.
nn MUST be an intrger power of 2(this is not checked for!). */
{
	/* 2023_0210_CY
	If the input data (input[]) length is n,
	have to create a new array (data[]) which data length is 2n,
	index = 0,2,4,6......2n-1 are real parts,
	index = 1,3,5,7......2n are imaginary parts.
	ex : four1(data - 1, n, 1)   fourier transform
	ex : four1(data - 1, n, -1)  inverse fourier transform, each value have to over n. */
	unsigned long n, mmax, m, j, istep, i;
	RealType wtemp, wr, wpr, wpi, wi, theta;					// Double precision for the trigonometric recurrences.
	RealType tempr, tempi;

	n = nn << 1;
	j = 1;
	for (i = 1; i < n; i += 2) {							// This is the bit-reversal section of the routine.
		if (j > i) {
			SWAP(data[j], data[i]);							// Exchange the two complex numbers.
			SWAP(data[j + 1], data[i + 1]);
		}
		m = n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	// Here begins the Danielson-Lanczos section of the routine.
	mmax = 2;
	while (n > mmax) {										// Outer loop executed log_2(nn) times.
		istep = mmax << 1;
		theta = isign * (6.28318530717959 / mmax);			// Initialize the trigonometric recurrence.
		wtemp = sin(0.5 * theta);
		wpr = -2.0 * wtemp * wtemp;
		wpi = sin(theta);
		wr = 1.0;
		wi = 0.0;
		for (m = 1; m < mmax; m += 2) {						// Here are the two nested inner loops.
			for (i = m; i <= n; i += istep) {
				j = i + mmax;								// This is the Danielson-Lanczos formula : 
				tempr = wr * data[j] - wi * data[j + 1];
				tempi = wr * data[j + 1] + wi * data[j];
				data[j] = data[i] - tempr;
				data[j + 1] = data[i + 1] - tempi;
				data[i] += tempr;
				data[i + 1] += tempi;
			}
			wr = (wtemp = wr) * wpr - wi * wpi + wr;		// Trigonometric recurrence.
			wi = wi * wpr + wtemp * wpi + wi;
		}
		mmax = istep;
	}
}

void autocorr_time_domain(int tau, int n, RealType *sig, RealType *acf, int dinc) {

	int i, j, k;
	RealType sum;
	for (i = 0, j = 0; i <= tau; i++, j += dinc) {
		sum = 0;
		for (k = 0; j + k < n; k += dinc) {
			sum += sig[k] * sig[j + k];
		}
		acf[i] = sum;
	}
}