
#include <string.h>
#include <fftw3.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>


fftw_plan fftplanin;
fftw_plan fftplanout;
fftw_complex * fft_from_input;
fftw_complex * fft_to_output;
double * raw_from_input;
double * raw_to_output;

void prepare_fftw(int inlen, int outlen) {

	raw_to_output = fftw_malloc(sizeof(double) * outlen);
	raw_from_input = fftw_malloc(sizeof(double) * inlen);

	fft_from_input = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * inlen);
	fft_to_output= (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * outlen);


	fftplanin = fftw_plan_dft_r2c_1d(inlen, raw_from_input, fft_from_input, FFTW_MEASURE);
	fftplanout = fftw_plan_dft_c2r_1d(outlen, fft_to_output, raw_to_output,  FFTW_MEASURE);
}

void do_fftw(int * input, int * output, int inlen, int outlen) {

	int i;
	double factor = (double)(outlen)/(double)(inlen);

	for(i=0; i < inlen; i++) {
		raw_from_input[i] = (double) input[i];
	//	raw_from_input[i] = ((i/10)%2 * 100);
	//	printf("%f\n", raw_from_input[i]);
	}

	fftw_execute(fftplanin);


	for(i = 0; i < outlen; i++) {
		fft_to_output[i][0] = 0;
		fft_to_output[i][1] = 0;
	}

	for (i = 0; i < inlen; i++) {
		int offset;
		offset = (int) (i * factor);
		assert(offset < outlen);
		fft_to_output[offset][0] += fft_from_input[i][0];
		fft_to_output[offset][1] += fft_from_input[i][1];
	}

	fftw_execute(fftplanout);

	for(i = 0; i < outlen; i++) {
		double c;
		//printf("%f\n", raw_to_output[i]/inlen);

		c = raw_to_output[i]/inlen;
		if (c > 0x7fffffff) {
			output[i] = 0x7fffffff;
		} else if (c < -0x7fffffff) {
			output[i] = -0x7fffffff;
		} else {
			output[i] = c;
		}
	}
}

int main(int argc, char ** argv) {

	int il = 3072;
	int ol = 768;
	int atoncei = 1024;
	int atonceo=256;
	int i;
	int inarray[il];
	int oarray[ol];
	int rv;
	prepare_fftw(il, ol);

	while (1) {
		memmove(inarray, inarray+atoncei, sizeof(int)*(il-atoncei));
		for (i = il-atoncei; i < il; i++) {
			do {
				rv = read(0, inarray + i, 4);
				rv = read(0, inarray + i, 4);
			} while (rv == 0);
			
			//printf("%i\n", inarray[i]);
		}

		do_fftw(inarray, oarray, il, ol);

		for(i = atonceo; i < (2*atonceo); i++) {
			write(1, oarray + i, 4);
			write(1, oarray + i, 4);
			write(1, oarray + i, 4);
			write(1, oarray + i, 4);
			write(1, oarray + i, 4);
			write(1, oarray + i, 4);
			write(1, oarray + i, 4);
			write(1, oarray + i, 4);
		//	printf("%i\n", oarray[i]);
		}
	}
}

