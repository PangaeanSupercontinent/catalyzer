
#include <string.h>
#include <fftw3.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>


/*
 * TODO: move all of this stuff into a struct that we pass around to 'helper' functions
 */
fftw_plan fftplanin;
fftw_plan fftplanout;
fftw_complex * fft_from_input;
fftw_complex * fft_to_output;
double * raw_from_input;
double * raw_to_output;

/*
 * Prepare FFTW
 */
void prepare_fftw(int inlen, int outlen) {

	raw_to_output = fftw_malloc(sizeof(double) * outlen);
	raw_from_input = fftw_malloc(sizeof(double) * inlen);

	fft_from_input = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * inlen);
	fft_to_output= (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * outlen);


	fftplanin = fftw_plan_dft_r2c_1d(inlen, raw_from_input, fft_from_input, FFTW_MEASURE);
	fftplanout = fftw_plan_dft_c2r_1d(outlen, fft_to_output, raw_to_output,  FFTW_MEASURE);
}

/*
 * Convert inlen samples from input to outlen samples to output, scaling frequences
 * as needed.
 */
void do_fftw(int * input, int * output, int inlen, int outlen) {

	int i;
	double factor = (double)(outlen)/(double)(inlen);

	// STEP 1: FFTW operates on doubles, so convert to a double.
	for(i=0; i < inlen; i++) {
		raw_from_input[i] = (double) input[i];
	}

	// STEP 2: Compute the forward FFT
	fftw_execute(fftplanin);


	// STEP 3: Clear the output space.
	for(i = 0; i < outlen; i++) {
		fft_to_output[i][0] = 0;
		fft_to_output[i][1] = 0;
	}

	// STEP 4: scale frequencies. 
	for (i = 0; i < inlen; i++) {
		int offset;
		offset = (int) (i * factor);
		assert(offset < outlen);
		fft_to_output[offset][0] += fft_from_input[i][0];
		fft_to_output[offset][1] += fft_from_input[i][1];
	}

	// STEP 5: inverse FFT
	fftw_execute(fftplanout);

	// STEP 6: normalize output/
	for(i = 0; i < outlen; i++) {
		double c;

		// FIXME: not quite the right normalization constant here...
		c = raw_to_output[i]/inlen;

		// Prevent overflow.
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

	/* FIXME: This is a horrible way to get data! */
	while (1) {
		/* IDEA: To prevent noise from sharp transations, let each "chunk" of 
 		 * data to convert overlap a bit with the previous chunk. We'll split into
		 * 3 parts and do the FFT over the whole thing but only extract data from
		 * the center. 
		 * hmmmmm. In practise this doesn't seem to work so well.
		 */

		// Slide data we have over.
		memmove(inarray, inarray+atoncei, sizeof(int)*(il-atoncei));

		// Read more data from STDIN
		for (i = il-atoncei; i < il; i++) {
			do {
				// Right Right and Left... Throwaway Right
				rv = read(0, inarray + i, 4);
				rv = read(0, inarray + i, 4);
			} while (rv == 0);
			
			//printf("%i\n", inarray[i]);
		}

		// Convert the data.
		do_fftw(inarray, oarray, il, ol);

		// Output the data. Since we use a 1:4 ratio, we need to output
		// four samples for every sample we recieve since teh input and output
		// clocks have the same frequency. 
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

