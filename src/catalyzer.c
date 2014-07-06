#include <assert.h>
#include <fftw3.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define MAX_INT 0x7fffffff

#define FFT_OPTION_HANN 1
typedef struct {
	/* What file descriptor to use */
	int cio_fd;

	/* How many times does each (ints) worth of data
	 * have to be repeated.
	 */
	int cio_duplication;

	/* Should we give up if you run into a problem with IO?
	 * This should be false when working with files but true
	 * when working with audio devices because data might not
	 * be ready yet.
	 */
	int cio_stop_on_eof;
} catalyzer_io;

typedef struct fftw_holder {
	fftw_plan fftplanin;
	fftw_plan fftplanout;
	fftw_complex * fft_from_input;
	fftw_complex * fft_to_output;
	double * raw_from_input;
	double * raw_to_output;
} catalyzer_fftw_state;




/*
 * readwrite: This manages low-level read and write stuff. 
 * cio is a catalyzer_io structure that specfies how to read or write. 
 * data is a pointer to an integer array where the data will go
 * len is thenumber of ints to copy
 * shouldwrite is 1 if you are writing data and 0 if you are reading data
 * */
int readwrite(catalyzer_io * cio, int * data, int len, int shouldwrite) {

	int i,j;

	for (i = 0; i < len; i++) {
		int * dptr;
		dptr = data + i;
		for (j = 0; j < cio->cio_duplication; j++) {
			int rv;
			/* Keep trying until we read/write the data we want */
			do {
				if (shouldwrite) {
					rv = write(cio->cio_fd, dptr, 4);
				} else {
					rv = read(cio->cio_fd, dptr, 4);
				}
			} while (rv == 0 && !cio->cio_stop_on_eof);

			if (rv < 1 && cio->cio_stop_on_eof) {
				fprintf(stderr, "Giving up on fd %i", cio->cio_fd);
				return 1;
			}

		}
	}
	return 0;
}

/*
 * Prepare FFTW
 */
catalyzer_fftw_state * prepare_fftw(int inlen, int outlen) {

	catalyzer_fftw_state * mydata = NULL;

	mydata = (catalyzer_fftw_state *) malloc(sizeof(catalyzer_fftw_state));
	mydata->raw_to_output = fftw_malloc(sizeof(double) * outlen);
	mydata->raw_from_input = fftw_malloc(sizeof(double) * inlen);

	mydata->fft_from_input = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * inlen);
	mydata->fft_to_output = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * outlen);

	mydata->fftplanin = fftw_plan_dft_r2c_1d(inlen, mydata->raw_from_input, mydata->fft_from_input, FFTW_MEASURE);
	mydata->fftplanout = fftw_plan_dft_c2r_1d(outlen, mydata->fft_to_output, mydata->raw_to_output,  FFTW_MEASURE);

	return mydata;
}


/*
 * Convert inlen samples from input to outlen samples to output, scaling frequences
 * as needed.
 */
void do_fftw(catalyzer_fftw_state * mydata, int * input, int * output, int inlen, int outlen, int options) {

	int i;
	double factor = (double)(outlen)/(double)(inlen);

	// STEP 1: FFTW operates on doubles, so convert to a double.
	for(i=0; i < inlen; i++) {
		mydata->raw_from_input[i] = (double) input[i];
	}

	if (options & FFT_OPTION_HANN) {
		// Hann window
		for (i = 0; i < inlen; i++) {
			double scale;
			scale = 0.5 * (1 - cos((2 * M_PI * i)/(inlen - 1)));
			mydata->raw_from_input[i] *= scale;
		}
	}


	// STEP 2: Compute the forward FFT
	fftw_execute(mydata->fftplanin);


	// STEP 3: Clear the output space.
	for(i = 0; i < outlen; i++) {
		mydata->fft_to_output[i][0] = 0;
		mydata->fft_to_output[i][1] = 0;
	}

	// STEP 4: scale frequencies. 
	for (i = 0; i < inlen; i++) {
		int offset;
		offset = (int) (i * factor);
		assert(offset < outlen);
		mydata->fft_to_output[offset][0] += mydata->fft_from_input[i][0];
		mydata->fft_to_output[offset][1] += mydata->fft_from_input[i][1];
	}


	// STEP 5: inverse FFT
	fftw_execute(mydata->fftplanout);

	// STEP 6: normalize output/
	for(i = 0; i < outlen; i++) {
		double c;

		// FIXME: not quite the right normalization constant here...
		c = mydata->raw_to_output[i]/inlen;

		// Prevent overflow.
		if (c > MAX_INT) {
			output[i] = MAX_INT;
		} else if (c < -MAX_INT) {
			output[i] = -MAX_INT;
		} else {
			output[i] = c;
		}
	}
}

/*
 * Main loop continuously reads/writes from input/output, downsampleing as we go.
 * input: data inputinformation
 * output: data output information
 * downsample_factor how much to downsample by (usually 4)
 * overlap: how many samples to overlap between subsequent samples
 */
int main_loop(catalyzer_io * input, 
		catalyzer_io * output, 
		int downsample_factor,
		int sample_size,
		int overlap) {

	if (sample_size % downsample_factor) {
		fprintf(stderr,"downsample factor does not evenly divide sample size.\n");
		return -1;
	}

	int input_overlap = overlap * downsample_factor * 2;
	int output_overlap = overlap * 2;
	int input_buffer_size = sample_size + input_overlap;
	int output_buffer_size = sample_size / downsample_factor  + output_overlap;
	int * input_buffer = calloc(input_buffer_size, sizeof(int));
	int * output_buffer = calloc(output_buffer_size, sizeof(int));

	catalyzer_fftw_state * fftw_state = prepare_fftw(input_buffer_size, output_buffer_size);

	int ok = 1;

	while (ok) {
		/* Copy some data from the tail of the input buffer to the begining, 
		 * as per overlap requirements 
		 */
		memmove(input_buffer, 
			input_buffer+ (input_buffer_size - input_overlap) * sizeof(int),
			input_overlap * sizeof(int));

		/* Read input data */
		readwrite(input, input_buffer + input_overlap, sample_size, 0);
		
		/* Do the converting */
		do_fftw(fftw_state, input_buffer, output_buffer, input_buffer_size, output_buffer_size,0 );

		/* Write output data */
		readwrite(output, output_buffer + overlap, sample_size / downsample_factor, 1);
	}
}


int main(int argc, char ** argv) {

	catalyzer_io input, output;
	int downsample = 4;

	input.cio_fd=0;
	input.cio_duplication=2;
	input.cio_stop_on_eof=1;

	output.cio_fd=1;
	output.cio_duplication=input.cio_duplication * downsample;
	output.cio_stop_on_eof=1;


	main_loop(&input, &output, downsample, 8192, 0);
}

