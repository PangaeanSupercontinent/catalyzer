#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#define MAX_INT 0x7fffffff

int main(int argc, char ** argv) {

	if (argc < 4) {
		printf("supply <frequency>, <sample-rate>, and <timeout> \n");
		return -1;
	}

	int freeq = atoi(argv[1]);
	double sample_rate = strtod(argv[2], NULL);
	int timeout = atoi(argv[3]);
	int index;
	double increment;
	double current;
	int value;

	alarm(timeout);
	
	value = 0;
	current = 0;
	increment = 2.0 * M_PI * (double)freeq / sample_rate;

	while (1) {
		value = (sin(current) + 1) * MAX_INT/2;  // scale to int
		fprintf(stderr, "%i \n", value);
		write(1, &value, 4);
		current += increment;
		if (current> 2*M_PI) {
			current-= 2*M_PI;
		}
	}
}
