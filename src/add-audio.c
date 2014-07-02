#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#define MAX_INT 0x7fffffff

int main(int argc, char ** argv) {

	int mine, input1, input2;
	FILE *file1, *file2;

	if (argc < 3) {
		printf("supply <file1> <file2> \n");
		return -1;
	}


	file1 = fopen(argv[1], "r");
	file2 = fopen(argv[2], "r");

	while ((fread(&input1, 4, 1, file1) != 0) && (fread(&input2, 4, 1, file1) != 0)) {
		mine = 0.25 * (input1 + input2);
		write(1, &mine, 4);
	}

}
