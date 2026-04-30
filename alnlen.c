#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// Define default values
#define MAX_LINE 1024

// Function prototypes
void help(int exitvalue);

//declare variables
int     i, file_index;

// main function
// This program calculates the length of the first sequence in a FASTA file
// and prints the length to the standard output.
int main(int argc, char *argv[]) {
    if (argc < 2) {
        help(1);
    }

    // ====================== Parse args ==================================
    for (i=1; i<argc; i++){
    if (!strcmp(argv[i], "-h")) help(0);
     else if (!strcmp(argv[i], "-i")) {
       if (++i > argc) help(1);
       file_index = i;
   }
}
        FILE *fp = fopen(argv[file_index], "r");
        if (!fp) {
            perror("ERROR: Not posible to open input aln");
            return 1;
        }

    char line[MAX_LINE];
    int in_sequence = 0;
    size_t length = 0;

    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '>') {
            if (in_sequence) {
                // We've already read the first sequence; stop now
                break;
            }
            in_sequence = 1;
        } else if (in_sequence) {
            // Remove newline characters
            line[strcspn(line, "\r\n")] = 0;
            length += strlen(line);
        }
    }

    fclose(fp);

    printf("%zu\n", length);
    return 0;
}

void help(int exitvalue) {
    printf("Usage: alnlen [options] \n");
    printf("Options:\n");
    printf("-h print this usage screen and exit\n");
    printf("-i \t Input aln to count length\n");
    exit(exitvalue);
  }