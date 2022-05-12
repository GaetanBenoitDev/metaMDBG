#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "abpoa.h"
#include "abpoa_graph.h"
#include "abpoa_align.h"
#include "abpoa_seq.h"
#include "utils.h"

char NAME[20] = "abPOA";
char PROG[20] = "abpoa";
#define _ba BOLD UNDERLINE "a" NONE
#define _bb BOLD UNDERLINE "b" NONE
#define _bP BOLD UNDERLINE "P" NONE
#define _bO BOLD UNDERLINE "O" NONE
#define _bA BOLD UNDERLINE "A" NONE
char DESCRIPTION[100] = _ba "daptive " _bb "anded " _bP "artial " _bO "rder " _bA "lignment";
char VERSION[20] = "1.4.0";
char CONTACT[30] = "gaoy286@mail.sysu.edu.cn";

const struct option abpoa_long_opt [] = {
    { "align-mode", 1, NULL, 'm' },

    { "match", 1, NULL, 'M' },
    { "mismatch", 1, NULL, 'X' },
    { "matrix", 1, NULL, 't' },
    { "gap-open", 1, NULL, 'O' },
    { "gap-ext", 1, NULL, 'E' },

    { "extra-b", 1, NULL, 'b' },
    { "extra-f", 1, NULL, 'f' },
    { "zdrop", 1, NULL, 'z' },
    { "bouns", 1, NULL, 'e' },

    { "seeding", 0, NULL, 'S'},
    { "k-mer", 1, NULL, 'k' },
    { "window", 1, NULL, 'w' },
    { "min-poa-win", 1, NULL, 'n' },
    { "progressive", 0, NULL, 'p'},

    { "amino-acid", 0, NULL, 'c'},
    { "in-list", 0, NULL, 'l' },
    { "increment", 1, NULL, 'i' },
    

    { "amb-strand", 0, NULL, 's' },
    { "output", 1, NULL, 'o' },
    { "result", 1, NULL, 'r' },
    { "out-pog", 1, NULL, 'g' },
    { "max-num-cons", 1, NULL, 'd', },
    { "min-freq", 1, NULL, 'q', },

    { "help", 0, NULL, 'h' },
    { "version", 0, NULL, 'v' },

    { 0, 0, 0, 0}
};

int abpoa_usage(void)
{
    err_printf("\n");
    err_printf("%s: %s \n\n", PROG, DESCRIPTION);
    err_printf("Version: %s\t", VERSION);
	err_printf("Contact: %s\n\n", CONTACT);
    err_printf("Usage: %s [options] <in.fa/fq> > cons.fa/msa.out/abpoa.gfa\n\n", PROG);
    err_printf("Options:\n");
    err_printf("  Alignment:\n");
    err_printf("    -m --aln-mode   INT     alignment mode [%d]\n", ABPOA_GLOBAL_MODE);
    err_printf("                              %d: global, %d: local, %d: extension\n", ABPOA_GLOBAL_MODE, ABPOA_LOCAL_MODE, ABPOA_EXTEND_MODE);
    err_printf("    -M --match      INT     match score [%d]\n", ABPOA_MATCH);
    err_printf("    -X --mismatch   INT     mismatch penalty [%d]\n", ABPOA_MISMATCH);
    err_printf("    -t --matrix    FILE     scoring matrix file, \'-M\' and \'-X\' are not used when \'-t\' is used [Null]\n");
    err_printf("                            e.g., \'HOXD70.mtx, BLOSUM62.mtx\'\n");
    err_printf("    -O --gap-open INT(,INT) gap opening penalty (O1,O2) [%d,%d]\n", ABPOA_GAP_OPEN1, ABPOA_GAP_OPEN2);
    err_printf("    -E --gap-ext  INT(,INT) gap extension penalty (E1,E2) [%d,%d]\n", ABPOA_GAP_EXT1, ABPOA_GAP_EXT2);
    err_printf("                            %s provides three gap penalty modes, cost of a g-long gap:\n", NAME);
    err_printf("                            - convex (default): min{O1+g*E1, O2+g*E2}\n");
    err_printf("                            - affine (set O2 as 0): O1+g*E1\n");
    err_printf("                            - linear (set O1 as 0): g*E1\n");
    err_printf("    -s --amb-strand         ambiguous strand mode [False]\n");
    err_printf("                            for each input sequence, try the reverse complement if the current\n");
    err_printf("                            alignment score is too low, and pick the strand with a higher score\n");
    err_printf("  Adaptive banded DP:\n");
    err_printf("    -b --extra-b    INT     first adaptive banding parameter [%d]\n", ABPOA_EXTRA_B);
    err_printf("                            set b as < 0 to disable adaptive banded DP\n");
    err_printf("    -f --extra-f  FLOAT     second adaptive banding parameter [%.2f]\n", ABPOA_EXTRA_F);
    err_printf("                            the number of extra bases added on both sites of the band is\n");
    err_printf("                            b+f*L, where L is the length of the aligned sequence\n");
    // err_printf("    -z --zdrop    INT       Z-drop score in extension alignment [-1]\n");
    // err_printf("                            set as <= 0 to disable Z-drop extension\n");
    // err_printf("    -e --bonus    INT       end bonus score in extension alignment [-1]\n");
    // err_printf("                            set as <= 0 to disable end bounus\n");
    err_printf("  Minimizer-based seeding and partition (only effective in global alignment mode):\n");
    err_printf("    -S --seeding            enable minimizer-based seeding and anchoring [False]\n");
    err_printf("    -k --k-mer       INT    minimizer k-mer size [%d]\n", ABPOA_MMK);
    err_printf("    -w --window      INT    minimizer window size [%d]\n", ABPOA_MMW);
    err_printf("    -n --min-poa-win INT    min. size of window to perform POA [%d]\n", ABPOA_MIN_POA_WIN);
    err_printf("    -p --progressive        build guide tree and perform progressive partial order alignment [False]\n");
    // err_printf("    -n --par-size           minimal partition size [%d]\n", ABPOA_W);

    err_printf("  Input/Output:\n");
    err_printf("    -c --amino-acid         input sequences are amino acid (default is nucleotide) [False]\n");
    err_printf("    -l --in-list            input file is a list of sequence file names [False]\n");
    err_printf("                            each line is one sequence file containing a set of sequences\n");
    err_printf("                            which will be aligned by abPOA to generate a consensus sequence\n");
    err_printf("    -i --incrmnt    FILE    incrementally align sequences to an existing graph/MSA [Null]\n");
    err_printf("                            graph could be in GFA or MSA format generated by abPOA\n");
    err_printf("    -o --output     FILE    ouput to FILE [stdout]\n");
    err_printf("    -r --result      INT    output result mode [%d]\n", ABPOA_OUT_CONS);
    err_printf("                            - %d: consensus in FASTA format\n", ABPOA_OUT_CONS);
    err_printf("                            - %d: MSA in PIR format\n", ABPOA_OUT_MSA);
    err_printf("                            - %d: both 0 & 1\n", ABPOA_OUT_CONS_MSA);
    err_printf("                            - %d: graph in GFA format\n", ABPOA_OUT_GFA);
    err_printf("                            - %d: graph with consensus path in GFA format\n", ABPOA_OUT_CONS_GFA);
    err_printf("                            - %d: consensus in FASTQ format\n", ABPOA_OUT_CONS_FQ);
    err_printf("    -d --maxnum-cons INT    max. number of consensus sequence to generate [1]\n");
    err_printf("    -q --min-freq  FLOAT    min. frequency of each consensus sequence (only effective when -d/--num-cons > 1) [%.2f]\n", MULTIP_MIN_FREQ);
    err_printf("    -g --out-pog    FILE    dump final alignment graph to FILE (.pdf/.png) [Null]\n\n");

    err_printf("    -h --help               print this help usage information\n");
    err_printf("    -v --version            show version number\n");


    err_printf("\n");
    return 1;
}

int abpoa_main(char *file_fn, int is_list, abpoa_para_t *abpt){
    double realtime0 = realtime();
    // TODO abpoa_init for each input file ???
    abpoa_t *ab = abpoa_init();
    if (is_list) { // input file list
        FILE *list_fp = fopen(file_fn, "r"); char read_fn[1024];
        while (fgets(read_fn, sizeof(read_fn), list_fp)) {
            read_fn[strlen(read_fn)-1] = '\0';
            abpoa_msa1(ab, abpt, read_fn, stdout);
        }
        fclose(list_fp);
    } else // input file
        abpoa_msa1(ab, abpt, file_fn, stdout);

    abpoa_free(ab);
	err_func_printf(__func__, "Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB.", realtime() - realtime0, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);
    return 0;
}


