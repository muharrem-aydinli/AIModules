// ========================================================================
// Computational Biology and Informatics Laboratory (CBIL) (www.cbil.upenn.edu)
// University of Pennsylvania Center for Bioinformatics (www.pcbi.upenn.edu)
// ========================================================================

// ========================================================================
// The CBIL Software License, Version 1.0
// ========================================================================

/* Copyright (c) 2013 The Computational Biology and Informatics Laboratory (CBIL) 
   in the University of Pennsylvania.  All rights reserved.

   Redistribution and use in source and binary forms, with or without modification, 
   are permitted provided that the following conditions are met:

      1. Redistributions of source code must retain the above copyright  notice, 
         this list of conditions and the following disclaimer.

      2. Redistributions in binary form must reproduce the above copyright  notice, 
         this list of conditions and the following disclaimer in the documentation
         and/or other materials provided with the distribution.

      3. The end-user documentation included with the redistribution, if any, must 
         include the following acknowledgment: 

 "This product includes software developed by CBIL at the Center for Bioinformatics at the University of Pennsylvania." 

         Alternately, this acknowledgment may appear in the software itself, if 
	 and wherever such third-party acknowledgments normally appear.

      4. The names "CBIL" and "Penn Center for Bioinformatics" must not be used to
         endorse or promote products derived from this software without prior 
	 written permission.

   THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESSED OR IMPLIED 
   WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
   OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
   DISCLAIMED.  IN NO EVENT SHALL CBIL OR THE PENN CENTER FOR 
   BIOINFORMATICS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF 
   USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND 
   ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
   OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
   OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF 
   SUCH DAMAGE. 
   ====================================================================

   This license is based on the open source license from the Apache Software 
   Foundation  http://www.apache.org.

*/

// ========================================================================
// tessWms
// ========================================================================

/* 

Transcription Element Search System: Partial Weight Matrices (tessWms) 
is a program to read (selected) PWMs from a file and predict
binding sites on DNA sequences read from another file.  Hits that have
a good enough score are reported in one of four different formats.

As the program is farily simple and has a limited set of requirements,
it is written in a fairly simple fashion that includes using a large
set of global variables that contain the user's arguments, the current
sequence and a list of the PWMs to be processed.  The goal was to keep
the program as simple as possible to preserve good performance and
maintainability across multiple platforms.  The also tries to minimize
the amount of exra work that is done to keep the speed fairly high.

*/

// ========================================================================
// ------------------------------- Includes -------------------------------
// ========================================================================

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>

// ========================================================================
// -------------------------- Macro Definitions ---------------------------
// ========================================================================

char *G_version = "1.1";

// #define __DEBUG__ 1

// ---------------------------- Storage Sizes -----------------------------

// length of longest sequence or matrix title
#define MAX_TTL_LENGTH 32767

// length of longest DNA sequence
#define MAX_SEQ_LENGTH 250000000

// number of positions in longest PWM
#define MAX_MTX_LENGTH 4096

// maximum number of PWMs to store
#define MAX_MTX_COUNT  50000

// ----------------------------- DNA Encoding -----------------------------

/*

As sequences are read in they are stored in both a human readable form
and a binary internal form for rapid processing.  The program is
capable of handling and scoring the IUPAC ambiguous base codes.
Therefore we use a bit-based code.  The macros below define the bit
code for the un ambiguous bases.  Ambiguous bases are represented by
having multiple bits set in a value.

*/

#define DNA_ERROR -1
#define DNA_NONE   0
#define DNA_A      1
#define DNA_C      2
#define DNA_G      4
#define DNA_T      8
#define DNA_ALL    (DNA_A|DNA_C|DNA_G|DNA_T)
#define DNA_SIZE   (DNA_ALL+1)

// ----------------------------- Comparisons ------------------------------

#define MIN(a,b)   ((a) < (b) ? (a) : (b))
#define MAX(a,b)   ((a) > (b) ? (a) : (b))

// -------------------------------- Scores --------------------------------

#define VERY_BAD_SCORE -1e6

// ========================================================================
// ------------------------ Function Declarations -------------------------
// ========================================================================

int  main              ( int, char *[] );

// ------------------------ Command-line Arguments ------------------------

void cla               ( int, char *[] );

void claParseArgs      ( int, char *[] );

void claErrorCheck     ( void);

void claSetDependentValues ( void );

void usage             ( void );

// ------------------------- DNA Sequence Loading -------------------------

char encode            ( char );

char decode            ( char );

int  loadNextSequence  ( FILE * );

int  loadTitle         ( FILE * );

void loadSequence      ( FILE * );

void updateLogOddsDist ( ) ;

// ------------------- Matrix Scoring and Manipulation --------------------

void loadAllMatrices   ( void );

void loadNextMatrix    ( FILE * );

void scoreWithMatrices ( );

void scoreOneMatrix    ( );

void flipMatrix        ( );

void makeProbs         ( );

void expandMatrix      ( );

double maxScore        ( );

// ----------------------- Ambiguous Base Treatment -----------------------

void clearTwo        ( int , int , int );

void clearThree      ( int , int , int , int );

void clearFour       ( int , int , int , int, int );

void averageTwo        ( int , int , int );

void averageThree      ( int , int , int , int );

void averageFour       ( int , int , int , int, int );

void minimumTwo        ( int , int , int );

void minimumThree      ( int , int , int , int );

void minimumFour       ( int , int , int , int, int );

void maximumTwo        ( int , int , int );

void maximumThree      ( int , int , int , int );

void maximumFour       ( int , int , int , int, int );

// -------------------------- Outputing Results ---------------------------

void outputHead        ( );

void outputSeqHead     ( );

void outputMatHead     ( int counterForMatrix );

void outputMatTail     ( );

void outputSeqTail     ( );

void outputTail        ( );

// ========================================================================
// --------------------- Global Command Line Options ----------------------
// ========================================================================

/*

The values of these variables is controlled by the command line
arguments.  The code in claParseArgs() should provide away for the user
to manipulate all of these values.  The help in usage() should
accurately reflect the default values set here.

*/

// there was an error parsing the args
int    G_CLA_error_b       = 0;

// remain silent
int    G_CLA_Silent_b      = 0;

// read sequence from this file
char * G_CLA_Sequence_f    = "-";

// read matrices from this file
char * G_CLA_Matrix_f      = NULL;

// write output to this file
char * G_CLA_Output_f      = "-";

// select the mats that has this id
char * G_CLA_SelectId      = NULL;

// select mats that contain this string.
char * G_CLA_Select        = NULL;

// score threshold
double G_CLA_MinScore      = 6.0;

// score Lmax - Lobs
double G_CLA_MaxDeficit    = 100.0;

// start position
int G_CLA_Start            = 0;

// end position
int G_CLA_End              = 1000000000 ;

// do not consider sites on the reverse complement strand
int G_CLA_NoRc             = 0;

// do not normalize scores
int G_CLA_NoNormalize      = 0;

// total pseudocounts
double G_CLA_PseudoCounts  = 1.0;

// position offset
int G_CLA_PositionOffset   = 0;

// format
char G_CLA_Format          = 'w';

// GFF feature
char * G_CLA_GffFeature    = NULL;

// segmented output
char G_CLA_SegmentedOutput = 0;

// IUPAC abiguous base treatment
char G_CLA_IupacTreatment  = 'p';

// use sequence as background
char G_CLA_UseSequenceBackground = 0;

// log-odds background distribution
double G_CLA_LogOddsBkg[]  = { 0.25, 0.25, 0.25, 0.25 };

// ========================================================================
// --------------------------- Current Sequence ---------------------------
// ========================================================================

/*

As sequences are read in from the sequence file, the information is
placed in these variables which can be accessed by all functions.

*/

// sequence title
char  G_sequenceTitle [MAX_TTL_LENGTH+1];

// sequence id
char  G_sequenceId    [MAX_TTL_LENGTH+1];

// sequence (human)
char  G_sequence_c    [MAX_SEQ_LENGTH+1];

// sequence (binary)
char  G_sequence_b        [MAX_SEQ_LENGTH+1];

// sequence ordinal
int   G_sequence_o;

// sequence length
int   G_sequenceLength_n = 0;

// ========================================================================
// ---------------- Weight Matrices - Loading and Current -----------------
// ========================================================================

/* 

The matrix iteration is the inner loop so consequently we try to speed
this up by loading the matrices at the start of the program.  These
variables hold a list of PWMs as well as some temporary locations to
hold the PWMs as they are read in.  The temporary locations are
staging areas; the info is moved into the list ASAP.

One has to be careful to note which matrix manipulation functions use
the temporary store or a pointer to the current list element.

*/

// ------------------------------- PWM List -------------------------------

// PWM typedef
typedef struct {
	double ** matrix;
  double ** counts;
	char   *  matrixTitle;
	char   *  matrixId;
	int    matrix_n;
	int    matrix_o;
} matrix_T;

// array of matrices that have been read in.
matrix_T   G_matrices[MAX_MTX_COUNT];

// number of matrices read in
int        G_matrixCache_n = 0;

matrix_T * G_theMatrix;

// ------------------------ PWM Temporary Storage -------------------------

// matrix title
char  G_matrixTitle[MAX_TTL_LENGTH+1];

// matrix
double     G_matrix[MAX_MTX_LENGTH][DNA_SIZE];

// matrix id
char     G_matrixId[MAX_TTL_LENGTH+1];

// matrix length
int      G_matrix_n;

// matrix ordinal
int G_matrix_o = 0;

// ========================================================================
// ----------------------- Secondary Global Values ------------------------
// ========================================================================

// ......................................................................

// do we need to compute the oligo for each hit?  Depends on the
// output format.
int    G_computeOligo_b = 0;

// log-2 lob
double G_logBkg[4];

// write all non-error output to this file.
FILE * G_out_F = NULL;

// ========================================================================
// --------------------------- CODE STARTS HERE ---------------------------
// ========================================================================

// --------------------------------- cla ----------------------------------

void cla ( int    Arg_n,
					 char * Args[]
					 )
{

	claParseArgs(Arg_n, Args);
	claErrorCheck();
	claSetDependentValues();
}

// ----------------------------- claParseArgs -----------------------------

void claParseArgs ( int    Arg_n,
										char * Args[]
										)
{
	int i;

	for (i = 1; i < Arg_n; i++) {

		if (strcmp(Args[i], "-seq") == 0) {
			G_CLA_Sequence_f = Args[++i];
		}

		else if (strcmp(Args[i], "-mat") == 0) {
			G_CLA_Matrix_f = Args[++i];
		}

		else if (strcmp(Args[i], "-out") == 0) {
			G_CLA_Output_f = Args[++i];
		}

		else if (strcmp(Args[i], "-pat") == 0) {
			G_CLA_Select = Args[++i];
		}

		else if (strcmp(Args[i], "-id") == 0) {
			G_CLA_SelectId = Args[++i];
		}

		else if (strcmp(Args[i], "-mlo") == 0) {
			G_CLA_MinScore = atof(Args[++i]);
		}

		else if (strcmp(Args[i], "-mxd") == 0) {
			G_CLA_MaxDeficit = atof(Args[++i]);
		}

		else if (strcmp(Args[i], "-beg") == 0) {
			G_CLA_Start = atoi(Args[++i]);
		}

		else if (strcmp(Args[i], "-end") == 0) {
			G_CLA_End = atoi(Args[++i]);
		}

		else if (strcmp(Args[i], "-nrc") == 0) {
			G_CLA_NoRc = 1;
		}

		else if (strcmp(Args[i], "-nno") == 0) {
			G_CLA_NoNormalize = 1;
		}

		else if (strcmp(Args[i], "-tpc") == 0) {
			G_CLA_PseudoCounts = atof(Args[++i]);
		}

		else if (strcmp(Args[i], "-po") == 0) {
			G_CLA_PositionOffset = atoi(Args[++i]);
		}

		else if (strcmp(Args[i], "-abt") == 0) {
			G_CLA_IupacTreatment = Args[++i][0];
		}

		else if (strcmp(Args[i], "-fmt") == 0) {
			G_CLA_Format = Args[++i][0];
		}

    else if (strcmp(Args[i], "-gf") == 0) {
      G_CLA_GffFeature = Args[++i];
    }

    else if (strcmp(Args[i], "-so") == 0) {
      G_CLA_SegmentedOutput = 1;
    }

    else if (strcmp(Args[i], "-lob") == 0) {
      if (sscanf(Args[++i],
                 "%lf,%lf,%lf,%lf",
                 &G_CLA_LogOddsBkg[0],
                 &G_CLA_LogOddsBkg[1],
                 &G_CLA_LogOddsBkg[2],
                 &G_CLA_LogOddsBkg[3]
                 ) != 4) {
        fprintf(stderr, "-lob value format is float,float,float,float, i.e., 0.3,0.2,0.2,0.3\n");
        G_CLA_error_b = 1;
      }
    }

    else if (strcmp(Args[i], "-usb") == 0) {
      G_CLA_UseSequenceBackground = 1;
    }

    else if (strcmp(Args[i], "-q") == 0) {
      G_CLA_Silent_b = 1;
    }

		else {
			fprintf(stderr, "Unrecognized option '%s'.\n", Args[i]);
			G_CLA_error_b = 1;
		}
	}
}

// ---------------------------- claErrorCheck -----------------------------

void claErrorCheck ( void ) {

  // we need a matrix file.
	if (G_CLA_Matrix_f == NULL) {
		fprintf(stderr, "You must use the -mat option.\n");
		G_CLA_error_b = 1;
	}

	if(strcmp(G_CLA_Output_f,"-") == 0
		&& G_CLA_Format == 'j')
	{
		fprintf(stderr, "'-fmt j' is only possible with '-out'.\n");
		G_CLA_error_b = 1;
	}

  // legal formats
	switch (G_CLA_Format) {
	case 'b':	case 'g': case 'r':	case 't':	case 'w':	case 'j':	break;

	default:
		fprintf(stderr,
						"Illegal -fmt option '%s'; should be one of t (tab), b (bulk), g (GFF), r (ResultSet), j (json) or x (XML)\n",
						G_CLA_Format
						);
		G_CLA_error_b = 1;
		break;
	}


  // legal ABT values
	switch (G_CLA_IupacTreatment) {
	case 'o':	case 'p':	case 'a':	case 'n': break;

	default:
		fprintf(stderr,
						"Illegal -abt option '%s'; should be one of o (optimistic), p (pessimistic), a (average), n (none)\n",
						G_CLA_IupacTreatment
						);
		G_CLA_error_b = 1;
		break;
	}

  // -lob values must be strictly positive
  if ( G_CLA_LogOddsBkg[0] <= 0 ||
       G_CLA_LogOddsBkg[1] <= 0 ||
       G_CLA_LogOddsBkg[2] <= 0 ||
       G_CLA_LogOddsBkg[3] <= 0
       ) {
    fprintf(stderr,
            "-lob values must be greater than zero\n"
            );
    G_CLA_error_b = 1;
  }
       
  // bail out if there is any problem.
  if (G_CLA_error_b) usage();
}

// ------------------------ claSetDependentValues -------------------------

void claSetDependentValues ( void ) {

  int i;

  double lobTotal_d;

	// we only need to compute oligos for certain output formats.
	switch (G_CLA_Format) {
	case 't': G_computeOligo_b = 1;	break;
	case 'b':	G_computeOligo_b = 0;	break;
  case 'g': G_computeOligo_b = 0; break;
	case 'r':	G_computeOligo_b = 0;	break;
	case 'w':	G_computeOligo_b = 1;	break;
	case 'j':	G_computeOligo_b = 1;	break;
	}

  // normalize the log-odd background
  lobTotal_d =
    G_CLA_LogOddsBkg[0] +
    G_CLA_LogOddsBkg[1] +
    G_CLA_LogOddsBkg[2] +
    G_CLA_LogOddsBkg[3];

  for (i = 0; i < 4; i++) {
    G_CLA_LogOddsBkg[i] /= lobTotal_d;
    G_logBkg[i]          = log(G_CLA_LogOddsBkg[i]) / log(2.0);
  }
  
	// convert boundaries to 0-based for speed in checking.
	G_CLA_Start--;
	G_CLA_End--;

}

// -------------------------------- usage ---------------------------------

/*

Be careful to keep this text in sync with the code above!

*/

void usage ( void ) {

  printf("tessWms [OPTIONS]\n");
  printf("  -seq string : FASTA file with one or more sequences\n");
  printf("                default is '-', i.e., stdin.\n");
  printf("  -mat string : integer matrix file\n");
  printf("                required.\n");
  printf("  -out string : write output to this file.\n");
	printf("                default is -, i.e., stdout.\n");
	printf("  -pat string : only process matrices that contain this string in the defline.\n");
	printf("                default is all PWMs.\n");
  printf("  -id  string : only process matrices that have this ID (first word on defline.)\n");
	printf("                default is all PWMs.\n");
  printf("  -mlo float  : only report hits with a log-odds score better than this\n");
	printf("                default is 6.0.\n");
  printf("  -mxd float  : only report hits with a log-odd score not more that this amount less than best possible score for PWM\n");
  printf("                default is 100.0\n");
  printf("  -beg int    : report hits 3' of this 1-based position\n");
	printf("                default is all positions.\n");
  printf("  -end int    : report hits 5' of this 1-based position\n");
	printf("                default is all positions.\n");
  printf("  -nrc        : do not consider reverse complement strand\n");
	printf("                default is both orientations.\n");
  printf("  -nno        : do not normalize counts to make probabilities\n");
  printf("                default is to compute probabilities\n");
  printf("  -tpc float  : total pseudocounts added to counts when computing log-odd scores; negative is sqrt(n)\n");
	printf("                default is 1.0.\n");
  printf("  -lob p x 4  : probabilities of each base for log-odds calculation\n");
  printf("                default is uniform\n");
  printf("  -usb        : use sequence composition as background\n");
  printf("                default is use uniform or -lob\n");
  printf("  -po  int    : add this amount to reported positions\n");
	printf("                default is 0.\n");
  printf("  -abt char   : ambiguous base treatment: a (average), o (optimistic), p (pessimistic), or n (none)\n");
	printf("                default is pessimistic.\n");
	printf("  -fmt char   : output format: b (bulk), t (tab), g (GFF), r (ResultSet-XML), w (WMS-XML), j (json).\n");
	printf("                default is WMS-XML.\n");
  printf("  -so         : segment the output at each sequence.\n");
	printf("                default is no segmentation.\n");
  printf("  -q          : remain silent expect to report results.\n");
  printf("                default is not to be silent and to report progress on stderr.\n");
  printf("\n\n");

	//printf("Code Version: $Revision $\n");

  exit(0);
}

// ----------------------------------------------------------------------

int main ( int    Arg_n,
					 char * Args[]
           )
{

  // sequence file
  FILE * seq_F;

  // process the user's arguments
	cla(Arg_n, Args);

  // open the sequence file
  seq_F = strcmp(G_CLA_Sequence_f,"-") == 0 ? stdin : fopen(G_CLA_Sequence_f,"r");
  if (seq_F == NULL) {
    fprintf(stderr,"Could not open sequence file: '%s'\n", G_CLA_Sequence_f);
    exit(-1);
  }
	else if (!G_CLA_Silent_b) {
		fprintf(stderr,"Opened sequence file '%s'\n", G_CLA_Sequence_f);
	}

  // open the output file
  G_out_F   = strcmp(G_CLA_Output_f,"-") == 0 ? stdout : fopen(G_CLA_Output_f,"w");
  if (G_out_F == NULL) {
    fprintf(stderr,"Could not open output file: '%s'\n", G_CLA_Output_f);
    exit(-1);
  }
	else if (!G_CLA_Silent_b) {
		fprintf(stderr,"Opened output file '%s'\n", G_CLA_Output_f);
	}

  // load all of the PWMs of interest.
	loadAllMatrices();

  // ......................................................................

  // sequence length
  G_sequenceLength_n = 0;

  // sequence ordinal
  G_sequence_o = 0;

	outputHead();

  // process the sequences
  while (1) {
    loadNextSequence(seq_F);
    if (G_sequenceLength_n < 0 && G_CLA_Format == 'j') 
	{
		// remove ,\n from last sequence Block when outputting json
		int charsToDelete = 2;
    	fseeko(G_out_F,-charsToDelete,SEEK_END);
    	off_t position = ftello(G_out_F);
    	ftruncate(fileno(G_out_F), position);
		fprintf(G_out_F, "\n");
		break;
	}
	else if (G_sequenceLength_n < 0)
	{
		break;
	}
    G_sequence_o++;
    if (G_CLA_UseSequenceBackground) updateLogOddsDist();
		outputSeqHead();
    scoreWithMatrices();
		outputSeqTail();
  }

	outputTail();

  // fprintf(stderr, "Done\n");

  exit(0);
}

// ----------------------------------------------------------------------

char encode ( char Base ) {
  switch (tolower(Base)) {
  case 'a': return DNA_A;
  case 'b': return DNA_ALL - DNA_A;
  case 'c': return DNA_C;
  case 'd': return DNA_ALL - DNA_C;
  case 'g': return DNA_G;
  case 'h': return DNA_ALL - DNA_G;
  case 'k': return DNA_G | DNA_T;
  case 'm': return DNA_A | DNA_C;
  case 'n': return DNA_ALL;
  case 'r': return DNA_A | DNA_G;
  case 's': return DNA_C | DNA_G;
  case 't': return DNA_T;
  case 'v': return DNA_ALL - DNA_T;
  case 'w': return DNA_A | DNA_T;
  case 'x': return DNA_NONE;
  case 'y': return DNA_C | DNA_T;
  default : return DNA_ERROR;
  }
}

char decode ( char Code ) {
  switch (Code) {
  case DNA_A          : return 'a';
  case DNA_ALL - DNA_A: return 'b';
  case DNA_C          : return 'c';
  case DNA_ALL - DNA_C: return 'd';
  case DNA_G          : return 'g';
  case DNA_ALL - DNA_G: return 'h';
  case DNA_G | DNA_T  : return 'k';
  case DNA_A | DNA_C  : return 'm';
  case DNA_ALL        : return 'n';
  case DNA_A | DNA_G  : return 'r';
  case DNA_C | DNA_G  : return 's';
  case DNA_T          : return 't';
  case DNA_ALL - DNA_T: return 'v';
  case DNA_A | DNA_T  : return 'w';
  case DNA_C | DNA_T  : return 'y';
  case DNA_NONE       : return 'x';
  default             : return 'X';
  }
}

// ========================================================================
// --------------------- Reading FASTA Sequence Files ---------------------
// ========================================================================

// ----------------------------------------------------------------------

int loadNextSequence ( FILE * File ) {
  G_sequenceLength_n = -1;
  if (loadTitle(File)) loadSequence(File);
}

int loadTitle ( FILE * File ) {

  // index into title
  int i = 0;

  // read in to this character
  int c;
	
  // make sure we get a '>'.
  c = fgetc(File);
  if (c != '>') {
    return 0;
  }

  // read the rest of the title.
  while (c = fgetc(File)) {
    if (c == '\n') {
      G_sequenceTitle[i] = '\000';

			// set the id.
			for (i = 0; i < strlen(G_sequenceTitle); i++) {
				if  (G_sequenceTitle[i] != ' ') {
					G_sequenceId[i] = G_sequenceTitle[i];
				}
				else {
					break;
				}
			}
			G_sequenceId[i] = '\000';
      return 1;

    } else if (i < MAX_TTL_LENGTH) {
      G_sequenceTitle[i++] = c;
    }
  }

  return 0;
}

void loadSequence ( FILE * File ) {

  // index into sequence
  int i = 0;

  // read into this char
  int c;

  // binary code
  int b;

  while (! feof(File)) {

    c = fgetc(File);

    // it seems we must check after reading as it does not get set
    // upon reading the last char of the file only after we have
    // attempted to read past the last char.
    if (feof(File)) break;

    //fprintf(stderr, "%c", c);

    if (c != '>') {
      b = encode(c);
      if (b != DNA_ERROR) {
        if (i < MAX_SEQ_LENGTH) {
          G_sequence_c[i] = c;
          G_sequence_b[i] = b;
          i++;
        }

        // warn that sequence is too long.
        else {
          printf("Sequence overrun\n");
        }
      }

      // white space is legal; but we'll ignore it.
      else if (isspace(c)) {
      }

      // bad base
      else {
        fprintf(stderr, "Skipping bad base '%c'.\n", c);
      }
    }

    // hit s '>', put it back and bail out
    else {
      ungetc(c,File);
      break;
    }
  }

  // set the global length.
  G_sequenceLength_n = i;
}

// -------------------------- updateLogOddsDist ---------------------------

void updateLogOddsDist ( void ) {

  int a_n = 1;
  int c_n = 1;
  int g_n = 1;
  int t_n = 1;
  int o_n = 1;

  int i;

  int n;

  for (i = 0; i < G_sequenceLength_n; i++) {
    switch (toupper(G_sequence_c[i])) {
    case 'A': a_n++; break;
    case 'C': c_n++; break;
    case 'G': g_n++; break;
    case 'T': t_n++; break;
    default:  o_n++; break;
    }
  }

  n = a_n + c_n + g_n + t_n;

  G_CLA_LogOddsBkg[0] = (double)a_n / (double)n;
  G_CLA_LogOddsBkg[1] = (double)c_n / (double)n;
  G_CLA_LogOddsBkg[2] = (double)g_n / (double)n;
  G_CLA_LogOddsBkg[3] = (double)t_n / (double)n;

  for (i = 0; i < 4; i++) {
    G_logBkg[i] = log(G_CLA_LogOddsBkg[i]) / log(2.0);
  }

}


// ========================================================================
// ---------------------------- Matrix Scoring ----------------------------
// ========================================================================

/* ......................................................................

Matrices are read from a file.  The defline and counts are read into a
temporary count buffer which is reused for each matrix.  This data is
then copied into a new matrix structure in an array of matrices.  The
counts may be copied into the scoring section as follows.  If -nno is
set the counts are used immediately as is.  If -usb is set, then the
counts are converted to log-odds after each sequence is read.
Otherwise the counts are converted to log-odds scores using the
specified or uniform background.

*/

// ----------------------------------------------------------------------

void scoreWithMatrices ( void ) {

	int i;
	// for json this is the counter for matrix+number
	int counterForMatrixName = 1;

	for (i = 0; i < G_matrixCache_n; i++) {
    //fprintf(stderr,"swm\t%d\n", i);
		G_theMatrix = & G_matrices[i];
    if (G_CLA_UseSequenceBackground) makeProbs();
		expandMatrix();
		outputMatHead(counterForMatrixName++);
		scoreOneMatrix();
		outputMatTail();
	}

}

// --------------------------- Loadallmatrices ----------------------------

void loadAllMatrices ( void ) {

	FILE * matrix_F;

	int i;

	int b;

  // ............................................................

	matrix_F = fopen(G_CLA_Matrix_f,"r");
  if (matrix_F == NULL) {
    fprintf(stderr, "Could not open matrix file '%s'.\n", G_CLA_Matrix_f);
    exit(-1);
  }
	else if (!G_CLA_Silent_b) {
		fprintf(stderr, "Opened matrix file '%s'\n", G_CLA_Matrix_f);
	}

  // matrix ordinal
  G_matrix_o = 0;

  while (1) {


		loadNextMatrix(matrix_F);

    // skip the empty ones
    if (G_matrix_n < 0) {
			break;
		}

    // make sure the defline matches any selections going on.
		else if ((G_CLA_SelectId == NULL || strcmp(G_matrixId,    G_CLA_SelectId) == 0) &&
						 (G_CLA_Select   == NULL || strstr(G_matrixTitle, G_CLA_Select)   != NULL)
						 ) {

			G_matrix_o++;

			//fprintf(stderr,"%d\n", G_matrix_o);

			// get an easy pointer to the current matrix
			G_theMatrix = &G_matrices[G_matrix_o - 1];

			// copy the ordinal
			G_theMatrix->matrix_o = G_matrix_o;

			// copy the length
			G_theMatrix->matrix_n = G_matrix_n;

			// copy the matrix title
			G_theMatrix->matrixTitle = (char *)malloc(strlen(G_matrixTitle)+1);
			strcpy(G_theMatrix->matrixTitle,G_matrixTitle);

			// copy the matrix id
			G_theMatrix->matrixId    = (char *)malloc(strlen(G_matrixId)+1);
			strcpy(G_theMatrix->matrixId,G_matrixId);

			// copy the matrix itself.
			G_theMatrix->matrix = (double **)malloc(sizeof(double *) * G_matrix_n);
      G_theMatrix->counts = (double **)malloc(sizeof(double *) * G_matrix_n);

			for (i = 0; i < G_matrix_n; i++) {

				G_theMatrix->matrix[i] = (double *)malloc(sizeof(double) * DNA_SIZE);
				G_theMatrix->counts[i] = (double *)malloc(sizeof(double) * DNA_SIZE);

        G_theMatrix->counts[i][DNA_A] = G_matrix[i][DNA_A];
        G_theMatrix->counts[i][DNA_C] = G_matrix[i][DNA_C];
        G_theMatrix->counts[i][DNA_G] = G_matrix[i][DNA_G];
        G_theMatrix->counts[i][DNA_T] = G_matrix[i][DNA_T];
			}

      if (!G_CLA_UseSequenceBackground) makeProbs();

			//fprintf(stderr,"%d\t%s\n", G_theMatrix->matrix_o, G_theMatrix->matrixId);
		}
	}

  fclose(matrix_F);

	G_matrixCache_n = G_matrix_o;

	if (!G_CLA_Silent_b) fprintf(stderr,"Loaded %d matrices.\n", G_matrix_o);
}

// ----------------------------------------------------------------------

void loadNextMatrix ( FILE * Matrix_F ) {

  // column counter
  int column_i = 0;

  // line
  char buffer[MAX_TTL_LENGTH+1];

	// working index
	int i;

	// default condition is to fail.
	G_matrix_n = -1;

  // read until > line
  while (1) {
    if (fgets(buffer, MAX_TTL_LENGTH, Matrix_F)) {

      // compare to 
      if (strncmp(buffer, ">", 1) == 0) {
        buffer[strlen(buffer)-1] = '\000';
        strcpy(G_matrixTitle, &buffer[1]);

				// set the id.
        for (i = 0; i < strlen(G_matrixTitle); i++) {
          if (!isblank(G_matrixTitle[i])) {
            G_matrixId[i] = G_matrixTitle[i];
          }
          else {
            break;
          }
        }
        G_matrixId[i] = '\000';

        break;
      }
    } else {
      return;
    }
  }

	//fprintf(stderr,"%s\n",G_matrixId);

  while (1) {
    if (fscanf(Matrix_F,
               "%lf %lf %lf %lf",
               &G_matrix[column_i][DNA_A],
               &G_matrix[column_i][DNA_C],
               &G_matrix[column_i][DNA_G],
               &G_matrix[column_i][DNA_T]
               ) == 4) {
      column_i++;
    } else {
      break;
    }
  }

	G_matrix_n = column_i;
  
  //fprintf(stderr, "Matrix: %d (%d)\n", G_matrix_o, G_matrix_n);
  fgets(buffer, 127, Matrix_F);
	//fprintf(stderr,"EOM-%s", buffer);


}

// ----------------------------------------------------------------------

/*

To compute matches to the reverse complement of the current PWM we
compute the r-c of the PWM in this function.  To do this we walk
through the first half of the columns of the PWM exchanging data
between the complementary unambiguous bases and columns.  The middle
column, when there is one, is a special case.  There we only need to
flip the first two bases (since we've arranged them in this order).
Once the exchange is done, we expand the matrix to get the full
augmented PWM.

*/

void flipMatrix   ( void ) {

  int column_i;

  int antiColumn_i;

	// order is important, so that we only need to flip first two bases
	// in mid-matrix column.
  int bases[4]     = { DNA_A, DNA_C, DNA_G, DNA_T };
  int antibases[4] = { DNA_T, DNA_G, DNA_C, DNA_A };

	// working base
  int base_i;

	// maximum base to consider
	int lastBase_i;

	// temporary working area
  double tmp;


  for (column_i = 0; column_i < (G_theMatrix->matrix_n+1)/2; column_i++) {

    antiColumn_i = G_theMatrix->matrix_n - column_i - 1;

		lastBase_i   = column_i == antiColumn_i ? 2 : 4;

		for (base_i = 0; base_i < lastBase_i; base_i++) {
      tmp                                                  = G_theMatrix->matrix[column_i    ][bases    [base_i]];
      G_theMatrix->matrix[column_i    ][bases    [base_i]] = G_theMatrix->matrix[antiColumn_i][antibases[base_i]];
      G_theMatrix->matrix[antiColumn_i][antibases[base_i]] = tmp;
    }

  }

  expandMatrix();
}

// ----------------------------------------------------------------------

void makeProbs    ( void ) {

  // column checker
  int column_i;

  // total
  double total_d;

  // per base pseudocounts
  double pbpc_d = G_CLA_PseudoCounts / 4;
  double pc_d   = G_CLA_PseudoCounts;

  // numerator
  double normNum_d;

  // log base converter
  double toBase2_d = log((double)(2.0));

  // just copy as is
  if (G_CLA_NoNormalize) {

    for (column_i = 0; column_i < G_theMatrix->matrix_n; column_i++) {

      G_theMatrix->matrix[column_i][DNA_A] = G_theMatrix->counts[column_i][DNA_A];
      G_theMatrix->matrix[column_i][DNA_C] = G_theMatrix->counts[column_i][DNA_C];
      G_theMatrix->matrix[column_i][DNA_G] = G_theMatrix->counts[column_i][DNA_G];
      G_theMatrix->matrix[column_i][DNA_T] = G_theMatrix->counts[column_i][DNA_T];

    }
  }

  // convert to probabilities, then log-odds
  else {

    for (column_i = 0; column_i < G_theMatrix->matrix_n; column_i++) {
      total_d
        = G_theMatrix->counts[column_i][DNA_A]
        + G_theMatrix->counts[column_i][DNA_C] 
        + G_theMatrix->counts[column_i][DNA_G] 
        + G_theMatrix->counts[column_i][DNA_T];

      // determine zero
      if (G_CLA_PseudoCounts < 0) {
        pc_d   = sqrt(total_d);
        if (pc_d <= 0) pc_d = 1;
        pbpc_d = pc_d / 4;
      }
    
      normNum_d = total_d + pc_d;

      G_theMatrix->matrix[column_i][DNA_A] = log((double)((G_theMatrix->counts[column_i][DNA_A] + pbpc_d) / normNum_d)) / toBase2_d - G_logBkg[0];
      G_theMatrix->matrix[column_i][DNA_C] = log((double)((G_theMatrix->counts[column_i][DNA_C] + pbpc_d) / normNum_d)) / toBase2_d - G_logBkg[1];
      G_theMatrix->matrix[column_i][DNA_G] = log((double)((G_theMatrix->counts[column_i][DNA_G] + pbpc_d) / normNum_d)) / toBase2_d - G_logBkg[2];
      G_theMatrix->matrix[column_i][DNA_T] = log((double)((G_theMatrix->counts[column_i][DNA_T] + pbpc_d) / normNum_d)) / toBase2_d - G_logBkg[3];
    }
  }

#ifdef __DEBUG__
  for (column_i = 0; column_i < G_theMatrix->matrix_n; column_i++) {

    fprintf(stderr,
            "  <!-- %d / %lf %lf %lf %lf / %lf %lf %lf %lf-->\n",
            column_i,
          
            G_theMatrix->counts[column_i][DNA_A],
            G_theMatrix->counts[column_i][DNA_C],
            G_theMatrix->counts[column_i][DNA_G],
            G_theMatrix->counts[column_i][DNA_T],
          
            G_theMatrix->matrix[column_i][DNA_A],
            G_theMatrix->matrix[column_i][DNA_C],
            G_theMatrix->matrix[column_i][DNA_G],
            G_theMatrix->matrix[column_i][DNA_T]
            
            );
  }
#endif

}

// ----------------------------------------------------------------------

void expandMatrix ( void ) {

  // column checker
  int column_i;

	// working
	int i;

  for (column_i = 0; column_i < G_theMatrix->matrix_n; column_i++) {

    //fprintf(stderr,"Expanding column %d\n", column_i);

		switch (G_CLA_IupacTreatment) {

		case 'n':
			clearTwo(column_i, DNA_A, DNA_C);
			clearTwo(column_i, DNA_A, DNA_G);
			clearTwo(column_i, DNA_A, DNA_T);
			clearTwo(column_i, DNA_C, DNA_G);
			clearTwo(column_i, DNA_C, DNA_T);
			clearTwo(column_i, DNA_G, DNA_T);

			clearThree(column_i, DNA_A, DNA_C, DNA_G);
			clearThree(column_i, DNA_A, DNA_C, DNA_T);
			clearThree(column_i, DNA_A, DNA_G, DNA_T);
			clearThree(column_i, DNA_C, DNA_G, DNA_T);

			clearFour(column_i, DNA_A, DNA_C, DNA_G, DNA_T);

		case 'a':
			averageTwo(column_i, DNA_A, DNA_C);
			averageTwo(column_i, DNA_A, DNA_G);
			averageTwo(column_i, DNA_A, DNA_T);
			averageTwo(column_i, DNA_C, DNA_G);
			averageTwo(column_i, DNA_C, DNA_T);
			averageTwo(column_i, DNA_G, DNA_T);

			averageThree(column_i, DNA_A, DNA_C, DNA_G);
			averageThree(column_i, DNA_A, DNA_C, DNA_T);
			averageThree(column_i, DNA_A, DNA_G, DNA_T);
			averageThree(column_i, DNA_C, DNA_G, DNA_T);

			averageFour(column_i, DNA_A, DNA_C, DNA_G, DNA_T);

      break;

		case 'p':
			minimumTwo(column_i, DNA_A, DNA_C);
			minimumTwo(column_i, DNA_A, DNA_G);
			minimumTwo(column_i, DNA_A, DNA_T);
			minimumTwo(column_i, DNA_C, DNA_G);
			minimumTwo(column_i, DNA_C, DNA_T);
			minimumTwo(column_i, DNA_G, DNA_T);

			minimumThree(column_i, DNA_A, DNA_C, DNA_G);
			minimumThree(column_i, DNA_A, DNA_C, DNA_T);
			minimumThree(column_i, DNA_A, DNA_G, DNA_T);
			minimumThree(column_i, DNA_C, DNA_G, DNA_T);
			
			minimumFour(column_i, DNA_A, DNA_C, DNA_G, DNA_T);

			break;

		case 'o':
			maximumTwo(column_i, DNA_A, DNA_C);
			maximumTwo(column_i, DNA_A, DNA_G);
			maximumTwo(column_i, DNA_A, DNA_T);
			maximumTwo(column_i, DNA_C, DNA_G);
			maximumTwo(column_i, DNA_C, DNA_T);
			maximumTwo(column_i, DNA_G, DNA_T);

			maximumThree(column_i, DNA_A, DNA_C, DNA_G);
			maximumThree(column_i, DNA_A, DNA_C, DNA_T);
			maximumThree(column_i, DNA_A, DNA_G, DNA_T);
			maximumThree(column_i, DNA_C, DNA_G, DNA_T);

			maximumFour(column_i, DNA_A, DNA_C, DNA_G, DNA_T);
			
			break;
		}

		G_theMatrix->matrix[column_i][DNA_NONE] = VERY_BAD_SCORE;

#ifdef __DEBUG__
		fprintf(stderr, "%d", column_i);
		for (i = DNA_NONE; i <= DNA_ALL; i++) {
			fprintf(stderr, "\t%0.1f", G_theMatrix->matrix[column_i][i]);
		}
		fprintf(stderr, "\n");
#endif
	}

}

// ----------------------------------------------------------------------

void clearTwo ( int Column_i, int Base1_i, int Base2_i ) {

  G_theMatrix->matrix[Column_i][Base1_i|Base2_i] = VERY_BAD_SCORE;
}

// ----------------------------------------------------------------------

void clearThree ( int Column_i, int Base1_i, int Base2_i, int Base3_i ) {

  G_theMatrix->matrix[Column_i][Base1_i|Base2_i|Base3_i] = VERY_BAD_SCORE;
}


// ----------------------------------------------------------------------

void clearFour ( int Column_i, int Base1_i, int Base2_i, int Base3_i, int Base4_i ) {

  G_theMatrix->matrix[Column_i][Base1_i|Base2_i|Base3_i|Base4_i] = VERY_BAD_SCORE;
}

// ----------------------------------------------------------------------

void averageTwo ( int Column_i, int Base1_i, int Base2_i ) {

  G_theMatrix->matrix[Column_i][Base1_i|Base2_i] =
    ( G_theMatrix->matrix[Column_i][Base1_i] +
      G_theMatrix->matrix[Column_i][Base2_i] ) / 2;
}

// ----------------------------------------------------------------------

void averageThree ( int Column_i, int Base1_i, int Base2_i, int Base3_i ) {

  G_theMatrix->matrix[Column_i][Base1_i|Base2_i|Base3_i] =
    ( G_theMatrix->matrix[Column_i][Base1_i] +
      G_theMatrix->matrix[Column_i][Base2_i] +
      G_theMatrix->matrix[Column_i][Base3_i] ) / 3;
}


// ----------------------------------------------------------------------

void averageFour ( int Column_i, int Base1_i, int Base2_i, int Base3_i, int Base4_i ) {

  G_theMatrix->matrix[Column_i][Base1_i|Base2_i|Base3_i|Base4_i] =
    ( G_theMatrix->matrix[Column_i][Base1_i] +
      G_theMatrix->matrix[Column_i][Base2_i] +
      G_theMatrix->matrix[Column_i][Base3_i] +
      G_theMatrix->matrix[Column_i][Base3_i]
			) / 4;
}

// ----------------------------------------------------------------------

void minimumTwo ( int Column_i, int Base1_i, int Base2_i ) {

  G_theMatrix->matrix[Column_i][Base1_i|Base2_i] =
    MIN( G_theMatrix->matrix[Column_i][Base1_i],
         G_theMatrix->matrix[Column_i][Base2_i]
         );
}

// ----------------------------------------------------------------------

void minimumThree ( int Column_i, int Base1_i, int Base2_i, int Base3_i ) {

  G_theMatrix->matrix[Column_i][Base1_i|Base2_i|Base3_i] =
    MIN(MIN(G_theMatrix->matrix[Column_i][Base1_i],
            G_theMatrix->matrix[Column_i][Base2_i]),
        G_theMatrix->matrix[Column_i][Base3_i]
        );
}

// ----------------------------------------------------------------------

void minimumFour ( int Column_i, int Base1_i, int Base2_i, int Base3_i, int Base4_i ) {

  G_theMatrix->matrix[Column_i][Base1_i|Base2_i|Base3_i|Base4_i] =
    MIN(MIN(MIN(G_theMatrix->matrix[Column_i][Base1_i],
								G_theMatrix->matrix[Column_i][Base2_i]
								),
						G_theMatrix->matrix[Column_i][Base3_i]
						),
				G_theMatrix->matrix[Column_i][Base4_i]
        );
}

// ----------------------------------------------------------------------

void maximumTwo ( int Column_i, int Base1_i, int Base2_i ) {

  G_theMatrix->matrix[Column_i][Base1_i|Base2_i] =
    MAX( G_theMatrix->matrix[Column_i][Base1_i],
         G_theMatrix->matrix[Column_i][Base2_i]
         );
}

// ----------------------------------------------------------------------

void maximumThree ( int Column_i, int Base1_i, int Base2_i, int Base3_i ) {

  G_theMatrix->matrix[Column_i][Base1_i|Base2_i|Base3_i] =
    MAX(MAX(G_theMatrix->matrix[Column_i][Base1_i],
            G_theMatrix->matrix[Column_i][Base2_i]),
        G_theMatrix->matrix[Column_i][Base3_i]
        );
}

// ----------------------------------------------------------------------

void maximumFour ( int Column_i, int Base1_i, int Base2_i, int Base3_i, int Base4_i ) {

  G_theMatrix->matrix[Column_i][Base1_i|Base2_i|Base3_i|Base4_i] =
    MAX(MAX(MAX(G_theMatrix->matrix[Column_i][Base1_i],
								G_theMatrix->matrix[Column_i][Base2_i]
								),
						G_theMatrix->matrix[Column_i][Base3_i]
						),
				G_theMatrix->matrix[Column_i][Base4_i]
        );
}

// ----------------------------------------------------------------------

double maxScore ( void ) {

  double totalScore = 0;

  double maxScoreInColumn;

  int column_i;

  int base_i;

  int bases[4] = { DNA_A, DNA_C, DNA_G, DNA_T };

  for (column_i = 0; column_i < G_theMatrix->matrix_n; column_i++) {
    maxScoreInColumn = G_theMatrix->matrix[column_i][bases[0]];
    for (base_i = 1; base_i < 4; base_i++) {
			maxScoreInColumn = MAX(maxScoreInColumn,
														 G_theMatrix->matrix[column_i][bases[base_i]]
														 );
    }
    totalScore += maxScoreInColumn;
  }

  return totalScore;
}

// ----------------------------------------------------------------------

void scoreOneMatrix ( void ) {
	
  // position in sequence
  int i;

  // position matrix
  int j;

  // temporary for printing
  int k;

  // local start and end
  int localStart = G_CLA_Start <  0                  ? 0                      : G_CLA_Start;
  int localEnd   = G_CLA_End   >= G_sequenceLength_n ? G_sequenceLength_n - 1 : G_CLA_End;

  // how far to go.
  int last_i = localEnd - G_theMatrix->matrix_n + 1;

  // best score (Lm)
  double lm = maxScore();

  // score accumulation
  double score;

  // coded base
  int b;

  // oligo
  char oligo[ MAX_MTX_LENGTH + 1 ];

  // senses
  char * senses = "NR";

  // sense counter
  int sense_i;

	// hit counter
	int hit_n = 0;

  // coded base pointer
  char * base_b;

  // hit coords
  int   hitBeg_n;
  int   hitEnd_n;

  // For json this is the counter for Hit+number
  int counterForHitName = 1;

#ifdef __DEBUG__
  printf("    <!-- LocalStart='%d' LocalEnd='%d' MatrixLength='%d' Last='%d' -->\n", localStart, localEnd, G_theMatrix->matrix_n, last_i);
#endif

  // consider each sense (we'll bail out in case of nrc
  for (sense_i = 0; sense_i < 2; sense_i++) {

    // consider all starting positions
    for (i = localStart; i <= last_i; i++) {

      // consider all positions in matrix
      for (j = 0, base_b = &G_sequence_b[i], score = 0.0;
           j < G_theMatrix->matrix_n;
           j++,   base_b++
           ) {
        b = *base_b;
        score += G_theMatrix->matrix[j][b];
#ifdef __DEBUG__
        fprintf(stderr,"%d %d %d %d %c %f\n", sense_i, i, j, b, G_sequence_c[i+j], score);
#endif
        if (b == DNA_NONE) { i += j; break; }
      }

      // print the score.
      if (score >= G_CLA_MinScore && lm - score <= G_CLA_MaxDeficit) {

        // count the hits
				hit_n++;

        // hits coords
        hitBeg_n = 1 + i + G_CLA_PositionOffset;
        hitEnd_n = 1 + i + G_CLA_PositionOffset + G_theMatrix->matrix_n - 1;

        // sequence of hit
				if (G_computeOligo_b) {
					for (k = 0; k < G_theMatrix->matrix_n; k++) {
						oligo[k] = G_sequence_c[i+k];
					}
					oligo[k] = '\000';
        }

				// output the hit.
				switch (G_CLA_Format) {
				case 'b':
					fprintf(G_out_F,
									"%s,%s,%lf,%d,%d,%d\n",
									G_theMatrix->matrixId,
									G_sequenceId,
									score,
                  hitBeg_n, hitEnd_n,
									sense_i
									);
					break;

        case 'g':
          fprintf(G_out_F,
                  "%s\t%s\t%s\t%d\t%d\t%lf\t%c\t%c\tModel %s ; \n",
                  G_sequenceId,
                  "tessWms",
                  G_CLA_GffFeature != NULL ? G_CLA_GffFeature : "misc_binding",
                  hitBeg_n, hitEnd_n,
                  score,
                  sense_i ? '-' : '+',
                  '.',
                  G_theMatrix->matrixId
                  );
          break;

				case 't':
					fprintf(G_out_F,
									"%s\t%s\t%lf\t%d\t%d\t%d\t%lf\t%s\n",
									G_sequenceId,
									G_theMatrix->matrixId,
									score,
                  hitBeg_n, hitEnd_n,
									sense_i,
                  lm-score,
                  oligo
									);
					break;

				case 'r':
					fprintf(G_out_F,
									"        <XchNode score='%f' index='%d' leaf='1' parts='0' sense='%s' start='%d' end='%d' contents='0' type='BindingSite' model='%s' />\n",
									score,
									hit_n,
									sense_i == 0 ? "+1" : "-1",
                  hitBeg_n, hitEnd_n,
									G_theMatrix->matrixId
									);
					break;

				case 'w':
					fprintf(G_out_F,
									"      <Hit Sense='%c' Start='%d' Stop='%d' La='%lf' Lm='%lf' Ld='%lf' Oligo='%s' />\n",
									senses[sense_i],
                  hitBeg_n, hitEnd_n,
									score, lm, lm-score, oligo
									);
					break;

				case 'j':
					fprintf(G_out_F,
									"      \"Hit%d\" : {\"Sense\" : \"%c\", \"Start\" : \"%d\", \"Stop\" : \"%d\", \"La\" : \"%lf\", \"Lm\" : \"%lf\", \"Ld\" : \"%lf\", \"Oligo\" : \"%s\"},\n",
									counterForHitName++, senses[sense_i], hitBeg_n, hitEnd_n, score, lm, lm-score, oligo
									);
					break;
				} // eo switch
			} // eo good enough score

			if (G_CLA_Format == 'j' && i == last_i && sense_i == 1) 
			{
				// remove ,\n from last Hit Block when outputting json
				int charsToDelete = 2;
	    		fseeko(G_out_F,-charsToDelete,SEEK_END);
	    		off_t position = ftello(G_out_F);
	    		ftruncate(fileno(G_out_F), position);
				fprintf(G_out_F, "\n");
			}

		} // eo start scan	

		if (G_CLA_NoRc) break;
			
		// flip sense
		flipMatrix();
		
	} // eo senses

}

// ========================================================================
// -------------------------------- Output --------------------------------
// ========================================================================

// ------------------------------ outputHead ------------------------------

/* Stuff output at the beginning of the run. */

void outputHead ( void ) {
	switch (G_CLA_Format) {
	case 'b':
		break;

  case 'g':
    fprintf(G_out_F, "##source-version tessWms %s\n", G_version);
    // fprintf(G_out_F, "##date \n");
    break;

	case 't':
		break;

	case 'r':
		fprintf(G_out_F, "<opt>\n");
		fprintf(G_out_F, "  <ExtraInfo\n");
    fprintf(G_out_F, "    mlo = '%f'\n", G_CLA_MinScore);
    fprintf(G_out_F, "    beg = '%d'\n", G_CLA_Start+1);
    fprintf(G_out_F, "    end = '%d'\n", G_CLA_End+1);
    fprintf(G_out_F, "    nrc = '%d'\n", G_CLA_NoRc);
    fprintf(G_out_F, "    nno = '%d'\n", G_CLA_NoNormalize);
    fprintf(G_out_F, "    tpc = '%f'\n", G_CLA_PseudoCounts);
    fprintf(G_out_F, "    po  = '%d'\n", G_CLA_PositionOffset);
    fprintf(G_out_F, "    abt = '%c'\n", G_CLA_IupacTreatment);
    fprintf(G_out_F, "    lob = '%f,%f,%f,%f'\n",
            G_CLA_LogOddsBkg[0],G_CLA_LogOddsBkg[1],G_CLA_LogOddsBkg[2],G_CLA_LogOddsBkg[3]
            );
    fprintf(G_out_F, "  />\n");
    fprintf(G_out_F, "  <Grammar></Grammar>\n");
		fprintf(G_out_F, "  <Result>\n");
		break;

	case 'w':
		fprintf(G_out_F, "<WMS>\n");
		break;

	case 'j':
		fprintf(G_out_F, "{\"WMS\" : \n  {\n");
		break;

	}
}

// ---------------------------- outputSeqHead -----------------------------

/* Stuff output before each sequence */

void outputSeqHead ( void ) {

	switch (G_CLA_Format) {
	case 'b':
		if (G_sequence_o > 1 && G_CLA_SegmentedOutput) {
			fprintf(G_out_F, "#\n");
		}
		break;

  case 'g':
    if (G_CLA_SegmentedOutput) {
      fprintf(G_out_F, "##DNA %s\n", G_sequenceTitle);
    }
    break;

	case 't':
		if (G_sequence_o > 1 && G_CLA_SegmentedOutput) {
			fprintf(G_out_F,
              "#%s\t%s\t%s\t%s\t%s\t%s\n",
              "SeqId",
              "PwmId",
              "Lod",
              "Beg",
              "End",
              "Sense"
              );
		}
		break;

	case 'r':
		if (G_sequence_o > 1 && G_CLA_SegmentedOutput) {
			outputTail();
			fprintf(G_out_F, "//\n");
			outputHead();
		}

		fprintf(G_out_F,
						"    <Sequence SeqId='%s' Start='1' Stop='%d'>\n",
						G_sequenceId, G_sequenceLength_n
						);
		fprintf(G_out_F,
						"      <FeasibleIntervals Id='%s' Start='%d' End='%d' Sense='+'>\n",
						G_sequenceId, MAX(1,G_CLA_Start+1), MIN(G_sequenceLength_n,G_CLA_End+1)
						);
		break;

	case 'w':
		if (G_sequence_o > 1 && G_CLA_SegmentedOutput) {
			outputTail();
			fprintf(G_out_F, "//\n");
			outputHead();
		}
			
		fprintf(G_out_F,
						"  <Sequence Ordinal='%d' SAC='na' Length='%d' >\n    <Title>%s</Title>\n",
						G_sequence_o, G_sequenceLength_n, G_sequenceTitle
						);
		break;

	case 'j':
		if (G_sequence_o > 1 && G_CLA_SegmentedOutput) {
			outputTail();
			fprintf(G_out_F, "//\n");
			outputHead();
		}
			
		fprintf(G_out_F,
						"  \"Sequence%d\" : {\"Ordinal\" : \"%d\", \"SAC\" : \"na\", \"Length\" : \"%d\",\n    \"Title\" : \"%s\",\n",
						G_sequence_o, G_sequence_o, G_sequenceLength_n, G_sequenceTitle
						);
		break;
	}
}

// ---------------------------- outputMatHead -----------------------------

/* Stuff output before each matrix. */

void outputMatHead ( int counterForMatrixName ) {

	switch (G_CLA_Format) {
	case 'b':
		break;

	case 't':
		break;

  case 'g':
    break;

	case 'r':
		break;

	case 'w':
    fprintf(G_out_F,
						"    <Matrix MAC='%d' FID='%s' Length='%d'>\n",
						G_theMatrix->matrix_o, G_theMatrix->matrixId, G_theMatrix->matrix_n
						);
    fprintf(G_out_F,
						"      <Title>%s</Title>\n",
						G_theMatrix->matrixTitle
						);
		break;

	case 'j':
    fprintf(G_out_F,
						"    \"Matrix%d\" : {\"MAC\" : \"%d\", \"FID\" : \"%s\", \"Length\" : \"%d\",\n",
						counterForMatrixName, G_theMatrix->matrix_o, G_theMatrix->matrixId, G_theMatrix->matrix_n
						);
    fprintf(G_out_F,
						"      \"Title\" : \"%s\",\n",
						G_theMatrix->matrixTitle
						);
		break;
	}
}

// ---------------------------- outputMatTail -----------------------------

/* Stuff to output after each matrix. */

void outputMatTail ( void ) {

  // For json output: last Matrix Block without comma
  int i = G_matrix_o - G_theMatrix->matrix_o;
  char *comma;
  if (i == 0)
	comma = "";
  else
	comma = ",";

  switch (G_CLA_Format) {
  case 'b':
    break;

  case 't':
    break;

  case 'g':
    break;

	case 'r':
		break;

  case 'w':
		fprintf(G_out_F, "    </Matrix>\n");
		break;

  case 'j':
		fprintf(G_out_F, "    }%s\n", comma );
		break;
	}

}

// ---------------------------- outputSeqTail -----------------------------

/* Stuff to output after each sequence. */

void outputSeqTail ( void ) {

	switch (G_CLA_Format) {
	case 'b':
		break;

	case 't':
		break;

  case 'g':
    break;

	case 'r':
		fprintf(G_out_F, "      </FeasibleIntervals>\n");
		//fprintf(G_out_F, "      <ExtraInfo Title="%s" />\n", G_sequenceTitle);
		fprintf(G_out_F, "      <ExtraInfo>\n");
		fprintf(G_out_F, "        <Title>%s</Title>\n", G_sequenceTitle);
		fprintf(G_out_F, "      </ExtraInfo>\n");
		fprintf(G_out_F, "    </Sequence>\n");
		break;

	case 'w':
		fprintf(G_out_F, "  </Sequence>\n");
		break;

	case 'j':
		fprintf(G_out_F, "   },\n");
		break;
	}
}
	
// ------------------------------ outputTail ------------------------------

/* Stuff to output after all is said and done. */

void outputTail ( void ) {
	switch (G_CLA_Format) {
	case 'b':
		break;

	case 't':
		break;

  case 'g':
    break;

	case 'r':
    fprintf(G_out_F, "    <ExtraInfo></ExtraInfo>\n");
    fprintf(G_out_F, "  </Result>\n");
    fprintf(G_out_F, "</opt>\n");
		break;

	case 'w':
		fprintf(G_out_F, "</WMS>\n");
		break;

	case 'j':
		fprintf(G_out_F, "  }\n}");
		break;

	}
}
