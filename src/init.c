#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
/*extern SEXP _microseq_read_fasta(SEXP);
extern SEXP _microseq_read_fastq(SEXP);
extern SEXP _microseq_read_fastq_Sanger(SEXP);*/
extern SEXP _microseq_revComp(SEXP, SEXP);
extern SEXP _microseq_transl(SEXP, SEXP);
extern SEXP _microseq_translCodon(SEXP);
/*extern SEXP _microseq_write_fasta(SEXP, SEXP, SEXP, SEXP);
extern SEXP _microseq_write_fastq(SEXP, SEXP, SEXP, SEXP);*/
extern SEXP _microseq_extractSeq(SEXP, SEXP, SEXP, SEXP);
extern SEXP _microseq_ORF_index(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
/*    {"_microseq_read_fasta",        (DL_FUNC) &_microseq_read_fasta,        1},
    {"_microseq_read_fastq",        (DL_FUNC) &_microseq_read_fastq,        1},
    {"_microseq_read_fastq_Sanger", (DL_FUNC) &_microseq_read_fastq_Sanger, 1},*/
    {"_microseq_revComp",           (DL_FUNC) &_microseq_revComp,           2},
    {"_microseq_transl",            (DL_FUNC) &_microseq_transl,            2},
    {"_microseq_translCodon",       (DL_FUNC) &_microseq_translCodon,       1},
    /*    {"_microseq_write_fasta",       (DL_FUNC) &_microseq_write_fasta,       4},
    {"_microseq_write_fastq",       (DL_FUNC) &_microseq_write_fastq,       4},*/
    {"_microseq_extractSeq",        (DL_FUNC) &_microseq_extractSeq,        4},
    {"_microseq_ORF_index",         (DL_FUNC) &_microseq_ORF_index,         3},
    {NULL, NULL, 0}
};

void R_init_microseq(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
