/*******************************************************************************************
 *
 *  Filter Expression Parser & Evaluator
 *
 *  Author:  Gene Myers
 *  Date  :  Oct. 31, 2016
 *
 ********************************************************************************************/

#ifndef _FILTER_EXPR
#define _FILTER_EXPR

#include "sam.h"
#include "bax.h"

typedef void *Filter;

Filter *parse_filter(char *expr);

int evaluate_bam_filter(Filter *v, samRecord *s);
int evaluate_bax_filter(Filter *v, BaxData *b, SubRead *s);

#endif // _FILTER_EXPR
