#include "c_ncbi_blast_aux.h"

#include <ctype.h>

/** Declared in blast_def.h. */
const int kSegWindow = 12;
const double kSegLocut = 2.2;
const double kSegHicut = 2.5;

void
__sfree(void **x)
{
    free(*x);
    *x = NULL;
    return;
}

void * BlastMemDup (const void *orig, size_t size)
{
	void*	copy;

	if (orig == NULL || size == 0)
		return NULL;

	if ((copy = malloc (size)) == NULL)
		return NULL;

	memcpy(copy, orig, size);
    return copy;
}

/*****************************************************************************
*
*   ListNodeNew(vnp)
*      adds after last node in list if vnp not NULL
*
*****************************************************************************/
ListNode* ListNodeNew (ListNode* vnp)
{
	ListNode* newnode;

	newnode = (ListNode*) calloc(1, sizeof(ListNode));
	if (vnp != NULL)
    {
        while (vnp->next != NULL)
            vnp = vnp->next;
		vnp->next = newnode;
    }
	return newnode;
}

/*****************************************************************************
*
*   ListNodeAdd(head)
*      adds after last node in list if *head not NULL
*      If *head is NULL, sets it to the new ListNode
*      returns pointer to the NEW node added
*
*****************************************************************************/
ListNode* ListNodeAdd (ListNode** head)
{
	ListNode* newnode;

	if (head != NULL)
	{
		newnode = ListNodeNew(*head);
		if (*head == NULL)
			*head = newnode;
	}
	else
		newnode = ListNodeNew(NULL);

	return newnode;
}

/*   ListNodeAddPointer (head, choice, value)
*      adds like ListNodeAdd()
*      sets newnode->choice = choice (if choice does not matter, use 0)
*      sets newnode->ptr = value
*
*****************************************************************************/
ListNode* ListNodeAddPointer (ListNode** head, Uint1 choice, 
                                void *value)
{
	ListNode* newnode;

	newnode = ListNodeAdd(head);
	if (newnode != NULL)
	{
		newnode->choice = choice;
		newnode->ptr = value;
	}

	return newnode;
}

/*****************************************************************************
*
*   ListNodeCopyStr (head, choice, str)
*      adds like ListNodeAdd()
*      sets newnode->choice = choice (if choice does not matter, use 0)
*      sets newnode->ptr = str
*         makes a COPY of str
*      if str == NULL, does not add a ListNode
*
*****************************************************************************/
ListNode* ListNodeCopyStr (ListNode** head, Uint1 choice, const char* str)
{
	ListNode* newnode;

	if (str == NULL) return NULL;

	newnode = ListNodeAdd(head);
	if (newnode != NULL)
	{
		newnode->choice = choice;
		newnode->ptr = strdup(str);
	}

	return newnode;
}

/*****************************************************************************
*
*   ListNodeFree(vnp)
*   	frees whole chain of ListNodes
*       Does NOT free associated data pointers
*           see ListNodeFreeData()
*
*****************************************************************************/
ListNode* ListNodeFree (ListNode* vnp)
{
	ListNode* next;

	while (vnp != NULL)
	{
		next = vnp->next;
		sfree(vnp);
		vnp = next;
	}
	return NULL;
}

/*****************************************************************************
*
*   ListNodeFreeData(vnp)
*   	frees whole chain of ListNodes
*       frees associated data pointers - BEWARE of this if these are not
*           allocated single memory block structures.
*
*****************************************************************************/
ListNode* ListNodeFreeData (ListNode* vnp)
{
	ListNode* next;

	while (vnp != NULL)
	{
		sfree(vnp->ptr);
		next = vnp->next;
		sfree(vnp);
		vnp = next;
	}
	return NULL;
}

char* BLAST_StrToUpper(const char* string)
{
    char* retval = NULL;        /* the return value */
    char* p = NULL;             /* auxiliary pointer */

    if ( ! string ) {
        return NULL;
    }

    retval = strdup(string);
    if ( !retval ) {
        return NULL;
    }

    for (p = retval; *p != NULLB; p++) {
        *p = toupper((unsigned char)(*p));
    }
    return retval;
}

/****************************************************************************/

void**
_PSIAllocateMatrix(unsigned int ncols, unsigned int nrows, 
                   unsigned int data_type_sz)
{
    void** retval = NULL;
    unsigned int i = 0;

    retval = (void**) malloc(sizeof(void*) * ncols);
    if ( !retval ) {
        return NULL;
    }

    for (i = 0; i < ncols; i++) {
        retval[i] = (void*) calloc(nrows, data_type_sz);
        if ( !retval[i] ) {
            retval = _PSIDeallocateMatrix(retval, i);
            break;
        }
    }
    return retval;
}

void**
_PSIDeallocateMatrix(void** matrix, unsigned int ncols)
{
    unsigned int i = 0;

    if (!matrix)
        return NULL;

    for (i = 0; i < ncols; i++) {
        sfree(matrix[i]);
    }
    sfree(matrix);
    return NULL;
}

/** Implements the generic copy matrix functions. Prototypes must be defined
 * in the header file manually following the naming convention for 
 * _PSICopyMatrix_int
 */
#define DEFINE_COPY_MATRIX_FUNCTION(type)                           \
void _PSICopyMatrix_##type(type** dest, type** src,                 \
                          unsigned int ncols, unsigned int nrows)   \
{                                                                   \
    unsigned int i = 0;                                             \
    unsigned int j = 0;                                             \
                                                                    \
    ASSERT(dest);                                                   \
    ASSERT(src);                                                    \
                                                                    \
    for (i = 0; i < ncols; i++) {                                   \
        for (j = 0; j < nrows; j++) {                               \
            dest[i][j] = src[i][j];                                 \
        }                                                           \
    }                                                               \
}                                                                   \

DEFINE_COPY_MATRIX_FUNCTION(int)
DEFINE_COPY_MATRIX_FUNCTION(double)