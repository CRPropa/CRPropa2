#ifndef __cvector_ph__
#define __cvector_ph__

#ifdef __cplusplus
extern "C" {
#endif

void New_dCVector (dCVector* pVector, const int n);
void Delete_dCVector (dCVector* pVector);
void Initialize_dCVector (dCVector* pVector);
void New_iCMatrix (iCMatrix* pMatrix, const int n1, const int n2);
void Delete_iCMatrix (iCMatrix* pMatrix);
void Initialize_iCMatrix (iCMatrix* pMatrix);

#ifdef __cplusplus
}
#endif

#endif

