#ifndef __vector_ph__
#define __vector_ph__

#ifdef __cplusplus
extern "C" {
#endif

dVector New_dVector (const int n);
iVector New_iVector (const int n);
dMatrix New_dMatrix (const int n1, const int n2);
iMatrix New_iMatrix (const int n1, const int n2);
d3Tensor New_d3Tensor (const int n1, const int n2, const int n3);
i3Tensor New_i3Tensor (const int n1, const int n2, const int n3);
void Delete_dVector (dVector vector);
void Delete_iVector (iVector vector);
void Delete_dMatrix (dMatrix matrix);
void Delete_iMatrix (iMatrix matrix);
void Delete_d3Tensor (d3Tensor tensor);
void Delete_i3Tensor (i3Tensor tensor);

#ifdef __cplusplus
}
#endif

#endif

