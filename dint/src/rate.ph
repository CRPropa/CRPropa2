#ifndef __rate_ph__
#define __rate_ph__

#ifdef __cplusplus
extern "C" {
#endif

void NewRawTotalRate (RawTotalRate* pRate, const int num_main_bins,		     const int num_bg_bins);
void DeleteRawTotalRate (RawTotalRate* pRate);
void InitializeRawTotalRate (RawTotalRate* pRate);
void CopyRawTotalRate (RawTotalRate* pLRate, const RawTotalRate* pRRate);
void ClipRawTotalRate (RawTotalRate* pRate, const int newSize);
void EnlargeRawTotalRate (RawTotalRate* pRate, const int newSize);
void ModifyRawTotalRate (RawTotalRate* pRate, const int newSize);
void ReadRawTotalRate (RawTotalRate* pRate, const char* filename);
void NewRawDiffRate (RawDiffRate* pRate, const int num_main_bins, 		    const int num_bg_bins, const int num_elements);
void DeleteRawDiffRate (RawDiffRate* pRate);
void InitializeRawDiffRate (RawDiffRate* pRate);
void CheckRawDiffRate (RawDiffRate* pRate);
double RawDiffRateElement (const RawDiffRate* pRate, const int i, const int j,			  const int k);
void CopyRawDiffRate (RawDiffRate* pLRate, const RawDiffRate* pRRate);
void CopyRawDiffRateBound (RawDiffRate* pLRate, const RawDiffRate* pRRate);
void ClipRawDiffRate (RawDiffRate* pRate, const int newSize);
void EnlargeRawDiffRate (RawDiffRate* pRate, const int newSize);
void ModifyRawDiffRate (RawDiffRate* pRate, const int newSize);
void ReadRawDiffRate (RawDiffRate* pRate, const char* filename);
void NewTotalRate (TotalRate* pRate, const int num_main_bins);
void DeleteTotalRate (TotalRate* pRate);
void InitializeTotalRate (TotalRate* pRate);
void CopyTotalRate (TotalRate* pLRate, const TotalRate* pRRate);
void ClipTotalRate (TotalRate* pRate, const int newSize);
void EnlargeTotalRate (TotalRate* pRate, const int newSize);
void ModifyTotalRate (TotalRate* pRate, const int newSize);
void ReadTotalRate (TotalRate* pRate, const char* filename);
void NewDiffRate (DiffRate* pRate, const int num_main_bins);
void DeleteDiffRate (DiffRate* pRate);
void InitializeDiffRate (DiffRate* pRate);
void CopyDiffRate (DiffRate* pLRate, const DiffRate* pRRate);
void CopyDiffRateBound (DiffRate* pLRate, const DiffRate* pRRate);
void SetStandardDiffRateBound (DiffRate* pRate);
void ClipDiffRate (DiffRate* pRate, const int newSize);
void EnlargeDiffRate (DiffRate* pRate, const int newSize);
void ModifyDiffRate (DiffRate* pRate, const int newSize);
void ReadDiffRate (DiffRate* pRate, const char* filename);
void InvalidateDiffRateBound (DiffRate* pRate);

#ifdef __cplusplus
}
#endif

#endif

