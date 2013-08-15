#ifndef __check_ph__
#define __check_ph__

#ifdef __cplusplus
extern "C" {
#endif

void DumpArray (const dCVector* pVector);
void CheckIndex (const int lowerLimit, const int upperLimit, const int i,		const char* functionName);

#ifdef __cplusplus
}
#endif

#endif

