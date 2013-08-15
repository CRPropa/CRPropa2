#include <stdio.h>
#include <stdlib.h>
#include "error.h"

void Error(const char* errorMessage, const ErrorCode errorCode)
{
    printf("%s\n", errorMessage);
    exit (errorCode);
}
