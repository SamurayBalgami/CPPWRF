#include <stdio.h>

void wrf_debug_(int* level, const char* msg, int len) {}
void wrf_message_(const char* msg, int len) {}
void wrf_error_fatal3_(const char* file, int* line, const char* msg, int flen, int mlen) {
    fprintf(stderr, "FATAL WRF ERROR\n");
}
void __module_wrf_error_MOD_wrf_err_message(const char* msg, int len) {}
