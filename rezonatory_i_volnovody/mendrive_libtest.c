// test_load.c
#include <dlfcn.h>
#include <stdio.h>

int main() {
    void* h = dlopen("./mendrive.so", RTLD_NOW);
    if (!h) {
        printf("dlopen failed: %s\n", dlerror());
        return 1;
    }
    void* f = dlsym(h, "compute_det_contours");
    if (!f) {
        printf("symbol not found: %s\n", dlerror());
        return 1;
    }
    printf("Success!\n");
    dlclose(h);
    return 0;
}