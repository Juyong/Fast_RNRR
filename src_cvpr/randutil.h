#ifndef INCLUDED_RANDUTIL_MODULE
#define INCLUDED_RANDUTIL_MODULE

#ifdef __cplusplus
extern "C" {
#endif

extern int randbiased (double x);
extern size_t randslot (const double slots[/* n-1 */], size_t n);
extern size_t randrange (size_t n);
extern double drand (void);
extern unsigned int randreseed (void);

#ifdef __cplusplus
}
#endif

#endif
