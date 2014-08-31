#if !defined(PRNG_H_)
#define PRNG_H_

void init_urandom (void);
void free_urandom (void);
void uniform_random (size_t, double *);
double poisson (const double /* lambda */ ,
		const double * /* random, (0,1), array */ ,
		size_t * /* index into random */ ,
		const size_t /* length of random */ );
double exponential (const double /* lambda */ ,
		    const double /* random, (0,1) */ );

#endif /* PRNG_H_ */
