#include <math.h>

typedef struct sQ
{
	double x;
	double y;
	double z;
	double w;
}sQ;

typedef struct sV
{
	double x;
	double y;
	double z;
}sV;

double v_norm(sV *v);
double v_norm(sV v);
double q_norm(sQ *q);
void q_renorm(sQ *q);
void v_renorm(sV *q);
sQ q_getconj(sQ *q);
void q_makeconj(sQ *q);
sQ q_mult(sQ *q1, sQ *q2);
void rotate_v(sQ *q, sV *v);
sQ v_mult(sQ *q1, sQ *q2);
sV v_mult(sV *v1, sV *v2);
double v_dot(sV *v1, sV *v2);