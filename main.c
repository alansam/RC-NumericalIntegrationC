// @see: https://rosettacode.org/wiki/Numerical_integration#C

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double int_leftrect(double from, double to, double n, double (*func)()) {
  double h = (to - from) / n;
  double sum = 0.0, x;
  for (x = from; x <= (to - h); x += h) {
    sum += func(x);
  }
  return h * sum;
}

double int_rightrect(double from, double to, double n, double (*func)()) {
  double h = (to - from) / n;
  double sum = 0.0, x;
  for (x = from; x <= (to - h); x += h) {
    sum += func(x + h);
  }
  return h * sum;
}

double int_midrect(double from, double to, double n, double (*func)()) {
  double h = (to - from) / n;
  double sum = 0.0, x;
  for (x = from; x <= (to - h); x += h) {
    sum += func(x + h / 2.0);
  }
  return h * sum;
}

double int_trapezium(double from, double to, double n, double (*func)()) {
   double h = (to - from) / n;
   double sum = func(from) + func(to);
   int i;
   for (i = 1;i < n;i++) {
     sum += 2.0 * func(from + i * h);
   }
   return  h * sum / 2.0;
}

double int_simpson(double from, double to, double n, double (*func)()) {
  double h = (to - from) / n;
  double sum1 = 0.0;
  double sum2 = 0.0;

  // double x;

  for (int i = 0; i < n; i++) {
    sum1 += func(from + h * i + h / 2.0);
  }

  for (int i = 1; i < n; i++) {
    sum2 += func(from + h * i);
  }

  return h / 6.0 * (func(from) + func(to) + 4.0 * sum1 + 2.0 * sum2);
}

/* test */
double f3(double x) {
  return x;
}

double f3a(double x) {
  return x * x / 2.0;
}

double f2(double x) {
  return 1.0 / x;
}

double f2a(double x) {
  return log(x);
}

double f1(double x) {
  return x * x * x;
}

double f1a(double x) {
  return x * x * x * x / 4.0;
}

typedef double (*pfunc)(double, double, double, double (*)());
typedef double (*rfunc)(double);

#define INTG(F, A, B) (F((B)) - F((A)))

typedef struct arguments arguments;
struct arguments {
  double ival1;
  double ival2;
  double approx;
};

typedef struct int_f int_f;
struct int_f {
  pfunc const func;
  char const name[16];
};

typedef struct rfuncs rfuncs;
struct rfuncs {
  rfunc const rf;
  rfunc const If;
};

/*
 *  MARK:  main()
 */
int main(int argc, char const * argv[]) {
  double ic;
  int_f int_f_list[] = {
    { int_leftrect,  "leftrect",  },
    { int_rightrect, "rightrect", },
    { int_midrect,   "midrect",   },
    { int_trapezium, "trapezium", },
    { int_simpson,   "simpson",   },
  };

  rfuncs rfl[] = {
    { .rf = f1, .If = f1a, },
    { .rf = f2, .If = f2a, },
    { .rf = f3, .If = f3a, },
    { .rf = f3, .If = f3a, },
};

  arguments args[] = {
    { 0.0,    1.0,     100.0, },
    { 1.0,  100.0,    1000.0, },
    { 0.0, 5000.0, 5000000.0, },
    { 0.0, 6000.0, 6000000.0, },
  };

  //  loop over integration functions
  for (size_t j_ = 0ul; j_ < (sizeof(rfl) / sizeof(*rfl)); j_++) {
    for (size_t i_ = 0ul; i_ < (sizeof(int_f_list) / sizeof(*int_f_list)); i_++) {
      ic = (*int_f_list[i_].func)(args[j_].ival1,
                                  args[j_].ival2,
                                  args[j_].approx,
                                  rfl[j_].rf);

    printf("%10s [ 0,1] num: %+16.6lf, an: %16.6lf\n",
           int_f_list[i_].name,
           ic,
           INTG((*rfl[j_].If),
                args[j_].ival1,
                args[j_].ival2));
    }
    putchar('\n');
  }
}
