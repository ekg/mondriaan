#ifndef SOLVER_H_DONE

#define SOLVER_H_DONE

#include <stdio.h>
#include <string.h>
#include "options.h"
#include "matrix.h"
#include "solution.h"

/* Declaration for externally used function */
void solve(struct mat *a, struct solution *sol);

#endif
