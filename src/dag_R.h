#ifndef IDLBNS_DAG_R_H
#define IDLBNS_DAG_R_H

#include <Rinternals.h>
#include "dag.h"

/*
 * Validate an external pointer and return the DAG it holds, or raise an R
 * error. Shared by dag_R.c and nbhd.c so that every entry point checks the
 * tag and the NULL address the same way -- a saved and reloaded external
 * pointer comes back with a NULL address, and this turns what would be a
 * segfault into an error message.
 */
idl_dag *idlBNs_dag_from_extptr(SEXP x);

#endif /* IDLBNS_DAG_R_H */
