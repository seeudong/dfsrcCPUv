#pragma once

#include <stdio.h>

// constant fields - boundary
int *patch_type = NULL;

/**
 * non-constant fields - internal
*/

// intermediate fields
double *grad_u = NULL;

/**
 * non-constant fields - boundary
*/

// intermediate fields
double *boundary_grad_u = NULL;

// boundary coeff fields
double *value_internal_coeffs = NULL;
double *value_boundary_coeffs= NULL;
double *gradient_internal_coeffs= NULL;
double *gradient_boundary_coeffs= NULL;


// non-constant fields - ldu
double *lower = NULL;
double *upper = NULL;
double *diag = NULL;

double *source = NULL;
double *internal_coeffs = NULL;
double *boundary_coeffs = NULL;

void process();


