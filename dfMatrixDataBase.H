#pragma once

#include <stdio.h>

enum BoundaryConditions{
    zeroGradient,
    fixedValue,
    coupled,
    empty,
    gradientEnergy,
    calculated,
    cyclic,
    processor,
    extrapolated,
    fixedEnergy,
    processorCyclic
};

// constant values -- basic
int num_cells = 0;

int num_surfaces = 0;
int num_boundary_surfaces = 0;

int num_patches = 0;

int *patch_size = NULL;

double rdelta_t = 0;

// constant indexes
int *owner = NULL;
int *neighbor = NULL;

// constant fields - internal
double *weight = NULL;

//constant fields - internal
double *volume = NULL;

// non-constant fields - internal 
// fields solved by eqns - new
double *rho = NULL;
double *u = NULL;

// fields solved by eqns - old 
double *rho_old = NULL;

// other shared fields between eqns
double *phi = NULL;

// other shared fields between eqns
double *boundary_phi = NULL;
