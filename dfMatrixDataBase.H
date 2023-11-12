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

// communication info
int *neighbProcNo = NULL;

// cyclic info
int *cyclicNeighbor = NULL;
int *patchSizeOffset = NULL;

// constant values -- basic
int num_cells = 0;

int num_surfaces = 0;
int num_boundary_surfaces = 0;

int num_patches = 0;

int *patch_size = NULL;
int *patch_type_calculated = NULL;
int *patch_type_extropolated = NULL;

double rdelta_t = 0;

// constant indexes
int *owner = NULL;
int *neighbor = NULL;

// constant fields - internal
double *sf = NULL;
double *mag_sf = NULL;
double *weight = NULL;

double *delta_coeffs = NULL;
double *volume = NULL;

// constant fields - boundary
double *boundary_sf = NULL;
double *boundary_mag_sf = NULL;
double *boundary_weight = NULL;
double *boundary_delta_coeffs = NULL;

double *boundary_face_cell = NULL;

/**
 * non-constant fields - internal
 * */  

// fields solved by eqns - new
double *rho = NULL;
double *u = NULL;

double *p =NULL;

// fields solved by eqns - old 
double *rho_old = NULL;

// other shared fields between eqns
double *phi = NULL;

// other shared fields between eqns - thermophysical
double *mu = NULL;

// internal fields used between eqns
double *rAU = NULL;

/**
 * non-constant fields - boundary
 * */

// fields solved by eqns - new
double *boundary_u = NULL;

double *boundary_p = NULL;

// other shared fields between eqns
double *boundary_phi = NULL;

// other shared fields between eqns - thermophysical
double *boundary_mu = NULL;

// boundary fields used between eqns
double *boundary_rAU = NULL;

