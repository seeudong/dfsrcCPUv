#include <stdio.h>

#include "dfUEqn.H"
#include "dfMatrixDataBase.H"


/**
 * fmv_ddt_vector
 **/
void fvm_ddt_vector(int num_cells, double rDeltaT,
        const double *rho, const double *rho_old, const double *vf, const double *volume,
        double *diag, double *source, double sign){
    
    //fvm_ddt_vector_kernel
    for(int index = 0; index < num_cells; index++){

        double vol = volume[index];
        double rho_old_kernel = rho_old[index];

        diag[index] += rDeltaT * rho[index] * vol * sign;

        source[0 + 3 * index] += rDeltaT * rho_old_kernel * vf[0 + 3 * index] * vol * sign;
        source[1 + 3 * index] += rDeltaT * rho_old_kernel * vf[1 + 3 * index] * vol * sign;
        source[2 + 3 * index] += rDeltaT * rho_old_kernel * vf[2 + 3 * index] * vol * sign;
    }
}


/**
 * fvm_div_vector
*/
void fvm_div_vector(int num_surfaces, int num_boundary_surfaces, 
        const int *lowerAddr, const int *upperAddr,
        const double *phi, const double *weight,
        double *lower, double *upper, double *diag, // end for internal
        int num_patches, const int *patch_size, const int *patch_type,
        const double *boundary_phi, const double *value_internal_coeffs, const double *value_boundary_coeffs,
        double *internal_coeffs, double *boundary_coeffs, double sign){
    
    //fvm_div_vector_internal
    for(int index = 0; index < num_cells; index++){

        double w = weight[index];
        double f = phi[index];

        double lower_value = (-w) * f * sign;
        double upper_value = (1-w) * f *sign;
        lower[index] += lower_value;
        upper[index] += upper_value;

        int owner = lowerAddr[index];
        int neighbor = upperAddr[index];
        diag[owner] -= lower_value;
        diag[neighbor] -= upper_value;
    }

    int offset = 0;
    for (int i = 0; i < num_patches; i++){
        if (patch_size[i] == 0) continue;

        int num = patch_size[i];
        //fvm_div_vector_boundary
        for (int index = 0; index < num; index++){
            int start_index = offset + index;
            double boundary_f = boundary_phi[start_index];
            internal_coeffs[0 + start_index * 3] += boundary_f * value_internal_coeffs[0 + start_index * 3] * sign;
            internal_coeffs[1 + start_index * 3] += boundary_f * value_internal_coeffs[1 + start_index * 3] * sign;
            internal_coeffs[2 + start_index * 3] += boundary_f * value_internal_coeffs[2 + start_index * 3] * sign;
            boundary_coeffs[0 + start_index * 3] -= boundary_f * value_boundary_coeffs[0 + start_index * 3] * sign;
            boundary_coeffs[1 + start_index * 3] -= boundary_f * value_boundary_coeffs[1 + start_index * 3] * sign;
            boundary_coeffs[2 + start_index * 3] -= boundary_f * value_boundary_coeffs[2 + start_index * 3] * sign;
        }
        if(patch_type[i] == processor|| patch_type[i] == processorCyclic) 
            offset += 2 + patch_size[i];
        else
            offset += 1 + patch_size[i];

    }

    
}


void process(){

    fvm_ddt_vector(num_cells, rdelta_t, rho, rho_old, u, volume, diag, source, 1.);

    fvm_div_vector(num_surfaces, num_boundary_surfaces, owner, neighbor,
                phi, weight,lower, upper, diag, // end for internal
                num_patches, patch_size, patch_type,
                boundary_phi, value_internal_coeffs, value_boundary_coeffs,
                internal_coeffs, boundary_coeffs, 1.);

}


   

