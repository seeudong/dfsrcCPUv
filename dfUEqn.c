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
    for(int index = 0; index < num_surfaces; index++){

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
            offset += 2 * num;
        else
            offset += num;

    }

}

/**
 * fvm_laplacian_vector
*/
void fvm_laplacian_vector(int num_surfaces, int num_boundary_surfaces, 
        const int *lowerAddr, const int *upperAddr,
        const double *weight, const double *mag_sf, const double *delta_coeffs, const double *gamma,
        double *lower, double *upper, double *diag, // end for internal
        int num_patches, const int *patch_size, const int *patch_type,
        const double *boundary_mag_sf, const double *boundary_gamma,
        const double *gradient_internal_coeffs, const double *gradient_boundary_coeffs,
        double *internal_coeffs, double *boundary_coeffs, double sign){
    
    // fvm_laplacian_internal
    for(int index = 0; index < num_surfaces; index++){

        int owner = lowerAddr[index];
        int neighbor = upperAddr[index];

        double w = weight[index];
        double face_gamma = w * gamma[owner] + (1 - w) * gamma[neighbor];
    
        // for fvm::laplacian, lower = upper
        double upper_value = face_gamma * mag_sf[index] * delta_coeffs[index];
        double lower_value = upper_value;

        lower_value = lower_value * sign;
        upper_value = upper_value * sign;

        lower[index] += lower_value;
        upper[index] += upper_value;

        diag[owner] -= lower_value;
        diag[neighbor] -= upper_value;

    }

    int offset = 0;
    for (int i = 0; i < num_patches; i++){
        if (patch_size[i] == 0) continue;

        int num = patch_size[i];
        //fvm_laplacian_vector_boundary
        for (int index = 0; index < num; index++){
            int start_index = offset + index;
            double boundary_value = boundary_gamma[start_index] * boundary_mag_sf[start_index];
            internal_coeffs[0 + start_index * 3] += boundary_value * gradient_internal_coeffs[0 + start_index * 3] * sign;
            internal_coeffs[1 + start_index * 3] += boundary_value * gradient_internal_coeffs[1 + start_index * 3] * sign;
            internal_coeffs[2 + start_index * 3] += boundary_value * gradient_internal_coeffs[2 + start_index * 3] * sign;
            boundary_coeffs[0 + start_index * 3] -= boundary_value * gradient_internal_coeffs[0 + start_index * 3] * sign;
            boundary_coeffs[1 + start_index * 3] -= boundary_value * gradient_internal_coeffs[1 + start_index * 3] * sign;
            boundary_coeffs[2 + start_index * 3] -= boundary_value * gradient_internal_coeffs[2 + start_index * 3] * sign;
        }
        if(patch_type[i] == processor|| patch_type[i] == processorCyclic) 
            offset += 2 * num;
        else
            offset += num;

    }
}

/**
 * fvc_grad_vector
*/
void fvc_grad_vector(int num_cells, int num_surfaces, int num_boundary_surfaces,
        const int *neighbor_peer, const int *lowerAddr, const int *upperAddr, 
        const double *weight, const double *Sf, const double *vf, double *output, // end for internal
        int num_patches, const int *patch_size, const int *patch_type,
        const int *boundary_cell_face, const double *boundary_vf, const double *boundary_Sf, const double *boundary_weight, 
        const double *volume, const double *boundary_mag_Sf, double *boundary_output,
        const int *cyclicNeighbor, const int *patchSizeOffset,
        const double *boundary_deltaCoeffs, double sign){
            
    // fvc_grad_vector_internal
    for (int index = 0; index < num_surfaces; index++){

                double w = weight[index];
                double Sfx = Sf[0 + index * 3];
                double Sfy = Sf[1 + index * 3];
                double Sfz = Sf[2 + index * 3];

                int owner = lowerAddr[index];
                int neighbor = upperAddr[index];

                double ssfx = (w * (vf[0 + owner * 3] - vf[0 + neighbor * 3]) + vf[0 + neighbor * 3]);
                double ssfy = (w * (vf[1 + owner * 3] - vf[1 + neighbor * 3]) + vf[1 + neighbor * 3]);
                double ssfz = (w * (vf[2 + owner * 3] - vf[2 + neighbor * 3]) + vf[2 + neighbor * 3]);
                
                double grad_xx = Sfx * ssfx;
                double grad_xy = Sfx * ssfy;
                double grad_xz = Sfx * ssfz;
                double grad_yx = Sfy * ssfx;
                double grad_yy = Sfy * ssfy;
                double grad_yz = Sfy * ssfz;
                double grad_zx = Sfz * ssfx;
                double grad_zy = Sfz * ssfy;
                double grad_zz = Sfz * ssfz;

                output[0 + owner * 9] += grad_xx;
                output[0 + neighbor * 9] -= grad_xx;
                output[1 + owner * 9] += grad_xy;
                output[1 + neighbor * 9] -= grad_xy;
                output[2 + owner * 9] += grad_xz;
                output[2 + neighbor * 9] -= grad_xz;
                output[3 + owner * 9] += grad_yx;
                output[3 + neighbor * 9] -= grad_yx;
                output[4 + owner * 9] += grad_yy;
                output[4 + neighbor * 9] -= grad_yy;
                output[5 + owner * 9] += grad_yz;
                output[5 + neighbor * 9] -= grad_yz;
                output[6 + owner * 9] += grad_zx;
                output[6 + neighbor * 9] -= grad_zx;
                output[7 + owner * 9] += grad_zy;
                output[7 + neighbor * 9] -= grad_zy;
                output[8 + owner * 9] += grad_zz;
                output[8 + neighbor * 9] -= grad_zz;

    }

    int offset = 0;

    // finish conctruct grad field except dividing cell volume
    for (int i = 0; i < num_patches; i++) {
        if (patch_size[i] == 0) continue;
        if (patch_type[i] == zeroGradient || patch_type[i] == fixedValue
                || patch_type[i] == calculated || patch_type[i] == cyclic){
            
            // fvc_grad_vector_boundary_zeroGradient

            int num = patch_size[i];
            for (int index = 0; index < num; index++){
                int start_index = offset + index;

                double bouSfx = boundary_Sf[0 + start_index * 3];
                double bouSfy = boundary_Sf[1 + start_index * 3];
                double bouSfz = boundary_Sf[2 + start_index * 3];

                double boussfx = boundary_vf[0 + start_index * 3];
                double boussfy = boundary_vf[1 + start_index * 3];
                double boussfz = boundary_vf[2 + start_index * 3];

                int cellIndex = boundary_cell_face[start_index];

                double grad_xx = bouSfx * boussfx;
                double grad_xy = bouSfx * boussfy;
                double grad_xz = bouSfx * boussfz;
                double grad_yx = bouSfy * boussfx;
                double grad_yy = bouSfy * boussfy;
                double grad_yz = bouSfy * boussfz;
                double grad_zx = bouSfz * boussfx;
                double grad_zy = bouSfz * boussfy;
                double grad_zz = bouSfz * boussfz;

                output[0 + cellIndex * 9] += grad_xx;
                output[1 + cellIndex * 9] += grad_xy;
                output[2 + cellIndex * 9] += grad_xz;
                output[3 + cellIndex * 9] += grad_yx;
                output[4 + cellIndex * 9] += grad_yy;
                output[5 + cellIndex * 9] += grad_yz;
                output[6 + cellIndex * 9] += grad_zx;
                output[7 + cellIndex * 9] += grad_zy;
                output[8 + cellIndex * 9] += grad_zz;

            }
        } else if (patch_type[i] == processor || patch_type[i] == processorCyclic){

            // fvc_grad_vector_boundary_processor
            int num = patch_size[i];
            for (int index = 0; index < num; index++){

                int neighbor_start_index = offset + index;
                int internal_start_index = offset + num + index;

                double bouWeight = boundary_weight[neighbor_start_index];

                double bouSfx = boundary_Sf[0 + neighbor_start_index * 3];
                double bouSfy = boundary_Sf[1 + neighbor_start_index * 3];
                double bouSfz = boundary_Sf[2 + neighbor_start_index * 3];

                double boussfx = (1 - bouWeight) * boundary_vf[0 + neighbor_start_index * 3] + 
                        bouWeight * boundary_vf[0 + internal_start_index * 3];
                double boussfy = (1 - bouWeight) * boundary_vf[1 + neighbor_start_index * 3] + 
                        bouWeight * boundary_vf[1 + internal_start_index * 3];
                double boussfz = (1 - bouWeight) * boundary_vf[2 + neighbor_start_index * 3] + 
                        bouWeight * boundary_vf[2 + internal_start_index * 3];

                int cellIndex = boundary_cell_face[neighbor_start_index];

                double grad_xx = bouSfx * boussfx;
                double grad_xy = bouSfx * boussfy;
                double grad_xz = bouSfx * boussfz;
                double grad_yx = bouSfy * boussfx;
                double grad_yy = bouSfy * boussfy;
                double grad_yz = bouSfy * boussfz;
                double grad_zx = bouSfz * boussfx;
                double grad_zy = bouSfz * boussfy;
                double grad_zz = bouSfz * boussfz;

                output[0 + cellIndex * 9] += grad_xx;
                output[1 + cellIndex * 9] += grad_xy;
                output[2 + cellIndex * 9] += grad_xz;
                output[3 + cellIndex * 9] += grad_yx;
                output[4 + cellIndex * 9] += grad_yy;
                output[5 + cellIndex * 9] += grad_yz;
                output[6 + cellIndex * 9] += grad_zx;
                output[7 + cellIndex * 9] += grad_zy;
                output[8 + cellIndex * 9] += grad_zz;

            }
            offset += 2 * patch_size[i]; // patchNeighbourFields and patchInternalFields
            continue;
        } else {
            fprintf(stderr, "%s %d, boundaryConditions other than zeroGradient are not support yet!\n", __FILE__, __LINE__);
        }
        offset += patch_size[i];
    }
    
    // divide cell volume

    // divide_cell_volume_tsr
    for (int index = 0; index < num_cells; index++){

        double vol = volume[index];
        output[0 + index * 9] = output[0 + index * 9] / vol;
        output[1 + index * 9] = output[1 + index * 9] / vol;
        output[2 + index * 9] = output[2 + index * 9] / vol;
        output[3 + index * 9] = output[3 + index * 9] / vol;
        output[4 + index * 9] = output[4 + index * 9] / vol;
        output[5 + index * 9] = output[5 + index * 9] / vol;
        output[6 + index * 9] = output[6 + index * 9] / vol;
        output[7 + index * 9] = output[7 + index * 9] / vol;
        output[8 + index * 9] = output[8 + index * 9] / vol;

    }

    // correct boundary conditions
    offset = 0;
    for (int i = 0; i < num_patches; i++) {
        if (patch_size[i] == 0) continue;
        if (patch_type[i] == zeroGradient){
            
            // fvc_grad_vector_correctBC_zeroGradient
            int num = patch_size[i];
            for (int index = 0; index < num; index++){

                int start_index = offset + index;

                int cellIndex = boundary_cell_face[start_index];

                double grad_xx = output[0 + cellIndex * 9];
                double grad_xy = output[1 + cellIndex * 9];
                double grad_xz = output[2 + cellIndex * 9];
                double grad_yx = output[3 + cellIndex * 9];
                double grad_yy = output[4 + cellIndex * 9];
                double grad_yz = output[5 + cellIndex * 9];
                double grad_zx = output[6 + cellIndex * 9];
                double grad_zy = output[7 + cellIndex * 9];
                double grad_zz = output[8 + cellIndex * 9];

                double n_x = boundary_sf[0 + start_index * 3] / boundary_mag_sf[start_index];
                double n_y = boundary_sf[1 + start_index * 3] / boundary_mag_sf[start_index];
                double n_z = boundary_sf[2 + start_index * 3] / boundary_mag_sf[start_index];

                double grad_correction_x = - (n_x * grad_xx + n_y * grad_yx + n_z * grad_zx); // sn_grad_x = 0
                double grad_correction_y = - (n_x * grad_xy + n_y * grad_yy + n_z * grad_zy);
                double grad_correction_z = - (n_x * grad_xz + n_y * grad_yz + n_z * grad_zz);

                boundary_output[0 + start_index * 9] = grad_xx + n_x * grad_correction_x;
                boundary_output[1 + start_index * 9] = grad_xy + n_x * grad_correction_y;
                boundary_output[2 + start_index * 9] = grad_xz + n_x * grad_correction_z;
                boundary_output[3 + start_index * 9] = grad_yx + n_y * grad_correction_x;
                boundary_output[4 + start_index * 9] = grad_yy + n_y * grad_correction_y;
                boundary_output[5 + start_index * 9] = grad_yz + n_y * grad_correction_z;
                boundary_output[6 + start_index * 9] = grad_zx + n_z * grad_correction_x;
                boundary_output[7 + start_index * 9] = grad_zy + n_z * grad_correction_y;
                boundary_output[8 + start_index * 9] = grad_zz + n_z * grad_correction_z;

            }

        } else if (patch_type[i] == fixedValue){

            // fvc_grad_vector_correctBC_fixedValue
            int num = patch_size[i];
            for (int index = 0; index < num; index++){

                int start_index = offset + index;

                int cellIndex = boundary_cell_face[start_index];

                double grad_xx = output[0 + cellIndex * 9];
                double grad_xy = output[1 + cellIndex * 9];
                double grad_xz = output[2 + cellIndex * 9];
                double grad_yx = output[3 + cellIndex * 9];
                double grad_yy = output[4 + cellIndex * 9];
                double grad_yz = output[5 + cellIndex * 9];
                double grad_zx = output[6 + cellIndex * 9];
                double grad_zy = output[7 + cellIndex * 9];
                double grad_zz = output[8 + cellIndex * 9];

                double n_x = boundary_sf[0 + start_index * 3] / boundary_mag_sf[start_index];
                double n_y = boundary_sf[1 + start_index * 3] / boundary_mag_sf[start_index];
                double n_z = boundary_sf[2 + start_index * 3] / boundary_mag_sf[start_index];

                // sn_grad: solving according to fixedValue BC
                double sn_grad_x = boundary_deltaCoeffs[start_index] * (boundary_vf[0 + start_index * 3] - vf[0 + cellIndex * 3]);
                double sn_grad_y = boundary_deltaCoeffs[start_index] * (boundary_vf[1 + start_index * 3] - vf[1 + cellIndex * 3]);
                double sn_grad_z = boundary_deltaCoeffs[start_index] * (boundary_vf[2 + start_index * 3] - vf[2 + cellIndex * 3]);

                double grad_correction_x = sn_grad_x - (n_x * grad_xx + n_y * grad_yx + n_z * grad_zx); // sn_grad_x = 0
                double grad_correction_y = sn_grad_y - (n_x * grad_xy + n_y * grad_yy + n_z * grad_zy);
                double grad_correction_z = sn_grad_z - (n_x * grad_xz + n_y * grad_yz + n_z * grad_zz);

                boundary_output[0 + start_index * 9] = grad_xx + n_x * grad_correction_x;
                boundary_output[1 + start_index * 9] = grad_xy + n_x * grad_correction_y;
                boundary_output[2 + start_index * 9] = grad_xz + n_x * grad_correction_z;
                boundary_output[3 + start_index * 9] = grad_yx + n_y * grad_correction_x;
                boundary_output[4 + start_index * 9] = grad_yy + n_y * grad_correction_y;
                boundary_output[5 + start_index * 9] = grad_yz + n_y * grad_correction_z;
                boundary_output[6 + start_index * 9] = grad_zx + n_z * grad_correction_x;
                boundary_output[7 + start_index * 9] = grad_zy + n_z * grad_correction_y;
                boundary_output[8 + start_index * 9] = grad_zz + n_z * grad_correction_z;

            }
        } else if (patch_type[i] == processor || patch_type[i] == processorCyclic){

            // fvc_grad_vector_correctBC_processor
            int num = patch_size[i];

            int neighbor_start_index_tmp = offset;
            int internal_start_index_tmp = offset + num;

            // correct_internal_boundary_field_tensor
            for (int index = 0; index < num; index++){

                int neighbor_start_index = offset + index;
                int internal_start_index = offset + num + index;

                int cellIndex = boundary_cell_face[neighbor_start_index];

                boundary_output[0 + internal_start_index * 9] = output[0 + cellIndex * 9];
                boundary_output[1 + internal_start_index * 9] = output[1 + cellIndex * 9];
                boundary_output[2 + internal_start_index * 9] = output[2 + cellIndex * 9];
                boundary_output[3 + internal_start_index * 9] = output[3 + cellIndex * 9];
                boundary_output[4 + internal_start_index * 9] = output[4 + cellIndex * 9];
                boundary_output[5 + internal_start_index * 9] = output[5 + cellIndex * 9];
                boundary_output[6 + internal_start_index * 9] = output[6 + cellIndex * 9];
                boundary_output[7 + internal_start_index * 9] = output[7 + cellIndex * 9];
                boundary_output[8 + internal_start_index * 9] = output[8 + cellIndex * 9];
                
            }

            //TODO: nccl

            offset += 2 * patch_size[i]; // patchNeighbourFields and patchInternalFields
            continue;
        } else if (patch_type[i] == cyclic){

            // fvc_grad_vector_correctBC_cyclic
            int num = patch_size[i];
            for (int index = 0; index < num; index++){

                int internal_start_index = offset + index;
                int neighbor_start_index = patchSizeOffset[cyclicNeighbor[i]] + index;

                double weight = boundary_weight[internal_start_index];

                int internal_cellIndex = boundary_cell_face[internal_start_index];
                int neighbor_cellIndex = boundary_cell_face[neighbor_start_index];

                boundary_output[0 + internal_start_index * 9] = weight * output[0 + internal_cellIndex * 9] 
                        + (1 - weight) * output[0 + neighbor_cellIndex * 9];
                boundary_output[1 + internal_start_index * 9] = weight * output[1 + internal_cellIndex * 9]
                        + (1 - weight) * output[1 + neighbor_cellIndex * 9];
                boundary_output[2 + internal_start_index * 9] = weight * output[2 + internal_cellIndex * 9]
                        + (1 - weight) * output[2 + neighbor_cellIndex * 9];
                boundary_output[3 + internal_start_index * 9] = weight * output[3 + internal_cellIndex * 9]
                        + (1 - weight) * output[3 + neighbor_cellIndex * 9];
                boundary_output[4 + internal_start_index * 9] = weight * output[4 + internal_cellIndex * 9]
                        + (1 - weight) * output[4 + neighbor_cellIndex * 9];
                boundary_output[5 + internal_start_index * 9] = weight * output[5 + internal_cellIndex * 9]
                        + (1 - weight) * output[5 + neighbor_cellIndex * 9];
                boundary_output[6 + internal_start_index * 9] = weight * output[6 + internal_cellIndex * 9]
                        + (1 - weight) * output[6 + neighbor_cellIndex * 9];
                boundary_output[7 + internal_start_index * 9] = weight * output[7 + internal_cellIndex * 9]
                        + (1 - weight) * output[7 + neighbor_cellIndex * 9];
                boundary_output[8 + internal_start_index * 9] = weight * output[8 + internal_cellIndex * 9]
                        + (1 - weight) * output[8 + neighbor_cellIndex * 9];

            }

        } else {
            fprintf(stderr, "%s %d, boundaryConditions other than zeroGradient are not support yet!\n", __FILE__, __LINE__);
        }
        offset += patch_size[i];
    }

}

/**
 * scale_dev2T_tensor
*/
void scale_dev2T_tensor(int num_cells, const double *vf1, double *vf2,
        int num_boundary_surfaces, const double *boundary_vf1, double *boundary_vf2){
    
    // scale_dev2t_tensor_kernel
    for (int index = 0; index < num_cells; index++){

        double scale = vf1[index];
        double val_xx = vf2[0 + index * 9];
        double val_xy = vf2[1 + index * 9];
        double val_xz = vf2[2 + index * 9];
        double val_yx = vf2[3 + index * 9];
        double val_yy = vf2[4 + index * 9];
        double val_yz = vf2[5 + index * 9];
        double val_zx = vf2[6 + index * 9];
        double val_zy = vf2[7 + index * 9];
        double val_zz = vf2[8 + index * 9];
        double trace_coeff = (2. / 3.) * (val_xx + val_yy + val_zz);
        vf2[0 + index * 9] = scale * (val_xx - trace_coeff);
        vf2[1 + index * 9] = scale * val_yx;
        vf2[2 + index * 9] = scale * val_zx;
        vf2[3 + index * 9] = scale * val_xy;
        vf2[4 + index * 9] = scale * (val_yy - trace_coeff);
        vf2[5 + index * 9] = scale * val_zy;
        vf2[6 + index * 9] = scale * val_xz;
        vf2[7 + index * 9] = scale * val_yz;
        vf2[8 + index * 9] = scale * (val_zz - trace_coeff);
    }

    // scale_dev2t_tensor_kernel (boundary)
    for (int index = 0; index < num_boundary_surfaces; index++){

        double scale = boundary_vf1[index];
        double val_xx = boundary_vf2[0 + index * 9];
        double val_xy = boundary_vf2[1 + index * 9];
        double val_xz = boundary_vf2[2 + index * 9];
        double val_yx = boundary_vf2[3 + index * 9];
        double val_yy = boundary_vf2[4 + index * 9];
        double val_yz = boundary_vf2[5 + index * 9];
        double val_zx = boundary_vf2[6 + index * 9];
        double val_zy = boundary_vf2[7 + index * 9];
        double val_zz = boundary_vf2[8 + index * 9];
        double trace_coeff = (2. / 3.) * (val_xx + val_yy + val_zz);
        boundary_vf2[0 + index * 9] = scale * (val_xx - trace_coeff);
        boundary_vf2[1 + index * 9] = scale * val_yx;
        boundary_vf2[2 + index * 9] = scale * val_zx;
        boundary_vf2[3 + index * 9] = scale * val_xy;
        boundary_vf2[4 + index * 9] = scale * (val_yy - trace_coeff);
        boundary_vf2[5 + index * 9] = scale * val_zy;
        boundary_vf2[6 + index * 9] = scale * val_xz;
        boundary_vf2[7 + index * 9] = scale * val_yz;
        boundary_vf2[8 + index * 9] = scale * (val_zz - trace_coeff);
    }

}

/**
 * fvc_div_cell_tensor
*/
void fvc_div_cell_tensor(int num_cells, int num_surfaces, int num_boundary_surfaces, 
        const int *lowerAddr, const int *upperAddr,
        const double *weight, const double *Sf, const double *vf, double *output, // end for internal
        int num_patches, const int *patch_size, const int *patch_type, const double *boundary_weight, 
        const int *boundary_cell_face, const double *boundary_vf, const double *boundary_Sf,
        const double *volume, double sign){
    
    // fvc_div_cell_tensor_internal
    for (int index = 0; index < num_surfaces; index++){
        
        double w = weight[index];
        double Sfx = Sf[0 + index * 3];
        double Sfy = Sf[1 + index * 3];
        double Sfz = Sf[2 + index * 3];
        int owner = lowerAddr[index];
        int neighbor = upperAddr[index];

        double ssf_xx = (w * (vf[0 + owner * 9] - vf[0 + neighbor * 9]) + vf[0 + neighbor * 9]);
        double ssf_xy = (w * (vf[1 + owner * 9] - vf[1 + neighbor * 9]) + vf[1 + neighbor * 9]);
        double ssf_xz = (w * (vf[2 + owner * 9] - vf[2 + neighbor * 9]) + vf[2 + neighbor * 9]);
        double ssf_yx = (w * (vf[3 + owner * 9] - vf[3 + neighbor * 9]) + vf[3 + neighbor * 9]);
        double ssf_yy = (w * (vf[4 + owner * 9] - vf[4 + neighbor * 9]) + vf[4 + neighbor * 9]);
        double ssf_yz = (w * (vf[5 + owner * 9] - vf[5 + neighbor * 9]) + vf[5 + neighbor * 9]);
        double ssf_zx = (w * (vf[6 + owner * 9] - vf[6 + neighbor * 9]) + vf[6 + neighbor * 9]);
        double ssf_zy = (w * (vf[7 + owner * 9] - vf[7 + neighbor * 9]) + vf[7 + neighbor * 9]);
        double ssf_zz = (w * (vf[8 + owner * 9] - vf[8 + neighbor * 9]) + vf[8 + neighbor * 9]);
        double div_x = (Sfx * ssf_xx + Sfy * ssf_yx + Sfz * ssf_zx) * sign;
        double div_y = (Sfx * ssf_xy + Sfy * ssf_yy + Sfz * ssf_zy) * sign;
        double div_z = (Sfx * ssf_xz + Sfy * ssf_yz + Sfz * ssf_zz) * sign;

    }
    
    int offset = 0;
    for (int i = 0; i < num_patches; i++) {
        if (patch_size[i] == 0) continue;
        if (patch_type[i] == zeroGradient || patch_type[i] == fixedValue
                || patch_type[i] == calculated || patch_type[i] == cyclic) {
            
            // fvc_div_cell_tensor_boundary_zeroGradient
            for (int index = 0; index < patch_size[i]; index++){

                int start_index = offset + index;

                double bouSfx = boundary_Sf[0 + start_index * 3];
                double bouSfy = boundary_Sf[1 + start_index * 3];
                double bouSfz = boundary_Sf[2 + start_index * 3];

                double boussf_xx = boundary_vf[0 + start_index * 9];
                double boussf_xy = boundary_vf[1 + start_index * 9];
                double boussf_xz = boundary_vf[2 + start_index * 9];
                double boussf_yx = boundary_vf[3 + start_index * 9];
                double boussf_yy = boundary_vf[4 + start_index * 9];
                double boussf_yz = boundary_vf[5 + start_index * 9];
                double boussf_zx = boundary_vf[6 + start_index * 9];
                double boussf_zy = boundary_vf[7 + start_index * 9];
                double boussf_zz = boundary_vf[8 + start_index * 9];
                int cellIndex = boundary_cell_face[start_index];

                double bouDiv_x = (bouSfx * boussf_xx + bouSfy * boussf_yx + bouSfz * boussf_zx) * sign;
                double bouDiv_y = (bouSfx * boussf_xy + bouSfy * boussf_yy + bouSfz * boussf_zy) * sign;
                double bouDiv_z = (bouSfx * boussf_xz + bouSfy * boussf_yz + bouSfz * boussf_zz) * sign;

                output[0 + cellIndex * 3] += bouDiv_x;
                output[1 + cellIndex * 3] += bouDiv_y;
                output[2 + cellIndex * 3] += bouDiv_z;
            }
        } else if (patch_type[i] == processor || patch_type[i] == processorCyclic){

            //fvc_div_cell_tensor_boundary_processor
            for(int index = 0; index < patch_size[i]; index++){

                int neighbor_start_index = offset + index;
                int internal_start_index = offset + patch_size[i] + index;
                
                double bouWeight = boundary_weight[neighbor_start_index];

                double bouSfx = boundary_Sf[0 + neighbor_start_index * 3];
                double bouSfy = boundary_Sf[1 + neighbor_start_index * 3];
                double bouSfz = boundary_Sf[2 + neighbor_start_index * 3];

                double boussf_xx = (1 - bouWeight) * boundary_vf[0 + neighbor_start_index * 9] + 
                        bouWeight * boundary_vf[0 + internal_start_index * 9];
                double boussf_xy = (1 - bouWeight) * boundary_vf[1 + neighbor_start_index * 9] +
                        bouWeight * boundary_vf[1 + internal_start_index * 9];
                double boussf_xz = (1 - bouWeight) * boundary_vf[2 + neighbor_start_index * 9] +
                        bouWeight * boundary_vf[2 + internal_start_index * 9];
                double boussf_yx = (1 - bouWeight) * boundary_vf[3 + neighbor_start_index * 9] +
                        bouWeight * boundary_vf[3 + internal_start_index * 9];
                double boussf_yy = (1 - bouWeight) * boundary_vf[4 + neighbor_start_index * 9] +
                        bouWeight * boundary_vf[4 + internal_start_index * 9];
                double boussf_yz = (1 - bouWeight) * boundary_vf[5 + neighbor_start_index * 9] +
                        bouWeight * boundary_vf[5 + internal_start_index * 9];
                double boussf_zx = (1 - bouWeight) * boundary_vf[6 + neighbor_start_index * 9] +
                        bouWeight * boundary_vf[6 + internal_start_index * 9];
                double boussf_zy = (1 - bouWeight) * boundary_vf[7 + neighbor_start_index * 9] +
                        bouWeight * boundary_vf[7 + internal_start_index * 9];
                double boussf_zz = (1 - bouWeight) * boundary_vf[8 + neighbor_start_index * 9] + 
                        bouWeight * boundary_vf[8 + internal_start_index * 9];
                int cellIndex = boundary_cell_face[neighbor_start_index];

                double bouDiv_x = (bouSfx * boussf_xx + bouSfy * boussf_yx + bouSfz * boussf_zx) * sign;
                double bouDiv_y = (bouSfx * boussf_xy + bouSfy * boussf_yy + bouSfz * boussf_zy) * sign;
                double bouDiv_z = (bouSfx * boussf_xz + bouSfy * boussf_yz + bouSfz * boussf_zz) * sign;

                output[0 + cellIndex * 3] += bouDiv_x;
                output[1 + cellIndex * 3] += bouDiv_y;
                output[2 + cellIndex * 3] += bouDiv_z;

            }
            offset += 2 * patch_size[i];
            continue;
        } else {
            fprintf(stderr, "%s %d, boundaryConditions other than zeroGradient are not support yet!\n", __FILE__, __LINE__);
        }
        offset += patch_size[i];
        
    }

}

/**
 * fvc_grad_cell_scalar
*/
void fvc_grad_cell_scalar(int num_cells, int num_surfaces, int num_boundary_surfaces, 
        const int *lowerAddr, const int *upperAddr, 
        const double *weight, const double *Sf, const double *vf, double *output, // end for internal
        int num_patches, const int *patch_size, const int *patch_type, const double *boundary_weight,
        const int *boundary_cell_face, const double *boundary_vf, const double *boundary_Sf, const double *volume, double sign){
    
    // fvc_grad_scalar_internal
    for (int index = 0; index <num_surfaces; index++){

        double w = weight[index];
        double Sfx = Sf[num_surfaces * 0 + index];
        double Sfy = Sf[num_surfaces * 1 + index];
        double Sfz = Sf[num_surfaces * 2 + index];

        int owner = lowerAddr[index];
        int neighbor = upperAddr[index];

        double ssf = (w * (vf[owner] - vf[neighbor]) + vf[neighbor]);

        double grad_x = Sfx * ssf * sign;
        double grad_y = Sfy * ssf * sign;
        double grad_z = Sfz * ssf * sign;

        // owner
        output[0 + owner * 3] += grad_x;
        output[1 + owner * 3] += grad_y;
        output[2 + owner * 3] += grad_z;

        // neighbour
        output[0 + neighbor * 3] -= grad_x;
        output[1 + neighbor * 3] -= grad_y;
        output[2 + neighbor * 3] -= grad_z;

    }

    int offset = 0;
    for (int i = 0; i < num_patches; i++){
        if (patch_size[i] == 0) continue;
        if (patch_type[i] == zeroGradient || patch_type[i] == fixedValue
                || patch_type[i] == calculated || patch_type[i] == cyclic){
            
            // fvc_grad_scalar_boundary_zeroGradient
            for (int index = 0; index < patch_size[i]; index++){

                int start_index = offset + index;

                double bouvf = boundary_vf[start_index];
                double bouSfx = boundary_Sf[0 + start_index * 3];
                double bouSfy = boundary_Sf[1 + start_index * 3];
                double bouSfz = boundary_Sf[2 + start_index * 3];

                int cellIndex = boundary_cell_face[start_index];

                double grad_x = bouSfx * bouvf;
                double grad_y = bouSfy * bouvf;
                double grad_z = bouSfz * bouvf;

                output[0 + cellIndex * 3] += grad_x * sign;
                output[1 + cellIndex * 3] += grad_y * sign;
                output[2 + cellIndex * 3] += grad_z * sign;

            }
        } else if (patch_type[i] == processor || patch_type[i] == processorCyclic){

            // fvc_grad_scalar_boundary_processor
            for (int index = 0; index < patch_size[i]; index++){

                int neighbor_start_index = offset + index;
                int internal_start_index = offset + patch_size[i] + index;

                double bouWeight = boundary_weight[neighbor_start_index];

                double bouSfx = boundary_Sf[0 + neighbor_start_index * 3];
                double bouSfy = boundary_Sf[1 + neighbor_start_index * 3];
                double bouSfz = boundary_Sf[2 + neighbor_start_index * 3];

                double bouvf = (1 - bouWeight) * boundary_vf[neighbor_start_index] + bouWeight * boundary_vf[internal_start_index];

                int cellIndex = boundary_cell_face[neighbor_start_index];

                double grad_x = bouSfx * bouvf;
                double grad_y = bouSfy * bouvf;
                double grad_z = bouSfz * bouvf;

                output[0 + cellIndex * 3] += grad_x * sign;
                output[1 + cellIndex * 3] += grad_y * sign;
                output[2 + cellIndex * 3] += grad_z * sign;

            }
            offset += 2 * patch_size[i];
            continue;
        } else {
            fprintf(stderr, "%s %d, boundaryConditions other than zeroGradient are not support yet!\n", __FILE__, __LINE__);
        }
        offset += patch_size[i];
    }

}

/**
 * getrAU
*/
void getrAU(int num_cells, int num_surfaces, int num_boundary_surfaces, 
        const int *neighbor_peer, int num_patches, const int *patch_size, const int *patch_type,
        const int *boundary_cell_face, const double *boundary_delta_coeffs, const double *internal_coeffs, const double *volume, 
        const double *diag, double *rAU, double *boundary_rAU){
    
    // addAveInternaltoDiagUeqn
    for (int index = 0; index < num_boundary_surfaces; index++){

        int cellIndex = boundary_cell_face[index];

        double internal_x = internal_coeffs[0 + index * 3];
        double internal_y = internal_coeffs[1 + index * 3];
        double internal_z = internal_coeffs[2 + index * 3];

        double ave_internal = (internal_x + internal_y + internal_z) / 3;

        rAU[cellIndex] += ave_internal;

    }

    // divide_cell_volume_scalar_reverse
    for (int index = 0; index < num_cells; index++){

        double vol = volume[index];

        rAU[index] = 1/ (rAU[index] / vol);
    }

    // correct_boundary_conditions_scalar
    int offset = 0;
    int gradient_offset = 0;
    for (int i = 0; i < num_patches; i++) {
        
        if (patch_size[i] == 0) continue;
        if(patch_type[i] == zeroGradient || patch_type[i] == extrapolated){
            
            // correct_boundary_conditions_zeroGradient_scalar
            for (int index = 0; index < patch_size[i]; index++){

                int start_index = offset + index;

                int cellIndex = boundary_cell_face[start_index];
                boundary_rAU[start_index] = rAU[cellIndex];
            }

        } else if (patch_type[i] == fixedValue || patch_type[i] == calculated){
            
            // No operation needed in this condition
        } else if (patch_type[i] == processor || patch_type[i] == processorCyclic){

            // correct_boundary_conditions_processor_scalar
            int neighbor_start_index_tmp = offset;
            int internal_start_index_tmp = offset + patch_size[i];

            // correct_internal_boundary_field_scalar
            for (int index = 0; index = patch_size[i]; index++){
                
                int neighbor_start_index = offset + index;
                int internal_start_index = offset + patch_size[i] + index;

                int cellIndex = boundary_cell_face[neighbor_start_index];
                boundary_rAU[internal_start_index] = rAU[cellIndex];

            }

            //TODO: nccl

            offset += 2 * patch_size[i];
            continue;
        } else if (patch_size[i] == gradientEnergy){

            // correct_boundary_conditions_gradientEnergy_scalar
            for (int index = 0; index < patch_size[i]; index++){

                int bou_start_index = offset + index;
                int gradient_start_index = gradient_offset + index;
                int cellIndex = boundary_cell_face[bou_start_index];

                boundary_rAU[bou_start_index] = rAU[cellIndex] + patchSizeOffset[gradient_start_index] / boundary_delta_coeffs[bou_start_index];

            }

            gradient_offset += patch_size[i];
        } else if (patch_size[i] == fixedEnergy){

            // GPUThermo->calculateEnthalpyGPU

            // ?

        } else if (patch_size[i] == cyclic){

            // correct_boundary_conditions_cyclic_scalar
            for (int index = 0; index < patch_size[i]; index++){

                int internal_start_index = offset + index;
                int neighbor_start_index = patchSizeOffset[cyclicNeighbor[i]] + index;

                double weight = boundary_weight[internal_start_index];

                int internal_cellIndex = boundary_cell_face[internal_start_index];
                int neighbor_cellIndex = boundary_cell_face[neighbor_start_index];

                boundary_rAU[internal_start_index] = weight * rAU[internal_cellIndex] + (1 - weight) * rAU[neighbor_cellIndex];

            }
        } else {
            fprintf(stderr, "%s %d, boundaryConditions %d are not support yet!\n", __FILE__, __LINE__, patch_type[i]);
        }
        offset += patch_size[i];
    }

}

void process(){

    fvm_ddt_vector(num_cells, rdelta_t, rho, rho_old, u, volume, diag, source, 1.);

    fvm_div_vector(num_surfaces, num_boundary_surfaces, owner, neighbor,
                phi, weight,lower, upper, diag, // end for internal
                num_patches, patch_size, patch_type,
                boundary_phi, value_internal_coeffs, value_boundary_coeffs,
                internal_coeffs, boundary_coeffs, 1.);

    fvm_laplacian_vector(num_surfaces, num_boundary_surfaces, 
               owner, neighbor,
               weight, mag_sf, delta_coeffs, mu,
               lower, upper, diag, // end for internal
               num_patches, patch_size, patch_type,
               boundary_mag_sf, boundary_mu,
               gradient_internal_coeffs, gradient_boundary_coeffs,
               internal_coeffs, boundary_coeffs, -1.);

    fvc_grad_vector(num_cells, num_surfaces, num_boundary_surfaces,
                neighbProcNo, owner, neighbor,
                weight, sf, u, grad_u,
                num_patches, patch_size, patch_type,
                boundary_face_cell, boundary_u, boundary_sf, boundary_weight, 
                volume, boundary_mag_sf, boundary_grad_u, cyclicNeighbor,
                patchSizeOffset, boundary_delta_coeffs, 1.);
    
    scale_dev2T_tensor(num_cells, mu, grad_u, // end for internal
                num_boundary_surfaces, boundary_mu, boundary_grad_u);

    fvc_div_cell_tensor(num_cells, num_surfaces, num_boundary_surfaces,
                owner, neighbor, weight, sf, grad_u, source, // end for internal
                num_patches, patch_size, patch_type_calculated, boundary_weight,
                boundary_face_cell, boundary_grad_u, boundary_sf, volume, 1.);

    fvc_grad_cell_scalar(num_cells, num_surfaces, num_boundary_surfaces,
                owner, neighbor,weight, sf, p, source,
                num_patches, patch_size, patch_type, boundary_weight,
                boundary_face_cell, boundary_p, boundary_sf, volume, -1.);

    getrAU(num_cells, num_surfaces, num_boundary_surfaces, 
                neighbProcNo, num_patches, patch_size, patch_type_extropolated,
                boundary_face_cell, boundary_delta_coeffs, internal_coeffs, volume, diag, 
                rAU, boundary_rAU);

}


   

