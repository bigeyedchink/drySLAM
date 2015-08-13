//
//  MapManagement.m
//  drySLAM
//
//  Created by Anthony Rodriguez on 4/17/15.
//  Copyright (c) 2015 ECE473 573 Common. All rights reserved.
//

#include "EKFUpdateLiInliers.hpp"


bool ekf_update_li_inliers(Filter* filter, NSMutableArray* features_info)
{
    //TODO: Find a more efficient way to do this:
    int counter = 0;
    for (int i=0; i<[features_info count]; ++i){
        Feature* feat = (Feature*)features_info[i];
        if (feat->low_innovation_inlier==1)
            ++counter;
    }
    //z and h are a single row of x,y coordinates
    double* z = (double*)malloc(counter*2*sizeof(double));
    double* h = (double*)malloc(counter*2*sizeof(double));
    double* H = (double*)malloc(counter*2*filter->state_size*sizeof(double));
    
    int current_element = 0;
    for (int i=0; i<[features_info count]; ++i){
        Feature* feat = ((Feature*)features_info[i]);
        if ( feat->low_innovation_inlier==1 ){
            z[current_element] = feat->z[0];
            z[current_element+1] = feat->z[1];
            h[current_element] = feat->h[0];
            h[current_element+1] = feat->h[1];
            double* matPtr = feat->H->ptr<double>(0);
            for (int j=0; j<(2*filter->state_size); ++j){
                H[current_element*filter->state_size+j] = matPtr[j];
            }
            current_element = current_element+2;
        }
    }
    
    //R is an identity matrix, counterxcounter size
    cv::Mat* R_mat = new cv::Mat(counter*2, counter*2, CV_64FC1, 0.0);
    for (int i=0; i<counter*2; ++i){
        double* matPtr = R_mat->ptr<double>(i) + i;
        *matPtr = 1.0;
    }
    
    cv::Mat* H_mat = new cv::Mat(counter*2, filter->state_size, CV_64FC1, H, cv::Mat::AUTO_STEP);
    cv::Mat* z_mat = new cv::Mat(counter*2, 1, CV_64FC1, z, cv::Mat::AUTO_STEP);
    cv::Mat* h_mat = new cv::Mat(counter*2, 1, CV_64FC1, h, cv::Mat::AUTO_STEP);
    
    //get new values for x_k_k and p_k_k
    for (int i=0; i<filter->state_size; ++i){
        filter->x_k_k[i] = filter->x_k_km1[i];
    }
    for (int i=0; i<(filter->state_size*filter->state_size); ++i){
        filter->p_k_k[i] = filter->p_k_km1[i];
    }
    
    bool success = update(filter->x_k_k, filter->p_k_k, filter->state_size, H_mat, R_mat, z_mat, h_mat);
    
    free(z);
    free(h);
    free(H);
    z_mat->release();
    h_mat->release();
    H_mat->release();
    R_mat->release();
    
    return success;
}

/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
/***************************************************************************************************************************************/


bool update(double* x_km1_k, double* p_km1_k, int state_size, cv::Mat* H, cv::Mat* R, cv::Mat* z, cv::Mat* h)
{
    if (x_km1_k==NULL || p_km1_k==NULL || H==NULL || R==NULL || z==NULL || h==NULL || state_size<=0)
        return false;
    
    if (z->cols==0)
        return true;
    
    //Miscellaneous variables used throughout function:
    int r = 0;
    int c = 0;
    double* mxPtr = NULL;
    
    //Convert p into an opencv matrix
    cv::Mat p_km1_k_mat(state_size, state_size, CV_64FC1, p_km1_k, cv::Mat::AUTO_STEP);
    //multiply H and p
    if (H->cols!=state_size){
        NSLog(@"Error: Cannot multiply matrix H(%dx%d) and p(%dx%d)", H->rows, H->cols, state_size, state_size);
        return false;
    }
    
    //S = H*p_km1_k*H' + R
    cv::Mat S = (*H)*p_km1_k_mat;
    cv::Mat H_transpose;
    cv::transpose(*H, H_transpose);
    S = S * H_transpose;
    cv::add(S, *R, S);
    
    //K = (p_km1_k*H') / S
    cv::Mat S_inverse;
    cv::invert(S, S_inverse);
    cv::Mat p_times_H_transpose = p_km1_k_mat * H_transpose;
    cv::Mat K = p_times_H_transpose * S_inverse;
    
    //Convert x to an opencv matrix
    cv::Mat x_km1_k_mat(state_size, 1, CV_64FC1, x_km1_k, cv::Mat::AUTO_STEP);
    
    //x_k_k = x_km1_k + K*(z-h)
    cv::Mat x_k_k_mat = x_km1_k_mat + (K * (*z - *h));
    
    cv::Mat K_transpose;
    cv::transpose(K, K_transpose);
    //p_k_k = p_km1_k - K*S*K'
    cv::Mat p_k_k_mat = p_km1_k_mat - (K * S * K_transpose);
    
    cv::Mat p_k_k_mat_transpose;
    cv::transpose(p_k_k_mat, p_k_k_mat_transpose);
    
    //p_k_k = 0.5*p_k_k + 0.5*p_k_k'
    p_k_k_mat = 0.5*p_k_k_mat + 0.5*p_k_k_mat_transpose;
    
    
    K.release();
    S.release();
    H_transpose.release();
    S_inverse.release();
    p_times_H_transpose.release();
    x_km1_k_mat.release();
    p_k_k_mat_transpose.release();
    
    
    double* Jnorm = (double*)malloc(16*sizeof(double));
    double* q = (double*)malloc(4*sizeof(double));
    //Set mxPtr to point to the first element of x_k_k row 0:
    
    //q = x_k_k(4:7) :
    q[0] = *(x_k_k_mat.ptr<double>(3));
    q[1] = *(x_k_k_mat.ptr<double>(4));
    q[2] = *(x_k_k_mat.ptr<double>(5));
    q[3] = *(x_k_k_mat.ptr<double>(6));
    normJac(Jnorm, q);
    free(q);
    cv::Mat Jnorm_mat(4, 4, CV_64FC1, Jnorm, cv::Mat::AUTO_STEP);
    cv::Mat JNorm_mat_transpose;
    cv::transpose(Jnorm_mat, JNorm_mat_transpose);
    
    //*******************************************************************
    //BUILDING p_k_k:
    //p_k_k is constructed of 9 smaller matrices, indexed by row_col here
    
    //0_0: p_k_k(1:3,1:3)   the first 3x3 matrix inside p_k_k:
    double* p_k_k_0_0 = (double*)malloc(9*sizeof(double));
    mxPtr = p_k_k_mat.ptr<double>(0);
    p_k_k_0_0[0] = mxPtr[0];
    p_k_k_0_0[1] = mxPtr[1];
    p_k_k_0_0[2] = mxPtr[2];
    mxPtr = p_k_k_mat.ptr<double>(1);
    p_k_k_0_0[3] = mxPtr[0];
    p_k_k_0_0[4] = mxPtr[1];
    p_k_k_0_0[5] = mxPtr[2];
    mxPtr = p_k_k_mat.ptr<double>(2);
    p_k_k_0_0[6] = mxPtr[0];
    p_k_k_0_0[7] = mxPtr[1];
    p_k_k_0_0[8] = mxPtr[2];
    
    //0_1: p_k_k(1:3,4:7)*Jnorm'  3x4 matrix multiplied across jnorm'
    double* p_k_k_0_1 = (double*)malloc(12*sizeof(double));
    mxPtr = p_k_k_mat.ptr<double>(0) + 3;
    p_k_k_0_1[0] = mxPtr[0];
    p_k_k_0_1[1] = mxPtr[1];
    p_k_k_0_1[2] = mxPtr[2];
    p_k_k_0_1[3] = mxPtr[3];
    mxPtr = p_k_k_mat.ptr<double>(1) + 3;
    p_k_k_0_1[4] = mxPtr[0];
    p_k_k_0_1[5] = mxPtr[1];
    p_k_k_0_1[6] = mxPtr[2];
    p_k_k_0_1[7] = mxPtr[3];
    mxPtr = p_k_k_mat.ptr<double>(2) + 3;
    p_k_k_0_1[8] = mxPtr[0];
    p_k_k_0_1[9] = mxPtr[1];
    p_k_k_0_1[10] = mxPtr[2];
    p_k_k_0_1[11] = mxPtr[3];
    cv::Mat p_k_k_0_1_mat(3, 4, CV_64FC1, p_k_k_0_1, cv::Mat::AUTO_STEP);
    p_k_k_0_1_mat = p_k_k_0_1_mat * JNorm_mat_transpose;
    
    //0_2: p_k_k(1:3, 8:size_p_k_k)   3x(state_size-7) matrix
    int p_k_featureSubMatrix_numCols = (state_size-7);
    double* p_k_k_0_2 = (double*)malloc(3*p_k_featureSubMatrix_numCols*sizeof(double));
    mxPtr = p_k_k_mat.ptr<double>(0) + 7;
    for(int i=0; i<p_k_featureSubMatrix_numCols; ++i){
        p_k_k_0_2[i] = mxPtr[i];
    }
    mxPtr = p_k_k_mat.ptr<double>(1) + 7;
    for (int i=0; i<p_k_featureSubMatrix_numCols; ++i){
        p_k_k_0_2[p_k_featureSubMatrix_numCols + i] = mxPtr[i];
    }
    mxPtr = p_k_k_mat.ptr<double>(2) + 7;
    for (int i=0; i<p_k_featureSubMatrix_numCols; ++i){
        p_k_k_0_2[p_k_featureSubMatrix_numCols*2 + i] = mxPtr[i];
    }
    
    //1_0: Jnorm*p_k_k(4:7, 1:3)  4x3 mx
    double* p_k_k_1_0 = (double*)malloc(12*sizeof(double));
    mxPtr = p_k_k_mat.ptr<double>(3);
    p_k_k_1_0[0] = mxPtr[0];
    p_k_k_1_0[1] = mxPtr[1];
    p_k_k_1_0[2] = mxPtr[2];
    mxPtr = p_k_k_mat.ptr<double>(4);
    p_k_k_1_0[3] = mxPtr[0];
    p_k_k_1_0[4] = mxPtr[1];
    p_k_k_1_0[5] = mxPtr[2];
    mxPtr = p_k_k_mat.ptr<double>(5);
    p_k_k_1_0[6] = mxPtr[0];
    p_k_k_1_0[7] = mxPtr[1];
    p_k_k_1_0[8] = mxPtr[2];
    mxPtr = p_k_k_mat.ptr<double>(6);
    p_k_k_1_0[9] = mxPtr[0];
    p_k_k_1_0[10] = mxPtr[1];
    p_k_k_1_0[11] = mxPtr[2];
    cv::Mat p_k_k_1_0_mat(4, 3, CV_64FC1, p_k_k_1_0, cv::Mat::AUTO_STEP);
    p_k_k_1_0_mat = Jnorm_mat * p_k_k_1_0_mat;
    
    //1_1: Jnorm*p_k_k(4:7, 4:7)*Jnorm', 4x4 mx
    double* p_k_k_1_1 = (double*)malloc(16*sizeof(double));
    mxPtr = p_k_k_mat.ptr<double>(3) + 3;
    p_k_k_1_1[0] = mxPtr[0];
    p_k_k_1_1[1] = mxPtr[1];
    p_k_k_1_1[2] = mxPtr[2];
    p_k_k_1_1[3] = mxPtr[3];
    mxPtr = p_k_k_mat.ptr<double>(4) + 3;
    p_k_k_1_1[4] = mxPtr[0];
    p_k_k_1_1[5] = mxPtr[1];
    p_k_k_1_1[6] = mxPtr[2];
    p_k_k_1_1[7] = mxPtr[3];
    mxPtr = p_k_k_mat.ptr<double>(5) + 3;
    p_k_k_1_1[8] = mxPtr[0];
    p_k_k_1_1[9] = mxPtr[1];
    p_k_k_1_1[10] = mxPtr[2];
    p_k_k_1_1[11] = mxPtr[3];
    mxPtr = p_k_k_mat.ptr<double>(6) + 3;
    p_k_k_1_1[12] = mxPtr[0];
    p_k_k_1_1[13] = mxPtr[1];
    p_k_k_1_1[14] = mxPtr[2];
    p_k_k_1_1[15] = mxPtr[3];
    cv::Mat p_k_k_1_1_mat(4, 4, CV_64FC1, p_k_k_1_1, cv::Mat::AUTO_STEP);
    p_k_k_1_1_mat = Jnorm_mat * p_k_k_1_1_mat;
    p_k_k_1_1_mat = p_k_k_1_1_mat * JNorm_mat_transpose;
    
    //1_2: Jnorm * p_k_k_(4:7, 8:state_size),   4x(p_k_featureSubMatrix_numCols) mx
    double* p_k_k_1_2 = (double*)malloc(4*p_k_featureSubMatrix_numCols*sizeof(double));
    mxPtr = p_k_k_mat.ptr<double>(3) + 7;
    for (int i=0; i<p_k_featureSubMatrix_numCols; ++i){
        p_k_k_1_2[i] = mxPtr[i];
    }
    mxPtr = p_k_k_mat.ptr<double>(4) + 7;
    for (int i=0; i<p_k_featureSubMatrix_numCols; ++i){
        p_k_k_1_2[p_k_featureSubMatrix_numCols + i] = mxPtr[i];
    }
    mxPtr = p_k_k_mat.ptr<double>(5) + 7;
    for (int i=0; i<p_k_featureSubMatrix_numCols; ++i){
        p_k_k_1_2[p_k_featureSubMatrix_numCols*2 + i] = mxPtr[i];
    }
    mxPtr = p_k_k_mat.ptr<double>(6) + 7;
    for (int i=0; i<p_k_featureSubMatrix_numCols; ++i){
        p_k_k_1_2[p_k_featureSubMatrix_numCols*3 + i] = mxPtr[i];
    }
    cv::Mat p_k_K_1_2_mat(4, p_k_featureSubMatrix_numCols, CV_64FC1, p_k_k_1_2, cv::Mat::AUTO_STEP);
    p_k_k_1_1_mat = Jnorm_mat * p_k_K_1_2_mat;
    
    //2_0: p_k_k(8:state_size, 1:3),   (p_k_featureSubMatrix_numCols)x3 mx
    double* p_k_k_2_0 = (double*)malloc(3*p_k_featureSubMatrix_numCols*sizeof(double));
    for (r = 7; r<state_size; ++r){
        mxPtr = p_k_k_mat.ptr<double>(r);
        int currRow = (r-7)*3;
        p_k_k_2_0[currRow] = mxPtr[0];
        p_k_k_2_0[currRow+1] = mxPtr[1];
        p_k_k_2_0[currRow+2] = mxPtr[2];
    }
    
    //2_1: p_k_k(8:state_size, 4:7) * Jnorm',  (p_k_K_featureSubMatrix_numCols)x4 mx
    double* p_k_k_2_1 = (double*)malloc(p_k_featureSubMatrix_numCols*4*sizeof(double));
    for (r=7; r<state_size; ++r){
        mxPtr = p_k_k_mat.ptr<double>(r) + 3;
        int currRow = (r-7)*4;
        p_k_k_2_1[currRow] = mxPtr[0];
        p_k_k_2_1[currRow+1] = mxPtr[1];
        p_k_k_2_1[currRow+2] = mxPtr[2];
        p_k_k_2_1[currRow+3] = mxPtr[3];
    }
    cv::Mat p_k_k_2_1_mat(p_k_featureSubMatrix_numCols, 4, CV_64FC1, p_k_k_2_1, cv::Mat::AUTO_STEP);
    p_k_k_2_1_mat = p_k_k_2_1_mat * JNorm_mat_transpose;
    
    //2_2: p_k_k(8:state_size, 8:state_size)
    double* p_k_k_2_2 = (double*)malloc(p_k_featureSubMatrix_numCols*p_k_featureSubMatrix_numCols*sizeof(double));
    for (r=7; r<state_size; ++r){
        mxPtr = p_k_k_mat.ptr<double>(r) + 7;
        int currRow = (r-7)*p_k_featureSubMatrix_numCols;
        for (c=0; c<p_k_featureSubMatrix_numCols; ++c){
            p_k_k_2_2[currRow + c] = mxPtr[c];
        }
    }
    
    p_k_k_mat.release();
    
    //Construct p_k_k using the 9 sub-matrices:
    int row_offset = 0;
    for (r=0; r<3; ++r){
        row_offset = r*state_size;
        //0_0
        for(c=0; c<3; ++c){
            p_km1_k[row_offset + c] = p_k_k_0_0[r*3+c];
        }
        //0_1
        mxPtr = p_k_k_0_1_mat.ptr<double>(r);
        for (c=3; c<7; ++c){
            p_km1_k[row_offset + c] = mxPtr[c - 3];
        }
        //0_2
        for (c=7; c<state_size; ++c){
            p_km1_k[row_offset + c] = p_k_k_0_2[r*p_k_featureSubMatrix_numCols + c - 7];
        }
    }
    for (r=3; r<7; ++r){
        row_offset = r*state_size;
        //1_0
        mxPtr = p_k_k_1_0_mat.ptr<double>(r-3);
        for (c=0; c<3; ++c){
            p_km1_k[row_offset + c] = mxPtr[c];
        }
        //1_1
        mxPtr = p_k_k_1_1_mat.ptr<double>(r-3);
        for (c=3; c<7; ++c){
            p_km1_k[row_offset + c] = mxPtr[c-3];
        }
        //1_2
        mxPtr = p_k_K_1_2_mat.ptr<double>(r-3);
        for (c=7; c<state_size; ++c){
            p_km1_k[row_offset + c] = mxPtr[c-7];
        }
    }
    for (r=7; r<state_size; ++r){
        row_offset = r*state_size;
        //2_0
        for (c=0; c<3; ++c){
            p_km1_k[row_offset + c] = p_k_k_2_0[(r-7)*3 + c];
        }
        //2_1
        mxPtr = p_k_k_2_1_mat.ptr<double>(r-7);
        for (c=3; c<7; ++c){
            p_km1_k[row_offset + c] = mxPtr[c-3];
        }
        //2_2
        for (c=7; c<state_size; ++c){
            p_km1_k[row_offset + c] = p_k_k_2_2[(r-7)*p_k_featureSubMatrix_numCols + c-7];
        }
    }
    
    free(p_k_k_0_0);free(p_k_k_0_1);free(p_k_k_0_2);
    free(p_k_k_1_0);free(p_k_k_1_1);free(p_k_k_1_2);
    free(p_k_k_2_0);free(p_k_k_2_1);free(p_k_k_2_2);
    p_k_k_0_1_mat.release();
    p_k_k_1_0_mat.release();
    p_k_k_1_1_mat.release();
    p_k_K_1_2_mat.release();
    p_k_k_2_1_mat.release();
    
    //update(normalize) x_km1_k:      x_k_k = x_k_k(4:7)/norm(x_k_k(4:7))
    double norm = (  (*x_k_k_mat.ptr<double>(3))*(*x_k_k_mat.ptr<double>(3)) + (*x_k_k_mat.ptr<double>(4))*(*x_k_k_mat.ptr<double>(4))
                   + (*x_k_k_mat.ptr<double>(5))*(*x_k_k_mat.ptr<double>(5)) + (*x_k_k_mat.ptr<double>(6))*(*x_k_k_mat.ptr<double>(6)) );
    norm = sqrt(norm);
    
    mxPtr = x_k_k_mat.ptr<double>(0);
    for (int i=0; i<state_size; ++i){
        x_km1_k[i] = mxPtr[i];
    }
    x_km1_k[3] = x_km1_k[3]/norm;
    x_km1_k[4] = x_km1_k[4]/norm;
    x_km1_k[5] = x_km1_k[5]/norm;
    x_km1_k[6] = x_km1_k[6]/norm;
    
    return true;
}