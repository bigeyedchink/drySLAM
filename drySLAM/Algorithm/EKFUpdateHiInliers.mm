//
//  MapManagement.m
//  drySLAM
//
//  Created by Anthony Rodriguez on 4/17/15.
//  Copyright (c) 2015 ECE473 573 Common. All rights reserved.
//

#import <Foundation/Foundation.h>
#include "EKFUpdateLiInliers.hpp"

bool ekf_update_hi_inliers(Filter* filter, NSMutableArray* features_info)
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
        if ( feat->high_innovation_inlier==1 ){
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
    /*for (int i=0; i<filter->state_size; ++i){
        filter->x_k_k[i] = filter->x_k_km1[i];
    }
    for (int i=0; i<(filter->state_size*filter->state_size); ++i){
        filter->p_k_k[i] = filter->p_k_km1[i];
    }*/
    
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
