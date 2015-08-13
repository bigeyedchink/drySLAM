//
//  MapManagement.m
//  drySLAM
//
//  Created by Anthony Rodriguez on 4/17/15.
//  Copyright (c) 2015 ECE473 573 Common. All rights reserved.
//

#import <Foundation/Foundation.h>
#include "RescueHiInliers.hpp"


bool rescue_hi_inliers(Filter* filter, NSMutableArray* features_info, Camera* cam){
    bool success1 = predict_camera_measurements(features_info, filter, cam);
    calculateDerivatives(features_info, filter, cam);
    for (int i=0; i<[features_info count]; ++i){
        Feature* feat = (Feature*)features_info[i];
        if (feat->individually_compatible==1 && feat->low_innovation_inlier==0){
            double* nui = (double*)malloc(2*sizeof(double));
            nui[0] = feat->z[0] - feat->h[0];
            nui[1] = feat->z[1] - feat->h[1];
            cv::Mat H = feat->H->clone();
            cv::Mat H_transpose;
            cv::transpose(H, H_transpose);
            cv::Mat p_k_k_mat(filter->state_size, filter->state_size, CV_64FC1, filter->p_k_k, cv::Mat::AUTO_STEP);
            cv::Mat Si;
            Si = H*p_k_k_mat*H_transpose;
            cv::Mat Si_inverse;
            cv::invert(Si, Si_inverse);
            cv::Mat nui_mat(2, 1, CV_64FC1, nui, cv::Mat::AUTO_STEP);
            cv::Mat nui_transpose;
            cv::invert(nui_mat, nui_transpose);
            cv::Mat compareVal = nui_transpose*Si_inverse*nui_mat;
            if (*(compareVal.ptr<double>(0)) < CHI_2_INV_2_95)
                ((Feature*)features_info[i])->high_innovation_inlier=1;
            else
                ((Feature*)features_info[i])->high_innovation_inlier=0;
            free(nui);
            H.release();
            H_transpose.release();
            p_k_k_mat.release();
            Si.release();
            Si_inverse.release();
            nui_mat.release();
            nui_transpose.release();
            compareVal.release();
        }
    }
    return true;
}