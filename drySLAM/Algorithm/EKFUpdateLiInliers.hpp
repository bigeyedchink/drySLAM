//
//  EKFUpdateLiInliers.hpp
//  drySLAM
//
//  Created by Anthony Rodriguez on 4/17/15.
//  Copyright (c) 2015 ECE473 573 Common. All rights reserved.
//

#ifndef EKFUPDATELIINLIERS_HPP
#define EKFUPDATELIINLIERS_HPP

#include "SLAMModel.hpp"



//Update state and covariance based on low innovation inliers list
//Precondition: Low innovation inliers are identified in features_info
//Postcondition: x_k and p_k are updated though "update" method
bool ekf_update_li_inliers(Filter* filter, NSMutableArray* features_info);


//Linear algebra calculations for matrices assembled in ekf_update_li_inliers
//Precondition: all matrices are built and ready to use (all have variable sizes, so check)
//Postcondition: x_km1_k and p_km1_k reflect new state values
bool update(double* x_km1_k, double* p_km1_k, int state_size, cv::Mat* H, cv::Mat* R, cv::Mat* z, cv::Mat* h);


#endif