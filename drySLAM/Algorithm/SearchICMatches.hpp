//
//  SearchICMatches.hpp
//  drySLAM
//
//  Created by Anthony Rodriguez on 4/17/15.
//  Copyright (c) 2015 ECE473 573 Common. All rights reserved.
//

#ifndef SEARCHICMATCHES_HPP
#define SEARCHICMATCHES_HPP

#include "SLAMModel.hpp"

//main gateway:
bool search_IC_matches(Filter* filter, NSMutableArray* features_info, Camera* cam, cv::Mat* im);


//predict_camera_measurements() is located in MapManagement.hpp

//Calculate jacobian H (derivatives of h matrix - note: this is an opencv Mat matrix) for each feature
//Dependencies: calculate_hi_inverse_depth, calculate_hi_cartesian
//Precondition: H has old values for all features
//PostCondition: H updated for all features
void calculateDerivatives(NSMutableArray* features_info, Filter* filter,Camera* cam);


//Calculate jacobian H for feature i in features_info
//Dependency: dh_dxv, dh_dy
//Precondition: Hi allocated as a matrix, 2rows, 13+(numberOfCartesianFeatures)*3+(numberOfInverseDepthFeatures)*6 rows
//              cartesian_features_index and inverse_depth_features_index are [features_info count]x1 arrays (binary: 1 if feature is that type)
//              xv_km1_k is a 13x1 state vector
//              yi is a 6x1 vector (first 6 states I think)
//Postcondition: Hi is filled with appropriate derivative values
bool calculate_hi_inverse_depth(double* Hi, double* inverse_depth_features_index, double* cartesian_features_index,
                                double* Xv_km1_k, double* yi, Camera* cam,
                                int i, NSMutableArray* features_info);


//Calculate jacobian H for feature i in features_info
//NOTE: Pgm currently does not use cartesian coordianates, so this is not necessary (yet)
//Dependency: calculate_derivatives
//Precondition: Hi allocated as a matrix, 2rows, 13+(numberOfCartesianFeatures)*3+(numberOfInverseDepthFeatures)*6 rows
//              cartesian_features_index and inverse_depth_features_index are [features_info count]x1 arrays (binary: 1 if feature is that type)
//              xv_km1_k is a 13x1 state vector
//              yi is a 6x1 vector (first 6 states I think)
//Postcondition: Hi is filled with appropriate derivative values
bool calculate_hi_cartesian(double* Hi, double* inverse_depth_features_index, double* cartesian_features_index,
                            double* Xv_km1_k, double* yi, Camera* cam,
                            int i, NSMutableArray* features_info);


//Auxiliary calculations for calculating hi (jacobian with respect to state vector?)
//Dependencies: q2r, jacob_undistor_fm
//These functions are part of the file calculate_Hi_inverse_depth
//Precondition: Xv_km1_k 1x13, yi 6x1, zi 1x2
//              Hi_features: 2x13
//Postcondition: Hi_features filled out appropriately
bool dh_dxv(double* Hi_features, Camera* cam, double* Xv_km1_k, double* yi, double* zi);
//Same thing as dh_dxh but Hi_features is the y vector (these are assembled in calculate_Hi_inverse_depth
bool dh_dy(double* Hi_features, Camera* cam, double* Xv_km1_k, double* yi, double* zi);
//critical internal function:
bool dh_dhrl(double* a, Camera* cam, double* Xv_km1_k, double* yi, double* zi);
bool dhrl_dqwr(double* a, double* Xv_km1_k, double* yi);


//Quick conversion from inverse depth to cartesian for predicting patches
//Precondition: inverse_depth is a 6-element array of parameters
//              cartesian is a 3-element cartesian coordinate array
//Postcondition: cartesian holds equivalent conversion values for coordinates
bool inverseDepth2Cartesian(double* inverse_depth, double* cartesian);


//For each feature, try to build a matrix that predicts where the feature will be in next frame and what it will look like
//Dependencies: inverseDepth2Cartesian, pred_patch_fc
//Precondition: Self-explanatory
//Postcondition: Features updated appropriately
bool predict_features_appearance(NSMutableArray* features_info, Filter* filter, Camera* cam);


//Predict what the feature patch_p_f_ini's image patch will look like
//Dependencies: rotate_with_dist_fc_c2c1, rotate_with_dist_fc_c2c2, interp2(this is a Matlab function that we will need to approximate)
//Precondition: R_Wk is a 3x3 matrix, r_Wk is a 3-element vector, XYZ_w is a 3-element vector, patch_pred is a CV_8UC1 13x13 image
//Postcondition: patch_pred holds the predicted image patch for the feature
bool pred_patch_fc(Camera* cam, Feature* patch_p_f_ini, double* R_Wk, double* r_Wk, double* XYZ_w, cv::Mat* patch_pred);


//Rotate some values, add some distortion (not sure exactly)
//Dependencies: undistort_fm, distort_fm
//Precondition: uv_c1: 2-element vector
//              R_c2c1: 3x3 matrix (camera rotation matrix)
//              t_c2c1: 3-element vector
//              n: 3-element vector
//              uv_c2 allocated as a 2-element array
//Postcondition: uv_c2 has appropriate update
bool rotate_with_dist_fc_c2c1(Camera* cam, double* uv_c1, double* R_c2c1, double* t_c2c1, double* n, double d, double* uv_c2);


//Rotate some values, add some distortion (not sure exactly)
//Similar to function above...
//Dependencies: undistort_fm, distort_fm
//Precondition: uv_c2: uv_c2_lengthx2 matrix
//              R_c1c2: 3x3 matrix (camera rotation matrix)
//              t_c1c2: 3-element vector
//              n: 3-element vector
//              uv_c1 allocated as a 2-element array
//Postcondition: uv_c1 has appropriate update
bool rotate_with_dist_fc_c1c2(Camera* cam, double* uv_c2, int uv_c2_length, double* R_c1c2, double* t_c1c2, double* n, double d, double* uv_c1);


//Interpolate function patch_p_f (defined at points from 1 to 41) at query points Xq,Yq
//Precondition: patch_p_f is a 41x41 matrix
//              Xq: X sample points
//              Yq: Y sample points
//              patch_pred: 13x13 interpolated result
bool interpPatchPred(uchar* patch_p_f, double* Xq, double* Yq, double* patch_pred);


//Match feature predictions with actual image
//Note: this function requires correlation coefficent calculation:
//openCV implementation:  void matchTemplate(Mat image, Mat templ, Mat result, CV_TM_CCOEFF)
//Dependencies: None
//Precondition: Old info
//Postcondition: individually_compatible and z elements of Features in features_info updated
bool matching(cv::Mat* im, NSMutableArray* features_info, Camera* cam);



#endif








