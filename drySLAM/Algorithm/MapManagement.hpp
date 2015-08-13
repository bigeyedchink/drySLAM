//
//  MapManagement.hpp
//  drySLAM
//
//  Created by Anthony Rodriguez on 4/17/15.
//  Copyright (c) 2015 ECE473 573 Common. All rights reserved.
//

#ifndef MAPMANAGEMENT_HPP
#define MAPMANAGEMENT_HPP

#include "SLAMModel.hpp"


//Call all appropriate functions in mapManagement
bool map_management(cv::Mat* im, Filter* filter, Camera* cam, NSMutableArray* features_info, int minNumberOfFeaturesInImage, int step);


//Finds all features that need to be deleted and deletes them
bool deleteFeatures(NSMutableArray* features_info, Filter* filter);


//Deletes the feature at index featToDelete from the state and covariance matrices and returns new matrices
//PreCondition: x_k_k, p_k_k & sizes are all old values
//PostCondition: x_k_k, p_k_k updated, x_k_k and p_k_k reflect new sizes
bool deleteOneFeature(NSMutableArray* features_info,Filter* filter,int featToDelete);


//Resets feature info parameters in each iteration
//Precondition:
//Postcondition: features_info is updated for all features
bool updateFeaturesInfo(NSMutableArray* features_info);


//NOTE: This function is currently unnecessary for our implementation (just return true)
//Find features using inverse depth coordinates & change them to cartesian
//TODO: We need to figure out how to test this because matlab code only changes one feature
//Note: We don't need to use cartesian at all for the current release
bool inversedepth_2_cartesian(NSMutableArray* features_info, Filter* filter);


//Initialize features (basically just calls initialize_a_feature in a loop)
//Dependencies: initialize_a_feature
//Precondition: Nothing special
//Postcondition: Filter, features_info are updated for new features
bool initialize_features(int step, Camera* cam, Filter* filter, NSMutableArray* features_info, int numFeaturesToInitialize, cv::Mat* im);


//Use FAST corner detection to find coordinates for a new feature to add to filter
//FAST opencv method: void FAST(InputArray image, vector<KeyPoint>& keypoints, int threshold, bool nonmaxSuppression=true )
//Dependencies: add_feature_to_info_vector, add_features_inverse_depth, hinv, add_a_feature_covariance_inverse_depth
//NOTE: If no elements are placed in uvd vector, DO NOT call add_features_inverse_depth as in matlab code
//Precondition: uv (coordinates of the new feature) is allocated as a 2-element array
//Postcondition: uv has coordinates of new feature, filter and features_info are modified accordingly
bool initialize_a_feature(Filter* filter, NSMutableArray* features_info, int step, Camera* cam, cv::Mat* im_k, NSMutableArray* uv);


//Add the feature at coordinate uvd to X_RES, P_RES matrices and instantiate newFeature with 6-variable id vector
//Dependencies: addOneFeatureCovarianceInverseDepth
//Precondition: uvd is a 2-element matrix with frame coordinates of feature,
//              filter holds current state/covariance sizes.
//Postcondition: new feature is added to filter
bool add_features_inverse_depth(double* uvd, int nNewFeat, Filter *filter, Camera* cam, double std_pxl, double initial_rho, double std_rho, double* newFeature);


//Initialize the 6-element inverse depth feature vector for a new feature located at coordinates in uvu
//Dependencies: undistort_fm
//Precondition: uvd is a 2-element coordinate vector, newFeature is allocated as a 6-element array
//Postcondition: newFeature initialized appropriately
bool hinv(double* uvd, double* Xv, Camera* cam, double initialRho, double* newFeature);


//Add a single feature at coordinates uvd to original covariance matrix P to produce P_RES
//Dependencies: undistort_fm, dRq_times_a_by_dq, jacob_undistor_fm
//Precondition: P_RES is allocated as the length of P + 6 in each direction (for id coordinate 6-element vector)
//              uvd is a 2-element array containing coordinates of feature in image
//              Xv is 13-element original state vector
//Postcondition: P_RES initialized as appropriate
bool addOneFeatureCovarianceInverseDepth(Filter *filter, double* uvd, double *Xv, double std_pxl, double std_rho, Camera* cam, double *P_RES);


//Update h values across features
//Dependencies: hi_inverse_depth, hi_cartesian (hi cartesian not applicable to this release)
//Precondition: each feature h value is old
//Postcondition: h values for all features are filled out
bool predict_camera_measurements(NSMutableArray *features_info,Filter *filter,Camera *cam);


/*Find the measurement for a given inverse depth feature (stored in zi, 2x1 vector of coordinates in image for a feature - I think)
 *Note: features_info not used in matlab version
 *Precondition:  yinit: 6x1 vector (array)
 *               t_wc: 3x1 vector used for linear algebra
 *               r_wc: 3x3 matrix used for linear algebra
 *               r_cw: 3x3 matrix used for linear algebra
 *               zi: 2-element array
 *Postcondition: zi is either given two elements for x,y coordinates, or is returned empty (count=0) if it cannot be found
 */
bool hiInverseDepth(Camera* cam, double* yinit, double* t_wc, double* r_wc, double* zi);


/*Find the measurement for a given cartesian feature (stored in zi, 2x1 vector of coordinates in image for a feature - I think)
 *Precondition:  yinit: 6x1 vector
 *               t_wc: 3x1 vector
 *               r_wc: 3x3 matrix
 *               r_cw: 3x3 matrix
 *              zi: Allocated array with size 0
 *Postcondition: zi is either given two elements for x,y coordinates, or is returned empty (count=0) if it cannot be found
 */
//TODO: Get rid of cv Matrices
//int hiCartesian(NSMutableArray* features_info,Camera* cam,cv::Mat yi,cv::Mat t_wc,cv::Mat r_wc,cv::Mat zi);

//TODO: Function "m"


//Not sure exactly what this accomplishes yet
//Precondition: yi is a vector of length 3, uv_u is an array of length 2
//              Note: I don't know if uv_u should be variable size...Most times it seems to be 2x1, but use uv_width just in case...
//              number of columns in uv_u will be the same as yi...Usually 1 but maybe more sometimes?
//Postcondition: uv_u is filled out
bool hu(double *uv_u,double *yi,int yi_numColums,Camera *cam);


//Distort points found in function "hu"
//Precondition: uv_u is a 2x(uv_u_numColumns) matrix, uvd is a vector of length 2
//Postcondition: Appropriate coordinate values are set for uvd
bool distortFM(double *uvd, double *uv, int uv_ncols,Camera* cam);


//TODO:add_features_inverse_depth()


//TODO: hinv()

//Undistort coordinates at uvd
//Precondition: uvu, uvd are 2-element coordinate arrays,uvd is input, uvu output
//Postcondition: uvu holds undistorted coordinates
 bool undistortFM(double *uvd, int uvd_nCols, Camera *cam, double *uvu);



//Initialize all of the data in a new feature
//Precondition: uv: 2-element coordinates
//              im_k is a grayscale image (widthxheight) matrix of unsigned chars (8-bit values) - opencv type CV_8UC1
//              X_RES: 91-element array? TODO: Find out if this changes
//Postcondition: new feature added to features_info
//              Note: new feature doesn't reference im_k - im_k gets clone()d into the new feature
bool add_feature_to_info_vector(int* uv, cv::Mat* im_k, double* X_RES, int X_RES_size, NSMutableArray* features_info, int step, double *newFeature);


//Get matrix dgx_dqwr for method addOneFeatureCovarianceInverseDepth
//Precondition: q is a 4-element array, aMat is a 3x1 matrix, dRqTimesABydqRES is a 3x4 matrix of zeros (already allocated/initialized)
//Postcondition: dRqTimesABydqRES is filled with correct values
//Should return false for null inputs of if the matrix data is invalid 
bool dRqTimesABydq(double *q,cv::Mat* aMat,cv::Mat* dRqTimesABydqRES);


//TODO: Heading
void JacobUndistorFM(Camera* cam,double* uvd,cv::Mat* J_unidistor);


//Get R_wc for calculating inverse depth covariance
//Precondition: q is a 4-element array, R is a 9-element (3x3) array
//Postcondition: R is filled with approrpriate values
bool q2r(double *q, double *R);


//Calculate the Jacobian of q
//Precondition: q is a 4-element array, J is allocated as a 4x4 (16-element) array
//Postcondition: J is the correct jacobian of q
bool normJac(double* J, double* q);

//Basically just initialize a new Feature
bool addFeatureToInfoVector(double *uv, double im_k, double *X_RES, NSMutableArray *features_info, int step, double newFeature);

#endif

