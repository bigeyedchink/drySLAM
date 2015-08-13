//
//  RANSAC.h
//  RANSAC
//
//  Created by guest on 3/27/15.
//  Copyright (c) 2015 dryslam. All rights reserved.
//

#ifndef RANSAC_H
#define RANSAC_H

#ifdef __cplusplus
#import <opencv2/opencv.hpp>
#import <opencv2/highgui/highgui_c.h>
#import <opencv2/highgui/cap_ios.h>
#endif
#import <Foundation/Foundation.h>
#include "SLAMModel.hpp"
#include "EKFilter.hpp"

//Return values: true: success / false: error(quit program)

//Finds all features that need to be deleted and deletes them
bool deleteFeatures(NSMutableArray* features_info, Filter* filter);


//Deletes the feature at index featToDelete from the state and covariance matrices and returns new matrices
//PreCondition: x_k_k, p_k_k & sizes are all old values
//PostCondition: x_k_k, p_k_k updated, x_k_k and p_k_k reflect new sizes
bool deleteOneFeature(NSMutableArray* features_info, Filter* filter,
                      int featToDelete, double* x_k_k, double* p_k_k,int parToDelete);


//Resets feature info parameters in each iteration
//Precondition:
//Postcondition: features_info is updated for all features
bool updateFeaturesInfo(NSMutableArray* features_info);


//Find features using inverse depth coordinates & change them to cartesian
//TODO: We need to figure out how to test this because matlab code only changes one feature
bool inversedepth_2_cartesian(NSMutableArray* features_info,
                              double *x_k_k, double* p_k_k, int state_size);

bool hu(double *uv_u,double *yi,int yi_numColums,Camera *cam);


//Generate state vector pattern
//Precondition: z_id, z_euc and state_vector_pattern are allocated but not initialized
//              z_id/z_euc end up being 2xn matrices
//              state_vector_pattern: (state_size x 4) matrix
//Postcondition: z_id, z_euc and state_vector_pattern are updated and set to correct values
//NOTE: z_id and z_euc are lists of x,y values stored as follows:   x1 x2 x3 x4
//                                                                  y1 y2 y3 y4
//       So, In order to store this we must create a list
//       of x values and a list of y values,
// then z is a 1d array stored in a way that looks like:            x1 x2 x3 x4 y1 y2 y3 y4
//***z_id and z_euc are stored as ints
bool generate_state_vector_pattern(NSMutableArray* features_info, int state_size, double* state_vector_pattern,
                                   NSMutableArray* z_id, NSMutableArray* z_euc);


//Select potential match for high innovation inlier (I think)
//Precondition: zi, position and num_IC_matches are allocated.
//Postcondition: zi, position and num_IC_matches are allocated and set correctly
bool select_random_match(NSMutableArray* features_info,
                         double* zi, int* position, int* num_IC_matches);


//Predict measurements and count matches under a threshold
//Precondition: hypothesis_support, positions_li_inliers_sizes not initialized, positions_li_inliers is NULL
//              Positions_li_inliers arrays are boolean arrays (1s and 0s)
//Postcondition: hypothesis_support, positions_li_inliers_id, positions_li_inliers_id_size set correctly
//              If there are no id or euc inliers, keep that array set to NULL and set its size to 0
bool compute_hypothesis_support(double* xi, int state_size, Camera* cam, double* state_vector_pattern,
                                double* z_id, double* z_euc, int threshold,
                                int* hypothesis_support,
                                NSMutableArray* positions_li_inliers_id,
                                NSMutableArray* positions_li_inliers_euc);


//Set appropriate features to be inliers
//Precondition: All variables are allocated and initialized
//Postcondition: Appropriate features_info elements have correct "low_innovation_inlier" values
bool set_as_most_supported_hypothesis(NSMutableArray* features_info,
                                      NSMutableArray* positions_li_inliers_id,
                                      NSMutableArray* positions_li_inliers_euc);


//Update state and covariance based on low innovation inliers list
//Precondition: Low innovation inliers are identified in features_info
//Postcondition: x_k and p_k are updated though "update" method
bool ekf_update_li_inliers(Filter* filter, NSMutableArray* features_info);


//Linear algebra calculations for matrices assembled in ekf_update_li_inliers
//Precondition: all matrices are built and ready to use (all have variable sizes, so check)
//Postcondition: x_km1_k and p_km1_k reflect new state values
bool update(double* x_km1_k, double* p_km1_k, int state_size, cv::Mat* H, cv::Mat* R, cv::Mat* z, cv::Mat* h);


//Identify features that are high innovation inliers
//Precondition: High innovation inliers not classified
//Postcondition: High-innovation inlier features identified
bool rescue_hi_inliers(Filter* filter, NSMutableArray* features_info, Camera* cam);


//Calculate jacobian H for feature i in features_info
//Dependency: dh_dxv
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


//Auxiliary calculations for calculating hi (jacobian with respect to 13-var state vector?)
//Precondition: Hi_features 2x13, Xv_km1_k 1x13, yi 6x1, zi 1x2
//Postcondition: Hi_features filled out appropriately
bool dh_dxv(double* Hi_features, Camera* cam, double* Xv_km1_k, double* yi, double* zi);


//Update state and covariance based on high innovation inliers list
//Precondition: High innovation inliers are identified in features_info
//Postcondition: x_k and p_k are updated though "update" method
bool ekf_update_hi_inliers(Filter* filter, NSMutableArray* features_info);






















//@interface RANSAC : NSObject

////////////////////////////////////////////////////////////////////////////////////////////////////
/////I.mapManagement
//I.1)deleteFeatures
//Not needed: void deleteFeatures(struct FeaturesInfo featureInfo,struct Filter filter,int stepnumber);



//I.2)updateFeaturesInfo:
//Resets feature info parameters in each iteration
//void updateFeaturesInfo(struct featureData* features_info, int features_info_size);

//I.3)inverseDepth2Cartesian0(corrsponding to inversedepth_2_caretesian)


//inversedepth_2_cartesian:
//Find features using inverse depth coordinates & change them to cartesian
//void inversedepth_2_cartesian(struct featureData* features_info, double *x_k_k, double* p_k_k);


/*
 double m(double *a,double *b);
 void inverseDepth2Cartesian(double *inverseDepth,double *Cartesian);
 
 //I.4)initializeFeatures
 
 double get_std_z(double *filter);
 void q2r(double *q, double *R);
 double hu(double *yi,double *cam);
 double distortFM(double *uv,double *cam);
 double distortOnePoint(double *uvu,double *cam);
 void fastCornerDetect9(double *im, double threshold,double *coords);
 void fastNonMax(double *im,double barrier,double c,double ret);
 void addFeaturesInveerseDepth(double *uvd,double *X,double *P,double *cam,double *std_px1,double *initialRho,double *stdRho,double *XRES,double *PRES,double *newFeature);
 void hinv(double *uvd,double *Xv,double *cam,double *initialRho,double *newFeature);
 void addOneFeatureCovarianceInverseDepth(double *P,double *uvd,double *Xv,double *std_px1,double *std_rho,double *cam,double *PRES);
 void undistortFM(double *uv,double *cam);
 void undistortOnePoint(double *uvu,double *cam);
 double dRqTimesaBydq(int q,double aMat);
 void JacobUndistorFM(double *uvd,double *cam,double *J_unidistor);
 void addFeatureToInfoVector(struct FeaturesInfo *featureInfo,double *uv,double *im_k,double *X_RES,int step,double *newFeature);
 ////////////////////////////////////////////////////////////////////////////////////////////////////
 /////II.searchINMatches
 
 //II.2)calculateDerivatives
 void calculate_derivatives(struct FeaturesInfo *featureInfo,double x_k_km1,double *cam);
 //II.3)predictFeaturesAppearance
 void pred_pathch_fc(struct FeaturesInfo *featureInfo,double *R_Wk,double *r_Wk,double *XYZ_w,double *cam,double *patchPred);
 void rotateWithDistFCC2C1(double *cam, double *uv_c2,double *R_c1c2,double *t_c1c2,double *n,int d,double *uv_c1);
 void rotateWithDistFCC1C2(double *cam, double *uv_c1,double *R_c2c1,double *t_c2c1,double *n,int d,double *uv_c2);
 //II.4)matching
 void matching(struct FeaturesInfo *featureInfo,double *im,double *cam);
 ////////////////////////////////////////////////////////////////////////////////////////////////////
 /////III.RANSACHypotheses
 //III.1)generateStateVectorPattern
 void generateStateVectorPattern(struct FeaturesInfo featureInfo,struct Filter filter,int stepnumber);
 //III.2)selectRandomMatch
 double selectRandomMatch(struct FeaturesInfo *featureInfo,double *position);
 //III.3)get
 double get_x_k_km1(double *filter);
 double get_p_k_km1(double *filter);
 //III.4)computeHypothesisSupportFast
 int computeHypothesisSupportFast(double *xi,double *cam,double *stateVectorPattern,double *z_id,double *z_euc,double threshold,double *positionsLIInlinersID,double *positionsHIInlinersID);
 //III.5)setAsMostSupportedHypothesis
 void setAsMostSupportedHypothesis(struct FeaturesInfo *featureInfo,double *positionsLIInlinersID,double *positionsHIInlinersID);
 ////////////////////////////////////////////////////////////////////////////////////////////////////
 /////MainLoop
 //not needed: void mapManagement(struct FeaturesInfo *featureInfo,double *filter,double *cam,double *im,int minNUmberOfFeaturesInImage,int step);     //I
 void searchICMatches(struct FeaturesInfo *featureInfo,double *filter,double *cam,double *im);     //II
 void RANSACHypothesis(struct FeaturesInfo featureInfo,struct Filter filter,double *cam);    //III
 void EKFUpdateLIInlier();
 void RescueHIInliers();
 void EKFUpdateHIIlier();
 */
//@end

#endif
