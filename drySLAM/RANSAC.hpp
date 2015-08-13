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
















//Identify features that are high innovation inliers
//Precondition: High innovation inliers not classified
//Postcondition: High-innovation inlier features identified
bool rescue_hi_inliers(Filter* filter, NSMutableArray* features_info, Camera* cam);







//Update state and covariance based on high innovation inliers list
//Precondition: High innovation inliers are identified in features_info
//Postcondition: x_k and p_k are updated though "update" method
bool ekf_update_hi_inliers(Filter* filter, NSMutableArray* features_info);


//Warp patches according to predicted motion (Note: patch_when_matching is an opencv Mat matrix)
//Precondition: for each featureData object, patch_when_matching is outdated
//Postcondition: patch_when_matching is updated for all features
void predictFeaturesAppearance(NSMutableArray* features_info,Filter* filter,Camera* cam);


















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
 */
// bool q2r(double *q, double *R);
/*
 double hu(double *yi,double *cam);
 double distortFM(double *uv,double *cam);
 double distortOnePoint(double *uvu,double *cam);
 void fastCornerDetect9(double *im, double threshold,double *coords);
 void fastNonMax(double *im,double barrier,double c,double ret);
 void addFeaturesInveerseDepth(double *uvd,double *X,double *P,double *cam,double *std_px1,double *initialRho,double *stdRho,double *XRES,double *PRES,double *newFeature);
 void hinv(double *uvd,double *Xv,double *cam,double *initialRho,double *newFeature);
 void addOneFeatureCovarianceInverseDepth(double *P,double *uvd,double *Xv,double *std_px1,double *std_rho,double *cam,double *PRES);
 */
// bool undistortFM(double *uvd, int uvd_nCols, Camera *cam, double *uvu);
/*
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
