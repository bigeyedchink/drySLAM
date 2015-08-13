//
//  SLAMModelController.h
//  drySLAM
//
//  Created by Anthony Rodriguez on 3/28/15.
//  Copyright (c) 2015 ECE473 573 Common. All rights reserved.
//

//This is the main loop called from the UI. It is responsible for holding EKFMonoSLAM relevant data structures and calling
//methods from the EKFilter and RANSAC libraries

#ifndef SLAMMODEL_H
#define SLAMMODEL_H

#ifdef __cplusplus
//#nclude <opencv2/opencv.hpp>
//#import <opencv2/highgui/highgui_c.h>
//#import <opencv2/highgui/cap_ios.h>
#endif
#import "Camera.hpp"
#import "Filter.hpp"
#import "Feature.hpp"
#import <Foundation/Foundation.h>

#include "Algorithm/EKFPrediction.hpp"
#include "Algorithm/MapManagement.hpp"
#include "Algorithm/EKFUpdateHiInliers.hpp"
#include "Algorithm/EKFUpdateLiInliers.hpp"
#include "Algorithm/RansacHypothesis.hpp"
#include "Algorithm/RescueHiInliers.hpp"
#include "Algorithm/SearchICMatches.hpp"


//Stdev for linear acceleration
#define SIGMA_A 0.007f
//Stdev for angular acceleration
#define SIGMA_ALPHA 0.007f
//stdev for measurement noise
#define SIGMA_IMAGE_NOISE 1.f
//This is 25 in MATLAB
#define MIN_NUMBER_OF_FEATURES_IN_IMAGE 25
//Time step info
#define DELTA_T 1.0



/**********************************THE MODEL*********************************/
//Note: These only handle allocations/deallocations WITHIN the structure
//#ifdef __cplusplus
/*
bool initializeCameraData(struct cameraData* cam);
bool destroyCameraData(struct cameraData* cam);
bool initializeFeatureData(struct featureData* feature);
bool destroyFeatureData(struct featureData* feature);
bool initializeFilterData(struct filterData* filter);
bool destroyFilterData(struct filterData* filter);
//#endif*/

@interface SLAMModel : NSObject{
    int step;
    //Camera info
    Camera* cam;
    //Keeps most matrices for the filter
    Filter* filter;
    //Array of featureData: One for each feature
    NSMutableArray* features_info;
    
    // array size is 7 x numImages
    double* trajectory;
    //TODO figure out what this does
    //X Y Z theta phi lambda
    double randSphere6D[6][1000];
#ifdef __cplusplus
    cv::Mat currentImage;
#endif
}
/*
bool ekf_prediction(Filter*, NSMutableArray*);

bool ransac_hypothesis(Filter*, NSMutableArray*, Camera*);

bool ekf_update_hi_inliers(Filter*, NSMutableArray*);
bool ekf_update_li_inliers(Filter*, NSMutableArray*);
*/
/*
+(bool) initializeCameraData:(struct cameraData*) cam;
+(bool) destroyCameraData:(struct cameraData*) cam;
+(bool) initializeFeatureData:(struct featureData**) feature;
+(bool) destroyFeatureData:(struct featureData**) feature;
+(bool) initializeFilterData:(struct filterData*) filter;
+(bool) destroyFilterData:(struct filterData*) filter;
*/
 
//Main gateway to run code (processes video at filePath)
-(void) processVideoFile:(NSString*)filePath;
//TODO:
//#ifdef __cplusplus
//Process a single frame and update data within SLAMModelController
//This is public to facilitate possibility of real-time processing in the future
//-(void) processVideoFrame:(cv::Mat)currentFrame;
//#endif
//Print relevant data to file to compare against EKFMonoSLAM data(may not be necessary)
//-(void) exportStateDataToFile;

@end


#endif