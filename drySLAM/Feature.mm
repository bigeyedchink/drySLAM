//
//  Feature.m
//  drySLAM
//
//  Created by Anthony Rodriguez on 4/1/15.
//  Copyright (c) 2015 ECE473 573 Common. All rights reserved.
//

#import <Foundation/Foundation.h>
#ifdef __cplusplus
#import <opencv2/opencv.hpp>
#import <opencv2/highgui/highgui_c.h>
#import <opencv2/highgui/cap_ios.h>
#endif
#import "Feature.hpp"


@implementation Feature : NSObject

-(id)init{
    self = [super init];
    Feature* feature = self;
    feature->patch_when_initialized = new cv::Mat(41, 41, CV_8UC1, 255);
   // *feature->patch_when_initialized = cv::Mat(41, 41, CV_8UC1, 255); //Please check later.
    feature->patch_when_matching = new cv::Mat(13, 13, CV_64FC1, 255);
    //*feature->patch_when_matching = cv::Mat(13, 13, CV_8UC1, 255); //Please check later.
    feature->r_wc_when_initialized = NULL;
    feature->R_wc_when_initialized = NULL;
    feature->uv_when_initialized[0] = 0;
    feature->uv_when_initialized[1] = 0;
    feature->half_patch_size_when_initialized = 20;
    feature->half_patch_size_when_matching = 6;
    feature->times_predicted = 0;
    feature->times_measured = 0;
    feature->individually_compatible = 0;
    feature->low_innovation_inlier = 0;
    feature->high_innovation_inlier = 0;
    feature->init_frame = 0;
    feature->init_measurement[0] = 0;
    feature->init_measurement[1] = 0;
    feature->type = [NSString stringWithFormat:@"inversedepth"];
    feature->z = NULL;
    feature->H = NULL;
    feature->h = NULL;
    feature->S = NULL;
    return self;
};

-(void)destroy{
    if (patch_when_initialized!=NULL && patch_when_initialized->data!=NULL){
        patch_when_initialized->release();
    }
    if (patch_when_matching!=NULL && patch_when_matching!=NULL){
        patch_when_matching->release();
    }
    free(R_wc_when_initialized);
    if (H!=NULL)
        H->release();
    if (h!=NULL)
        free(h);
    if (z!=NULL)
        free(z);
    if (S==NULL)
        free(S);
}



@end




