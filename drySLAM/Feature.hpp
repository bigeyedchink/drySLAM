//
//  Feature.h
//  drySLAM
//
//  Created by Anthony Rodriguez on 4/1/15.
//  Copyright (c) 2015 ECE473 573 Common. All rights reserved.
//

#ifdef __cplusplus
#import <opencv2/opencv.hpp>
#import <opencv2/highgui/highgui_c.h>
#import <opencv2/highgui/cap_ios.h>
#endif

#ifndef drySLAM_Feature_h
#define drySLAM_Feature_h

@interface Feature : NSObject{
#ifdef __cplusplus
    //41x41 matrix of 8UC1(8bit gray) pixels
    @public cv::Mat* patch_when_initialized;
    //13x13 8UC1 matrix
    @public cv::Mat* patch_when_matching;
#endif
    //TODO Not sure what this does yet (size=3?) YES.1*3
    @public double* r_wc_when_initialized;
    //TODO Not sure what this does yet (size=3x3?) YES.3*3
    @public double* R_wc_when_initialized;
    //TODO  Not sure what this does yet (size=2?) I think so.
    @public double uv_when_initialized[2];
    //TODO Not sure what this does yet (almost always =20)
    @public int half_patch_size_when_initialized;
    //TODO Not sure what this does yet (almost always =6)
    @public int half_patch_size_when_matching;
    
    @public int times_predicted;
    @public int times_measured;
    
    //Boolean values (0,1)
    @public int individually_compatible;
    @public int low_innovation_inlier;
    @public int high_innovation_inlier;
    
    @public int init_frame;
    @public int init_measurement[2];
    //type is almost always "inversedepth"
    @public NSString* type;
    //TODO Not sure what these do:
    @public double yi[6];
    
    //Normally z will be of size 2, but sometimes it must be empty.
    //In these cases, z will be NULL
    @public int* z;
    //S is the variance/covariance of gaussian in which we search for z
    @public double* S;
    @public double* h;
#ifdef __cplusplus
    //H is the jacobian of measurement h (calculate_derivatives.m)
    //Almost always a sparse matrix?
    @public cv::Mat* H;
#endif
    
    @public int state_size;
    @public int measurement_size;
    //R: Gaussian noise
    @public int R[2][2];
}


-(id)init;
-(void)destroy;


@end
#endif

