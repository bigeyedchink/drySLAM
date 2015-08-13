//
//  Camera.h
//  drySLAM
//
//  Created by Anthony Rodriguez on 4/1/15.
//  Copyright (c) 2015 ECE473 573 Common. All rights reserved.
//

#ifndef drySLAM_Camera_h
#define drySLAM_Camera_h

#import <Foundation/Foundation.h>
//#include "RecordViewController.h"
//#ifdef __cplusplus
#import <opencv2/opencv.hpp>
//#import <opencv2/highgui/highgui_c.h>
//#import <opencv2/highgui/cap_ios.h>
//#endif

@interface Camera : NSObject{
    @public double k1;
    @public double k2;
    @public int nRows;
    @public int nCols;
    @public double Cx;
    @public double Cy;
    @public double f;
    @public double dx;
    @public double dy;
    @public NSString* model;
#ifdef __cplusplus
    @public cv::Mat* k; //sparse matrix
#endif
    @public double* K; //TODO: Figure out which data structure is more appropriate
}


-(id)init;
-(void)destroy;

@end

#endif
