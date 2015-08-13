//
//  Camera.m
//  drySLAM
//
//  Created by Anthony Rodriguez on 4/1/15.
//  Copyright (c) 2015 ECE473 573 Common. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "Camera.hpp"

@implementation Camera: NSObject


-(id)init{
    //These values don't ever seem to change
    Camera* camera = self;
    camera->k1 = 0.0633;
    camera->k2 = 0.0139;;
    camera->nRows = 240;
    camera->nCols = 320;
    camera->Cx = 0160.2232;
    camera->Cy = 128.8861;
    camera->f = 2.1735;
    camera->dx = 0.0112;
    camera->dy = 0.0112;
    camera->K = (double*)malloc(9*sizeof(double));
    camera->k = new cv::Mat;
    *camera->k = cv::Mat(3,3,CV_64FC1,0.0);
    camera->model = [NSString stringWithFormat:@"two_distortion_parameters"];
    return camera;
};

-(void)destroy{
    if (k!=NULL && k->data!=NULL)
        k->release();
    if (K!=NULL)
        free(K);
}




@end