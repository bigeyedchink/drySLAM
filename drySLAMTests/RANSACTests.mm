//
//  RANSACTests.m
//  drySLAM
//
//  Created by Anthony Rodriguez on 3/27/15.
//  Copyright (c) 2015 ECE473 573 Common. All rights reserved.
//


#include "SLAMModel.hpp"
#include "RANSAC.hpp"

#import <UIKit/UIKit.h>
#import <XCTest/XCTest.h>
#include <stdio.h>
#include <stdlib.h>




@interface RANSACTests : XCTestCase

@end

@implementation RANSACTests


- (void)setUp {
    [super setUp];
    // Put setup code here. This method is called before the invocation of each test method in the class.
}

- (void)tearDown {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
    [super tearDown];
}

/************************************************deleteFeatures****************************************************/




















/**********************************PREDICTCAMERAMEASUREMENTS*******************************/
//001: Should return false if any arguments are inappropriately NULL or negative
/*
-(void) testPredictCameraMeasurements{
    int state_size = 10;
    double*     x_km1_k = (double*)malloc(10*sizeof(double));
    Filter* filter = [[Filter alloc]init];
    filter->x_k
    Feature* feat1 = [[Feature alloc]init];
    NSMutableArray* features_info = [[NSMutableArray alloc]init];
    [features_info addObject:feat1];
    Camera* cam = [[Camera alloc]init];
    //x_km1_k==NULL
    bool success = predict_camera_measurements(features_info,filter,cam);
    XCTAssert(success==false, @"PredictCameraMeasurements should retun false when x_km1_k is NULL");
    //state_size==-1
    success = predict_camera_measurements(x_km1_k, -1, features_info, cam);
    XCTAssert(success==false, @"PredictCameraMeasurements should retun false when state_size is negative");
    //features_info==NULL
    success = predict_camera_measurements(x_km1_k, state_size, NULL, cam);
    XCTAssert(success==false, @"PredictCameraMeasurements should retun false when features_info is NULL");
    //cam==NULL
    success = predict_camera_measurements(x_km1_k, state_size, features_info, NULL);
    XCTAssert(success==false, @"PredictCameraMeasurements should retun false when cam is NULL");
    features_info[0] = nil;
    features_info = nil;
    cam = nil;
    free(x_km1_k);
}*/


/**************************************HIINVERSEDEPTH****************************************/
//001: Should immediately return false for any invalid arguments
/*
-(void) testHiInverseDepth001{
    NSMutableArray* zi = [[NSMutableArray alloc]init];
    double* yinit = (double*)malloc(6*sizeof(double));
    double* r_wc = (double*)malloc(3*3*sizeof(double));
    double* t_wc = (double*)malloc(3*sizeof(double));
    Camera* cam = [[Camera alloc]init];
    bool success = hi_inverse_depth(NULL, yinit, t_wc, r_wc, cam);
    XCTAssert(success==false, @"Should return false for null argument zi");
    success = hi_inverse_depth(zi, NULL, t_wc, r_wc, cam);
    XCTAssert(success==false, @"Should return false for null argument yinit");
    success = hi_inverse_depth(zi, yinit, NULL, r_wc, cam);
    XCTAssert(success==false, @"Should return false for null argument t_wc");
    success = hi_inverse_depth(zi, yinit, t_wc, NULL, cam);
    XCTAssert(success==false, @"Should return false for null argument r_wc");
    success = hi_inverse_depth(zi, yinit, t_wc, r_wc, NULL);
    XCTAssert(success==false, @"Should return false for null argument cam");
    free(yinit);
    free(r_wc);
    free(t_wc);
    zi = nil;
    cam = nil;
}
//002: Various inputs from pgm
-(void) testHiInverseDepth002{
    NSMutableArray* zi = [[NSMutableArray alloc]init];
    double* yinit = (double*)malloc(6*sizeof(double));
    double* r_wc = (double*)malloc(3*3*sizeof(double));
    double* t_wc = (double*)malloc(3*sizeof(double));
    Camera* cam = [[Camera alloc]init];
    cam->nRows = 240;
    cam->nCols = 320;
    cam->Cx = 160.2232142857143;
    cam->Cy = 128.8660714285714;
    cam->f = 2.1735000000000;
    cam->dx = 0.011200000000;
    cam->dy = 0.011200000000;
    cam->k1 = 0.063330000000;
    cam->k2 = 0.013900000000;
    yinit[0] = 0.0; yinit[1] = 0.0; yinit[2] = 0.0;
    yinit[3] = -0.556172386921866;
    yinit[4] =  0.008733111023816;
    yinit[5] =  0.991308209869086;
    t_wc[0]  =  0.065845939558637;
    t_wc[1]  = -0.023791341798149;
    t_wc[2]  = -0.023240229497544;
    NSString* path = [[NSBundle mainBundle] pathForResource:@"hid-002_r_wc"
                                                     ofType:@"txt"];
    FILE* fp_r_wc = fopen((char*)[path UTF8String], "r");
    XCTAssert(fp_r_wc!=NULL && path!=nil, @"Setup problem");
    for (int i=0; i<9; ++i){
        fscanf(fp_r_wc, "%lf ", (r_wc+i));
        NSLog(@"%f", r_wc[i]);
    }
    bool success = hi_inverse_depth(zi, yinit, t_wc, r_wc, cam);
    XCTAssert(success==true, @"Should return true for valid data");
    double diff;
    if ([zi count]!=2){
        XCTAssert(NO, @"zi returned an empty set or too many values");
        free(yinit);
        free(r_wc);
        free(t_wc);
        return;
    }
    NSNumber* val = (NSNumber*)zi[0];
    double val1d = [val doubleValue];
    diff = fabs(30.80383321208684-val1d);
    XCTAssert(diff<0.0001, @"zi[0] returned value error exceeds limits");
    val = (NSNumber*)zi[1];
    val1d = [val doubleValue];
    diff = fabs(139.3039860264888-val1d);
    XCTAssert(diff<0.0001, @"zi[1] returned value error exceeds limits");
    free(r_wc);
    free(t_wc);
    free(yinit);
}*/

/************************************************HU*********************************************/
//001: valid input from program
/*
-(void)testHu001{
    double* uv_u = (double*)malloc(2*sizeof(double));
    double* yi = (double*)malloc(3*sizeof(double));
    Camera* cam = [[Camera alloc]init];
    int yi_numColumns = 1;
    yi[0] =    -0.480289739712039;
    yi[1] =     0.187875239166189;
    yi[2] =     0.919080447748411;
    cam->Cx = 160.2232142857143;
    cam->Cy = 128.8660714285714;
    cam->f =    2.1735000000000;
    cam->dx =   0.0112000000000;
    cam->dy =   0.0112000000000;
    bool success = hu(uv_u, yi, yi_numColumns, cam);
    XCTAssert(success = true, @"Expect return of true for valid data");
    double diff = fabs(uv_u[0] - 58.8107341908468);
    XCTAssert(diff<100, @"Error for uv_u exceeds limits");
    diff = fabs(uv_u[1] - 1.685356549672139);
    XCTAssert(diff<0.000001, @"Error for uv_u exceeds limits");
    free(uv_u);
    free(yi);
    cam = nil;
}

/******************************************DISTORT_FM********************************************/
//001: Testing against input from program
/*
-(void)testDistortFm{
    double* uv_u = (double*)malloc(2*sizeof(double));
    double* uvd = (double*)malloc(2*sizeof(double));
    Camera* cam = [[Camera alloc]init];
    cam->Cx = 160.2232142857143;
    cam->Cy = 128.8660714285714;
    cam->k1 = 0.063330000000;
    cam->k2 = 0.013900000000;
    cam->dx = 0.0112000000000;
    cam->dy = 0.0112000000000;
    uv_u[0] = 206.0657952179708;
    uv_u[1] = 152.0438869747222;
    bool success = distort_fm(uvd, uv_u, 1, cam);
    XCTAssert(success==true, @"Should return true for valid data");
    double diff = fabs(uvd[0]-205.1014227508562);
    XCTAssert(diff<100, @"uvd error exceeds limits");
    diff = fabs(uvd[1] - 251.5563043232655);
    XCTAssert(diff<0.00001, @"uvd error exceeds limits");
    free(uvd);
    free(uv_u);
    cam = nil;
}*/

@end



