//
//  MapManagementTests.m
//  drySLAM
//
//  Created by Anthony Rodriguez on 3/27/15.
//  Copyright (c) 2015 ECE473 573 Common. All rights reserved.
//


#include "SLAMModel.hpp"
#include "MapManagement.hpp"

#import <UIKit/UIKit.h>
#import <XCTest/XCTest.h>
#include <stdio.h>
#include <stdlib.h>




@interface MapManagementTests : XCTestCase

@end

@implementation MapManagementTests



- (void)setUp {
    [super setUp];
    // Put setup code here. This method is called before the invocation of each test method in the class.
}

- (void)tearDown {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
    [super tearDown];
}

/************************************************MAP_MANAGEMENT**************************************************/

-(void)testMapManagement000{
    XCTAssert(YES, @"Function requires validation");
}

/************************************************DELETE_FEATURES**************************************************/
//TODO: MAKE ANOTHER TEST FOR THIS
/*
-(void)testDeleteFeatures004{
    NSMutableArray* features_info = [[NSMutableArray alloc]init];
    for (int i=0; i<31; ++i){
        Feature* feat = [[Feature alloc]init];
        feat->times_measured = 12;
        feat->times_measured = 12;
        feat->type = [NSString stringWithFormat:@"inversedepth"];
        [features_info addObject:feat];
    }
    ((Feature*)features_info[22])->times_measured = 4;
    ((Feature*)features_info[22])->times_predicted = 9;
    //int stateSize0 = 199-6;
    //int stateSIze1 = 199;
    int stateSize0 = 199;
    int stateSIze1 = 199-6;
    
    NSString* path = [[NSBundle mainBundle] pathForResource:@"df-004_x_k" ofType:@"txt"];
    FILE* fp_x0 = fopen((char*)[path UTF8String], "r");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"df-004_p_k" ofType:@"txt"];
    FILE* fp_p0 = fopen((char*)[path UTF8String], "r");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"df-004_x_k_after" ofType:@"txt"];
    FILE* fp_x1 = fopen((char*)[path UTF8String], "r");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"df-004_p_k_after" ofType:@"txt"];
    FILE* fp_p1 = fopen((char*)[path UTF8String], "r");
    path = nil;
    
    Filter* filter = [[Filter alloc]init];
    filter->state_size = stateSize0;
    filter->x_k_k = (double*)malloc(stateSize0*sizeof(double));
    for (int i=0; i<stateSize0; ++i){
        fscanf(fp_x0, "%lf ", filter->x_k_k+i);
    }
    fclose(fp_x0);
    double* x_expected = (double*)malloc(stateSIze1*sizeof(double));
    for (int i=0; i<stateSIze1; ++i){
        fscanf(fp_x1, "%lf ", x_expected+i);
    }
    filter->p_k_k = (double*)malloc(stateSize0*stateSize0*sizeof(double));
    for (int i=0; i<stateSize0*stateSize0; ++i){
        fscanf(fp_p0, "%lf ", filter->p_k_k);
    }
    fclose(fp_p0);
    double* p_expected = (double*)malloc(stateSIze1*stateSIze1*sizeof(double));
    for (int i=0; i<stateSIze1*stateSIze1; ++i){
        fscanf(fp_p1, "%lf ", p_expected+i);
    }
    fclose(fp_p1);
    
    bool success = deleteFeatures(features_info, filter);
    XCTAssert(success==true, @"expected return value of true for valid data");
    
    XCTAssert(filter->state_size==stateSIze1, @"State size was not reset as expected");
    XCTAssert([features_info count]==30, @"features_info not resized as expected");
    
    for (int i=0; i<stateSIze1; ++i){
        double actual = filter->x_k_k[i];
        double expected = x_expected[i];
        double diff = fabs(actual-expected);
        if (diff>0.00001){
            XCTAssert(NO, @"x error exceeds limits");
            break;
        }
    }
    for (int i=0; i<stateSIze1*stateSIze1; ++i){
        double actual = filter->p_k_k[i];
        double expected = p_expected[i];
        double diff = fabs(actual-expected);
        if (diff>0.00001){
            XCTAssert(NO, @"x error exceeds limits");
            break;
        }
    }
    
    [filter destroy];
    for (int i=0; i<[features_info count]; ++i){
        [((Feature*)features_info[i]) destroy];
    }
    free(x_expected);
    free(p_expected);
    
}*/
//NOTE: We need to write more tests for this - In matlab, this function works with the designated feature being ==23, which means it should be 22 in our test...
//However, deleting 23 (the 24th feature) works and 22 doesn't
-(void)testDeleteFeatures005{
    NSMutableArray* features_info = [[NSMutableArray alloc]init];
    for (int i=0; i<32; ++i){
        Feature* feat = [[Feature alloc]init];
        feat->times_measured = 11;
        feat->times_measured = 11;
        feat->type = [NSString stringWithFormat:@"inversedepth"];
        [features_info addObject:feat];
    }
    ((Feature*)features_info[23])->times_measured = 4;
    ((Feature*)features_info[23])->times_predicted = 9;

    int stateSize0 = 205;
    int stateSize1 = 199;
    
    NSString* path = [[NSBundle mainBundle] pathForResource:@"df-005_x_k_k_before" ofType:@"txt"];
    FILE* fp_x0 = fopen((char*)[path UTF8String], "r");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"df-005_p_k_k_before" ofType:@"txt"];
    FILE* fp_p0 = fopen((char*)[path UTF8String], "r");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"df-005_x_k_k_after" ofType:@"txt"];
    FILE* fp_x1 = fopen((char*)[path UTF8String], "r");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"df-005_p_k_k_after" ofType:@"txt"];
    FILE* fp_p1 = fopen((char*)[path UTF8String], "r");
    path = nil;
    
    Filter* filter = [[Filter alloc]init];
    filter->state_size = stateSize0;
    filter->x_k_k = (double*)malloc(stateSize0*sizeof(double));
    double val;
    for (int i=0; i<stateSize0; ++i){
        fscanf(fp_x0, "%lf ", &val);
        filter->x_k_k[i] = val;
    }
    fclose(fp_x0);
    double* x_expected = (double*)malloc(stateSize1*sizeof(double));
    for (int i=0; i<stateSize1; ++i){
        fscanf(fp_x1, "%lf ", x_expected+i);
    }
    fclose(fp_x1);
    filter->p_k_k = (double*)malloc(stateSize0*stateSize0*sizeof(double));
    for (int i=0; i<stateSize0*stateSize0; ++i){
        fscanf(fp_p0, "%lf ", &val);
        filter->p_k_k[i] = val;
    }
    fclose(fp_p0);
    double* p_expected = (double*)malloc(stateSize1*stateSize1*sizeof(double));
    for (int i=0; i<stateSize1*stateSize1; ++i){
        fscanf(fp_p1, "%lf ", p_expected+i);
    }
    fclose(fp_p1);
    
    bool success = deleteFeatures(features_info, filter);
    XCTAssert(success==true, @"expected return value of true for valid data");
    
    XCTAssert(filter->state_size==stateSize1, @"State size was not reset as expected");
    XCTAssert([features_info count]==31, @"features_info not resized as expected");
    
    for (int i=0; i<filter->state_size; ++i){
        double actual = filter->x_k_k[i];
        double expected = x_expected[i];
        double diff = fabs(actual-expected);
        if (diff>0.00001){
            XCTAssert(NO, @"x error exceeds limits");
            break;
        }
    }
    for (int i=0; i<filter->state_size*filter->state_size; ++i){
        double actual = filter->p_k_k[i];
        double expected = p_expected[i];
        double diff = fabs(actual-expected);
        if (diff>0.00001){
            XCTAssert(NO, @"p error exceeds limits");
            break;
        }
    }
    
    [filter destroy];
    for (int i=0; i<[features_info count]; ++i){
        [((Feature*)features_info[i]) destroy];
    }
    free(x_expected);
    free(p_expected);
    
}


/************************************************DELETEONEFEATURE**************************************************/
-(void)testDeleteOneFeature001{
    NSMutableArray* features_info = [[NSMutableArray alloc]init];
    for(int i=0; i<32; ++i){
        Feature* feat = [[Feature alloc]init];
        feat->type = [NSString stringWithFormat:@"inversedepth"];
        [features_info addObject:feat];
    }
    int idxToDelete = 23;
    Filter* filter = [[Filter alloc]init];
    filter->state_size = 205;
    
    filter->x_k_k = (double*)malloc(205*sizeof(double));
    filter->p_k_k = (double*)malloc(205*205*sizeof(double));
    double* x_k_expected = (double*)malloc(199*sizeof(double));
    double* p_k_expected = (double*)malloc(199*199*sizeof(double));
    int stateSizeAfter_expected = 199;
    
    NSString* path = [[NSBundle mainBundle] pathForResource:@"daf-001_x_before"
                                                     ofType:@"txt"];
    FILE* fp_x0 = fopen((char*)[path UTF8String], "r");
    XCTAssert(fp_x0!=NULL && path!=nil, @"Setup problem");
    path = [[NSBundle mainBundle] pathForResource:@"daf-001_p_before"
                                           ofType:@"txt"];
    FILE* fp_p0 = fopen((char*)[path UTF8String], "r");
    XCTAssert(fp_p0!=NULL && path!=nil, @"Setup problem");
    path = [[NSBundle mainBundle] pathForResource:@"daf-001_x_after"
                                           ofType:@"txt"];
    FILE* fp_x1 = fopen((char*)[path UTF8String], "r");
    XCTAssert(fp_x1!=NULL && path!=nil, @"Setup problem");
    path = [[NSBundle mainBundle] pathForResource:@"daf-001_p_after"
                                           ofType:@"txt"];
    FILE* fp_p1 = fopen((char*)[path UTF8String], "r");
    XCTAssert(fp_p1!=NULL && path!=nil, @"Setup problem");
    
    for(int i=0; i<filter->state_size; ++i){
        fscanf(fp_x0, "%lf ", filter->x_k_k+i);
    }
    for(int i=0; i<(filter->state_size*filter->state_size); ++i){
        fscanf(fp_p0, "%lf ", filter->p_k_k+i);
    }
    for (int i=0; i<stateSizeAfter_expected; ++i){
        fscanf(fp_x1, "%lf ", x_k_expected+i);
    }
    for (int i=0; i<stateSizeAfter_expected*stateSizeAfter_expected; ++i){
        fscanf(fp_p1, "%lf ", p_k_expected+i);
    }
    
    bool success = deleteOneFeature(features_info, filter, idxToDelete);
    XCTAssert(success==true, @"Should return true for valid data");
    
    if (filter->state_size!=stateSizeAfter_expected){
        XCTAssert(filter->state_size==stateSizeAfter_expected, @"State size not reset correctly!");
        [filter destroy];
        for (int i=0; i<[features_info count]; ++i){
            Feature* feat = features_info[i];
            [feat destroy];
        }
        free(x_k_expected);
        free(p_k_expected);
        return;
    }
    
    for (int i=0; i<filter->state_size; ++i){
        double actual = filter->x_k_k[i];
        double expected = x_k_expected[i];
        double diff = fabs(actual-expected);
        if (diff>0.00001){
            XCTAssert(NO, @"x_k error exceeds limits");
            break;
        }
    }
    for (int i=0; i<filter->state_size*filter->state_size; ++i){
        double actual = filter->p_k_k[i];
        double expected = p_k_expected[i];
        double diff = fabs(actual-expected);
        if (diff>0.00001){
            XCTAssert(NO, @"p_k error exceeds limits");
            break;
        }
    }
    [filter destroy];
    for (int i=0; i<[features_info count]; ++i)
    {
        Feature* feat = features_info[i];
        [feat destroy];
    }
    
    free(x_k_expected);
    free(p_k_expected);
}

/******************************************************UPDATE_FEATURES_INFO*****************************************************/

-(void)testUpdateFeaturesInfo001{
    NSMutableArray* features_info = [[NSMutableArray alloc]init];
    for (int i=0; i<12; ++i){
        Feature* feat = [[Feature alloc]init];
        feat->times_predicted = 0;
        feat->times_measured = 0;
        feat->low_innovation_inlier = 1;
        feat->h = (double*)malloc(2*sizeof(double));
        feat->h[0] = 99;
        feat->h[1] = 100;
        feat->z = (int*)malloc(2*sizeof(int));
        feat->z[0] = 99;
        feat->z[1] = 100;
        feat->H = new cv::Mat(2,85,CV_64FC1,0.0);
        feat->S = (double*)malloc(4*sizeof(double));
        feat->S[0] = 99;
        feat->S[1] = 100;
        feat->S[2] = 101;
        feat->S[3] = 102;
        feat->times_measured = 0;
        if (i==8){
            feat->low_innovation_inlier = 0;
            feat->high_innovation_inlier = 1;
        }
        [features_info addObject:feat];
    }
    bool success = updateFeaturesInfo(features_info);
    XCTAssert(success==true, @"Should return true for valid data");
    for (int i=0; i<12; ++i){
        Feature* feat = (Feature*)features_info[i];
        XCTAssert(feat->times_measured==1, @"times_measured should be 1");
        XCTAssert(feat->times_predicted==1, @"times_predicted should be 1");
        XCTAssert(feat->low_innovation_inlier==0, @"low_innovation_inlier should be 0");
        XCTAssert(feat->high_innovation_inlier==0, @"high_innovation_inlier should be 0");
        XCTAssert(feat->individually_compatible==0, @"individually_compatible should be 1");
        XCTAssert(feat->h==NULL, @"h should be NULL");
        XCTAssert(feat->z==NULL, @"z should be NULL");
        XCTAssert(feat->H==NULL, @"H should be NULL");
        XCTAssert(feat->S==NULL, @"S should be NULL");
        [feat destroy];
    }
}

/**************************************************INVERSEDEPTH_2_CARTESIAN**************************************************/

-(void)testInverseDepth_2_Cartesian001{
    NSMutableArray* features_info = [[NSMutableArray alloc]init];
    Filter* filter = [[Filter alloc]init];
    bool success = inversedepth_2_cartesian(features_info, filter);
    XCTAssert(success==true, @"Expected true return value");
    features_info = nil;
    [filter destroy];
}

/************************************************INITIALIZE_FEATURES**************************************************/

-(void)testInitializeFeatures000{
    //Implement as necessary
    XCTAssert(YES, @"Function requires validation");
}

/************************************************INITIALIZE_A_FEATURE**************************************************/
//Just some dummy data so that we can monitor function in action
-(void)testInitializeAFeature001{
    Camera* cam = [[Camera alloc]init];
    cam->k1 = 0.063330000000000;
    cam->k2 = 0.013900000000000;
    cam->nRows = 240;
    cam->nCols = 320;
    cam->Cx = 100*1.602232142857143;
    cam->Cy = 100*1.288660714285714;
    cam->dx = 0.011200000000000;
    cam->dy = 0.011200000000000;
    cam->f = 2.1735000000;
    cam->K = (double*)malloc(9*sizeof(double));
    cam->K[0] = 100*1.940625000000000;
    cam->K[1] = 100*0;
    cam->K[2] = 100*0;
    cam->K[3] = 100*1.940625000000000;
    cam->K[4] = 100*0;
    cam->K[5] = 100*0;
    cam->K[6] = 100*1.602232142857143;
    cam->K[7] = 100*1.288660714285715;
    cam->K[8] = 100*0.010000000000000;
    uchar* img = (uchar*)malloc(240*320*sizeof(uchar));
    NSString* path = [[NSBundle mainBundle] pathForResource:@"iaf-001_im" ofType:@"txt"];
    FILE* fp_im = fopen((char*)[path UTF8String], "r");
    path = nil;
    for (int i=0; i<(320*240); ++i){
        double val;
        fscanf(fp_im, "%lf ", &val);
        img[i] = (uchar)((int)val);
    }
    cv::Mat* im_k = new cv::Mat(240, 320, CV_8UC1, img, cv::Mat::AUTO_STEP);
    fclose(fp_im);
    Filter* filter = [[Filter alloc]init];
    filter->state_size = 85;
    filter->x_k_k = (double*)malloc(85*sizeof(double));
    filter->p_k_k = (double*)malloc(85*85*sizeof(double));
    path = [[NSBundle mainBundle] pathForResource:@"iaf-001_x_k" ofType:@"txt"];
    FILE* fp_x = fopen((char*)[path UTF8String], "r");
    path = nil;
    for (int i=0; i<filter->state_size; ++i){
        fscanf(fp_x, "%lf ", filter->x_k_k+i);
    }
    fclose(fp_x);
    path = [[NSBundle mainBundle] pathForResource:@"iaf-001_p_k" ofType:@"txt"];
    FILE* fp_p = fopen((char*)[path UTF8String], "r");
    path = nil;
    for (int i=0; i<filter->state_size*filter->state_size; ++i){
        fscanf(fp_p, "%lf ", filter->p_k_k+i);
    }
    fclose(fp_p);
    NSMutableArray* features_info = [[NSMutableArray alloc]init];
    for (int i=0; i<18; ++i){
        Feature* feat = [[Feature alloc]init];
        feat->h = (double*)malloc(2*sizeof(double));
        feat->h[0] = i;
        feat->h[1] = 18 - i;
        feat->type = [NSString stringWithFormat:@"inversedepth"];
        [features_info addObject:feat];
    }
    int step = 93;
    NSMutableArray* uv = [[NSMutableArray alloc]init];
    bool success = initialize_a_feature(filter, features_info, step, cam, im_k, uv);
    free(img);
    im_k->release();
    [filter destroy];
    [cam destroy];
    for (int i=0; i<[features_info count]; ++i){
        [(Feature*)features_info[i] destroy];
    }
    XCTAssert(success==true);
}

/************************************************ADD_FEATURES_INVERSE_DEPTH**************************************************/

-(void)testAddFeaturesInverseDepth001{
    double* uvd = (double*)malloc(2*sizeof(double));
    uvd[0] = 266;
    uvd[1] = 35;
    int stateSizeBefore = 127;
    int stateSizeAfter = 133;
    double* x0 = (double*)malloc(stateSizeBefore*sizeof(double));
    double* x1_expected = (double*)malloc(stateSizeAfter*sizeof(double));

    double* p0 = (double*)malloc(stateSizeBefore*stateSizeBefore*sizeof(double));
    double* p1_expected = (double*)malloc(stateSizeAfter*stateSizeAfter*sizeof(double));

    double* newFeature_expected = (double*)malloc(6*sizeof(double));
    double* newFeature_actual = (double*)malloc(6*sizeof(double));
    newFeature_expected[0] = 0.029077829480646;
    newFeature_expected[1] = -0.006490103047387;
    newFeature_expected[2] = -0.000412392877165;
    newFeature_expected[3] = 0.631739664650969;
    newFeature_expected[4] = 0.459982779575005;
    newFeature_expected[5] = 1.000000000000000;
    
    Camera* cam = [[Camera alloc]init];
    cam->k1 = 0.063330000000000;
    cam->k2 = 0.013900000000000;
    cam->nRows = 240;
    cam->nCols = 320;
    cam->Cx = 100*1.602232142857143;
    cam->Cy = 100*1.288660714285714;
    cam->dx = 0.011200000000000;
    cam->dy = 0.011200000000000;
    cam->f = 2.1735000000;
    cam->K = (double*)malloc(9*sizeof(double));
    cam->K[0] = 100*1.940625000000000;
    cam->K[1] = 100*0;
    cam->K[2] = 100*0;
    cam->K[3] = 100*1.940625000000000;
    cam->K[4] = 100*0;
    cam->K[5] = 100*0;
    cam->K[6] = 100*1.602232142857143;
    cam->K[7] = 100*1.288660714285715;
    cam->K[8] = 100*0.010000000000000;
    
    NSString* path = [[NSBundle mainBundle] pathForResource:@"afid-001_X_before" ofType:@"txt"];
    FILE* fp_x0 = fopen((char*)[path UTF8String], "r");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"afid-001_X_after" ofType:@"txt"];
    FILE* fp_x1 = fopen((char*)[path UTF8String], "r");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"afid-001_P_before" ofType:@"txt"];
    FILE* fp_p0 = fopen((char*)[path UTF8String], "r");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"afid-001_P_after" ofType:@"txt"];
    FILE* fp_p1 = fopen((char*)[path UTF8String], "r");
    path = nil;
    
    for (int i=0; i<stateSizeBefore; ++i){
        fscanf(fp_x0, "%lf ", x0+i);
    }
    fclose(fp_x0);
    for (int i=0; i<stateSizeAfter; ++i){
        fscanf(fp_x1, "%lf ", x1_expected+i);
    }
    fclose(fp_x1);
    for (int i=0; i<stateSizeBefore*stateSizeBefore; ++i){
        fscanf(fp_p0, "%lf ", p0+i);
    }
    fclose(fp_p0);
    for (int i=0; i<stateSizeAfter*stateSizeAfter; ++i){
        fscanf(fp_p1, "%lf ", p1_expected+i);
    }
    fclose(fp_p1);
    
    Filter* filter = [[Filter alloc]init];
    filter->x_k_k = x0;
    filter->p_k_k = p0;
    filter->state_size = stateSizeBefore;
    
    bool success = add_features_inverse_depth(uvd, 1, filter, cam, 1.0, 1.0, 1.0, newFeature_actual);
    XCTAssert(success==true, @"Expected true return for valid data");
    
    for (int i=0; i<6; ++i){
        double actual = newFeature_actual[i];
        double expected = newFeature_expected[i];
        double diff = fabs(actual-expected);
        if (diff>0.00001){
            XCTAssert(NO, @"newFeature error exceeds limits");
            break;
        }
    }
    XCTAssert(filter->state_size==stateSizeAfter, @"State size in filter was never updated to correct value");
    for (int i=0; i<stateSizeAfter; ++i){
        double actual = filter->x_k_k[i];
        double expected = x1_expected[i];
        double diff = fabs(actual-expected);
        if (diff>0.00001){
            XCTAssert(NO, @"x_k_k error exceeds limits");
            break;
        }
    }
    for (int i=0; i<stateSizeAfter*stateSizeAfter; ++i){
        double actual = filter->p_k_k[i];
        double expected = p1_expected[i];
        double diff = fabs(actual-expected);
        if (diff>0.00001){
            XCTAssert(NO, @"p_k_k error exceeds limits");
            break;
        }
    }
    
    [filter destroy];
    [cam destroy];
    free(uvd);
    free(x1_expected);
    free(p1_expected);
    free(newFeature_expected);
    free(newFeature_actual);
}

/************************************************HINV**************************************************/

-(void)testHinv001{
    double* uvd = (double*)malloc(2*sizeof(double));
    uvd[0] = 100;
    uvd[1] = 165;
    
    double* xv = (double*)malloc(13*sizeof(double));
    xv[0] = 0.070312527690209;
    xv[1] = -0.027708469946418;
    xv[2] = -0.024897159877980;
    xv[3] = 0.999071851232890;
    xv[4] = 0.013917621295685;
    xv[5] = 0.038894493895543;
    xv[6] = 0.012204680912096;
    xv[7] = 0.003291561833459;
    xv[8] = -0.003954620955759;
    xv[9] = -0.001914400443389;
    xv[10] = 0.002122152156030;
    xv[11] = 0.001150856858162;
    xv[12] = -0.002426112234365;
    
    Camera* cam = [[Camera alloc]init];
    cam->k1 = 0.063330000000000;
    cam->k2 = 0.013900000000000;
    cam->nRows = 240;
    cam->nCols = 320;
    cam->Cx = 100*1.602232142857143;
    cam->Cy = 100*1.288660714285714;
    cam->dx = 0.011200000000000;
    cam->dy = 0.011200000000000;
    cam->f = 2.1735000000;
    cam->K = (double*)malloc(9*sizeof(double));
    cam->K[0] = 100*1.940625000000000;
    cam->K[1] = 100*0;
    cam->K[2] = 100*0;
    cam->K[3] = 100*1.940625000000000;
    cam->K[4] = 100*0;
    cam->K[5] = 100*0;
    cam->K[6] = 100*1.602232142857143;
    cam->K[7] = 100*1.288660714285715;
    cam->K[8] = 100*0.010000000000000;
    
    double* newFeature_expected = (double*)malloc(6*sizeof(double));
    double* newFeature_actual = (double*)malloc(6*sizeof(double));
    newFeature_expected[0] = 0.070312527690209;
    newFeature_expected[1] = -0.027708469946418;
    newFeature_expected[2] = -0.024897159877980;
    newFeature_expected[3] = -0.238300950432964;
    newFeature_expected[4] = -0.149505074503326;
    newFeature_expected[5] = 1.000000000000000;
    
    bool success = hinv(uvd, xv, cam, 1.0, newFeature_actual);
    XCTAssert(success==true, @"expected true return for valid data");
    for (int i=0; i<6; ++i){
        double actual = newFeature_actual[i];
        double expected = newFeature_expected[i];
        double diff = fabs(actual-expected);
        if (diff>0.00001){
            XCTAssert(NO, @"newFeature error exceeds limits");
            break;
        }
    }
    
    
    
    free(newFeature_expected);
    free(newFeature_actual);
    free(uvd);
    [cam destroy];
}

/************************************************ADDONEFEATURECOVARIANCEID**************************************************/

-(void)testAddOneFeatureCovarianceInverseDepth001{
    double* uvd = (double*)malloc(2*sizeof(double));
    uvd[0] = 44.0;
    uvd[1] = 89.0;
    double* xv = (double*)malloc(13*sizeof(double));
    xv[0] = 0.065845939558637;
    xv[1] = -0.023791341798149;
    xv[2] = -0.023240229497544;
    xv[3] = 0.999103798298479;
    xv[4] = 0.012692065225640;
    xv[5] = 0.038139953756084;
    xv[6] = 0.013261057022112;
    xv[7] = 0.004777281162448;
    xv[8] = -0.003589470132243;
    xv[9] = -0.004429777630095;
    xv[10] = 0.002614243593750;
    xv[11] = 0.003954661511372;
    xv[12] = -0.000404077038371;
    double stdPxl = 1.0;
    double stdRho = 1.0;
    Camera* cam = [[Camera alloc]init];
    cam->k1 = 0.063330000000000;
    cam->k2 = 0.013900000000000;
    cam->nRows = 240;
    cam->nCols = 320;
    cam->Cx = 100*1.602232142857143;
    cam->Cy = 100*1.288660714285714;
    cam->dx = 0.011200000000000;
    cam->dy = 0.011200000000000;
    cam->f = 2.1735000000;
    cam->K = (double*)malloc(9*sizeof(double));
    cam->K[0] = 100*1.940625000000000;
    cam->K[1] = 100*0;
    cam->K[2] = 100*0;
    cam->K[3] = 100*1.940625000000000;
    cam->K[4] = 100*0;
    cam->K[5] = 100*0;
    cam->K[6] = 100*1.602232142857143;
    cam->K[7] = 100*1.288660714285715;
    cam->K[8] = 100*0.010000000000000;
    cam->k = new cv::Mat(3, 3, CV_64FC1, cam->K, cv::Mat::AUTO_STEP);
    
    int stateSize = 169;
    int stateSizeAfter = 169+6;
    NSString* path = [[NSBundle mainBundle] pathForResource:@"aafciv-001_P"
                                                     ofType:@"txt"];
    FILE* fp_P = fopen((char*)[path UTF8String], "r");
    path = nil;
    double* P = (double*)malloc(169*169*sizeof(double));
    for (int i=0; i<stateSize; ++i){
        fscanf(fp_P, "%lf ", P+i);
    }
    fclose(fp_P);
    path = [[NSBundle mainBundle] pathForResource:@"aafciv-001_P_RES"
                                           ofType:@"txt"];
    fp_P = fopen((char*)[path UTF8String], "r");
    path = nil;
    double* P_RES_actual = (double*)malloc(stateSizeAfter*stateSizeAfter*sizeof(double));
    double* P_RES_expected = (double*)malloc(stateSizeAfter*stateSizeAfter*sizeof(double));
    for (int i=0; i<stateSizeAfter; ++i){
        fscanf(fp_P, "%lf ", P_RES_expected+i);
        P_RES_actual[i] = 0.0;
    }
    Filter* filter = [[Filter alloc]init];
    filter->p_k_k = P;
    filter->state_size = stateSize;
    bool success = addOneFeatureCovarianceInverseDepth(filter, uvd, xv, stdPxl, stdRho, cam, P_RES_actual);
    
    XCTAssert(success==true, @"Expected return value of true for valid data");
    int i,j;
    for (i=0; i<13; ++i)
    {
        for (j=0; j<13; ++j)
        {
            double actual = P_RES_actual[i];
            double expected = P_RES_expected[i];
            double diff = fabs(actual-expected);
            if (diff>0.00001)
            {
                XCTAssert(NO, @"P_RES error exceeds limits");
                break;
            }
        }
    }
    
    for (i=0; i<13; ++i)
    {
        for (j=13; j<stateSize; ++j)
        {
            double actual = P_RES_actual[i];
            double expected = P_RES_expected[i];
            double diff = fabs(actual-expected);
            if (diff>0.00001)
            {
                XCTAssert(NO, @"P_RES error exceeds limits");
                break;
            }
        }
    }
    
    for (i=0; i<13; ++i)
    {
        for (j=stateSize; j<stateSizeAfter; ++j)
        {
            double actual = P_RES_actual[i];
            double expected = P_RES_expected[i];
            double diff = fabs(actual-expected);
            if (diff>0.00001)
            {
                XCTAssert(NO, @"P_RES error exceeds limits");
                break;
            }
        }
    }
    for (i=13; i<stateSize; ++i)
    {
        for (j=0; j<13; ++j)
        {
            double actual = P_RES_actual[i];
            double expected = P_RES_expected[i];
            double diff = fabs(actual-expected);
            if (diff>0.00001)
            {
                XCTAssert(NO, @"P_RES error exceeds limits");
                break;
            }
        }
    }
    
    for (i=13; i<stateSize; ++i)
    {
        for (j=13; j<stateSize; ++j)
        {
            double actual = P_RES_actual[i];
            double expected = P_RES_expected[i];
            double diff = fabs(actual-expected);
            if (diff>0.00001)
            {
                XCTAssert(NO, @"P_RES error exceeds limits");
                break;
            }
        }
    }
    
    for (i=13; i<stateSize; ++i)
    {
        for (j=stateSize; j<stateSizeAfter; ++j)
        {
            double actual = P_RES_actual[i];
            double expected = P_RES_expected[i];
            double diff = fabs(actual-expected);
            if (diff>0.00001)
            {
                XCTAssert(NO, @"P_RES error exceeds limits");
                break;
            }
        }
    }
    
    for (i=stateSize; i<stateSizeAfter; ++i)
    {
        for (j=0; j<13; ++j)
        {
            double actual = P_RES_actual[i];
            double expected = P_RES_expected[i];
            double diff = fabs(actual-expected);
            if (diff>0.00001)
            {
                XCTAssert(NO, @"P_RES error exceeds limits");
                break;
            }
        }
    }
    
    for (i=stateSize; i<stateSizeAfter; ++i)
    {
        for (j=13; j<stateSize; ++j)
        {
            double actual = P_RES_actual[i];
            double expected = P_RES_expected[i];
            double diff = fabs(actual-expected);
            if (diff>0.00001)
            {
                XCTAssert(NO, @"P_RES error exceeds limits");
                break;
            }
        }
    }
    
    for (i=stateSize; i<stateSizeAfter; ++i)
    {
        for (j=stateSize; j<stateSizeAfter; ++j)
        {
            double actual = P_RES_actual[i];
            double expected = P_RES_expected[i];
            double diff = fabs(actual-expected);
            if (diff>0.00001)
            {
                XCTAssert(NO, @"P_RES error exceeds limits");
                break;
            }
        }
    }


    [filter destroy];
    free(P_RES_actual);
    free(P_RES_expected);
    [cam destroy];
    free(xv);
    free(uvd);
    
}

/**************************************************PREDICT_CAMERA_MEASUREMENTS**************************************************/

-(void)testPredictCameraMeasurements001{
    Camera* cam = [[Camera alloc]init];
    cam->k1 = 0.063330000000000;
    cam->k2 = 0.013900000000000;
    cam->nRows = 240;
    cam->nCols = 320;
    cam->Cx = 100*1.602232142857143;
    cam->Cy = 100*1.288660714285714;
    cam->dx = 0.011200000000000;
    cam->dy = 0.011200000000000;
    int stateSize = 121;
    double* x = (double*)malloc(stateSize*sizeof(double));
    NSMutableArray* features_info = [[NSMutableArray alloc]init];
    for (int i=0; i<18; ++i){
        Feature* feat = [[Feature alloc]init];
        [features_info addObject:feat];
    }
    
    NSString* path = [[NSBundle mainBundle] pathForResource:@"pcm-001_x"
                                                     ofType:@"txt"];
    FILE* fp_x = fopen((char*)[path UTF8String], "r");
    XCTAssert(fp_x!=NULL && path!=nil, @"Setup problem");
    for (int i=0; i<stateSize; ++i){
        fscanf(fp_x, "%lf ", x+i);
    }
    Filter* filt = [[Filter alloc]init];
    filt->x_k_k = x;
    filt->state_size = stateSize;
    bool success = predict_camera_measurements(features_info, filt, cam);
    
    double actual1 = ((Feature*)(features_info[0]))->h[0];
    double actual2 = ((Feature*)(features_info[0]))->h[1];
    double expected1 = 100*0.460152134362705;
    double expected2 = 100*1.293963086347413;
    XCTAssert(fabs(actual1-expected1)<0.00001, @"predict error");
    XCTAssert(fabs(actual2-expected2)<0.00001, @"predict error");
    
    actual1 = ((Feature*)(features_info[1]))->h[0];
    actual2 = ((Feature*)(features_info[1]))->h[1];
    expected1 = 100*0.859945102394409;
    expected2 = 100*1.549544527602405;
    XCTAssert(fabs(actual1-expected1)<0.00001, @"predict error");
    XCTAssert(fabs(actual2-expected2)<0.00001, @"predict error");
    
    actual1 = ((Feature*)(features_info[2]))->h[0];
    actual2 = ((Feature*)(features_info[2]))->h[1];
    expected1 = 100*2.280206057997468;
    expected2 = 100*1.431463510555592;
    XCTAssert(fabs(actual1-expected1)<0.00001, @"predict error");
    XCTAssert(fabs(actual2-expected2)<0.00001, @"predict error");
    
    actual1 = ((Feature*)(features_info[3]))->h[0];
    actual2 = ((Feature*)(features_info[3]))->h[1];
    expected1 = 100*0.268688391464448;
    expected2 = 100*1.638788869467983;
    XCTAssert(fabs(actual1-expected1)<0.00001, @"predict error");
    XCTAssert(fabs(actual2-expected2)<0.00001, @"predict error");
    
    actual1 = ((Feature*)(features_info[4]))->h[0];
    actual2 = ((Feature*)(features_info[4]))->h[1];
    expected1 = 68.006633078392227;
    expected2 = 44.952324240599211;
    XCTAssert(fabs(actual1-expected1)<0.00001, @"predict error");
    XCTAssert(fabs(actual2-expected2)<0.00001, @"predict error");
    
    actual1 = ((Feature*)(features_info[5]))->h[0];
    actual2 = ((Feature*)(features_info[5]))->h[1];
    expected1 = 76.994047750542023;
    expected2 = 63.918399520153400;
    XCTAssert(fabs(actual1-expected1)<0.00001, @"predict error");
    XCTAssert(fabs(actual2-expected2)<0.00001, @"predict error");
    
    actual1 = ((Feature*)(features_info[6]))->h[0];
    actual2 = ((Feature*)(features_info[6]))->h[1];
    expected1 = 100*0.710388690352725;
    expected2 = 100*1.331345943604877;
    XCTAssert(fabs(actual1-expected1)<0.00001, @"predict error");
    XCTAssert(fabs(actual2-expected2)<0.00001, @"predict error");
    
    actual1 = ((Feature*)(features_info[7]))->h[0];
    actual2 = ((Feature*)(features_info[7]))->h[1];
    expected1 = 100*1.200223388760884;
    expected2 = 100*0.979941942164523;
    XCTAssert(fabs(actual1-expected1)<0.00001, @"predict error");
    XCTAssert(fabs(actual2-expected2)<0.00001, @"predict error");
    
    actual1 = ((Feature*)(features_info[8]))->h[0];
    actual2 = ((Feature*)(features_info[8]))->h[1];
    expected1 = 100*1.589678983718885;
    expected2 = 100*1.039136105791849;
    XCTAssert(fabs(actual1-expected1)<0.00001, @"predict error");
    XCTAssert(fabs(actual2-expected2)<0.00001, @"predict error");

    actual1 = ((Feature*)(features_info[9]))->h[0];
    actual2 = ((Feature*)(features_info[9]))->h[1];
    expected1 = 100*2.768546402863326;
    expected2 = 100*0.321476737756434;
    XCTAssert(fabs(actual1-expected1)<0.00001, @"predict error");
    XCTAssert(fabs(actual2-expected2)<0.00001, @"predict error");
    
    actual1 = ((Feature*)(features_info[10]))->h[0];
    actual2 = ((Feature*)(features_info[10]))->h[1];
    expected1 = 91.001656238051112;
    expected2 = 56.869044789351747;
    XCTAssert(fabs(actual1-expected1)<0.00001, @"predict error");
    XCTAssert(fabs(actual2-expected2)<0.00001, @"predict error");
    //Skip some
    actual1 = ((Feature*)(features_info[16]))->h[0];
    actual2 = ((Feature*)(features_info[16]))->h[1];
    expected1 = 100*1.470694167420400;
    expected2 = 100*1.159773358026855;
    XCTAssert(fabs(actual1-expected1)<0.00001, @"predict error");
    XCTAssert(fabs(actual2-expected2)<0.00001, @"predict error");
    
    actual1 = ((Feature*)(features_info[17]))->h[0];
    actual2 = ((Feature*)(features_info[17]))->h[1];
    expected1 = 100*1.130861521303999;
    expected2 = 100*1.460379501887305;
    XCTAssert(fabs(actual1-expected1)<0.00001, @"predict error");
    XCTAssert(fabs(actual2-expected2)<0.00001, @"predict error");

    
    for (int i=0; i<[features_info count]; ++i){
        [(Feature*)features_info[i] destroy];
    }
    [filt destroy];
    [cam destroy];
    
}


/*******************************************************HI_INVERSE_DEPTH*******************************************************/

-(void)testHiInverseDepth001{
    Camera* cam = [[Camera alloc]init];
    cam->k1 = 0.063330000000000;
    cam->k2 = 0.013900000000000;
    cam->nRows = 240;
    cam->nCols = 320;
    cam->Cx = 100*1.602232142857143;
    cam->Cy = 100*1.288660714285714;
    cam->dx = 0.011200000000000;
    cam->dy = 0.011200000000000;
    double* yinit = (double*)malloc(6*sizeof(double));
    yinit[0] = 0;
    yinit[1] = 0;
    yinit[2] = 0;
    yinit[3] = -0.476444812646242;
    yinit[4] = 0.429516507868721;
    yinit[5] = 1.000000000000000;
    double* t_wc = (double*)malloc(3*sizeof(double));
    t_wc[0] = 0.009139450566591;
    t_wc[1] = 0.001725790855408;
    t_wc[2] = 0.000681597720107;
    double* r_wc = (double*)malloc(9*sizeof(double));
    r_wc[0] = 0.999945094340855;
    r_wc[1] = -0.000944152863288;
    r_wc[2] = 0.010436324977187;
    r_wc[3] = 0.001010630421979;
    r_wc[4] = 0.999979223654465;
    r_wc[5] = -0.006366387167216;
    r_wc[6] = -0.010430097305820;
    r_wc[7] = 0.006376584884048;
    r_wc[8] = 0.999925273325666;
    double* zi_actual = (double*)malloc(2*sizeof(double));
    double* zi_expected = (double*)malloc(2*sizeof(double));
    zi_expected[0] = 70.759105228688412;
    zi_expected[1] = 43.994098116953637;
    bool success = hiInverseDepth(cam, yinit, t_wc, r_wc, zi_actual);
    XCTAssert(success==true, @"expect true for valid data");
    double diff = fabs(zi_expected[0]-zi_actual[0]);
    XCTAssert(diff<0.00001, @"zi error exceeds limits");
    diff = fabs(zi_expected[1]-zi_actual[1]);
    XCTAssert(diff<0.00001, @"zi error exceeds limits");
    free(yinit);
    free(zi_actual);
    free(zi_expected);
    free(t_wc);
    free(r_wc);
    [cam destroy];
}

/*************************************************************HU***************************************************************/
//001: Data values
-(void)testHu001{
    Camera* cam = [[Camera alloc]init];
    cam->k1 = 0.06333;
    cam->k2 = 0.0139;
    cam->nRows = 240;
    cam->nCols = 320;
    cam->Cx = 160.2232142857143;
    cam->Cy = 128.8660714285714;
    cam->f = 2.1735;
    cam->dx = 0.0112;
    cam->dy = 0.0112;
    double* y = (double*)malloc(3*sizeof(double));
    y[0] = -0.172991252307490;
    y[1]= -0.163681444674160;
    y[2] = 0.850672908321281;
    double* uv_u_actual = (double*)malloc(2*sizeof(double));
    double* uv_u_expected = (double*)malloc(2*sizeof(double));
    uv_u_expected[0] = 100*1.207590270845824;
    uv_u_expected[1] = 100*0.915257141110178;
    bool success = hu(uv_u_actual,y,1,cam);
    XCTAssert(success==true, @"Should return true");
    double expected = uv_u_expected[0];
    double actual = uv_u_actual[0];
    double diff = fabs(expected-actual);
    XCTAssert(diff<0.00001, @"uv_u error exceeds limits");
    expected = uv_u_expected[1];
    actual = uv_u_actual[1];
    diff = fabs(expected-actual);
    XCTAssert(diff<0.00001, @"uv_u error exceeds limits");
}

/*********************************************DISTORT_FM************************************************/
//001: Values from program
-(void)testDistortFM001{
    double* uvd_actual = (double*)malloc(2*sizeof(double));
    double* uvd_expected = (double*)malloc(2*sizeof(double));
    double* uvd_out = (double*)malloc(2*sizeof(double));
    Camera* cam = [[Camera alloc]init];
    cam->k1 = 0.06333;
    cam->k2 = 0.0139;
    cam->nRows = 240;
    cam->nCols = 320;
    cam->Cx = 160.2232142857143;
    cam->Cy = 128.8660714285714;
    cam->f = 2.1735;
    cam->dx = 0.0112;
    cam->dy = 0.0112;
    uvd_actual[0] = 100*1.454732458030142;
    uvd_actual[1] = 100*0.902063477809479;
    uvd_expected[0] = 100*1.456745909312229;
    uvd_expected[1] = 100*0.907340741298874;
    
    bool success = distortFM(uvd_out,uvd_actual,1,cam);
    
    XCTAssert(success==true, @"should return true for valid data");
    double actual = uvd_out[0];
    double expected = uvd_expected[0];
    double diff = fabs(actual-expected);
    XCTAssert(diff<0.00001, @"uvd error exceeds limits");
    actual = uvd_out[1];
    expected = uvd_expected[1];
    diff = fabs(actual-expected);
    XCTAssert(diff<0.00001, @"uvd error exceeds limits");
    free(uvd_actual);
    free(uvd_expected);
    free(uvd_out);
    [cam destroy];
    cam = nil;
}

-(void)testDistortFM002{
    Camera* cam = [[Camera alloc]init];
    cam->k1 = 0.063330000000000;
    cam->k2 = 0.013900000000000;
    cam->nRows = 240;
    cam->nCols = 320;
    cam->Cx = 100*1.602232142857143;
    cam->Cy = 100*1.288660714285714;
    cam->dx = 0.011200000000000;
    cam->dy = 0.011200000000000;
    cam->f = 2.1735000000;
    cam->K = (double*)malloc(9*sizeof(double));
    cam->K[0] = 100*1.940625000000000;
    cam->K[1] = 100*0;
    cam->K[2] = 100*1.602232142857143;
    cam->K[3] = 100*0;
    cam->K[4] = 100*1.940625000000000;
    cam->K[5] = 100*1.288660714285715;
    cam->K[6] = 100*0;
    cam->K[7] = 100*0;
    cam->K[8] = 100*0.010000000000000;
    int cols = 169;
    double* uv = (double*)malloc(2*169*sizeof(double));
    double* uvd_actual = (double*)malloc(2*169*sizeof(double));
    double* uvd_expected = (double*)malloc(2*169*sizeof(double));
    
    NSString* path = [[NSBundle mainBundle] pathForResource:@"dfm-002_uv" ofType:@"txt"];
    FILE* fp_uv = fopen((char*)[path UTF8String], "r");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"dfm-002_uvd" ofType:@"txt"];
    FILE* fp_uvd = fopen((char*)[path UTF8String], "r");
    path = nil;
    for (int i=0; i<cols*2; ++i){
        fscanf(fp_uv, "%lf ", uv+i);
        fscanf(fp_uvd, "%lf ", uvd_expected+i);
        uvd_actual[i] = 0.0;
    }
    fclose(fp_uvd);
    fclose(fp_uv);
    
    bool success = distortFM(uvd_actual, uv, cols, cam);
    
    for (int i=0; i<cols*2; ++i){
        double actual = uvd_actual[i];
        double expected = uvd_expected[i];
        double diff = fabs(actual-expected);
        if (diff>0.00001){
            XCTAssert(NO, @"uvd error exceeds limits");
            break;
        }
    }
    free(uvd_expected);
    free(uvd_actual);
    free(uv);
    
}

/******************************************ADD_FEATURE_TO_INFO_VECTOR******************************************/

-(void)testAddFeatureToInfoVector001{
    NSMutableArray* features_info = [[NSMutableArray alloc]init];
    for (int i=0; i<12; ++i)
    {
        Feature* feat = [[Feature alloc]init];
        [features_info addObject:feat];
    }
    int X_RES_size = 91;
    int* uv = (int*)malloc(2*sizeof(int));
    uv[0] = 94;
    uv[1] = 117;
    int step = 92;
    double* newFeature = (double*)malloc(6*sizeof(double));
    newFeature[0] = 0.009139450566591;
    newFeature[1] = 0.001725790855408;
    newFeature[2] = 0.000681597720107;
    newFeature[3] = -0.330804250616683;
    newFeature[4] = 0.066217626541015;
    newFeature[5] = 1.000000000000000;
    
    NSString* path = [[NSBundle mainBundle] pathForResource:@"aftiv-001_im_k"
                                                     ofType:@"txt"];
    FILE* fp_im = fopen((char*)[path UTF8String], "r");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"aftiv-001_x"
                                           ofType:@"txt"];
    FILE* fp_x = fopen((char*)[path UTF8String], "r");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"aftiv-001_pwi"
                                           ofType:@"txt"];
    FILE* fp_pwi = fopen((char*)[path UTF8String], "r");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"aftiv-001_pwm"
                                           ofType:@"txt"];
    FILE* fp_pwm = fopen((char*)[path UTF8String], "r");
    path = nil;
    
    unsigned char* im_data = (unsigned char*)malloc(240*320*sizeof(unsigned char));
    unsigned char* patchWhenInit_expected = (unsigned char*)malloc(41*41*sizeof(unsigned char));
    unsigned char* patchWhenMatching_expected = (unsigned char*)malloc(13*13*sizeof(unsigned char));
    double* X_RES = (double*)malloc(X_RES_size*sizeof(double));
    
    for (int i=0; i<(240*320); ++i){
        double val;
        fscanf(fp_im, "%lf ", &val);
        im_data[i] = (unsigned char)((int)val);
    }
    fclose(fp_im);
    for (int i=0; i<(13*13); ++i){
        double val;
        fscanf(fp_pwm, "%lf ", &val);
        patchWhenMatching_expected[i] = (unsigned char)((int)val);
    }
    fclose(fp_pwm);
    for (int i=0; i<(41*41); ++i){
        double val;
        fscanf(fp_pwi, "%lf ", &val);
        patchWhenInit_expected[i] = (unsigned char)((int)val);
    }
    fclose(fp_pwi);
    for (int i=0; i<X_RES_size; ++i){
        fscanf(fp_x, "%lf ", X_RES+i);
    }
    fclose(fp_x);
    
    cv::Mat* im_k = new cv::Mat(240, 320, CV_8UC1, im_data, cv::Mat::AUTO_STEP);
    
    bool success = add_feature_to_info_vector(uv, im_k, X_RES, X_RES_size, features_info, step, newFeature);
    
    XCTAssert(success==true, @"init error");
    XCTAssert([features_info count]==13, @"init error");
    Feature* newFeat = (Feature*)features_info[12];
    XCTAssert(newFeat->patch_when_initialized, @"init error");
    XCTAssert(newFeat->patch_when_initialized->data, @"init error");
    XCTAssert(newFeat->patch_when_initialized->rows==41, @"init error");
    XCTAssert(newFeat->patch_when_initialized->cols==41, @"init error");
    unsigned char* mxPtr = newFeat->patch_when_initialized->ptr<uchar>(0);
    for (int i=0; i<(41*41); ++i){
        uchar actual = mxPtr[i];
        uchar expected = patchWhenInit_expected[i];
        if (actual!=expected){
            XCTAssert(NO, @"patchWhenInit has error data");
            break;
        }
    }
    XCTAssert(newFeat->patch_when_matching, @"init error");
    XCTAssert(newFeat->patch_when_matching->data, @"init error");
    XCTAssert(newFeat->patch_when_matching->rows==13, @"init error");
    XCTAssert(newFeat->patch_when_matching->cols==13, @"init error");
    mxPtr = newFeat->patch_when_matching->ptr<uchar>(0);
    for (int i=0; i<(13*13); ++i){
        uchar actual = mxPtr[i];
        uchar expected = patchWhenMatching_expected[i];
        if (actual!=expected){
            XCTAssert(NO, @"patchWhenMatching has error data");
            break;
        }
    }
    XCTAssert(newFeat->r_wc_when_initialized!=NULL, @"init error");
    double actual = newFeat->r_wc_when_initialized[0];
    double expected = 0.009139450566591;
    XCTAssert(fabs(actual-expected)<0.00001, @"init error");
    actual = newFeat->r_wc_when_initialized[1];
    expected = 0.001725790855408;
    XCTAssert(fabs(actual-expected)<0.00001, @"init error");
    actual = newFeat->r_wc_when_initialized[2];
    expected = 0.000681597720107;
    XCTAssert(fabs(actual-expected)<0.00001, @"init error");
    //R_wc
    XCTAssert(newFeat->R_wc_when_initialized!=NULL, @"init error");
    if (newFeat->R_wc_when_initialized==NULL)
        return;
    actual = newFeat->R_wc_when_initialized[0];
    expected = 0.999945094340855;
    XCTAssert(fabs(actual-expected)<0.0001, @"init error");
    actual = newFeat->R_wc_when_initialized[1];
    expected = -0.000944152863288;
    XCTAssert(fabs(actual-expected)<0.0001, @"init error");
    actual = newFeat->R_wc_when_initialized[2];
    expected = 0.010436324977187;
    XCTAssert(fabs(actual-expected)<0.0001, @"init error");
    actual = newFeat->R_wc_when_initialized[3];
    expected = 0.001010630421979;
    XCTAssert(fabs(actual-expected)<0.0001, @"init error");
    actual = newFeat->R_wc_when_initialized[4];
    expected = 0.999979223654465;
    XCTAssert(fabs(actual-expected)<0.0001, @"init error");
    actual = newFeat->R_wc_when_initialized[5];
    expected = -0.006366387167216;
    XCTAssert(fabs(actual-expected)<0.0001, @"init error");
    actual = newFeat->R_wc_when_initialized[6];
    expected = -0.010430097305820;
    XCTAssert(fabs(actual-expected)<0.0001, @"init error");
    actual = newFeat->R_wc_when_initialized[7];
    expected = 0.006376584884048 ;
    XCTAssert(fabs(actual-expected)<0.0001, @"init error");
    actual = newFeat->R_wc_when_initialized[8];
    expected = 0.999925273325666;
    XCTAssert(fabs(actual-expected)<0.0001, @"init error");
    XCTAssert(newFeat->uv_when_initialized[0]==94, @"init error");
    XCTAssert(newFeat->uv_when_initialized[1]==117, @"init error");
    XCTAssert(newFeat->half_patch_size_when_initialized==20, @"init error");
    XCTAssert(newFeat->half_patch_size_when_matching==6, @"init error");
    XCTAssert(newFeat->times_measured==0, @"init error");
    XCTAssert(newFeat->times_predicted==0, @"init error");
    XCTAssert(newFeat->init_frame==92, @"init error");
    XCTAssert(newFeat->init_measurement[0]==94, @"init error");
    XCTAssert(newFeat->init_measurement[1]==117, @"init error");
    XCTAssert([newFeat->type isEqualToString:@"inversedepth"], @"init error");
    actual = newFeat->yi[0];
    expected = newFeature[0];
    XCTAssert(fabs(actual-expected)<0.00001, @"init error");
    actual = newFeat->yi[1];
    expected = newFeature[1];
    XCTAssert(fabs(actual-expected)<0.00001, @"init error");
    actual = newFeat->yi[2];
    expected = newFeature[2];
    XCTAssert(fabs(actual-expected)<0.00001, @"init error");
    actual = newFeat->yi[3];
    expected = newFeature[3];
    XCTAssert(fabs(actual-expected)<0.00001, @"init error");
    actual = newFeat->yi[4];
    expected = newFeature[4];
    XCTAssert(fabs(actual-expected)<0.00001, @"init error");
    actual = newFeat->yi[5];
    expected = newFeature[5];
    XCTAssert(fabs(actual-expected)<0.00001, @"init error");
    XCTAssert(newFeat->individually_compatible==0, @"init error");
    XCTAssert(newFeat->low_innovation_inlier==0, @"init error");
    XCTAssert(newFeat->high_innovation_inlier==0, @"init error");
    XCTAssert(newFeat->z==NULL, @"init error");
    XCTAssert(newFeat->h==NULL, @"init error");
    XCTAssert(newFeat->H==NULL, @"init error");
    XCTAssert(newFeat->S==NULL, @"init error");
    XCTAssert(newFeat->state_size==6, @"init error");
    XCTAssert(newFeat->measurement_size==2, @"init error");
    XCTAssert(newFeat->R[0][0]==1, @"init error");
    XCTAssert(newFeat->R[1][0]==0, @"init error");
    XCTAssert(newFeat->R[0][1]==0, @"init error");
    XCTAssert(newFeat->R[1][1]==1, @"init error");
    
    for (int i=0; i<[features_info count]; ++i){
        [(Feature*)features_info[i] destroy];
    }
    
    free(uv);
    free(newFeature);
    free(im_data);
    free(patchWhenInit_expected);
    free(patchWhenMatching_expected);
    im_k->release();
}

/*********************************************UNDISTORT_FM************************************************/
//001: Pgm inputs
-(void)testUndistortFM001{
    double* uvd = (double*)malloc(2*sizeof(double));
    double* uvu_actual = (double*)malloc(2*sizeof(double));
    double* uvu_expected = (double*)malloc(2*sizeof(double));
    Camera* cam = [[Camera alloc]init];
    cam->k1 = 0.06333;
    cam->k2 = 0.0139;
    cam->nRows = 240;
    cam->nCols = 320;
    cam->Cx = 160.2232142857143;
    cam->Cy = 128.8660714285714;
    cam->f = 2.1735;
    cam->dx = 0.0112;
    cam->dy = 0.0112;
    uvd[0] = 163;
    uvd[1] = 30;
    uvu_expected[0] = 100*1.632739038288059;
    uvu_expected[1] = 100*0.202477906869224;
    bool success = undistortFM(uvd, 1, cam, uvu_actual);
    XCTAssert(success==true, @"Should return true for valid data");
    double expected = uvu_expected[0];
    double actual = uvu_actual[0];
    double diff = fabs(actual-expected);
    XCTAssert(diff<0.00001, @"uvu error exceeds limits");
    expected = uvu_expected[1];
    actual = uvu_actual[1];
    diff = fabs(actual-expected);
    XCTAssert(diff<0.00001, @"uvu error exceeds limits");
    free(uvd);
    free(uvu_expected);
    free(uvu_actual);
    [cam destroy];
}


//call from pred_patch_fm
-(void)testUndistortFM002{
    Camera* cam = [[Camera alloc]init];
    cam->k1 = 0.063330000000000;
    cam->k2 = 0.013900000000000;
    cam->nRows = 240;
    cam->nCols = 320;
    cam->Cx = 100*1.602232142857143;
    cam->Cy = 100*1.288660714285714;
    cam->dx = 0.011200000000000;
    cam->dy = 0.011200000000000;
    cam->f = 2.1735000000;
    cam->K = (double*)malloc(9*sizeof(double));
    cam->K[0] = 100*1.940625000000000;
    cam->K[1] = 100*0;
    cam->K[2] = 100*1.602232142857143;
    cam->K[3] = 100*0;
    cam->K[4] = 100*1.940625000000000;
    cam->K[5] = 100*1.288660714285715;
    cam->K[6] = 100*0;
    cam->K[7] = 100*0;
    cam->K[8] = 100*0.010000000000000;
    
    int cols = 169;
    NSString* path = [[NSBundle mainBundle] pathForResource:@"dfm-001_uvd" ofType:@"txt"];
    FILE* fp_uvd = fopen((char*)[path UTF8String], "r");
    path = nil;
    double* uvd = (double*)malloc(cols*2*sizeof(double));
    for (int i=0; i<(cols*2); ++i){
        fscanf(fp_uvd, "%lf ", uvd+i);
    }
    fclose(fp_uvd);
    path = [[NSBundle mainBundle] pathForResource:@"dfm-001_uvu" ofType:@"txt"];
    FILE* fp_uvu = fopen((char*)[path UTF8String], "r");
    path = nil;
    double* uvu_expected = (double*)malloc(2*cols*sizeof(double));
    double* uvu_actual = (double*)malloc(2*cols*sizeof(double));
    for (int i=0; i<(2*cols); ++i){
        fscanf(fp_uvu, "%lf ", uvu_expected+i);
    }
    fclose(fp_uvu);
    
    bool success =  undistortFM(uvd, cols, cam, uvu_actual);
    XCTAssert(success==true, @"expected return value of true");
    for (int i=0; i<2*cols; ++i){
        double actual = uvu_actual[i];
        double expected = uvu_expected[i];
        double diff = abs(actual-expected);
        if (diff>0.00001){
            XCTAssert(NO, @"uvu error exceeds limits");
            break;
        }
    }
    
    free(uvu_expected);
    free(uvu_actual);
    free(uvd);
    
}

//call from rotate_with_dist
-(void)testUndistortFM003{
    Camera* cam = [[Camera alloc]init];
    cam->k1 = 0.063330000000000;
    cam->k2 = 0.013900000000000;
    cam->nRows = 240;
    cam->nCols = 320;
    cam->Cx = 100*1.602232142857143;
    cam->Cy = 100*1.288660714285714;
    cam->dx = 0.011200000000000;
    cam->dy = 0.011200000000000;
    cam->f = 2.1735000000;
    cam->K = (double*)malloc(9*sizeof(double));
    cam->K[0] = 100*1.940625000000000;
    cam->K[1] = 100*0;
    cam->K[2] = 100*1.602232142857143;
    cam->K[3] = 100*0;
    cam->K[4] = 100*1.940625000000000;
    cam->K[5] = 100*1.288660714285715;
    cam->K[6] = 100*0;
    cam->K[7] = 100*0;
    cam->K[8] = 100*0.010000000000000;
    
    int cols = 1;
    double* uvd = (double*)malloc(cols*2*sizeof(double));
    uvd[0] = 53;
    uvd[1] = 127;

    double* uvu_expected = (double*)malloc(2*cols*sizeof(double));
    double* uvu_actual = (double*)malloc(2*cols*sizeof(double));
    uvu_expected[0] = 100*0.401024595594034;
    uvu_expected[1] = 100*1.267755361851874;
    uvu_actual[0] = -99;
    uvu_actual[1] = -99;
    
    bool success =  undistortFM(uvd, cols, cam, uvu_actual);
    XCTAssert(success==true, @"expected return value of true");
    for (int i=0; i<2*cols; ++i){
        double actual = uvu_actual[i];
        double expected = uvu_expected[i];
        double diff = abs(actual-expected);
        if (diff>0.00001){
            XCTAssert(NO, @"uvu error exceeds limits");
            break;
        }
    }
    
    free(uvu_expected);
    free(uvu_actual);
    free(uvd);
    
}

/******************************************dRqTimesABydq*************************************************/
//Test for various invalid inputs

-(void)testdRqTimesAByDq001{
    double* q = (double*)malloc(4*sizeof(double));
    cv::Mat* aMat = new cv::Mat;
    *aMat = cv::Mat(3,1,CV_64FC1,0);
    cv::Mat* dRqTimesABydqRES = new cv::Mat;
    *dRqTimesABydqRES = cv::Mat(3,1,CV_64FC1,0.0);
    cv::Mat* dummyMatrix = NULL;
    //Should return false for NULL Matrix:
    XCTAssert(!dummyMatrix, @"Problems...");
    bool success = dRqTimesABydq(q, dummyMatrix, dRqTimesABydqRES);
    XCTAssert(success==false, @"Should return false for unallocated Matrix (!matrix)");
    success = dRqTimesABydq(q, dRqTimesABydqRES, dummyMatrix);
    XCTAssert(success==false, @"Should return false for unallocated Matrix (!matrix)");
    //Should return false if matrix exists but internal data wasn't allocated:
    XCTAssert(!aMat->data, @"Data works?");
    success = dRqTimesABydq(q, aMat, dRqTimesABydqRES);
    XCTAssert(success==false, @"Should return false for unallocated Matrix data (!matrix->data)");
    *aMat = cv::Mat(3,1,CV_64FC1,0.0);
     success = dRqTimesABydq(q, aMat, dRqTimesABydqRES);
    XCTAssert(success==true, @"Should return true for valid data");
    aMat->release();
    free(q);
    dRqTimesABydqRES->release();
}

//Input from program
-(void)testdRqTimesAByDq002{
    double* q = (double*)malloc(3*sizeof(double));
    double* aMatData = (double*)malloc(3*sizeof(double));
    double* RES_expected = (double*)malloc(12*sizeof(double));
    aMatData[0] = -0.401821923129474;
    aMatData[1] = -0.325073879355352;
    aMatData[2] = 0.870710892613979;
    q[0] = 0.999521087165285;
    q[1] = -0.005511800542970;
    q[2] = -0.028417909460399;
    q[3] = -0.010937951755297;
    cv::Mat* aMat = new cv::Mat(3, 1, CV_64FC1, aMatData, cv::Mat::AUTO_STEP);
    RES_expected[0] = -0.859857822350002;
    RES_expected[1] = 0.003857777258639;
    RES_expected[2] = 1.721333402692439;
    RES_expected[3] =  0.631447807444415;
    RES_expected[4] = -0.631447807444415;
    RES_expected[5] = -1.721333402692439;
    RES_expected[6] = 0.003857777258639;
    RES_expected[7] = -0.859857822350002;
    RES_expected[8] = 1.721333402692439;
    RES_expected[9] = -0.631447807444415;
    RES_expected[10] = 0.859857822350002;
    RES_expected[11] = 0.003857777258639;
    cv::Mat* RES_actual = new cv::Mat(3,4,CV_64FC1,0.0);
    bool success = dRqTimesABydq(q, aMat, RES_actual);
    XCTAssert(success==true, @"Expected return value of true for valid data");
    double* mxPtr = RES_actual->ptr<double>(0);
    for (int i=0; i<12; ++i){
        double actual = mxPtr[i];
        double expected = RES_expected[i];
        double diff = fabs(actual-expected);
        if (diff>0.00001){
            XCTAssert(NO, @"dRq_times_a_by_dqRES matrix data error exceeds limits");
            break;
        }
    }
    RES_actual->release();
    aMat->release();
    free(RES_expected);
    free(q);
    free(aMatData);
}

-(void)testdRqTimesAByDq003{
    double* q = (double*)malloc(4*sizeof(double));
    q[0] = 1.000000000000000;
    q[1] = -0.000000000000001;
    q[2] = -0.000000000000001;
    q[3] = -0.000000000000001;
    double* aMat = (double*)malloc(3*sizeof(double));
    aMat[0] = 0.573044609938538;
    aMat[1] = -0.452109735526377;
    aMat[2] = 0.683532487935034;
    double* ans_expected = (double*)malloc(12*sizeof(double));
    double* ans_actual = (double*)malloc(12*sizeof(double));
    ans_expected[0] = 1.146089219877076;
    ans_expected[1] = -0.000000000000001;
    ans_expected[2] = 1.367064975870070;
    ans_expected[3] = 0.904219471052753;
    ans_expected[4] = -0.904219471052753;
    ans_expected[5] = -1.367064975870070;
    ans_expected[6] = -0.000000000000001;
    ans_expected[7] = 1.146089219877076;
    ans_expected[8] = 1.367064975870070;
    ans_expected[9] = -0.904219471052753;
    ans_expected[10] = -1.146089219877076;
    ans_expected[11] = -0.000000000000001;
    
    cv::Mat* aMat_mat = new cv::Mat(3, 1, CV_64FC1, aMat, cv::Mat::AUTO_STEP);
    cv::Mat* ans_actual_mat = new cv::Mat(3, 4, CV_64FC1, ans_actual, cv::Mat::AUTO_STEP);
    bool success = dRqTimesABydq(q, aMat_mat, ans_actual_mat);
    XCTAssert(success==true, @"expected true for valid data");
    double* ansPtr = ans_actual_mat->ptr<double>(0);
    for (int i=0; i<12; ++i){
        double actual = ansPtr[i];
        double expected = ans_expected[i];
        double diff = fabs(actual-expected);
        if (diff>0.00001){
            XCTAssert(NO, @"dRq error exceeds limits");
            break;
        }
    }
    free(q);
    free(aMat);
    free(ans_expected);
    free(ans_actual);
}

/*****************************************JACOB_UNDISTOR_FM**********************************************/
//TODO: Invalid input tests
//001: Program inputs:
-(void)testJacobUndistorFm001{
    Camera* cam = [[Camera alloc]init];
    cam->k1 = 0.06333;
    cam->k2 = 0.0139;
    cam->nRows = 240;
    cam->nCols = 320;
    cam->Cx = 160.2232142857143;
    cam->Cy = 128.8660714285714;
    cam->f = 2.1735;
    cam->dx = 0.0112;
    cam->dy = 0.0112;
    double* uvd_data = (double*)malloc(2*sizeof(double));
    double* jUndistor_expected = (double*)malloc(4*sizeof(double));
    cv::Mat* jUndistor_mat;
    cv::Mat* uvd_mat;
    uvd_data[0] = 100*0.612147429734590;
    uvd_data[1] = 100*1.370766495277684;
    jUndistor_expected[0] = 1.340110904688409;
    jUndistor_expected[1] = -0.019935437716640;
    jUndistor_expected[2] = -0.019935437716640;
    jUndistor_expected[3] = 1.101369684905456;
    jUndistor_mat = new cv::Mat(2,2,CV_64FC1,0.0);
    JacobUndistorFM(cam, uvd_data, jUndistor_mat);
    double* mxPtr = jUndistor_mat->ptr<double>(0);
    for (int i=0; i<4; ++i){
        double expected = jUndistor_expected[i];
        double actual = mxPtr[i];
        double diff = fabs(actual-expected);
        if (diff>0.00001){
            XCTAssert(NO, @"jUndistor error exceeds limits");
            break;
        }
    }
    free(uvd_data);
    free(jUndistor_expected);
    jUndistor_mat->release();
}

-(void)testJacobUndistorFm002{
    Camera* cam = [[Camera alloc]init];
    cam->k1 = 0.063330000000000;
    cam->k2 = 0.013900000000000;
    cam->nRows = 240;
    cam->nCols = 320;
    cam->Cx = 100*1.602232142857143;
    cam->Cy = 100*1.288660714285714;
    cam->dx = 0.011200000000000;
    cam->dy = 0.011200000000000;
    double* uvd_data = (double*)malloc(2*sizeof(double));
    double* jUndistor_expected = (double*)malloc(4*sizeof(double));
    cv::Mat* jUndistor_mat;
    cv::Mat* uvd_mat;
    uvd_data[0] = 119;
    uvd_data[1] = 197;
    jUndistor_expected[0] = 1.095602144760301;
    jUndistor_expected[1] = -0.060208227870147;
    jUndistor_expected[2] = -0.060208227870147;
    jUndistor_expected[3] = 1.158686684992629;
    jUndistor_mat = new cv::Mat(2,2,CV_64FC1,0.0);
    JacobUndistorFM(cam, uvd_data, jUndistor_mat);
    double* mxPtr = jUndistor_mat->ptr<double>(0);
    for (int i=0; i<4; ++i){
        double expected = jUndistor_expected[i];
        double actual = mxPtr[i];
        double diff = fabs(actual-expected);
        if (diff>0.00001){
            XCTAssert(NO, @"jUndistor error exceeds limits");
            //break;
        }
    }
    free(uvd_data);
    free(jUndistor_expected);
    jUndistor_mat->release();
}

-(void)testJacobUndistorFm003{
    Camera* cam = [[Camera alloc]init];
    cam->k1 = 0.063330000000000;
    cam->k2 = 0.013900000000000;
    cam->nRows = 240;
    cam->nCols = 320;
    cam->Cx = 100*1.602232142857143;
    cam->Cy = 100*1.288660714285714;
    cam->dx = 0.011200000000000;
    cam->dy = 0.011200000000000;
    double* uvd_data = (double*)malloc(2*sizeof(double));
    double* jUndistor_expected = (double*)malloc(4*sizeof(double));
    cv::Mat* jUndistor_mat;
    cv::Mat* uvd_mat;
    uvd_data[0] = 82.999999999999744;
    uvd_data[1] = 62.000000000000256;
    jUndistor_expected[0] = 1.255894471960513;
    jUndistor_expected[1] = 0.129178770487826;
    jUndistor_expected[2] = 0.129178770487826;
    jUndistor_expected[3] = 1.218560076849706;
    
    
    jUndistor_mat = new cv::Mat(2,2,CV_64FC1,0.0);
    JacobUndistorFM(cam, uvd_data, jUndistor_mat);
    double* mxPtr = jUndistor_mat->ptr<double>(0);
    for (int i=0; i<4; ++i){
        double expected = jUndistor_expected[i];
        double actual = mxPtr[i];
        double diff = fabs(actual-expected);
        if (diff>0.00001){
            XCTAssert(NO, @"jUndistor error exceeds limits");
            break;
        }
    }
    free(uvd_data);
    free(jUndistor_expected);
    jUndistor_mat->release();
}


/***********************************************Q2R******************************************************/
//Test for null inputs
-(void)testQ2R001{
    double* q = (double*)malloc(4*sizeof(double));
    double* R = (double*)malloc(9*sizeof(double));
    bool success = q2r(NULL, R);
    XCTAssert(success==false, @"Should return false for invalid data");
    success = q2r(q, NULL);
    XCTAssert(success==false, @"Should return false for invalid data");
    free(q);
    free(R);
}
//Test against pgm inputs
-(void)testQ2R002{
    double* q = (double*)malloc(4*sizeof(double));
    double* R = (double*)malloc(9*sizeof(double));
    double* R_act = (double*)malloc(9*sizeof(double));
    q[0] = 0.999703781918051;
    q[1] = 0.004106594044519;
    q[2] = 0.022265672336388;
    q[3] = 0.008928837522818;
    R[0] = 0.998849031391799;
    R[1] = -0.017669513124558;
    R[2] = 0.044591487905261;
    R[3] = 0.018035257434214;
    R[4] = 0.999806823491689;
    R[5] = -0.007813142052961;
    R[6] = -0.044444819461279;
    R[7] = 0.008608368335472;
    R[8] = 0.998974751441524;
    
    bool success = q2r(q, R_act);
    XCTAssert(success==true, @"Should return false for invalid data");
    
    for (int i=0; i<9; ++i){
        double actual = R_act[i];
        double expected = R[i];
        double diff = fabs(actual-expected);
        if (diff>0.00001){
            XCTAssert(NO, @"R error exceeds limits");
            break;
        }
    }
    free(q);
    free(R);
    free(R_act);
}

/*****************************************NORMJAC**************************************************/
//001: Testing for invalid inputs
-(void)testNormJac001{
    double* q = (double*)malloc(4*sizeof(double));
    double* jac = (double*)malloc(4*4*sizeof(double));
    bool success = normJac(NULL, q);
    XCTAssert(success==false, @"Should return false for invalid data");
    success = normJac(q,NULL);
    XCTAssert(success==false, @"Should return false for invalid data");
    free(q);
    free(jac);
}
//002: Data from program
-(void)testNormJac002{
    double* q = (double*)malloc(4*sizeof(double));
    double* jac = (double*)malloc(4*4*sizeof(double));
    double* jac_expected = (double*)malloc(4*4*sizeof(double));
    *q = 1.000000000000000;
    *(q+1) = 0.003192918697552;
    *(q+2) = 0.005333653154188;
    *(q+3) = 0.000596290674040;
    NSString* path = [[NSBundle mainBundle] pathForResource:@"nj-002_j"
                                                     ofType:@"txt"];
    FILE* fp_j = fopen((char*)[path UTF8String], "r");
    XCTAssert(fp_j!=NULL && path!=nil, @"Setup problem");
    for (int i=0; i<16; ++i){
        fscanf(fp_j, "%lf ", (jac_expected+i));
    }
    bool success = normJac(jac, q);
    XCTAssert(success==true, @"Should return true for valid data");
    
    for(int i=0; i<16; ++i){
        double expected = jac_expected[i];
        double actual = jac[i];
        double diff = fabs(expected-actual);
        if (diff>0.0001){
            XCTAssert(diff<0.0001, @"jacobian error exceeds limits");
            break;
        }
    }
    
    free(q);
    free(jac);
}


@end

