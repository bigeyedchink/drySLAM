//
//  EKFUpdateHiLnliersTests.m
//  drySLAM
//
//  Created by Anthony Rodriguez on 3/27/15.
//  Copyright (c) 2015 ECE473 573 Common. All rights reserved.
//


#include "SLAMModel.hpp"
#include "EKFUpdateLiInliers.hpp"

#import <UIKit/UIKit.h>
#import <XCTest/XCTest.h>
#include <stdio.h>
#include <stdlib.h>




@interface EKFUpdateLiInliers : XCTestCase

@end

@implementation EKFUpdateLiInliers


- (void)setUp {
    [super setUp];
    // Put setup code here. This method is called before the invocation of each test method in the class.
}

- (void)tearDown {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
    [super tearDown];
}

/***************************************UPDATE*****************************************/
//001: Test for NULL inputs
-(void)testUpdate001{
    double* x_vec = (double*)malloc(85*sizeof(double));
    double* p_vec = (double*)malloc(85*85*sizeof(double));
    int state_size = 85;
    double* H_data = (double*)malloc(2*85*sizeof(double));
    double* R_data = (double*)malloc(2*2*sizeof(double));
    double* z_data = (double*)malloc(2*sizeof(double));
    double* h_data = (double*)malloc(2*sizeof(double));
    
    cv::Mat H(2, 85, CV_64FC1, H_data, cv::Mat::AUTO_STEP);
    cv::Mat R(2, 2, CV_64FC1, R_data, cv::Mat::AUTO_STEP);
    cv::Mat z(2, 1, CV_64FC1, z_data, cv::Mat::AUTO_STEP);
    cv::Mat h(2, 1, CV_64FC1, h_data, cv::Mat::AUTO_STEP);
    
    
    bool success = update(NULL, p_vec, state_size, &H, &R, &z, &h);
    XCTAssert(success==false, @"Should return false for null pointer arguments");
    success = update(x_vec, NULL, state_size, &H, &R, &z, &h);
    XCTAssert(success==false, @"Should return false for null pointer arguments");
    success = update(x_vec, p_vec, -1, &H, &R, &z, &h);
    XCTAssert(success==false, @"Should return false for invalid arguments");
    success = update(x_vec, p_vec, state_size, NULL, &R, &z, &h);
    XCTAssert(success==false, @"Should return false for null pointer arguments");
    success = update(x_vec, p_vec, state_size, &H, NULL, &z, &h);
    XCTAssert(success==false, @"Should return false for null pointer arguments");
    success = update(x_vec, p_vec, state_size, &H, &R, NULL, &h);
    XCTAssert(success==false, @"Should return false for null pointer arguments");
    success = update(x_vec, p_vec, state_size, &H, &R, &z, NULL);
    XCTAssert(success==false, @"Should return false for null pointer arguments");
    H.release();
    R.release();
    z.release();
    h.release();
    free(x_vec);
    free(p_vec);
    free(H_data);
    free(R_data);
    free(z_data);
    free(h_data);
}
//002: Test against input from MATLAB
-(void)testUpdate002{
    int state_size = 127;
    double* x_vec = (double*)malloc(state_size*sizeof(double));
    double* p_vec = (double*)malloc(state_size*state_size*sizeof(double));
    double* x_vec_expected = (double*)malloc(state_size*sizeof(double));
    double* p_vec_expected = (double*)malloc(state_size*state_size*sizeof(double));
    double* H_data = (double*)malloc(34*state_size*sizeof(double));
    double* R_data = (double*)malloc(34*34*sizeof(double));
    double* z_data = (double*)malloc(34*1*sizeof(double));
    double* h_data = (double*)malloc(34*1*sizeof(double));
    
    NSString* path = [[NSBundle mainBundle] pathForResource:@"u-002_x_before"
                                                     ofType:@"txt"];
    FILE* fp_x_before = fopen((char*)[path UTF8String], "r");
    XCTAssert(fp_x_before!=NULL && path!=nil, @"Setup problem");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"u-002_p_before"
                                           ofType:@"txt"];
    FILE* fp_p_before = fopen((char*)[path UTF8String], "r");
    XCTAssert(fp_p_before!=NULL && path!=nil, @"Setup problem");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"u-002_R"
                                           ofType:@"txt"];
    FILE* fp_R = fopen((char*)[path UTF8String], "r");
    XCTAssert(fp_R!=NULL && path!=nil, @"Setup problem");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"u-002_h_upper"
                                           ofType:@"txt"];
    FILE* fp_H = fopen((char*)[path UTF8String], "r");
    XCTAssert(fp_H!=NULL && path!=nil, @"Setup problem");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"u-002_x_after"
                                           ofType:@"txt"];
    FILE* fp_x_after = fopen((char*)[path UTF8String], "r");
    XCTAssert(fp_x_after!=NULL && path!=nil, @"Setup problem");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"u-002_p_after"
                                           ofType:@"txt"];
    FILE* fp_p_after = fopen((char*)[path UTF8String], "r");
    XCTAssert(fp_p_after!=NULL && path!=nil, @"Setup problem");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"u-002_z"
                                           ofType:@"txt"];
    FILE* fp_z = fopen((char*)[path UTF8String], "r");
    XCTAssert(fp_z!=NULL && path!=nil, @"Setup problem");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"u-002_h_lower"
                                           ofType:@"txt"];
    FILE* fp_h_lower = fopen((char*)[path UTF8String], "r");
    XCTAssert(fp_h_lower!=NULL && path!=nil, @"Setup problem");
    path = nil;
    
    for (int i=0; i<state_size; ++i){
        fscanf(fp_x_before, "%lf ", (x_vec+i));
        fscanf(fp_x_after, "%lf ", (x_vec_expected+i));
    }
    fclose(fp_x_before);
    fclose(fp_x_after);
    for (int i=0; i<(state_size*state_size); ++i){
        fscanf(fp_p_before, "%lf ", (p_vec+i));
        fscanf(fp_p_after, "%lf ", (p_vec_expected+i));
    }
    fclose(fp_p_before);
    fclose(fp_p_after);
    for (int i=0; i<(34*state_size); ++i){
        fscanf(fp_H, "%lf ", (H_data+i));
    }
    fclose(fp_H);
    for (int i=0; i<(34*34); ++i){
        fscanf(fp_R, "%lf ", (R_data+i));
    }
    fclose(fp_R);
    for (int i=0; i<(34); ++i){
        fscanf(fp_z, "%lf ", (z_data+i));
    }
    fclose(fp_z);
    for (int i=0; i<(34); ++i){
        fscanf(fp_h_lower, "%lf ", (h_data+i));
    }
    fclose(fp_h_lower);
    
    cv::Mat H(34, 127, CV_64FC1, H_data, cv::Mat::AUTO_STEP);
    cv::Mat R(34, 34, CV_64FC1, R_data, cv::Mat::AUTO_STEP);
    cv::Mat z(34, 1, CV_64FC1, z_data, cv::Mat::AUTO_STEP);
    cv::Mat h(34, 1, CV_64FC1, h_data, cv::Mat::AUTO_STEP);
    
    bool success = update(x_vec, p_vec, state_size, &H, &R, &z, &h);
    
    XCTAssert(success==true, @"Expected return value of true for valid data");
    
    for (int i=0; i<state_size; ++i){
        double diff = fabs(x_vec[i]-x_vec_expected[i]);
        if (diff>0.00001){
            XCTAssert(diff<0.00001, @"x error exceeds limits");
            break;
        }
    }
    for (int i=0; i<(state_size*state_size); ++i){
        double expected = p_vec_expected[i];
        double actual = p_vec[i];
        double diff = fabs(expected-actual);
        if (diff>0.0001){
            XCTAssert(diff<0.0001, @"p error exceeds limits");
            break;
        }
    }
    
    H.release();
    R.release();
    z.release();
    h.release();
    free(x_vec_expected);
    free(x_vec);
    free(p_vec);
    free(p_vec_expected);
    free(z_data);
    free(R_data);
    free(H_data);
    free(h_data);
}
//003: Test against input from MATLAB
-(void)testUpdate003{
    int state_size = 85;
    double* x_vec = (double*)malloc(state_size*sizeof(double));
    double* p_vec = (double*)malloc(state_size*state_size*sizeof(double));
    double* x_vec_expected = (double*)malloc(state_size*sizeof(double));
    double* p_vec_expected = (double*)malloc(state_size*state_size*sizeof(double));
    double* H_data = (double*)malloc(22*state_size*sizeof(double));
    double* R_data = (double*)malloc(22*22*sizeof(double));
    double* z_data = (double*)malloc(22*1*sizeof(double));
    double* h_data = (double*)malloc(22*1*sizeof(double));
    
    NSString* path = [[NSBundle mainBundle] pathForResource:@"u-003_x_before"
                                                     ofType:@"txt"];
    FILE* fp_x_before = fopen((char*)[path UTF8String], "r");
    XCTAssert(fp_x_before!=NULL && path!=nil, @"Setup problem");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"u-003_p_before"
                                           ofType:@"txt"];
    FILE* fp_p_before = fopen((char*)[path UTF8String], "r");
    XCTAssert(fp_p_before!=NULL && path!=nil, @"Setup problem");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"u-003_R"
                                           ofType:@"txt"];
    FILE* fp_R = fopen((char*)[path UTF8String], "r");
    XCTAssert(fp_R!=NULL && path!=nil, @"Setup problem");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"u-003_h_upper"
                                           ofType:@"txt"];
    FILE* fp_H = fopen((char*)[path UTF8String], "r");
    XCTAssert(fp_H!=NULL && path!=nil, @"Setup problem");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"u-003_x_after"
                                           ofType:@"txt"];
    FILE* fp_x_after = fopen((char*)[path UTF8String], "r");
    XCTAssert(fp_x_after!=NULL && path!=nil, @"Setup problem");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"u-003_p_after"
                                           ofType:@"txt"];
    FILE* fp_p_after = fopen((char*)[path UTF8String], "r");
    XCTAssert(fp_p_after!=NULL && path!=nil, @"Setup problem");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"u-003_z"
                                           ofType:@"txt"];
    FILE* fp_z = fopen((char*)[path UTF8String], "r");
    XCTAssert(fp_z!=NULL && path!=nil, @"Setup problem");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"u-003_h_lower"
                                           ofType:@"txt"];
    FILE* fp_h_lower = fopen((char*)[path UTF8String], "r");
    XCTAssert(fp_h_lower!=NULL && path!=nil, @"Setup problem");
    path = nil;
    
    for (int i=0; i<state_size; ++i){
        fscanf(fp_x_before, "%lf ", (x_vec+i));
        fscanf(fp_x_after, "%lf ", (x_vec_expected+i));
    }
    fclose(fp_x_before);
    fclose(fp_x_after);
    for (int i=0; i<(state_size*state_size); ++i){
        fscanf(fp_p_before, "%lf ", (p_vec+i));
        fscanf(fp_p_after, "%lf ", (p_vec_expected+i));
    }
    fclose(fp_p_before);
    fclose(fp_p_after);
    for (int i=0; i<(22*state_size); ++i){
        fscanf(fp_H, "%lf ", (H_data+i));
    }
    fclose(fp_H);
    for (int i=0; i<(22*22); ++i){
        fscanf(fp_R, "%lf ", (R_data+i));
    }
    fclose(fp_R);
    for (int i=0; i<(22); ++i){
        fscanf(fp_z, "%lf ", (z_data+i));
    }
    fclose(fp_z);
    for (int i=0; i<(22); ++i){
        fscanf(fp_h_lower, "%lf ", (h_data+i));
    }
    fclose(fp_h_lower);
    
    cv::Mat H(22, 85, CV_64FC1, H_data, cv::Mat::AUTO_STEP);
    cv::Mat R(22, 22, CV_64FC1, R_data, cv::Mat::AUTO_STEP);
    cv::Mat z(22, 1, CV_64FC1, z_data, cv::Mat::AUTO_STEP);
    cv::Mat h(22, 1, CV_64FC1, h_data, cv::Mat::AUTO_STEP);
    
    bool success = update(x_vec, p_vec, state_size, &H, &R, &z, &h);
    
    XCTAssert(success==true, @"Expected return value of true for valid data");
    
    for (int i=0; i<state_size; ++i){
        double diff = fabs(x_vec[i]-x_vec_expected[i]);
        if (diff>0.00001){
            XCTAssert(diff<0.00001, @"x error exceeds limits");
            break;
        }
    }
    for (int i=0; i<(state_size*state_size); ++i){
        double expected = p_vec_expected[i];
        double actual = p_vec[i];
        double diff = fabs(expected-actual);
        if (diff>0.0001){
            XCTAssert(diff<0.0001, @"p error exceeds limits");
            break;
        }
    }
    H.release();
    R.release();
    z.release();
    h.release();
    free(x_vec_expected);
    free(x_vec);
    free(p_vec);
    free(p_vec_expected);
    free(z_data);
    free(R_data);
    free(H_data);
    free(h_data);
}

/*****************************************EKF_UPDATE_LI_INLIERS*********************************************/
//001: input from MATLAB
-(void)testEKFUpdateLiInliers001{
    //Each H is 2xstate_size
    int state_size = 85;
    Feature* feat1 = [[Feature alloc]init];
    Feature* feat2 = [[Feature alloc]init];
    Feature* feat3 = [[Feature alloc]init];
    Feature* feat4 = [[Feature alloc]init];
    Feature* feat5 = [[Feature alloc]init];
    Feature* feat6 = [[Feature alloc]init];
    Feature* feat7 = [[Feature alloc]init];
    Feature* feat8 = [[Feature alloc]init];
    Feature* feat9 = [[Feature alloc]init];
    Feature* feat10 =[[Feature alloc]init];
    Feature* feat11= [[Feature alloc]init];
    Feature* feat12= [[Feature alloc]init];
    feat1->z = (int*)malloc(2*sizeof(int));
    feat2->z = (int*)malloc(2*sizeof(int));
    feat3->z = (int*)malloc(2*sizeof(int));
    feat4->z = (int*)malloc(2*sizeof(int));
    feat5->z = (int*)malloc(2*sizeof(int));
    feat6->z = (int*)malloc(2*sizeof(int));
    feat7->z = (int*)malloc(2*sizeof(int));
    feat8->z = (int*)malloc(2*sizeof(int));
    feat9->z = (int*)malloc(2*sizeof(int));
    feat10->z = (int*)malloc(2*sizeof(int));
    feat11->z = (int*)malloc(2*sizeof(int));
    feat12->z = (int*)malloc(2*sizeof(int));
    feat1->z[0] = 49;       feat1->z[1] = 128;
    feat2->z[0] = 89;       feat2->z[1] = 154;
    feat3->z[0] = 232;      feat3->z[1] = 142;
    feat4->z[0] = 29;       feat4->z[1] = 163;
    feat5->z[0] = 71;       feat5->z[1] = 44;
    feat6->z[0] = 80;       feat6->z[1] = 63;
    feat7->z[0] = 75;       feat7->z[1] = 132;
    feat8->z[0] = 124;      feat8->z[1] = 97;
    feat9->z[0] = 163;      feat9->z[1] = 103;
    feat10->z[0] = 280;     feat10->z[1] = 32;
    feat11->z[0] = 94;      feat11->z[1] = 56;
    feat12->z[0] = 72;      feat12->z[1] = 196;
    feat1->h = (double*)malloc(2*sizeof(double));
    feat2->h = (double*)malloc(2*sizeof(double));
    feat3->h = (double*)malloc(2*sizeof(double));
    feat4->h = (double*)malloc(2*sizeof(double));
    feat5->h = (double*)malloc(2*sizeof(double));
    feat6->h = (double*)malloc(2*sizeof(double));
    feat7->h = (double*)malloc(2*sizeof(double));
    feat8->h = (double*)malloc(2*sizeof(double));
    feat9->h = (double*)malloc(2*sizeof(double));
    feat10->h = (double*)malloc(2*sizeof(double));
    feat11->h = (double*)malloc(2*sizeof(double));
    feat12->h = (double*)malloc(2*sizeof(double));
    feat1->h[0] = 53;       feat1->h[1] = 127;
    feat2->h[0] = 93;       feat2->h[1] = 153;
    feat3->h[0] = 236;      feat3->h[1] = 142;
    feat4->h[0] = 32;       feat4->h[1] = 162;
    feat5->h[0] = 74;       feat5->h[1] = 43;
    feat6->h[0] = 83;       feat6->h[1] = 62;
    feat7->h[0] = 79;       feat7->h[1] = 131;
    feat8->h[0] = 128;      feat8->h[1] = 96;
    feat9->h[0] = 168;      feat9->h[1] = 102;
    feat10->h[0] = 283;     feat10->h[1] = 32;
    feat11->h[0] = 97;      feat11->h[1] = 55;
    feat12->h[0] = 76;      feat12->h[1] = 195;
    
    feat1->low_innovation_inlier = 1; feat2->low_innovation_inlier = 1; feat3->low_innovation_inlier = 1;
    feat4->low_innovation_inlier = 1; feat5->low_innovation_inlier = 1; feat6->low_innovation_inlier = 1;
    feat7->low_innovation_inlier = 1; feat8->low_innovation_inlier = 1; feat9->low_innovation_inlier = 0;
    feat10->low_innovation_inlier = 1; feat11->low_innovation_inlier = 1; feat12->low_innovation_inlier = 1;
    
    NSString* path = [[NSBundle mainBundle] pathForResource:@"euli-001_H1"
                                                     ofType:@"txt"];
    FILE* fp_h1 = fopen((char*)[path UTF8String], "r");
    XCTAssert(fp_h1!=NULL && path!=nil, @"Setup problem");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"euli-001_H2"
                                           ofType:@"txt"];
    FILE* fp_h2 = fopen((char*)[path UTF8String], "r");
    XCTAssert(fp_h2!=NULL && path!=nil, @"Setup problem");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"euli-001_H3"
                                           ofType:@"txt"];
    FILE* fp_h3 = fopen((char*)[path UTF8String], "r");
    XCTAssert(fp_h3!=NULL && path!=nil, @"Setup problem");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"euli-001_H4"
                                           ofType:@"txt"];
    FILE* fp_h4 = fopen((char*)[path UTF8String], "r");
    XCTAssert(fp_h4!=NULL && path!=nil, @"Setup problem");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"euli-001_H5"
                                           ofType:@"txt"];
    FILE* fp_h5 = fopen((char*)[path UTF8String], "r");
    XCTAssert(fp_h5!=NULL && path!=nil, @"Setup problem");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"euli-001_H6"
                                           ofType:@"txt"];
    FILE* fp_h6 = fopen((char*)[path UTF8String], "r");
    XCTAssert(fp_h6!=NULL && path!=nil, @"Setup problem");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"euli-001_H7"
                                           ofType:@"txt"];
    FILE* fp_h7 = fopen((char*)[path UTF8String], "r");
    XCTAssert(fp_h7!=NULL && path!=nil, @"Setup problem");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"euli-001_H8"
                                           ofType:@"txt"];
    FILE* fp_h8 = fopen((char*)[path UTF8String], "r");
    XCTAssert(fp_h8!=NULL && path!=nil, @"Setup problem");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"euli-001_H9"
                                           ofType:@"txt"];
    FILE* fp_h9 = fopen((char*)[path UTF8String], "r");
    XCTAssert(fp_h9!=NULL && path!=nil, @"Setup problem");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"euli-001_H10"
                                           ofType:@"txt"];
    FILE* fp_h10 = fopen((char*)[path UTF8String], "r");
    XCTAssert(fp_h10!=NULL && path!=nil, @"Setup problem");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"euli-001_H11"
                                           ofType:@"txt"];
    FILE* fp_h11 = fopen((char*)[path UTF8String], "r");
    XCTAssert(fp_h11!=NULL && path!=nil, @"Setup problem");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"euli-001_H12"
                                           ofType:@"txt"];
    FILE* fp_h12 = fopen((char*)[path UTF8String], "r");
    XCTAssert(fp_h12!=NULL && path!=nil, @"Setup problem");
    path = nil;
    
    path = [[NSBundle mainBundle] pathForResource:@"euli-001_x_before"
                                           ofType:@"txt"];
    FILE* fp_x_before = fopen((char*)[path UTF8String], "r");
    XCTAssert(fp_x_before!=NULL && path!=nil, @"Setup problem");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"euli-001_x_after"
                                           ofType:@"txt"];
    FILE* fp_x_after = fopen((char*)[path UTF8String], "r");
    XCTAssert(fp_x_after!=NULL && path!=nil, @"Setup problem");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"euli-001_p_before"
                                           ofType:@"txt"];
    FILE* fp_p_before = fopen((char*)[path UTF8String], "r");
    XCTAssert(fp_p_before!=NULL && path!=nil, @"Setup problem");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"euli-001_p_after"
                                           ofType:@"txt"];
    FILE* fp_p_after = fopen((char*)[path UTF8String], "r");
    XCTAssert(fp_p_after!=NULL && path!=nil, @"Setup problem");
    path = nil;
    
    double* h1 = (double*)malloc(2*85*sizeof(double)); double* h2 = (double*)malloc(2*85*sizeof(double)); double* h3 = (double*)malloc(2*85*sizeof(double));
    double* h4 = (double*)malloc(2*85*sizeof(double)); double* h5 = (double*)malloc(2*85*sizeof(double)); double* h6 = (double*)malloc(2*85*sizeof(double));
    double* h7 = (double*)malloc(2*85*sizeof(double)); double* h8 = (double*)malloc(2*85*sizeof(double)); double* h9 = (double*)malloc(2*85*sizeof(double));
    double* h10 = (double*)malloc(2*85*sizeof(double)); double* h11 = (double*)malloc(2*85*sizeof(double)); double* h12 = (double*)malloc(2*85*sizeof(double));
    for (int i=0; i<2*85; ++i){
        fscanf(fp_h1, "%lf ", (h1+i)); fscanf(fp_h2, "%lf ", (h2+i)); fscanf(fp_h3, "%lf ", (h3+i));
        fscanf(fp_h4, "%lf ", (h4+i)); fscanf(fp_h5, "%lf ", (h5+i)); fscanf(fp_h6, "%lf ", (h6+i));
        fscanf(fp_h7, "%lf ", (h7+i)); fscanf(fp_h8, "%lf ", (h8+i)); fscanf(fp_h9, "%lf ", (h9+i));
        fscanf(fp_h10, "%lf ", (h10+i)); fscanf(fp_h11, "%lf ", (h11+i)); fscanf(fp_h12, "%lf ", (h12+i));
    }
    fclose(fp_h1); fclose(fp_h2); fclose(fp_h3); fclose(fp_h4); fclose(fp_h5); fclose(fp_h6);
    fclose(fp_h7); fclose(fp_h8); fclose(fp_h9); fclose(fp_h10); fclose(fp_h11); fclose(fp_h12);
    
    double* x_k_actual = (double*)malloc(state_size*sizeof(double));
    double* x_k_expected = (double*)malloc(state_size*sizeof(double));
    double* p_k_actual = (double*)malloc(state_size*state_size*sizeof(double));
    double* p_k_expected = (double*)malloc(state_size*state_size*sizeof(double));
    
    for (int i=0; i<state_size; ++i){
        fscanf(fp_x_after, "%lf ", (x_k_expected+i));
        fscanf(fp_x_before, "%lf ", (x_k_actual+i));
    }
    for (int i=0; i<state_size*state_size; ++i){
        fscanf(fp_p_after, "%lf ", (p_k_expected+i));
        fscanf(fp_p_before, "%lf ", (p_k_actual+i));
    }
    
    Filter* filter = [[Filter alloc]init];
    filter->p_k_km1 = p_k_actual;
    filter->x_k_km1 = x_k_actual;
    filter->p_k_k = (double*)malloc(state_size*state_size*sizeof(double));
    filter->x_k_k = (double*)malloc(state_size*sizeof(double));
    filter->state_size = state_size;
    
    cv::Mat* H1 = new cv::Mat(2,85,CV_64FC1,h1,cv::Mat::AUTO_STEP); cv::Mat* H2 = new cv::Mat(2,85,CV_64FC1,h2,cv::Mat::AUTO_STEP);
    cv::Mat* H3 = new cv::Mat(2,85,CV_64FC1,h3,cv::Mat::AUTO_STEP);
    cv::Mat* H4 = new cv::Mat(2,85,CV_64FC1,h4,cv::Mat::AUTO_STEP); cv::Mat* H5 = new cv::Mat(2,85,CV_64FC1,h5,cv::Mat::AUTO_STEP);
    cv::Mat* H6 = new cv::Mat(2,85,CV_64FC1,h6,cv::Mat::AUTO_STEP);
    cv::Mat* H7 = new cv::Mat(2,85,CV_64FC1,h7,cv::Mat::AUTO_STEP); cv::Mat* H8 = new cv::Mat(2,85,CV_64FC1,h8,cv::Mat::AUTO_STEP);
    cv::Mat* H9 = new cv::Mat(2,85,CV_64FC1,h9,cv::Mat::AUTO_STEP);
    cv::Mat* H10 = new cv::Mat(2,85,CV_64FC1,h10,cv::Mat::AUTO_STEP); cv::Mat* H11 = new cv::Mat(2,85,CV_64FC1,h11,cv::Mat::AUTO_STEP);
    cv::Mat* H12 = new cv::Mat(2,85,CV_64FC1,h12,cv::Mat::AUTO_STEP);
    
    feat1->H = H1; feat2->H = H2; feat3->H = H3; feat4->H = H4; feat5->H = H5; feat6->H = H6;
    feat7->H = H7; feat8->H = H8; feat9->H = H9; feat10->H = H10; feat11->H = H11; feat12->H = H12;
    
    NSMutableArray* features_info = [[NSMutableArray alloc]init];
    [features_info addObject: feat1]; [features_info addObject: feat2]; [features_info addObject: feat3];
    [features_info addObject: feat4]; [features_info addObject: feat5]; [features_info addObject: feat6];
    [features_info addObject: feat7];[features_info addObject: feat8]; [features_info addObject: feat9];
    [features_info addObject: feat10]; [features_info addObject: feat11]; [features_info addObject: feat12];
    
    bool success = ekf_update_li_inliers(filter, features_info);
    XCTAssert(success==true, @"Should return true for valid data");
    
    
    for(int i=0; i<state_size; ++i){
        double actual = filter->x_k_k[i];
        double expected = x_k_expected[i];
        double diff = fabs(actual-expected);
        if (diff>0.00001){
            XCTAssert(NO, @"x error exceeds limits");
            break;
        }
    }
    free(x_k_expected);
    for(int i=0; i<state_size*state_size; ++i){
        double actual = filter->p_k_k[i];
        double expected = p_k_expected[i];
        double diff = fabs(actual-expected);
        if (diff>0.00008){ //The accuracy of this matrix is not very good
            XCTAssert(NO, @"p error exceeds limits");
            break;
        }
    }
    free(p_k_expected);
    [filter destroy];
    for (int i=0; i<[features_info count]; ++i){
        [(Feature*)(features_info[i]) destroy];
    }
    
    free(h1);free(h2);free(h3);free(h4);free(h5);free(h6);free(h7);free(h8);free(h9);free(h10);free(h11);free(h12);
    
}


@end


