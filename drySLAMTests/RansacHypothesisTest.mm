//
//  RansacHypothesisTests.m
//  drySLAM
//
//  Created by Anthony Rodriguez on 3/27/15.
//  Copyright (c) 2015 ECE473 573 Common. All rights reserved.
//


#include "SLAMModel.hpp"
#include "RansacHypothesis.hpp"

#import <UIKit/UIKit.h>
#import <XCTest/XCTest.h>
#include <stdio.h>
#include <stdlib.h>




@interface RansacHypothesisTests : XCTestCase

@end

@implementation RansacHypothesisTests

- (void)setUp {
    [super setUp];
    // Put setup code here. This method is called before the invocation of each test method in the class.
}

- (void)tearDown {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
    [super tearDown];
}

/************************************************RANSAC_HYPOTHESIS**************************************************/


-(void)testRansacHypothesis000{
    //Implement as necessary
    XCTAssert(YES, @"Function requires validation");
}


/******************************************GENERATE_STATE_VECTOR_PATTERN********************************************/
//001: Test for return value of false when inputs are invalid
-(void)testGenerateStateVectorPattern001{
    NSMutableArray* features_info = [[NSMutableArray alloc]init];
    int state_size = 22;
    double* state_vector_pattern = (double*)malloc(22*4*sizeof(double));
    NSMutableArray* z_id = [[NSMutableArray alloc]init];
    NSMutableArray* z_euc = [[NSMutableArray alloc]init];
    for (int i=0; i<state_size; ++i){
        Feature* feat = [[Feature alloc]init];
        [features_info addObject:feat];
    }
    bool success = generate_state_vector_pattern(NULL, state_size,
                                                 state_vector_pattern, z_id, z_euc);
    XCTAssert(success==false, @"generate_shate_vector_pattern should return false when features_info is NULL");
    success = generate_state_vector_pattern(features_info, -1,
                                            state_vector_pattern, z_id, z_euc);
    XCTAssert(success==false, @"generate_shate_vector_pattern should return false when state_size<0");
    success = generate_state_vector_pattern(features_info, state_size,
                                            NULL, z_id, z_euc);
    XCTAssert(success==false, @"generate_shate_vector_pattern should return false when state_vector_pattern arg is NULL");
    success = generate_state_vector_pattern(features_info, state_size,
                                            state_vector_pattern, NULL, z_euc);
    XCTAssert(success==false, @"generate_shate_vector_pattern should return false when z_id is NULL");
    success = generate_state_vector_pattern(features_info, state_size,
                                            state_vector_pattern, z_id, NULL);
    XCTAssert(success==false, @"generate_shate_vector_pattern should return false when z_eud is NULL");
    for (int i=0; i<[features_info count]; ++i){
        Feature* currFeat = features_info[i];
        [currFeat destroy];
        currFeat= nil;
    }
    features_info = nil;
    free(state_vector_pattern);
    z_euc = nil;
    z_id = nil;
}
//002: From program
-(void)testGenerateStateVectorPattern002{
    int state_size = 121;
    double* state_vector_pattern_actual = (double*)malloc(121*4*sizeof(double));
    double* state_vector_pattern_expected = (double*)malloc(121*4*sizeof(double));
    NSMutableArray* z_id = [[NSMutableArray alloc]init];
    NSMutableArray* z_euc = [[NSMutableArray alloc]init];
    int* z_id_expected = (int*)malloc(2*18*sizeof(int));
    NSMutableArray* features_info = [[NSMutableArray alloc]init];
    Feature* feat1 = [[Feature alloc]init];
    feat1->z = (int*)malloc(2*sizeof(int));
    feat1->z[0] = 46; feat1->z[1] = 130;
    [features_info addObject:feat1];
    Feature* feat2 = [[Feature alloc]init];
    feat2->z = (int*)malloc(2*sizeof(int));
    feat2->z[0] = 86; feat2->z[1] = 155;
    [features_info addObject:feat2];
    Feature* feat3 = [[Feature alloc]init];
    feat3->z = (int*)malloc(2*sizeof(int));
    feat3->z[0] = 228; feat3->z[1] = 143;
    [features_info addObject:feat3];
    Feature* feat4 = [[Feature alloc]init];
    feat4->z = (int*)malloc(2*sizeof(int));
    feat4->z[0] = 27; feat4->z[1] = 164;
    [features_info addObject:feat4];
    Feature* feat5 = [[Feature alloc]init];
    feat5->z = (int*)malloc(2*sizeof(int));
    feat5->z[0] = 68; feat5->z[1] = 45;
    [features_info addObject:feat5];
    Feature* feat6 = [[Feature alloc]init];
    feat6->z = (int*)malloc(2*sizeof(int));
    feat6->z[0] = 77; feat6->z[1] = 64;
    [features_info addObject:feat6];
    Feature* feat7 = [[Feature alloc]init];
    feat7->z = (int*)malloc(2*sizeof(int));
    feat7->z[0] = 71; feat7->z[1] = 133;
    [features_info addObject:feat7];
    Feature* feat8 = [[Feature alloc]init];
    feat8->z = (int*)malloc(2*sizeof(int));
    feat8->z[0] = 120; feat8->z[1] = 98;
    [features_info addObject:feat8];
    Feature* feat9 = [[Feature alloc]init];
    feat9->z = (int*)malloc(2*sizeof(int));
    feat9->z[0] = 159; feat9->z[1] = 104;
    [features_info addObject:feat9];
    Feature* feat10 = [[Feature alloc]init];
    feat10->z = (int*)malloc(2*sizeof(int));
    feat10->z[0] = 277; feat10->z[1] = 32;
    [features_info addObject:feat10];
    Feature* feat11 = [[Feature alloc]init];
    feat11->z = (int*)malloc(2*sizeof(int));
    feat11->z[0] = 91; feat11->z[1] = 57;
    [features_info addObject:feat11];
    Feature* feat12 = [[Feature alloc]init];
    feat12->z = (int*)malloc(2*sizeof(int));
    feat12->z[0] = 69; feat12->z[1] = 197;
    [features_info addObject:feat12];
    Feature* feat13 = [[Feature alloc]init];
    feat13->z = (int*)malloc(2*sizeof(int));
    feat13->z[0] = 90; feat13->z[1] = 118;
    [features_info addObject:feat13];
    Feature* feat14 = [[Feature alloc]init];
    feat14->z = (int*)malloc(2*sizeof(int));
    feat14->z[0] = 229; feat14->z[1] = 157;
    [features_info addObject:feat14];
    Feature* feat15 = [[Feature alloc]init];
    feat15->z = (int*)malloc(2*sizeof(int));
    feat15->z[0] = 166; feat15->z[1] = 186;
    [features_info addObject:feat15];
    Feature* feat16 = [[Feature alloc]init];
    feat16->z = (int*)malloc(2*sizeof(int));
    feat16->z[0] = 33; feat16->z[1] = 154;
    [features_info addObject:feat16];
    Feature* feat17 = [[Feature alloc]init];
    feat17->z = (int*)malloc(2*sizeof(int));
    feat17->z[0] = 147; feat17->z[1] = 116;
    [features_info addObject:feat17];
    Feature* feat18 = [[Feature alloc]init];
    feat18->z = (int*)malloc(2*sizeof(int));
    feat18->z[0] = 113; feat18->z[1] = 146;
    [features_info addObject:feat18];
    NSString* path = [[NSBundle mainBundle] pathForResource:@"gsvp-002_state_vector_pattern"
                                                     ofType:@"txt"];
    FILE* fp_svp = fopen((char*)[path UTF8String], "r");
    XCTAssert(fp_svp!=NULL && path!=nil, @"Setup problem");
    for (int i=0; i<(121*4); ++i){
        fscanf(fp_svp, "%lf ", (state_vector_pattern_expected+i));
        state_vector_pattern_actual[i] = -1.0;
    }
    fclose(fp_svp);
    path = [[NSBundle mainBundle] pathForResource:@"gsvp-002_z_id"
                                           ofType:@"txt"];
    FILE* fp_z_id = fopen((char*)[path UTF8String], "r");
    XCTAssert(fp_z_id!=NULL && path!=nil, @"Setup problem");
    for (int i=0; i<(2*18); ++i){
        double num;
        fscanf(fp_z_id, "%lf ", &num);
        z_id_expected[i] = (int)num;
    }
    fclose(fp_z_id);
    
    bool success = generate_state_vector_pattern(features_info, state_size, state_vector_pattern_actual, z_id, z_euc);
    
    XCTAssert(success==true, @"Expect function to return true for valid data");
    
    if ([z_euc count]!=0){
        XCTAssert(NO,@"z_euc should be an empty set");
        free(state_vector_pattern_actual);
        free(state_vector_pattern_expected);
        for (int i=0; i<18; ++i){
            Feature* currFeat = (Feature*)features_info[i];
            [currFeat destroy];
            currFeat = nil;
        }
        features_info = nil;
        return;
    }
    if ([z_id count]!=(2*18)){
        XCTAssert(NO,@"z_id should be a 2x18 matrix (36 elements)");
        free(state_vector_pattern_actual);
        free(state_vector_pattern_expected);
        for (int i=0; i<18; ++i){
            Feature* currFeat = (Feature*)features_info[i];
            [currFeat destroy];
            currFeat = nil;
        }
        features_info = nil;
        return;
    }
    
    for (int i=0; i<(2*18); ++i){
        int diff = abs(z_id_expected[i] - (int)[z_id[i] integerValue]);
        if (diff!=0){
            XCTAssert(NO, @"Error in z_id exceeds limits");
            break;
        }
    }
    for (int i=0; i<(121*4); ++i){
        double actual = state_vector_pattern_actual[i];
        double expected = state_vector_pattern_expected[i];
        double diffD = fabs(actual-expected);
        if (diffD>=0.00001){
            XCTAssert(NO, @"Error in state_vector_pattern exceeds limits");
            break;
        }
    }
    free(state_vector_pattern_actual);
    free(state_vector_pattern_expected);
    free(z_id_expected);
    for (int i=0; i<18; ++i){
        Feature* currFeat = (Feature*)features_info[i];
        [currFeat destroy];
        currFeat = nil;
    }
    features_info = nil;
}

/***************************************************SELECT_RANDOM_MATCH**************************************************/


-(void)testSelectRandomMatch001{
    NSMutableArray* features_info = [[NSMutableArray alloc]init];
    for (int i=0; i<18; ++i){
        Feature* feat = [[Feature alloc]init];
        feat->individually_compatible = 1;
        if (feat->z==NULL)
            feat->z = (int*)malloc(2*sizeof(int));
        feat->z[0] = (int)(10*i*cos(rand()));
        feat->z[1] = (int)(10*i*sin(rand()));
        [features_info addObject:feat];
    }
    
    int* zi = (int*)malloc(2*sizeof(int));
    int numMatches = 0;
    int position = -1;
    
    bool success = select_random_match(features_info, zi, &position, &numMatches);
    XCTAssert(success==true, @"expected true for valid data");
    if (position<0 || position>=[features_info count]){
        XCTAssert(NO, @"position not valid");
        return;
    }
    Feature* feat = features_info[position];
    XCTAssert(feat->z[0]==zi[0], @"z values don't match");
    XCTAssert(feat->z[1]==zi[1], @"z values don't match");
    XCTAssert(numMatches==18, @"num matches wasn't set correctly");
    
    for (int i=0; i<[features_info count]; ++i){
        [(Feature*)features_info[i] destroy];
    }
    features_info = nil;
    
}

/************************************************COMPUTE_HYPOTHESIS_SUPPORT**************************************************/


-(void)testComputeHypothesisSupport001{
    int stateSize = 121;
    double* stateVectorPattern = (double*)malloc(4*stateSize*sizeof(double));
    double* xi = (double*)malloc(stateSize*sizeof(double));
    double* z_id = (double*)malloc(2*18*sizeof(double));
    double* z_euc = NULL;
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
    
    NSString* path = [[NSBundle mainBundle] pathForResource:@"chs-001_hi"
                                                     ofType:@"txt"];
    FILE* fp_x = fopen((char*)[path UTF8String], "r");
    path = nil;
    for (int i=0; i<stateSize; ++i){
        fscanf(fp_x, "%lf ", xi+i);
    }
    fclose(fp_x);
    
    path = [[NSBundle mainBundle] pathForResource:@"chs-001_svp" ofType:@"txt"];
    FILE* fp_svp = fopen((char*)[path UTF8String], "r");
    path = nil;
    for (int i=0; i<(4*stateSize); ++i){
        fscanf(fp_svp, "%lf ", stateVectorPattern+i);
    }
    fclose(fp_svp);
    
    path = [[NSBundle mainBundle] pathForResource:@"chs-001_zid"
                                           ofType:@"txt"];
    FILE* fp_zid = fopen((char*)[path UTF8String], "r");
    path = nil;
    for (int i=0; i<(2*18); ++i){
        fscanf(fp_zid, "%lf ", z_id+i);
    }
    fclose(fp_zid);
    
    int* pos_li_inliers_id_expected = (int*)malloc(18*sizeof(int));
    path = [[NSBundle mainBundle] pathForResource:@"chs-001_pos_li_id"
                                           ofType:@"txt"];
    FILE* fp_pli = fopen((char*)[path UTF8String], "r");
    path = nil;
    for (int i=0; i<(18); ++i){
        double val;
        fscanf(fp_pli, "%lf ", &val);
        pos_li_inliers_id_expected[i] = (int)val;
    }
    fclose(fp_pli);
    
    NSMutableArray* positions_li_inliers_id = [[NSMutableArray alloc]init];
    NSMutableArray* positions_li_inliers_euc = [[NSMutableArray alloc]init];
    int threshold = 1;
    int* hypothesis_support = (int*)malloc(sizeof(int));
    
    bool success = compute_hypothesis_support(xi, stateSize, cam, stateVectorPattern, z_id, 18, z_euc, threshold, hypothesis_support, positions_li_inliers_id, positions_li_inliers_euc);
    
    XCTAssert(success==true, @"Expected true for valid data");
    XCTAssert(*hypothesis_support==14, @"hypothesis support incorrect");
    XCTAssert([positions_li_inliers_euc count]==0, @"Should have empty euclidian vector");
    if ([positions_li_inliers_id count]!=18){
        XCTAssert(NO, @"positions_li_inliers_id should have 18 elements here");
        return;
    }
    for (int i=0; i<[positions_li_inliers_id count]; ++i){
        int actual = (int)[(NSNumber*)positions_li_inliers_id[i] doubleValue];
        int expected = pos_li_inliers_id_expected[i];
        if (actual!=expected){
            XCTAssert(NO, @"positions_li_inliers error exceeds limits");
            break;
        }
    }
    
    free(xi);
    free(z_id);
    free(stateVectorPattern);
    positions_li_inliers_id = nil;
    
}


/*************************SET_AS_MOST_SUPPORTED_HYPOTHESIS***********************/
//features_info length=12, positions_li_inliers euc=zeros, id=zeros
//Test 001: Function should return false and not crash if features_info==NULL
-(void)testSetAsMostSupportedHypothesis001{
    bool success = false;
    NSMutableArray* features_info = [[NSMutableArray alloc]initWithCapacity:1];
    Feature* feat1 = [[Feature alloc]init];
    [features_info addObject:feat1];
    NSMutableArray* positions_li_inliers_euc = [[NSMutableArray alloc]initWithCapacity:12];
    NSMutableArray* positions_li_inliers_id = [[NSMutableArray alloc]initWithCapacity:12];
    for(int i=0; i<12; ++i){
        [positions_li_inliers_euc addObject:[NSNumber numberWithInt:0]];
        [positions_li_inliers_id addObject:[NSNumber numberWithInt:0]];
    }
    success = set_as_most_supported_hypothesis(NULL, positions_li_inliers_id, positions_li_inliers_euc);
    XCTAssert(success==false, @"set_as_most_supported_hypothesis should return false and not crash if features_info is NULL");
    positions_li_inliers_euc = nil;
    positions_li_inliers_id = nil;
}

//features_info length=12, positions_li_inliers euc=NULL, id=NULL
//Test 002: Function should return true if all arrays are empty
-(void)testSetAsMostSupportedHypothesis002{
    bool success = false;
    NSMutableArray* features_info = [[NSMutableArray alloc]initWithCapacity:1];
    Feature* feat1 = [[Feature alloc]init];
    [features_info addObject:feat1];
    
    success = set_as_most_supported_hypothesis(features_info, NULL, NULL);
    XCTAssert(success==false, @"set_as_most_supported_hypothesis should return false and not crash if inlier arrays are NULL");
}


@end
