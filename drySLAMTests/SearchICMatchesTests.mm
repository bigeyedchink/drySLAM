//
//  SearchICMatchesTests.m
//  drySLAM
//
//  Created by Anthony Rodriguez on 3/27/15.
//  Copyright (c) 2015 ECE473 573 Common. All rights reserved.
//


#include "SLAMModel.hpp"
#include "SearchICMatches.hpp"

#import <UIKit/UIKit.h>
#import <XCTest/XCTest.h>
#include <stdio.h>
#include <stdlib.h>

//Thresholds used for testing accuracy of matching methods:
#define MATCHING_ERROR_THRESHOLD 10         //pixel intensity difference at which a prediction is considered an error
#define MATCHING_ACCURACY_THRESHOLD 0.82f   //this many pixels must be matched to be considered a success


@interface SearchICMatchesTests: XCTestCase

@end

@implementation SearchICMatchesTests

- (void)setUp {
    [super setUp];
    // Put setup code here. This method is called before the invocation of each test method in the class.
}

- (void)tearDown {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
    [super tearDown];
}



/*************************************SEARCH_IC_MATCHES********************************************/

-(void)testSearchICMatches000{
    //Implement as necessary
    XCTAssert(YES, @"Function requires validation");
}

/************************************CALCULATE_HI_INVERSE_DEPTH*****************************************/

-(void)testCalculateHiInverseDepth001{
    NSMutableArray* features_info = [[NSMutableArray alloc]init];
    for (int i=0; i<19; ++i){
        Feature* feat = [[Feature alloc]init];
        feat->h = (double*)malloc(2*sizeof(double));
        feat->type = [NSString stringWithFormat:@"inversedepth"];
        [features_info addObject:feat];
    }
    ((Feature*)features_info[0])->h[0] = 100*0.429427804012463;
    ((Feature*)features_info[0])->h[1] = 100*1.308839303854355;
    ((Feature*)features_info[1])->h[0] = 100*0.822513926104724;
    ((Feature*)features_info[1])->h[1] = 100*1.561597910402890;
    ((Feature*)features_info[2])->h[0] = 100*2.241159282550573;
    ((Feature*)features_info[2])->h[1] = 100*1.435015824234898;
    ((Feature*)features_info[3])->h[0] = 100*0.243981927633748;
    ((Feature*)features_info[3])->h[1] = 100*1.651445465805020;
    ((Feature*)features_info[4])->h[0] = 64.974168231795005;
    ((Feature*)features_info[4])->h[1] = 46.023644563616074;
    ((Feature*)features_info[5])->h[0] = 73.954718431010420;
    ((Feature*)features_info[5])->h[1] = 64.964188943574641;
    ((Feature*)features_info[6])->h[0] = 100*0.677703250825545;
    ((Feature*)features_info[6])->h[1] = 100*1.346059961255739;
    ((Feature*)features_info[7])->h[0] = 100*1.167676807691521;
    ((Feature*)features_info[7])->h[1] = 100*0.989237751389494;
    ((Feature*)features_info[8])->h[0] = 100*1.548898773407871;
    ((Feature*)features_info[8])->h[1] = 100*1.048418620936615;
    ((Feature*)features_info[9])->h[0] = 100*2.738749086630378;
    ((Feature*)features_info[9])->h[1] = 100*0.318905419580472;
    ((Feature*)features_info[10])->h[0] = 87.795521807369525;
    ((Feature*)features_info[10])->h[1] = 57.847845801414564;
    ((Feature*)features_info[11])->h[0] = 100*0.653642274463169;
    ((Feature*)features_info[11])->h[1] = 100*1.980867120847255;
    ((Feature*)features_info[12])->h[0] = 100*0.860104624324212;
    ((Feature*)features_info[12])->h[1] = 100*1.196873629314743;
    ((Feature*)features_info[13])->h[0] = 100*2.250580205529030;
    ((Feature*)features_info[13])->h[1] = 100*1.580435691862814;
    ((Feature*)features_info[14])->h[0] = 100*1.620438098428692;
    ((Feature*)features_info[14])->h[1] = 100*1.869535272062445;
    ((Feature*)features_info[15])->h[0] = 100*0.308063055303360;
    ((Feature*)features_info[15])->h[1] = 100*1.549606647956691;
    ((Feature*)features_info[16])->h[0] = 100*1.430338884294521;
    ((Feature*)features_info[16])->h[1] = 100*1.169878306665135;
    ((Feature*)features_info[17])->h[0] = 100*1.090804326966315;
    ((Feature*)features_info[17])->h[1] = 100*1.471692623147241;
    ((Feature*)features_info[18])->h[0] = 100*1.121091110146125;
    ((Feature*)features_info[18])->h[1] = 100*1.740652752322901;
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
    int index = 17; //in matlab this is 18
    double* cartesianFeaturesIndex = (double*)malloc(19*sizeof(double));
    double* inverseDepthFeaturesIndex = (double*)malloc(19*sizeof(double));
    for (int i=0; i<19; ++i){
        cartesianFeaturesIndex[i] = 0;
        inverseDepthFeaturesIndex[i] = 1;
    }
    NSString* path = [[NSBundle mainBundle] pathForResource:@"chid-001_cam_K"
                                                     ofType:@"txt"];
    FILE* fp_cam_k = fopen((char*)[path UTF8String], "r");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"chid-001_Hi_after"
                                           ofType:@"txt"];
    FILE* fp_Hi_after = fopen((char*)[path UTF8String], "r");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"chid-001_xv_km1_k"
                                           ofType:@"txt"];
    FILE* fp_xv = fopen((char*)[path UTF8String], "r");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"chid-001_yi"
                                           ofType:@"txt"];
    FILE* fp_yi = fopen((char*)[path UTF8String], "r");
    path = nil;
    
    double* xv_km1_k = (double*)malloc(13*sizeof(double));
    for (int i=0; i<13; ++i){
        fscanf(fp_xv, "%lf ", xv_km1_k+i);
    }
    fclose(fp_xv);
    
    double* Hi_actual = (double*)malloc(2*127*sizeof(double));
    double* Hi_expected = (double*)malloc(2*127*sizeof(double));
    for (int i=0; i<2*127; ++i){
        fscanf(fp_Hi_after, "%lf ", Hi_expected+i);
        Hi_actual[i] = -99;
    }
    fclose(fp_Hi_after);
    
    double* yi= (double*)malloc(6*sizeof(double));
    for (int i=0; i<6; ++i){
        fscanf(fp_yi, "%lf ", yi+i);
    }
    fclose(fp_yi);
    
    double* camK= (double*)malloc(9*sizeof(double));
    for (int i=0; i<9; ++i){
        fscanf(fp_cam_k, "%lf ", camK+i);
    }
    fclose(fp_cam_k);
    
    bool success = calculate_hi_inverse_depth(Hi_actual, inverseDepthFeaturesIndex, cartesianFeaturesIndex, xv_km1_k, yi,
                                              cam, index, features_info);
    
    XCTAssert(success==true, @"Should return true for valid data");
    
    for (int i=0; i<2*127; ++i){
        double actual = Hi_actual[i];
        double expected = Hi_expected[i];
        double diff = fabs(actual-expected);
        if (diff>0.00001){
            XCTAssert(NO, @"Hi error exceeds limits");
            break;
        }
    }
    
    for (int i=0; i<[features_info count]; ++i){
        [(Feature*)features_info[i] destroy];
    }
    
    free(xv_km1_k);
    free(Hi_expected);
    free(Hi_actual);
    free(inverseDepthFeaturesIndex);
    free(cartesianFeaturesIndex);
    free(yi);
    [cam destroy];
    
    
}

/**********************************************************DH_DXV*************************************************************/

-(void)testDhDxv001{
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
    
    double* Xv_km1_k = (double*)malloc(13*sizeof(double));
    double* yi = (double*)malloc(6*sizeof(double));
    double* zi = (double*)malloc(2*sizeof(double));
    double* Hi_expected = (double*)malloc(2*13*sizeof(double));
    double* Hi_actual = (double*)malloc(2*13*sizeof(double));
    
    Xv_km1_k[0] = 0.008468930540656;
    Xv_km1_k[1] = 0.001899369789428;
    Xv_km1_k[2] = 0.001070000974077;
    Xv_km1_k[3] = 0.999980501496129;
    Xv_km1_k[4] = 0.003192856440414;
    Xv_km1_k[5] = 0.005333549155932;
    Xv_km1_k[6] = 0.000596279047264;
    Xv_km1_k[7] = 0.008468930540656;
    Xv_km1_k[8] = 0.001899369789428;
    Xv_km1_k[9] = 0.001070000974077;
    Xv_km1_k[10] = 0.006385837395103;
    Xv_km1_k[11] = 0.010667306308377;
    Xv_km1_k[12] = 0.001192581348079;
    yi[0] = 0.0;
    yi[1] = 0.0;
    yi[2] = 0.0;
    yi[3] = -0.555613960977581;
    yi[4] = 0.008786183032201;
    yi[5] =  1.000000000000000;
    zi[0] = 100*0.492691391184923;
    zi[1] = 100*1.279211391030913;
    
    NSString* path = [[NSBundle mainBundle] pathForResource:@"dhdhv-001_Hi1"
                                                     ofType:@"txt"];
    FILE* fp_hi = fopen((char*)[path UTF8String], "r");
    path = nil;
    
    for (int i=0; i<26; ++i){
        Hi_actual[i] = 0;
        fscanf(fp_hi, "%lf ", Hi_expected + i);
    }
    fclose(fp_hi);
    
    bool success = dh_dxv(Hi_actual, cam, Xv_km1_k, yi, zi);
    XCTAssert(success==true, @"Expected return value of true for valid data");
    
    for (int i=0; i<26; ++i){
        double actual = Hi_actual[i];
        double expected = Hi_expected[i];
        double diff = fabs(actual-expected);
        if (diff>0.00001){
            XCTAssert(diff<0.00001, @"Hi error exceeds limits");
            break;
        }
    }
    free(Xv_km1_k);
    free(yi);
    free(zi);
    free(Hi_expected);
    free(Hi_actual);
    
}

-(void) testDhDxv002{
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
    double* Xv_km1_k = (double*)malloc(13*sizeof(double));
    double* yi = (double*)malloc(6*sizeof(double));
    double* zi = (double*)malloc(2*sizeof(double));
    double* Hi_actual = (double*)malloc(26*sizeof(double));
    double* Hi_expected = (double*)malloc(26*sizeof(double));
    Xv_km1_k[0] = 0;
    Xv_km1_k[1] = 0;
    Xv_km1_k[2] = 0;
    Xv_km1_k[3] = 1.000000000000000;
    Xv_km1_k[4] = 0.000000000000001;
    Xv_km1_k[5] = 0.000000000000001;
    Xv_km1_k[6] = 0.000000000000001;
    Xv_km1_k[7] = 0;
    Xv_km1_k[8] = 0;
    Xv_km1_k[9] = 0;
    Xv_km1_k[10] = 0.000000000000001;
    Xv_km1_k[11] = 0.000000000000001;
    Xv_km1_k[12] = 0.000000000000001;

    yi[0] = 0;
    yi[1] = 0;
    yi[2] = 0;
    yi[3] = -0.418255484817285;
    yi[4] = -0.010670574181629;
    yi[5] = 1.000000000000000;

    zi[0] = 100*0.789999999999998;
    zi[1] = 100*1.310000000000003;
    
    NSString* path = [[NSBundle mainBundle] pathForResource:@"dhdxv-002_Hi1"
                                                     ofType:@"txt"];
    FILE* fp_h = fopen((char*)[path UTF8String], "r");
    path = nil;
    
    for (int i=0; i<26; ++i){
        Hi_actual[i] = 0.0;
        fscanf(fp_h, "%lf ", Hi_expected+i);
    }
    fclose(fp_h);
    bool success =  dh_dxv(Hi_actual, cam, Xv_km1_k, yi, zi);
    XCTAssert(success==true, @"Should return true for valid data");
    for (int i=0; i<26; ++i){
        double actual = Hi_actual[i];
        double expected = Hi_expected[i];
        double diff = fabs(actual-expected);
        if (diff>0.00001){
            XCTAssert(diff<0.00001, @"Hi error exceeds limits");
            break;
        }
    }
    free(Xv_km1_k);
    free(yi);
    free(zi);
    free(Hi_expected);
    free(Hi_actual);
    [cam destroy];
}

/*************************************************************DH_DY******************************************************************/

-(void)testDhDy001{
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
    
    double* Xv_km1_k = (double*)malloc(13*sizeof(double));
    double* yi = (double*)malloc(6*sizeof(double));
    double* zi = (double*)malloc(2*sizeof(double));
    double* Hi_expected = (double*)malloc(12*sizeof(double));
    double* Hi_actual = (double*)malloc(12*sizeof(double));
    
    Xv_km1_k[0] = 0.008468930540656;
    Xv_km1_k[1] = 0.001899369789428;
    Xv_km1_k[2] = 0.001070000974077;
    Xv_km1_k[3] = 0.999980501496129;
    Xv_km1_k[4] = 0.003192856440414;
    Xv_km1_k[5] = 0.005333549155932;
    Xv_km1_k[6] = 0.000596279047264;
    Xv_km1_k[7] = 0.008468930540656;
    Xv_km1_k[8] = 0.001899369789428;
    Xv_km1_k[9] = 0.001070000974077;
    Xv_km1_k[10] = 0.006385837395103;
    Xv_km1_k[11] = 0.010667306308377;
    Xv_km1_k[12] = 0.001192581348079;
    yi[0] = 0.0;
    yi[1] = 0.0;
    yi[2] = 0.0;
    yi[3] = 0.389224640622712;
    yi[4] = -0.064068769482697;
    yi[5] = 1.000000000000000;
    zi[0] = 100* 2.322566479278837;
    zi[1] = 100*1.423737141884500;
    
    
    NSString* path = [[NSBundle mainBundle] pathForResource:@"dhdy-001_Hii"
                                                     ofType:@"txt"];
    FILE* fp_hi = fopen((char*)[path UTF8String], "r");
    path = nil;
    
    for (int i=0; i<12; ++i){
        Hi_actual[i] = 0;
        fscanf(fp_hi, "%lf ", Hi_expected + i);
    }
    fclose(fp_hi);
    
    bool success = dh_dy(Hi_actual, cam, Xv_km1_k, yi, zi);
    XCTAssert(success==true, @"Expected return value of true for valid data");
    
    for (int i=0; i<12; ++i){
        double actual = Hi_actual[i];
        double expected = Hi_expected[i];
        double diff = fabs(actual-expected);
        if (diff>0.00001){
            XCTAssert(diff<0.00001, @"Hi error exceeds limits");
            break;
        }
    }
    free(Xv_km1_k);
    free(yi);
    free(zi);
    free(Hi_expected);
    free(Hi_actual);
}

/***************************************************************DH_DHRL****************************************************************/

-(void)testDhDhrl001{
    double* a_expected = (double*)malloc(6*sizeof(double));
    double* a_actual = (double*)malloc(6*sizeof(double));
    double* Xv_km1_k = (double*)malloc(13*sizeof(double));
    double* yi = (double*)malloc(13*sizeof(double));
    double* zi = (double*)malloc(2*sizeof(double));
    
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
    
    a_expected[0] = 100*1.784569751516351;
    a_expected[1] = 100*-0.276497630445984;
    a_expected[2] = 100*0.781590981160553;
    a_expected[3] = 100*-0.276497630445984;
    a_expected[4] = 100*1.786865061852187;
    a_expected[5] = 100*0.778353574176348;
    
    Xv_km1_k[0] = 0;
    Xv_km1_k[1] = 0;
    Xv_km1_k[2] = 0;
    Xv_km1_k[3] = 1.000000000000000;
    Xv_km1_k[4] = 0.000000000000001;
    Xv_km1_k[5] = 0.000000000000001;
    Xv_km1_k[6] = 0.000000000000001;
    Xv_km1_k[7] = 0;
    Xv_km1_k[8] = 0;
    Xv_km1_k[9] = 0;
    Xv_km1_k[10] = 0.000000000000001;
    Xv_km1_k[11] = 0.000000000000001;
    Xv_km1_k[12] = 0.000000000000001;
    
    yi[0] = 0;
    yi[1] = 0;
    yi[2] = 0;
    yi[3] = -0.477847760223840;
    yi[4] = 0.429457243275458;
    yi[5] = 1.000000000000000;
    
    zi[0] = 73.999999999999730;
    zi[1] = 43.000000000000270;
    
    bool success = dh_dhrl(a_actual, cam, Xv_km1_k, yi, zi);
    XCTAssert(success==true, @"Expected true return for valid data");
    for (int i=0; i<6; ++i){
        double actual = a_actual[i];
        double expected = a_expected[i];
        double diff = fabs(actual-expected);
        if (diff>0.00001){
            XCTAssert(NO, @"dh_dhrl error exceeds limits");
            break;
        }
    }
    
    free(a_expected);
    free(a_actual);
    free(Xv_km1_k);
    free(yi);
    free(zi);
    [cam destroy];
    
}

/***************************************************************DHRL_DQWR****************************************************************/

-(void)testDhrlDqwr001{
    double* Xv_km1_k = (double*)malloc(13*sizeof(double));
    double* yi = (double*)malloc(6*sizeof(double));
    double* a_expected = (double*)malloc(12*sizeof(double));
    double* a_actual = (double*)malloc(12*sizeof(double));
    
    Xv_km1_k[0] = 0;
    Xv_km1_k[1] = 0;
    Xv_km1_k[2] = 0;
    Xv_km1_k[3] = 1.000000000000000;
    Xv_km1_k[4] = 0.000000000000001;
    Xv_km1_k[5] = 0.000000000000001;
    Xv_km1_k[6] = 0.000000000000001;
    Xv_km1_k[7] = 0;
    Xv_km1_k[8] = 0;
    Xv_km1_k[9] = 0;
    Xv_km1_k[10] = 0.000000000000001;
    Xv_km1_k[11] = 0.000000000000001;
    Xv_km1_k[12] = 0.000000000000001;
    
    yi[0] = 0;
    yi[1] = 0;
    yi[2] = 0;
    yi[3] = -0.418255484817285;
    yi[4] = -0.010670574181629;
    yi[5] = 1.000000000000000;
    
    a_expected[0] = -0.812287625548368;
    a_expected[1] = 0.000000000000001;
    a_expected[2] = -1.827493744463462;
    a_expected[3] = 0.021340743377605;
    a_expected[4] = 0.021340743377605;
    a_expected[5] = 1.827493744463462;
    a_expected[6] = 0.000000000000001;
    a_expected[7] = 0.812287625548368;
    a_expected[8] = 1.827493744463462;
    a_expected[9] = -0.021340743377605;
    a_expected[10] = -0.812287625548368;
    a_expected[11] = 0.000000000000001;
    
    bool success = dhrl_dqwr(a_actual, Xv_km1_k, yi);
    XCTAssert(success==true, @"Expected return value of true for valid data");
    
    for (int i=0; i<12; ++i){
        double actual = a_actual[i];
        double expected = a_expected[i];
        double diff = fabs(actual-expected);
        if (diff>0.00001){
            XCTAssert(NO, @"a error exceeds limits");
            break;
        }
    }
    
    free(Xv_km1_k);
    free(yi);
    free(a_expected);
    free(a_actual);
}

/*********************************************************DRQ_TIMES_A_BY_DQ************************************************************/
//Test function which is located in mapManagement


/***************************************************INVERSEDEPTH_2_CARTESIAN*********************************************************/

-(void)testInverseDepth2Cartesian001{
    double* inverse_depth = (double*)malloc(6*sizeof(double));
    inverse_depth[0] = 0;
    inverse_depth[1] = 0;
    inverse_depth[2] = 0;
    inverse_depth[3] = -0.555548750902877;
    inverse_depth[4] = 0.008813479020505;
    inverse_depth[5] = 1.000000000000000;
    double* cartesian_actual = (double*)malloc(3*sizeof(double));
    double* cartesian_expected = (double*)malloc(3*sizeof(double));
    cartesian_expected[0] = -0.527389120701184;
    cartesian_expected[1] = -0.008813364919574;
    cartesian_expected[2] = 0.849578154124049;
    bool success = inverseDepth2Cartesian(inverse_depth, cartesian_actual);
    XCTAssert(success==true, @"expect return value of true for valid data");
    double diff = fabs(cartesian_actual[0]-cartesian_expected[0]);
    XCTAssert(diff<0.00001, @"coordinate vector exceeds limits");
    diff = fabs(cartesian_actual[1]-cartesian_expected[1]);
    XCTAssert(diff<0.00001, @"coordinate vector exceeds limits");
    diff = fabs(cartesian_actual[2]-cartesian_expected[2]);
    XCTAssert(diff<0.00001, @"coordinate vector exceeds limits");
    free(inverse_depth);
    free(cartesian_expected);
    free(cartesian_actual);
}

/*************************************PREDICT_FEATURES_APPEARANCE********************************************/


//Basically, we just pick a random feature from features_info and make sure its' patch_when_matching value
//is updated to the expected matrix value
//This method should likely pass if both inversedepth2cartesian and pred_match_fc are working
-(void)testPredictFeaturesAppearance001{
    NSMutableArray* features_info = [[NSMutableArray alloc]init];
    for (int i=0; i<19; ++i){
        Feature* feat = [[Feature alloc]init];
        feat->type = [NSString stringWithFormat:@"inversedepth"];
        [features_info addObject:feat];
    }
    Feature* testFeat = (Feature*)features_info[4];
    if (testFeat->h==NULL)
        testFeat->h = (double*)malloc(2*sizeof(double));
    testFeat->h[0] = 64.877343094374879;
    testFeat->h[1] = 46.049408739718260;
    testFeat->half_patch_size_when_matching = 6;
    testFeat->half_patch_size_when_initialized = 20;
    testFeat->uv_when_initialized[0] = 74;
    testFeat->uv_when_initialized[1] = 43;
    if (testFeat->R_wc_when_initialized==NULL)
        testFeat->R_wc_when_initialized = (double*)malloc(9*sizeof(double));
    testFeat->R_wc_when_initialized[0] = 1.0;
    testFeat->R_wc_when_initialized[1] = 0.0;
    testFeat->R_wc_when_initialized[2] = 0.0;
    testFeat->R_wc_when_initialized[3] = 0.0;
    testFeat->R_wc_when_initialized[4] = 1.0;
    testFeat->R_wc_when_initialized[5] = 0.0;
    testFeat->R_wc_when_initialized[6] = 0.0;
    testFeat->R_wc_when_initialized[7] = 0.0;
    testFeat->R_wc_when_initialized[8] = 1.0;
    if (testFeat->r_wc_when_initialized==NULL)
        testFeat->r_wc_when_initialized = (double*)malloc(3*sizeof(double));
    testFeat->r_wc_when_initialized[0] = 0.0;
    testFeat->r_wc_when_initialized[1] = 0.0;
    testFeat->r_wc_when_initialized[2] = 0.0;
    
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
    
    double* patchVal = (double*)malloc(sizeof(double));
    uchar* initializedPtr;
    double* matchingPtr;
    
    NSString* path = [[NSBundle mainBundle] pathForResource:@"pfa-001_pwi"
                                                     ofType:@"txt"];
    FILE* fp_pwi = fopen((char*)[path UTF8String], "r");
    initializedPtr = testFeat->patch_when_initialized->ptr<uchar>(0);
    for (int i=0; i<(41*41); ++i){
        fscanf(fp_pwi, "%lf ", patchVal);
        initializedPtr[i] = (uchar)((int)(*patchVal));
    }
    fclose(fp_pwi);
    
    path = [[NSBundle mainBundle] pathForResource:@"pfa-001_pwm_before"
                                           ofType:@"txt"];
    FILE* fp_pwm = fopen((char*)[path UTF8String], "r");
    matchingPtr = testFeat->patch_when_matching->ptr<double>(0);
    for (int i=0; i<(13*13); ++i){
        fscanf(fp_pwm, "%lf ", matchingPtr+i);
        //matchingPtr[i] = (uchar)((int)round(*patchVal));
    }
    fclose(fp_pwm);
    
    double* patch_when_matching_expected = (double*)malloc(13*13*sizeof(double));
    path = [[NSBundle mainBundle] pathForResource:@"pfa-001_pwm_after"
                                           ofType:@"txt"];
    fp_pwm = fopen((char*)[path UTF8String], "r");
    for (int i=0; i<(13*13); ++i){
        fscanf(fp_pwm, "%lf ", patch_when_matching_expected+i);
        //patch_when_matching_expected[i] = (uchar)((int)round(*patchVal));
    }
    fclose(fp_pwm);
    
    int stateSize = 127;
    Filter* filter = [[Filter alloc]init];
    filter->state_size = stateSize;
    filter->x_k_k = NULL;
    filter->x_k_km1 = (double*)malloc(stateSize*sizeof(double));
    
    path = [[NSBundle mainBundle] pathForResource:@"pfa-001_x"
                                           ofType:@"txt"];
    FILE* fp_x = fopen((char*)[path UTF8String], "r");
    for (int i=0; i<(stateSize); ++i){
        fscanf(fp_x, "%lf ", patchVal);
        filter->x_k_km1[i] = *patchVal;
    }
    fclose(fp_x);
    
    bool success = predict_features_appearance(features_info, filter, cam);
    XCTAssert(success==true, @"expected return value of true for valid data");
    
    matchingPtr = ((Feature*)features_info[4])->patch_when_matching->ptr<double>(0);
    double numErrors = 0;
    for (int i=0; i<(13*13); ++i){
        double actual = matchingPtr[i];
        double expected = patch_when_matching_expected[i];
        double diff = fabs(actual-expected);
        if (diff>MATCHING_ERROR_THRESHOLD){
            numErrors = numErrors + 1.0;
        }
    }
    
    double accuracy = 1-numErrors/(13*13);
    if (accuracy<MATCHING_ACCURACY_THRESHOLD){
        XCTAssert(NO, @"Accuracy is not above threshold");
    }
    NSLog(@"---testPredictFeaturesAppearance001 accuracy %lf---", accuracy);
    
    free(patch_when_matching_expected);
    for (int i=0; i<[features_info count]; ++i){
        [(Feature*)features_info[i] destroy];
    }
    [filter destroy];
    [cam destroy];
    
}
/*****************/
//Second one
-(void)testPredictFeaturesAppearance002{
    NSMutableArray* features_info = [[NSMutableArray alloc]init];
    int numFeatures = 25;
    int featIdx = 4;
    for (int i=0; i<numFeatures; ++i){
        Feature* feat = [[Feature alloc]init];
        feat->type = [NSString stringWithFormat:@"inversedepth"];
        [features_info addObject:feat];
    }
    Feature* testFeat = (Feature*)features_info[featIdx];
    if (testFeat->h==NULL)
        testFeat->h = (double*)malloc(2*sizeof(double));
    testFeat->h[0] = 59.370934808257715;
    testFeat->h[1] = 47.947155033786501;
    testFeat->half_patch_size_when_matching = 6;
    testFeat->half_patch_size_when_initialized = 20;
    
    testFeat->uv_when_initialized[0] = 74;
    testFeat->uv_when_initialized[1] = 43;
    if (testFeat->R_wc_when_initialized==NULL)
        testFeat->R_wc_when_initialized = (double*)malloc(9*sizeof(double));
    testFeat->R_wc_when_initialized[0] = 1.0;
    testFeat->R_wc_when_initialized[1] = 0.0;
    testFeat->R_wc_when_initialized[2] = 0.0;
    testFeat->R_wc_when_initialized[3] = 0.0;
    testFeat->R_wc_when_initialized[4] = 1.0;
    testFeat->R_wc_when_initialized[5] = 0.0;
    testFeat->R_wc_when_initialized[6] = 0.0;
    testFeat->R_wc_when_initialized[7] = 0.0;
    testFeat->R_wc_when_initialized[8] = 1.0;
    if (testFeat->r_wc_when_initialized==NULL)
        testFeat->r_wc_when_initialized = (double*)malloc(3*sizeof(double));
    testFeat->r_wc_when_initialized[0] = 0.0;
    testFeat->r_wc_when_initialized[1] = 0.0;
    testFeat->r_wc_when_initialized[2] = 0.0;
    
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
    
    double* patchVal = (double*)malloc(sizeof(double));
    uchar* initializedPtr;
    double* matchingPtr;
    
    NSString* path = [[NSBundle mainBundle] pathForResource:@"pfa-002_pwi"
                                                     ofType:@"txt"];
    FILE* fp_pwi = fopen((char*)[path UTF8String], "r");
    initializedPtr = testFeat->patch_when_initialized->ptr<uchar>(0);
    for (int i=0; i<(41*41); ++i){
        fscanf(fp_pwi, "%lf ", patchVal);
        initializedPtr[i] = (uchar)((int)(*patchVal));
    }
    fclose(fp_pwi);
    
    path = [[NSBundle mainBundle] pathForResource:@"pfa-002_pwm_before"
                                           ofType:@"txt"];
    FILE* fp_pwm = fopen((char*)[path UTF8String], "r");
    matchingPtr = testFeat->patch_when_matching->ptr<double>(0);
    for (int i=0; i<(13*13); ++i){
        fscanf(fp_pwm, "%lf ", matchingPtr+i);
        //mxPtr2[i] = (uchar)((int)round(*patchVal));
    }
    fclose(fp_pwm);
    
    double* patch_when_matching_expected = (double*)malloc(13*13*sizeof(double));
    path = [[NSBundle mainBundle] pathForResource:@"pfa-002_pwm_after"
                                           ofType:@"txt"];
    fp_pwm = fopen((char*)[path UTF8String], "r");
    for (int i=0; i<(13*13); ++i){
        fscanf(fp_pwm, "%lf ", patch_when_matching_expected+i);
        //patch_when_matching_expected[i] = (uchar)((int)round(*patchVal));
    }
    fclose(fp_pwm);
    
    int stateSize = 163;
    Filter* filter = [[Filter alloc]init];
    filter->state_size = filter->state_km1_size = stateSize;
    filter->x_k_k = NULL;
    filter->x_k_km1 = (double*)malloc(stateSize*sizeof(double));
    
    path = [[NSBundle mainBundle] pathForResource:@"pfa-002_x"
                                           ofType:@"txt"];
    FILE* fp_x = fopen((char*)[path UTF8String], "r");
    for (int i=0; i<(stateSize); ++i){
        fscanf(fp_x, "%lf ", filter->x_k_km1+i);
        //filter->x_k_k[i] = *patchVal;
    }
    fclose(fp_x);
    
    bool success = predict_features_appearance(features_info, filter, cam);
    XCTAssert(success==true, @"expected return value of true for valid data");
    
    matchingPtr = ((Feature*)features_info[featIdx])->patch_when_matching->ptr<double>(0);
    double numErrors = 0;
    for (int i=0; i<(13*13); ++i){
        double actual = matchingPtr[i];
        double expected = patch_when_matching_expected[i];
        double diff = fabs(actual-expected);
        if (diff>MATCHING_ERROR_THRESHOLD){
            numErrors = numErrors + 1.0;
        }
    }
    
    double accuracy = 1-numErrors/(13*13);
    if (accuracy<MATCHING_ACCURACY_THRESHOLD){
        XCTAssert(NO, @"Accuracy is not above threshold");
    }
    NSLog(@"---testPredictFeaturesAppearance002 accuracy %lf---", accuracy);
    
    free(patch_when_matching_expected);
    for (int i=0; i<[features_info count]; ++i){
        [(Feature*)features_info[i] destroy];
    }
    [filter destroy];
    [cam destroy];
    
}

/*************************************PRED_PATCH_FC********************************************/


-(void)testPredPatchFC001{
    Feature* feat = [[Feature alloc]init];
    if (feat->h==NULL)
        feat->h = (double*)malloc(2*sizeof(double));
    feat->h[0] = 32.0;
    feat->h[1] = 162.0;
    feat->half_patch_size_when_matching = 6;
    feat->half_patch_size_when_initialized = 20;
    feat->uv_when_initialized[0] = 32;
    feat->uv_when_initialized[1] = 162;
    if (feat->R_wc_when_initialized==NULL)
        feat->R_wc_when_initialized = (double*)malloc(9*sizeof(double));
    feat->R_wc_when_initialized[0] = 1.0;
    feat->R_wc_when_initialized[1] = 0.0;
    feat->R_wc_when_initialized[2] = 0.0;
    feat->R_wc_when_initialized[3] = 0.0;
    feat->R_wc_when_initialized[4] = 1.0;
    feat->R_wc_when_initialized[5] = 0.0;
    feat->R_wc_when_initialized[6] = 0.0;
    feat->R_wc_when_initialized[7] = 0.0;
    feat->R_wc_when_initialized[8] = 1.0;
    if (feat->r_wc_when_initialized==NULL)
        feat->r_wc_when_initialized = (double*)malloc(3*sizeof(double));
    feat->r_wc_when_initialized[0] = 0.0;
    feat->r_wc_when_initialized[1] = 0.0;
    feat->r_wc_when_initialized[2] = 0.0;
    
    double* patchVal = (double*)malloc(sizeof(double));
    uchar* mxPtr;
    
    NSString* path = [[NSBundle mainBundle] pathForResource:@"ppfc-001_pwi"
                                           ofType:@"txt"];
    FILE* fp_pwi = fopen((char*)[path UTF8String], "r");
    mxPtr = feat->patch_when_initialized->ptr<uchar>(0);
    for (int i=0; i<(41*41); ++i){
        fscanf(fp_pwi, "%lf ", patchVal);
        mxPtr[i] = (uchar)((int)round(*patchVal));
    }
    fclose(fp_pwi);
    
    path = [[NSBundle mainBundle] pathForResource:@"ppfc-001_pwm"
                                           ofType:@"txt"];
    FILE* fp_pwm = fopen((char*)[path UTF8String], "r");
    double* mxPtr2 = feat->patch_when_matching->ptr<double>(0);
    for (int i=0; i<(13*13); ++i){
        fscanf(fp_pwm, "%lf ", mxPtr2+i);
    }
    fclose(fp_pwm);
    
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
    
    double* R_Wk = (double*)malloc(9*sizeof(double));
    R_Wk[0] = 1.000000000000000;
    R_Wk[1] = -0.000000000000001;
    R_Wk[2] = 0.000000000000001;
    R_Wk[3] = 0.000000000000001;
    R_Wk[4] = 1.000000000000000;
    R_Wk[5] = -0.000000000000001;
    R_Wk[6] = -0.000000000000001;
    R_Wk[7] = 0.000000000000001;
    R_Wk[8] = 1.000000000000000;
    
    double* r_Wk = (double*)malloc(3*sizeof(double));
    r_Wk[0] = 0.0;
    r_Wk[1] = 0.0;
    r_Wk[2] = 0.0;
    
    double* XYZ_w = (double*)malloc(3*sizeof(double));
    XYZ_w[0] = -0.615448111983574;
    XYZ_w[1] = 0.159036831945620;
    XYZ_w[2] = 0.771965612926220;
    
    double* patch_pred_expected = (double*)malloc(13*13*sizeof(double));
    path = [[NSBundle mainBundle] pathForResource:@"ppfc-001_patch_pred"
                                           ofType:@"txt"];
    FILE* fp_pp = fopen((char*)[path UTF8String], "r");
    for (int i=0; i<(13*13); ++i){
        fscanf(fp_pp, "%lf ", patch_pred_expected+i);
        //patch_pred_expected[i] = (uchar)((int)(*patchVal));
    }
    fclose(fp_pp);
    
    cv::Mat* patch_pred_actual = new cv::Mat(13, 13, CV_64FC1, 255);
    
    bool success = pred_patch_fc(cam, feat, R_Wk, r_Wk, XYZ_w, patch_pred_actual);
    
    XCTAssert(success==true, @"expected return value true for valid data");
    mxPtr2 = patch_pred_actual->ptr<double>(0);
    
    double errorCounter = 0.0;
    for (int i=0; i<(13*13); ++i){
        double expected = (patch_pred_expected[i]);
        double actual = (mxPtr2[i]);
        double diff = fabs(expected-actual);
        if (diff>MATCHING_ERROR_THRESHOLD){
            errorCounter = errorCounter+1;
        }
    }
    double successRate = 1-errorCounter/(13*13);
    NSLog(@"---testPredPatchFC001 accuracy %lf---", successRate);
    if (successRate<MATCHING_ACCURACY_THRESHOLD){
        XCTAssert(NO, @"Error rate too high");
    }
    
    patch_pred_actual->release();
    free(patch_pred_expected);
    free(r_Wk);
    free(R_Wk);
    [cam destroy];
    [feat destroy];
}

//Better case:
-(void)testPredPatchFC002{
    Feature* feat = [[Feature alloc]init];
    if (feat->h==NULL)
        feat->h = (double*)malloc(2*sizeof(double));
    feat->h[0] = 100*0.349723184912363;
    feat->h[1] = 100*1.349633651567365;
    feat->half_patch_size_when_matching = 6;
    feat->half_patch_size_when_initialized = 20;
    feat->uv_when_initialized[0] = 53;
    feat->uv_when_initialized[1] = 127;
    if (feat->R_wc_when_initialized==NULL)
        feat->R_wc_when_initialized = (double*)malloc(9*sizeof(double));
    feat->R_wc_when_initialized[0] = 1.0;
    feat->R_wc_when_initialized[1] = 0.0;
    feat->R_wc_when_initialized[2] = 0.0;
    feat->R_wc_when_initialized[3] = 0.0;
    feat->R_wc_when_initialized[4] = 1.0;
    feat->R_wc_when_initialized[5] = 0.0;
    feat->R_wc_when_initialized[6] = 0.0;
    feat->R_wc_when_initialized[7] = 0.0;
    feat->R_wc_when_initialized[8] = 1.0;
    if (feat->r_wc_when_initialized==NULL)
        feat->r_wc_when_initialized = (double*)malloc(3*sizeof(double));
    feat->r_wc_when_initialized[0] = 0.0;
    feat->r_wc_when_initialized[1] = 0.0;
    feat->r_wc_when_initialized[2] = 0.0;
    
    double* patchVal = (double*)malloc(sizeof(double));
    
    NSString* path = [[NSBundle mainBundle] pathForResource:@"ppfc-002_pwi"
                                                     ofType:@"txt"];
    FILE* fp_pwi = fopen((char*)[path UTF8String], "r");
    uchar* mxPtr = feat->patch_when_initialized->ptr<uchar>(0);
    for (int i=0; i<(41*41); ++i){
        fscanf(fp_pwi, "%lf ", patchVal);
        mxPtr[i] = (uchar)((int)(*patchVal));
    }
    fclose(fp_pwi);
    
    path = [[NSBundle mainBundle] pathForResource:@"ppfc-002_pwm"
                                           ofType:@"txt"];
    FILE* fp_pwm = fopen((char*)[path UTF8String], "r");
    double* mxPtr2 = feat->patch_when_matching->ptr<double>(0);
    for (int i=0; i<(13*13); ++i){
        fscanf(fp_pwm, "%lf ", mxPtr2+i);
        //mxPtr[i] = (uchar)((int)(*patchVal));
    }
    fclose(fp_pwm);
    
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
    
    double* R_Wk = (double*)malloc(9*sizeof(double));
    R_Wk[0] = 0.998150681654144;
    R_Wk[1] = -0.022307661159424;
    R_Wk[2] = 0.056547192387989;
    R_Wk[3] = 0.023040376655367;
    R_Wk[4] = 0.999658388621707;
    R_Wk[5] = -0.012338845238982;
    R_Wk[6] = -0.056252624444969;
    R_Wk[7] = 0.013618895397537;
    R_Wk[8] = 0.998323678939453;
    
    double* r_Wk = (double*)malloc(3*sizeof(double));
    r_Wk[0] = 0.051969751686999;
    r_Wk[1] = -0.013018922325676;
    r_Wk[2] = -0.013114667021295;
    
    double* XYZ_w = (double*)malloc(3*sizeof(double));
    XYZ_w[0] = -0.486895351219034;
    XYZ_w[1] = -0.007537476386865;
    XYZ_w[2] = 0.784390873223130;
    
    double* patch_pred_expected = (double*)malloc(13*13*sizeof(double));
    path = [[NSBundle mainBundle] pathForResource:@"ppfc-002_patch_pred"
                                           ofType:@"txt"];
    FILE* fp_pp = fopen((char*)[path UTF8String], "r");
    for (int i=0; i<(13*13); ++i){
        fscanf(fp_pp, "%lf ", patch_pred_expected+i);
       // patch_pred_expected[i] = (uchar)((int)(*patchVal));
    }
    fclose(fp_pp);
    
    cv::Mat* patch_pred_actual = new cv::Mat(13, 13, CV_64FC1, 255);
    
    bool success = pred_patch_fc(cam, feat, R_Wk, r_Wk, XYZ_w, patch_pred_actual);
    
    XCTAssert(success==true, @"expected return value true for valid data");
    mxPtr2 = patch_pred_actual->ptr<double>(0);
    
    double errorCounter = 0.0;
    for (int i=0; i<(13*13); ++i){
        double expected = (patch_pred_expected[i]);
        double actual = (mxPtr2[i]);
        int diff = abs(expected-actual);
        if (diff>MATCHING_ERROR_THRESHOLD){
            errorCounter = errorCounter+1;
        }
    }
    double successRate = 1-errorCounter/(13*13);
    NSLog(@"---testPredPatchFC002 accuracy %lf---", successRate);
    if (successRate<MATCHING_ACCURACY_THRESHOLD){
        XCTAssert(NO, @"Error rate too high");
    }
    
    patch_pred_actual->release();
    free(patch_pred_expected);
    free(r_Wk);
    free(R_Wk);
    [cam destroy];
    [feat destroy];
}

/*************************************ROTATE_WITH_DIST_c2c1********************************************/

-(void)testRotateWithDistc2c1001{
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
    
    double* uv_c1 = (double*)malloc(2*sizeof(double));
    uv_c1[0] = 93;
    uv_c1[1] = 153;
    
    double* R_c2c1 = (double*)malloc(9*sizeof(double));
    R_c2c1[0] = 1.000000000000000;
    R_c2c1[1] = -0.000000000000001;
    R_c2c1[2] = 0.000000000000001;
    R_c2c1[3] = 0.000000000000001;
    R_c2c1[4] = 1.000000000000000;
    R_c2c1[5] = -0.000000000000001;
    R_c2c1[6] = -0.000000000000001;
    R_c2c1[7] = 0.000000000000001;
    R_c2c1[8] = 1.000000000000000;
    
    double* t_c2c1 = (double*)malloc(3*sizeof(double));
    t_c2c1[0] = 0.0;
    t_c2c1[1] = 0.0;
    t_c2c1[2] = 0.0;
    
    double* n = (double*)malloc(3*sizeof(double));
    n[0] = -0.325081302808470;
    n[1] = 0.116708030480980;
    n[2] = -0.938456915465808;
    
    double d = 0.751660634838213;
    
    double* uv_c2_expected = (double*)malloc(2*sizeof(double));
    double* uv_c2_actual = (double*)malloc(2*sizeof(double));
    uv_c2_expected[0] = 100*0.929999999999998;
    uv_c2_expected[1] = 100*1.530000000000003;
    
    bool success = rotate_with_dist_fc_c2c1(cam, uv_c1, R_c2c1, t_c2c1, n, d, uv_c2_actual);
    
    double actual = uv_c2_actual[0];
    double expected = uv_c2_expected[0];
    double diff = fabs(actual-expected);
    XCTAssert(diff<0.00001, @"uv_c2 error exceeds limits");
    actual = uv_c2_actual[1];
    expected = uv_c2_expected[1];
    diff = fabs(actual-expected);
    XCTAssert(diff<0.00001, @"uv_c2 error exceeds limits");
    
    free(uv_c2_actual);
    free(uv_c2_expected);
    [cam destroy];
    free(R_c2c1);
    free(t_c2c1);
    free(n);
    free(uv_c1);
}

-(void)testRotateWithDistc2c1002{
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
    
    double* uv_c1 = (double*)malloc(2*sizeof(double));
    uv_c1[0] = 53;
    uv_c1[1] = 127;
    
    double* R_c2c1 = (double*)malloc(9*sizeof(double));
    R_c2c1[0] = 0.999780382872853;
    R_c2c1[1] = -0.001821696588341;
    R_c2c1[2] = 0.020877438639632;
    R_c2c1[3] = 0.002087600154716;
    R_c2c1[4] = 0.999916896701652;
    R_c2c1[5] = -0.012721698633884;
    R_c2c1[6] = -0.020852528580621;
    R_c2c1[7] = 0.012762488475112;
    R_c2c1[8] = 0.999701100799492;
    
    double* t_c2c1 = (double*)malloc(3*sizeof(double));
    t_c2c1[0] = 0.018297059067958;
    t_c2c1[1] = 0.003472111748255;
    t_c2c1[2] = 0.001025677445689;
    
    double* n = (double*)malloc(3*sizeof(double));
    n[0] =  -0.504853723024900;
    n[1] =  0.000548454960802;
    n[2] =  -0.863204736748504;
    
    double d = 0.467110359542819;
    
    double* uv_c2_expected = (double*)malloc(2*sizeof(double));
    double* uv_c2_actual = (double*)malloc(2*sizeof(double));
    uv_c2_expected[0] = 100*0.462293777351376;
    uv_c2_expected[1] = 100*1.287443264232931;
    
    bool success = rotate_with_dist_fc_c2c1(cam, uv_c1, R_c2c1, t_c2c1, n, d, uv_c2_actual);
    XCTAssert(success==true, @"No success");
    double actual = uv_c2_actual[0];
    double expected = uv_c2_expected[0];
    double diff = fabs(actual-expected);
    XCTAssert(diff<0.00000001, @"uv_c2 error exceeds limits");
    actual = uv_c2_actual[1];
    expected = uv_c2_expected[1];
    diff = fabs(actual-expected);
    XCTAssert(diff<0.00000001, @"uv_c2 error exceeds limits");
    
    free(uv_c2_actual);
    free(uv_c2_expected);
    [cam destroy];
    free(R_c2c1);
    free(t_c2c1);
    free(n);
    free(uv_c1);
}



/*************************************ROTATE_WITH_DIST_c1c2********************************************/

-(void)testRotateWithDistc1c2001{
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
    
    int uv_c2_length = 169;
    
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
    
    double* R_c1c2 = (double*)malloc(9*sizeof(double));
    R_c1c2[0] = 1.000000000000000;
    R_c1c2[1] = -0.000000000000001;
    R_c1c2[2] = 0.000000000000001;
    R_c1c2[3] = 1.000000000000000;
    R_c1c2[4] = -0.000000000000001;
    R_c1c2[5] = 0.000000000000001;
    R_c1c2[6] =  0.000000000000001;
    R_c1c2[7] = -0.000000000000001;
    R_c1c2[8] = 1.000000000000000;
    
    double* t_c1c2 = (double*)malloc(3*sizeof(double));
    t_c1c2[0] = 0.f;
    t_c1c2[1] = 0.f;
    t_c1c2[2] = 0.f;
    
    double* n = (double*)malloc(3*sizeof(double));
    n[0] = -0.483593482113592;
    n[1] = -0.008416274274438;
    n[2] = -0.875252255286777;
    
    double d = 0.489599986784764;
    
    double* uv_c2 = (double*)malloc(uv_c2_length*2*sizeof(double));
    NSString* path = [[NSBundle mainBundle] pathForResource:@"rwdc1c2-001_uv_c2"
                                                     ofType:@"txt"];
    FILE* fp_uv = fopen((char*)[path UTF8String], "r");
    path = nil;
    
    for (int i=0; i<(2*uv_c2_length); ++i){
        fscanf(fp_uv, "%lf ", uv_c2+i);
    }
    fclose(fp_uv);
    
    double* uv_c1_actual = (double*)malloc(uv_c2_length*2*sizeof(double));
    double* uv_c1_expected = (double*)malloc(uv_c2_length*2*sizeof(double));
    path = [[NSBundle mainBundle] pathForResource:@"rwdc1c2-001_uv_c1"
                                           ofType:@"txt"];
    fp_uv = fopen((char*)[path UTF8String], "r");
    for (int i=0; i<(2*uv_c2_length); ++i){
        uv_c1_actual[i] = -1.0;
        fscanf(fp_uv, "%lf ", uv_c1_expected+i);
    }
    fclose(fp_uv);
    
    bool success  = rotate_with_dist_fc_c1c2(cam, uv_c2, uv_c2_length, R_c1c2, t_c1c2, n, d, uv_c1_actual);
    XCTAssert(success==true, @"Expected return value of true for valid data");
    
    for (int i=0; i<(2*uv_c2_length); ++i){
        double actual = uv_c1_actual[i];
        double expected = uv_c1_expected[i];
        double diff = fabs(actual-expected);
        if (diff>0.00001){
            XCTAssert(NO, @"uv_c1 error exceeds limits");
            break;
        }
    }
    
    free(uv_c1_actual);
    free(uv_c1_expected);
    free(uv_c2);
    free(R_c1c2);
    free(n);
    free(t_c1c2);
}

//Nonzero t values
-(void)testRotateWithDistc1c2002{
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
    
    int uv_c2_length = 169;
    
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
    
    double* R_c1c2 = (double*)malloc(9*sizeof(double));
    R_c1c2[0] = 0.999780382872853;
    R_c1c2[1] = -0.001821696588341;
    R_c1c2[2] = 0.020877438639632;
    R_c1c2[3] = 0.002087600154716;
    R_c1c2[4] = 0.999916896701652;
    R_c1c2[5] = -0.012721698633884;
    R_c1c2[6] = -0.020852528580621;
    R_c1c2[7] = 0.012762488475112;
    R_c1c2[8] = 0.999701100799492;
    
    double* t_c1c2 = (double*)malloc(3*sizeof(double));
    t_c1c2[0] = 0.018297059067958;
    t_c1c2[1] = 0.003472111748255;
    t_c1c2[2] = 0.001025677445689;
    
    double* n = (double*)malloc(3*sizeof(double));
    n[0] = -0.504853723024900;
    n[1] = 0.000548454960802;
    n[2] = -0.863204736748504;
    
    double d = 0.467110359542819;
    
    double* uv_c2 = (double*)malloc(uv_c2_length*2*sizeof(double));
    NSString* path = [[NSBundle mainBundle] pathForResource:@"rwdc1c2-002_uv_c2"
                                                     ofType:@"txt"];
    FILE* fp_uv = fopen((char*)[path UTF8String], "r");
    path = nil;
    
    for (int i=0; i<(2*uv_c2_length); ++i){
        fscanf(fp_uv, "%lf ", uv_c2+i);
    }
    fclose(fp_uv);
    
    double* uv_c1_actual = (double*)malloc(uv_c2_length*2*sizeof(double));
    double* uv_c1_expected = (double*)malloc(uv_c2_length*2*sizeof(double));
    path = [[NSBundle mainBundle] pathForResource:@"rwdc1c2-002_uv_c1"
                                           ofType:@"txt"];
    fp_uv = fopen((char*)[path UTF8String], "r");
    for (int i=0; i<(2*uv_c2_length); ++i){
        uv_c1_actual[i] = -1.0;
        fscanf(fp_uv, "%lf ", uv_c1_expected+i);
    }
    fclose(fp_uv);
    
    bool success  = rotate_with_dist_fc_c1c2(cam, uv_c2, uv_c2_length, R_c1c2, t_c1c2, n, d, uv_c1_actual);
    XCTAssert(success==true, @"Expected return value of true for valid data");
    
    for (int i=0; i<(2*uv_c2_length); ++i){
        double actual = uv_c1_actual[i];
        double expected = uv_c1_expected[i];
        double diff = fabs(actual-expected);
        if (diff>0.00001){
            XCTAssert(NO, @"uv_c1 error exceeds limits");
            break;
        }
    }
    
    free(uv_c1_actual);
    free(uv_c1_expected);
    free(uv_c2);
    free(R_c1c2);
    free(n);
    free(t_c1c2);
}

/*********************************INTERP_PATCH_PRED***************************************/

//The best we can do is test for a certain accuracy
//Here, we test that 90% of interpolated pixels are within 10 intensity of matlab value
-(void)testInterpPatchPred001{
    uchar* patch_p_f = (uchar*)malloc(41*41*sizeof(uchar));
    double* u_pred_imak_dist = (double*)malloc(13*13*sizeof(double));
    double* v_pred_imak_dist = (double*)malloc(13*13*sizeof(double));
    double* patch_pred_actual = (double*)malloc(13*13*sizeof(double));
    double* patch_pred_expected = (double*)malloc(13*13*sizeof(double));
    
    NSString* path = [[NSBundle mainBundle] pathForResource:@"ipp-001_patch_p_f" ofType:@"txt"];
    FILE* fp_patch_pf = fopen((char*)[path UTF8String], "r");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"ipp-001_patch_pred" ofType:@"txt"];
    FILE* fp_patch_pred = fopen((char*)[path UTF8String], "r");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"ipp-001_u_pred_imak_dist" ofType:@"txt"];
    FILE* fp_up = fopen((char*)[path UTF8String], "r");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"ipp-001_v_pred_imak_dist" ofType:@"txt"];
    FILE* fp_vp = fopen((char*)[path UTF8String], "r");
    path = nil;
    
    double* temp = (double*)malloc(sizeof(double));
    for (int i=0; i<(13*13); ++i){
        fscanf(fp_up, "%lf ", u_pred_imak_dist+i);
        fscanf(fp_vp, "%lf ", v_pred_imak_dist+i);
        fscanf(fp_patch_pred, "%lf ", patch_pred_expected+i);
        patch_pred_actual[i] = 0;
    }
    for (int i=0; i<(41*41); ++i){
        fscanf(fp_patch_pf, "%lf ", temp);
        patch_p_f[i] = (unsigned char)((int)*temp);
    }
    fclose(fp_patch_pf);
    fclose(fp_patch_pred);
    fclose(fp_vp);
    fclose(fp_up);
    
    bool success = interpPatchPred(patch_p_f, u_pred_imak_dist, v_pred_imak_dist, patch_pred_actual);
    XCTAssert(success==true, @"Expected return value of true for valid data");
    
    double missCounter = 0;
    
    int maxDist = 10;
    
    for (int i=0; i<(13*13); ++i){
        double actual = patch_pred_actual[i];
        double expected = patch_pred_expected[i];
        double diff = abs( actual-expected );
        if (diff>MATCHING_ERROR_THRESHOLD){
            missCounter = missCounter+1.0;
        }
    }
    double percentAccuracy = 1-missCounter/(13*13);
    XCTAssert(percentAccuracy>MATCHING_ACCURACY_THRESHOLD, @"Interpolation producing too many errors");
    NSLog(@"--InterpPatchPred001 threshold %d accuracy %lf--", maxDist, percentAccuracy);
    free(patch_pred_expected);
    free(patch_pred_actual);
    free(patch_p_f);
    free(u_pred_imak_dist);
    free(v_pred_imak_dist);
    
}

-(void)testInterpPatchPred002{
    uchar* patch_p_f = (uchar*)malloc(41*41*sizeof(uchar));
    double* u_pred_imak_dist = (double*)malloc(13*13*sizeof(double));
    double* v_pred_imak_dist = (double*)malloc(13*13*sizeof(double));
    double* patch_pred_actual = (double*)malloc(13*13*sizeof(double));
    double* patch_pred_expected = (double*)malloc(13*13*sizeof(double));
    
    NSString* path = [[NSBundle mainBundle] pathForResource:@"ipp-002_patch_p_f" ofType:@"txt"];
    FILE* fp_patch_pf = fopen((char*)[path UTF8String], "r");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"ipp-002_patch_pred" ofType:@"txt"];
    FILE* fp_patch_pred = fopen((char*)[path UTF8String], "r");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"ipp-002_u_pred_imak_dist" ofType:@"txt"];
    FILE* fp_up = fopen((char*)[path UTF8String], "r");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"ipp-002_v_pred_imak_dist" ofType:@"txt"];
    FILE* fp_vp = fopen((char*)[path UTF8String], "r");
    path = nil;
    
    double* temp = (double*)malloc(sizeof(double));
    for (int i=0; i<(13*13); ++i){
        fscanf(fp_up, "%lf ", u_pred_imak_dist+i);
        fscanf(fp_vp, "%lf ", v_pred_imak_dist+i);
        fscanf(fp_patch_pred, "%lf ", patch_pred_expected+i);
        patch_pred_actual[i] = 0;
    }
    for (int i=0; i<(41*41); ++i){
        fscanf(fp_patch_pf, "%lf ", temp);
        patch_p_f[i] = (unsigned char)((int)*temp);
    }
    fclose(fp_patch_pf);
    fclose(fp_patch_pred);
    fclose(fp_vp);
    fclose(fp_up);
    
    bool success = interpPatchPred(patch_p_f, u_pred_imak_dist, v_pred_imak_dist, patch_pred_actual);
    XCTAssert(success==true, @"Expected return value of true for valid data");
    
    double missCounter = 0;
    
    int maxDist = 10;
    
    for (int i=0; i<(13*13); ++i){
        double actual = patch_pred_actual[i];
        double expected = patch_pred_expected[i];
        double diff = abs( actual-expected );
        if (diff>MATCHING_ERROR_THRESHOLD){
            missCounter = missCounter+1.0;
        }
    }
    double percentAccuracy = 1-missCounter/(13*13);
    XCTAssert(percentAccuracy>MATCHING_ACCURACY_THRESHOLD, @"Interpolation producing too many errors");
    NSLog(@"--InterpPatchPred002 threshold %d accuracy %lf--", maxDist, percentAccuracy);
    free(patch_pred_expected);
    free(patch_pred_actual);
    free(patch_p_f);
    free(u_pred_imak_dist);
    free(v_pred_imak_dist);
    
}


/*************************************MATCHING********************************************/

//Condition which creates null z vector
//TODO: Test the case in which the z vector is updated
-(void)testMatching001{
    NSMutableArray* features_info = [[NSMutableArray alloc]init];
    for (int i=0; i<19; ++i){
        Feature* currFeat = [[Feature alloc]init];
        currFeat->h = NULL;
        [features_info addObject:currFeat];
    }
    Feature* testFeat = (Feature*)features_info[3];
    if (testFeat->h==NULL){
        testFeat->h = (double*)malloc(2*sizeof(double));
    }

    testFeat->h[0] = 100*0.243349126967757;
    testFeat->h[1] = 100*1.649166495910018;
    
    if (testFeat->S==NULL){
        testFeat->S = (double*)malloc(4*sizeof(double));
    }
    testFeat->S[0] = 6.923750801999334;
    testFeat->S[1] = -0.221264224538376;
    testFeat->S[2] = -0.221264224538375;
    testFeat->S[3] = 5.987312609501809;
    testFeat->half_patch_size_when_matching = 6;
    
    NSString* path = [[NSBundle mainBundle] pathForResource:@"mat-001_pwm"
                                                     ofType:@"txt"];
    FILE* fp_pwm = fopen((char*)[path UTF8String], "r");
    double* pwmPtr = testFeat->patch_when_matching->ptr<double>(0);
    double* patchVal = (double*)malloc(sizeof(double));
    for (int i=0; i<(13*13); ++i){
        fscanf(fp_pwm, "%lf ", pwmPtr+i);
        //mxPtr[i] = (uchar)((int)(*patchVal));
    }
    fclose(fp_pwm);
    
    cv::Mat* im = new cv::Mat(240, 320, CV_8UC1, 0.0);
    uchar* imPtr = im->ptr<uchar>(0);
    path = [[NSBundle mainBundle] pathForResource:@"mat-001_im"
                                           ofType:@"txt"];
    FILE* fp_im = fopen((char*)[path UTF8String], "r");
    for (int i=0; i<(240*320); ++i){
        fscanf(fp_im, "%lf ", patchVal);
        imPtr[i] = (uchar)((int)(*patchVal));
    }
    fclose(fp_im);
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
    bool success = matching(im, features_info, cam);
    XCTAssert(success==true, @"expected true return");
    XCTAssert(((Feature*)features_info[3])->z==NULL, @"Expected this z to be empty");
    
    
    for (int i=0; i<[features_info count]; ++i){
        [(Feature*)features_info[i] destroy];
    }
    [cam destroy];
    im->release();
}


@end




