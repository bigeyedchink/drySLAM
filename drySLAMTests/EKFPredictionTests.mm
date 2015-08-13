//
//  EKFPredictionTests.m
//  drySLAM
//
//  Created by Anthony Rodriguez on 3/27/15.
//  Copyright (c) 2015 ECE473 573 Common. All rights reserved.
//


#include "SLAMModel.hpp"
#include "EKFPrediction.hpp"

#import <UIKit/UIKit.h>
#import <XCTest/XCTest.h>
#include <stdio.h>
#include <stdlib.h>




@interface EKFPredictionTests : XCTestCase

@end

@implementation EKFPredictionTests


- (void)setUp {
    [super setUp];
    // Put setup code here. This method is called before the invocation of each test method in the class.
}

- (void)tearDown {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
    [super tearDown];
}

/************************************************EKF_PREDICTION**************************************************/

-(void)testEKFPrediction000{
    //Implement a test here if we suspect that EKF prediction works improperly
    XCTAssert(YES, @"Function requires validation");
}

/******************************PREDICT_STATE_AND_COVARIANCE**********************************/
//test001 deleted

//Test 002: should return false if any values are null
-(void) testPredictStateAndCovariance002{
    int state_size = 85;
    double* x_k = NULL;
    double* p_k = (double*)malloc(state_size*state_size*sizeof(double));
    double* x_km1_k = (double*)malloc(state_size*sizeof(double));
    double* p_km1_k = (double*)malloc(state_size*state_size*sizeof(double));
    //x_k==NULL
    bool success = predict_state_and_covariance(NULL, p_k, 121, x_km1_k,  p_km1_k,
                                                @"constant_velocity", 0.0070, 0.0070);
    XCTAssert(success==false, @"Should return false when x_k==NULL");
    x_k = (double*)malloc(state_size*sizeof(double));
    free(p_k);
    //p_k==NULL
    success = predict_state_and_covariance(x_k, NULL, 121, x_km1_k,  p_km1_k,
                                           @"constant_velocity", 0.0070, 0.0070);
    XCTAssert(success==false, @"Should return false when p_k==NULL");
    p_k = (double*)malloc(state_size*state_size*sizeof(double));
    free(x_km1_k);
    //x_km1_k==NULL
    success = predict_state_and_covariance(x_k, p_k, 121, NULL,  p_km1_k,
                                           @"constant_velocity", 0.0070, 0.0070);
    XCTAssert(success==false, @"Should return false when x_km1_k==NULL");
    x_km1_k = (double*)malloc(state_size*sizeof(double));
    free(p_km1_k);
    //p_km1_k==NULL
    success = predict_state_and_covariance(x_k, p_k, 121, x_km1_k,  NULL,
                                           @"constant_velocity", 0.0070, 0.0070);
    XCTAssert(success==false, @"Should return false when p_km1_k==NULL");
    free(x_k);
    free(x_km1_k);
    free(p_k);
}

//Test 003: Another test from program data:
-(void)testPredictStateAndCovariance003{
    int state_size = 85;
    double* X_k = (double*)malloc(state_size*sizeof(double));
    double* P_k = (double*)malloc(state_size*state_size*sizeof(double));
    double* X_km1_k_actual = (double*)malloc(state_size*sizeof(double));
    double* P_km1_k_actual = (double*)malloc(state_size*state_size*sizeof(double));
    double* X_km1_k_expected = (double*)malloc(state_size*sizeof(double));
    double* P_km1_k_expected = (double*)malloc(state_size*state_size*sizeof(double));
    NSString* type = [NSString stringWithFormat:@"constant_velocity"];
    double filter_params = 0.0070;
    
    NSString* path = [[NSBundle mainBundle] pathForResource:@"psac-003_x_k"
                                                     ofType:@"txt"];
    FILE* fp_x_k = fopen((char*)[path UTF8String], "r");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"psac-003_p_k"
                                           ofType:@"txt"];
    FILE* fp_p_k = fopen((char*)[path UTF8String], "r");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"psac-003_x_km1_k_expected"
                                           ofType:@"txt"];
    FILE* fp_x_km1 = fopen((char*)[path UTF8String], "r");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"psac-003_p_km1_k_expected"
                                           ofType:@"txt"];
    FILE* fp_p_km1 = fopen((char*)[path UTF8String], "r");
    path = nil;
    
    for (int i=0; i<state_size; ++i){
        fscanf(fp_x_k, "%lf ", X_k+i);
        fscanf(fp_x_km1, "%lf ", X_km1_k_expected+i);
        X_km1_k_actual[i] = 0;
    }
    for (int i=0; i<state_size*state_size; ++i){
        fscanf(fp_p_k, "%lf ", P_k+i);
        fscanf(fp_p_km1, "%lf ", P_km1_k_expected+i);
        P_km1_k_actual[i] = 0;
    }
    fclose(fp_x_k);
    fclose(fp_x_km1);
    fclose(fp_p_k);
    fclose(fp_p_km1);
    
    bool success = predict_state_and_covariance(X_k, P_k, state_size, X_km1_k_actual, P_km1_k_actual, type, filter_params, filter_params);
    
    XCTAssert(success==true, @"Should return true for valid params");
    
    for (int i=0; i<state_size; ++i){
        double actual = X_km1_k_actual[i];
        double expected = X_km1_k_expected[i];
        double diff = fabs(actual-expected);
        if (diff>0.00001){
            XCTAssert(NO, @"X_km1_k error exceeds limits");
            break;
        }
    }
    for (int i=0; i<state_size*state_size; ++i){
        double actual = P_km1_k_actual[i];
        double expected = P_km1_k_expected[i];
        double diff = fabs(actual-expected);
        if (diff>0.00001){
            XCTAssert(NO, @"P_km1_k error exceeds limits");
            // break;
        }
    }
    
    free(X_km1_k_expected); free(X_km1_k_actual); free(P_km1_k_actual); free(P_km1_k_expected);
    free(X_k); free(P_k);
    
}

//Test 004: Another test from program data:
-(void)testPredictStateAndCovariance004{
    int state_size = 205;
    double* X_k = (double*)malloc(state_size*sizeof(double));
    double* P_k = (double*)malloc(state_size*state_size*sizeof(double));
    double* X_km1_k_actual = (double*)malloc(state_size*sizeof(double));
    double* P_km1_k_actual = (double*)malloc(state_size*state_size*sizeof(double));
    double* X_km1_k_expected = (double*)malloc(state_size*sizeof(double));
    double* P_km1_k_expected = (double*)malloc(state_size*state_size*sizeof(double));
    NSString* type = [NSString stringWithFormat:@"constant_velocity"];
    double filter_params = 0.0070;
    
    NSString* path = [[NSBundle mainBundle] pathForResource:@"psac-004_x_k"
                                                     ofType:@"txt"];
    FILE* fp_x_k = fopen((char*)[path UTF8String], "r");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"psac-004_p_k"
                                           ofType:@"txt"];
    FILE* fp_p_k = fopen((char*)[path UTF8String], "r");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"psac-004_x_km1_k_expected"
                                           ofType:@"txt"];
    FILE* fp_x_km1 = fopen((char*)[path UTF8String], "r");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"psac-004_p_km1_k_expected"
                                           ofType:@"txt"];
    FILE* fp_p_km1 = fopen((char*)[path UTF8String], "r");
    path = nil;
    
    for (int i=0; i<state_size; ++i){
        fscanf(fp_x_k, "%lf ", X_k+i);
        fscanf(fp_x_km1, "%lf ", X_km1_k_expected+i);
        X_km1_k_actual[i] = 0;
    }
    for (int i=0; i<state_size*state_size; ++i){
        fscanf(fp_p_k, "%lf ", P_k+i);
        fscanf(fp_p_km1, "%lf ", P_km1_k_expected+i);
        P_km1_k_actual[i] = 0;
    }
    fclose(fp_x_k);
    fclose(fp_x_km1);
    fclose(fp_p_k);
    fclose(fp_p_km1);
    
    bool success = predict_state_and_covariance(X_k, P_k, state_size, X_km1_k_actual, P_km1_k_actual, type, filter_params, filter_params);
    
    XCTAssert(success==true, @"Should return true for valid params");
    
    for (int i=0; i<state_size; ++i){
        double actual = X_km1_k_actual[i];
        double expected = X_km1_k_expected[i];
        double diff = fabs(actual-expected);
        if (diff>0.00001){
            XCTAssert(NO, @"X_km1_k error exceeds limits");
            break;
        }
    }
    for (int i=0; i<state_size*state_size; ++i){
        double actual = P_km1_k_actual[i];
        double expected = P_km1_k_expected[i];
        double diff = fabs(actual-expected);
        if (diff>0.00001){
            XCTAssert(NO, @"P_km1_k error exceeds limits");
            //break;
        }
    }
    
    free(X_km1_k_expected); free(X_km1_k_actual); free(P_km1_k_actual); free(P_km1_k_expected);
    free(X_k); free(P_k);
    
}

/*******************************************FV*********************************************/
//001: matlab inputs
-(void)testFv{
    double* xv = (double*)malloc(13*sizeof(double));
    double* x_km1_k_expected = (double*)malloc(13*sizeof(double));
    double* x_km1_k_actual = (double*)malloc(13*sizeof(double));
    xv[0] = 0.009139450566591;
    xv[1] = 0.001725790855408;
    xv[2] = 0.000681597720107;
    xv[3] = 0.999981198738380;
    xv[4] = 0.003185802909930;
    xv[5] = 0.005216703651362;
    xv[6] = 0.000488705009587;
    xv[7] = 0.009139450566591;
    xv[8] = 0.001725790855407;
    xv[9] = 0.000681597720107;
    xv[10] = 0.006371725606322;
    xv[11] = 0.010433603301728;
    xv[12] = 0.000977428243326;
    x_km1_k_expected[0] = 0.018278901133183;
    x_km1_k_expected[1] = 0.003451581710815;
    x_km1_k_expected[2] = 0.001363195440213;
    x_km1_k_expected[3] = 0.999924794718833;
    x_km1_k_expected[4] = 0.006371525949650;
    x_km1_k_expected[5] = 0.010433276442552;
    x_km1_k_expected[6] = 0.000977397691232;
    x_km1_k_expected[7] = 0.009139450566591;
    x_km1_k_expected[8] = 0.001725790855407;
    x_km1_k_expected[9] = 0.000681597720107;
    x_km1_k_expected[10] = 0.006371725606322;
    x_km1_k_expected[11] = 0.010433603301728;
    x_km1_k_expected[12] = 0.000977428243326;
    
    fv(x_km1_k_actual, xv, 1, [NSString stringWithFormat:@"constant_velocity"], 0.0070, 0.0070);
    
    for (int i=0; i<13; ++i){
        double actual = x_km1_k_actual[i];
        double expected = x_km1_k_expected[i];
        double diff = fabs(actual-expected);
        if (diff>0.00001){
            XCTAssert(NO, @"x error exceeds limits");
        }
    }
    free(xv);
    free(x_km1_k_expected);
    free(x_km1_k_actual);
}

/****************************************FUNC_Q*******************************************/
//001: Test from matlab state
-(void)testFuncQ001{
    //double delta_t = 1.0;
    NSString* type = [NSString stringWithFormat:@"constant_velocity"];
    double* u = (double*)malloc(6*sizeof(double));
    double* Pn = (double*)malloc(6*6*sizeof(double));
    double* Xv = (double*)malloc(13*sizeof(double));
    //double* Q_actual = (double*)malloc(13*13*sizeof(double));
    double* Q_expected = (double*)malloc(13*13*sizeof(double));
    
    NSString* path = [[NSBundle mainBundle] pathForResource:@"fq-001_q"
                                                     ofType:@"txt"];
    FILE* fp_q = fopen((char*)[path UTF8String], "r");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"fq-001_xv"
                                           ofType:@"txt"];
    FILE* fp_xv = fopen((char*)[path UTF8String], "r");
    path = nil;
    
    for (int i=0; i<6; ++i){
        u[i] = 0.0;
    }
    for (int i=0; i<6; ++i){
        for (int j=0; j<6; ++j){
            if (i==j)
                Pn[i*6+j] = 0.49 * 0.0001;
            else
                Pn[i*6+j] = 0.0;
        }
    }
    for (int i=0; i<13; ++i){
        fscanf(fp_xv, "%lf ", Xv+i);
    }
    fclose(fp_xv);
    for (int i=0; i<13*13; ++i){
        //Q_actual[i] = 0;
        fscanf(fp_q, "%lf ", Q_expected + i);
    }
    fclose(fp_q);
    
    double* Q_actual = func_Q(Xv, NULL, Pn, DELTA_T, type, NULL);
    
    XCTAssert(Q_actual!=NULL, @"funq_Q returned NULL!");
    
    for (int i=0; i<13*13; ++i){
        double actual = Q_actual[i];
        double expected = Q_expected[i];
        double diff = fabs(actual-expected);
        if (diff>0.00001){
            XCTAssert(NO, @"Q error exceeds limits");
            break;
        }
    }
    
    free(Pn);
    free(Xv);
    free(u);
    free(Q_actual);
    free(Q_expected);
}


/******************************************QPROD******************************************/
//001: Outputs from pgm
-(void)testQProd001{
    double* p = (double*)malloc(4*sizeof(double));
    double* q = (double*)malloc(4*sizeof(double));
    double* qp_expected = (double*)malloc(4*sizeof(double));
    double* qp_actual = (double*)malloc(4*sizeof(double));
    
    p[0] = 0.999987820425606;
    p[1] = 0.001707669656903;
    p[2] = 0.004383886762825;
    p[3] = 0.001491442804760;
    q[0] = 0.999949394169134;
    q[1] = 0.003610774945383;
    q[2] = 0.009110643525261;
    q[3] = 0.002273231099308;
    qp_expected[0] = 0.999887718776568;
    qp_expected[1] = 0.005321936622550;
    qp_expected[2] = 0.013492694137553;
    qp_expected[3] = 0.003764842000361;
    
    qp_actual = qprod(q, p, qp_actual);
    
    for (int i=0; i<4; ++i){
        double actual = qp_actual[i];
        double expected = qp_expected[i];
        double diff = fabs(actual-expected);
        if (diff>0.00001){
            XCTAssert(NO, @"qp error exceeds limits");
            break;
        }
    }
    free(qp_expected);
    free(qp_actual);
    free(q);
    free(p);
    
}


/***************************************DFV_BY_DXV*****************************************/
//001: Should return false if inputs aren't valid
-(void) testDfvByDxv001{
    double* dfvByDxv = (double*)malloc(13*13*sizeof(double));
    double* xv = (double*)malloc(13*sizeof(double));
    double* u = (double*)malloc(6*sizeof(double));
    bool success = dfv_by_dxv(NULL, xv, u, 1, @"constant_velocity");
    XCTAssert(success==false, @"dfvbydxv should immediately return false when dfvbydxv is NULL");
    success = dfv_by_dxv(dfvByDxv, NULL, u, 1, @"constant_velocity");
    XCTAssert(success==false, @"dfvbydxv should immediately return false when xv is NULL");
    success = dfv_by_dxv(dfvByDxv,  xv, NULL, 1, @"constant_velocity");
    XCTAssert(success==false, @"dfvbydxv should immediately return false when u is NULL");
    success = dfv_by_dxv(dfvByDxv,  xv, u, 1, @"It's a trap!");
    XCTAssert(success==false, @"dfvbydxv should immediately return false when type is invalid");
    free(dfvByDxv);
    free(xv);
    free(u);
}
//002: Step 104 of demo
-(void)testDfvByDxv002{
    double* dfvbydxv = (double*)malloc(13*13*sizeof(double));
    double* dfvbydxv_expected = (double*)malloc(13*13*sizeof(double));
    double* u = (double*)malloc(6*sizeof(double));
    double dt = 1.0;
    double* xv = (double*)malloc(13*sizeof(double));
    NSString* type = [NSString stringWithFormat:@"constant_velocity"];
    NSString* path = [[NSBundle mainBundle] pathForResource:@"dfvbdxv-002_xv"
                                                     ofType:@"txt"];
    FILE* fp_xv = fopen((char*)[path UTF8String], "r");
    XCTAssert(fp_xv!=NULL && path!=nil, @"Setup problem");
    path = nil;
    path = [[NSBundle mainBundle] pathForResource:@"dfvbdxv-002_dfv_by_dxvRES"
                                           ofType:@"txt"];
    FILE* fp_dfv_by_dxv_expected = fopen((char*)[path UTF8String], "r");
    XCTAssert(fp_dfv_by_dxv_expected!=NULL && path!=nil, @"Setup problem");
    path = nil;
    for (int i=0; i<6; ++i)
        u[i] = 0.0;
    for (int i=0; i<13; ++i){
        fscanf(fp_xv, "%lf ", (xv+i));
    }
    for (int i=0; i<(13*13); ++i){
        fscanf(fp_dfv_by_dxv_expected, "%lf ", (dfvbydxv_expected+i));
    }
    bool success = dfv_by_dxv(dfvbydxv, xv, u, dt, type);
    XCTAssert(success==true, @"dfv_by_dxv should return true for valid data");
    XCTAssert(dfvbydxv!=NULL, @"dfvbydxf memory lost");
    XCTAssert(xv!=NULL, @"xv memory lost");
    XCTAssert(u!=NULL, @"u memory lost");
    
    for(int i=0; i<(13*13); ++i){
        double actual = dfvbydxv[i];
        double expected = dfvbydxv_expected[i];
        double diff = fabs(actual-expected);
        if (diff>0.000001){
            XCTAssert(diff<0.000001, @"dfv_by_dxvRES error exceeds limits");
            // break;
        }
    }
    
    free(dfvbydxv);
    free(dfvbydxv_expected);
    free(u);
    free(xv);
}


/*************************************DQOMEGADT_BY_DOMEGA**********************************/
//001: Should return false if input arguments are incorrect (null,etc)
-(void) testDqomegaByDomega001{
    bool success = dqomegadt_by_domega(NULL, 0.01, 0.01, 0.01, 0.01);
    XCTAssert(success==false, @"Should return false when dqomegadt_by_domega is NULL");
}
//002: Data from demo step 91
-(void) testDqOmegaByDomega002{
    double* dqomegadt_by_domegaVar = (double*)malloc(12*sizeof(double));
    double* expected_answer = (double*)malloc(12*sizeof(double));
    expected_answer[0] = -0.0015929214182781;
    expected_answer[1] = -0.0026083844779267;
    expected_answer[2] = -0.0002443555293842;
    expected_answer[3] = 0.4999951747600563;
    expected_answer[4] = -0.0000027699919723;
    expected_answer[5] = -0.0000002594950478;
    expected_answer[6] = -0.0000027699919723;
    expected_answer[7] = 0.4999923305544742;
    expected_answer[8] = -0.0000004249191749;
    expected_answer[9] = -0.0000002594950478;
    expected_answer[10] = -0.0000004249191749;
    expected_answer[11] = 0.4999968265672295;
    bool success = dqomegadt_by_domega(dqomegadt_by_domegaVar, 0.006371725606322, 0.010433603301728,
                                       0.000977428243326, 1);
    
    XCTAssert(success==true, @"dqomegadt_by_domega should return true for real data!");
    for (int i=0; i<12; ++i){
        double diff = fabs(expected_answer[i]-dqomegadt_by_domegaVar[i]);
        if (diff>0.00000001){
            XCTAssert(NO, @"dqomegadt_by_domega error out of bounds!");
            break;
        }
    }
    free(dqomegadt_by_domegaVar);
    free(expected_answer);
}

@end



