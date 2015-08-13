//
//  EKFilterTests.m
//  drySLAM
//
//  Created by Anthony Rodriguez on 3/27/15.
//  Copyright (c) 2015 ECE473 573 Common. All rights reserved.
//


#include "EKFilter.hpp"
#import <UIKit/UIKit.h>
#import <XCTest/XCTest.h>


@interface EKFilterTests : XCTestCase

@end

@implementation EKFilterTests


- (void)setUp {
    [super setUp];
    // Put setup code here. This method is called before the invocation of each test method in the class.
}

- (void)tearDown {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
    [super tearDown];
}




/***********************************GENERATERANDOM6DSPHERE*********************************/
//001: Should return false if input arguments are NULL or negative
-(void) testGenerateRandom6DSphere001{
    bool success = true;
    success = generate_random_6D_sphere(NULL, 100);
    XCTAssert(success==false, @"Error: generateRandom6DSphere should return false if array is NULL");
    double* sphere = (double*)malloc(10*sizeof(double));
    success = generate_random_6D_sphere(sphere, -1);
    XCTAssert(success == false, @"Error: generateRandom6DSphere should return false if size<0");
    free(sphere);
}




@end







