//
//  EKFilter.h
//  EKFMonoSLAM
//
//  Created by Shreenivaas Devarajanon 3/12/15.
//  Copyright (c) 2015 dryslam. All rights reserved.
//

#ifndef EKFILTER_H
#define EKFILTER_H

#include "SLAMModel.hpp"

//Return values: true: success / false: error(quit program)




//Generate 1000 random points in a sphere
//PreCondition: randSphere6D allocated but not initialized
//PostCondition: randSphere6D filled with values [X,Y,Z,theta,phi,lambda]
double* generate_random_6D_sphere(double* randSphere6D, int nPointsRand);























/*
//initialize_x_and_p
- (void) initialize_x_and_p; // Method to call methods to set the values for the matrix for State and Covariance

- (double *) initialize_x; // Method which initializes the value for x_k_k matrix
- (double *) initialize_p; // Method which initializes the values for p_k_k matrix

// Initialize the EKF

- (int) ekf_filter:(NSArray *) varargin, ...; // Initializes the variables for State and and Measurement Vectors

// Methods for updating the EKF

-(double) EKF_PredictionState:(NSArray *) xk_kcap
              Probabilitycond: (NSArray *) pk_kcap; // Code to predict the State vector to initialize the EKF

// Method to update the state and covariance

- (double *) predict_state_and_covariance: (NSArray *)fx_k_k
                               Covariance: (NSArray *)fp_k_K
                                     Type: (NSString *) type
                                     stda: (int) std_a
                            StandardAlpha: (int) std_alpha;
- (double *) fv: (double * ) x_k_k
         DeltaT: (int) delta_t
           Type:(NSString *) type
           Stda: (double) std_a
  StandardAlpha: (double) std_alpha;

- (double *) func_Q: (double *) x_k_k
 ZeroMatrixSIxByOn1: (double *) zeroessixbyone
           PnMatrix: (double *) pn
             DeltaT: (int) delta_t
               Type:(NSString *) type;

-(double *) dq_by_deuler: (double *) eulerangles;








- (double) EKF_Measurement; // Code to return the vector of measurements

- (double) EKF_getmeasurementvector;    // Method to return the Measurement Vector

- (double) EKF_StateUpdateMeasurement:(double *) z_i
                                State:(double *)xk_kcap; // Method to update the State Value of the EKF

- (double) EKF_Update;      // Updates the EKF with the Low and High Innovation Inliers at the end

- (void) update; //Update method to update the values of Low and High Innovation Inliers



- (void) dfv_by_dxv;

- (void) dq3_by_dq2;

- (void) dq3_by_dq1;

- (void) dqomegadt_by_domega;

- (void) dq0_by_domegaA;

- (void) normJac;




@end*/

#endif
