//
//  EKFPrediction.hpp
//  drySLAM
//
//  Created by Anthony Rodriguez on 4/17/15.
//  Copyright (c) 2015 ECE473 573 Common. All rights reserved.
//

#ifndef EKFPREDICTION_HPP
#define EKFPREDICTION_HPP

#include "SLAMModel.hpp"


//Basically just call predict_state_and_covariance and update the appropriate matrices in filter
//Dependencies: predict_state_and_covariance
//Precondition: old state
//Postcondition: filter and features_info are updated with new values
bool ekf_prediction(Filter* filter, NSMutableArray* features_info);


//Updates state and covariance based on prediction, updates sizes appropriately
//PreCondition: x_km1_k, p_km1_k allocated but not set (sizes before and after are the same)
//PostCondition: x_km1_k, p_km1_k populated with values, size values set to sizes of data
bool predict_state_and_covariance(double* x_k, double* p_k, int state_size, double* x_km1_k, double* p_km1_k,
                                     NSString* type, double SD_A_component_filter, double SD_alpha_component_filter);


//Calculate portions of new state vector
//Precondition: xv is our 13-var state vector, xv_km1_k is our returned vector
//Postcondition: xv_km1_k is initialized properly
double* fv (double *xv_km1_k, double* xv, double delta_t, NSString* type, double std_a, double std_alphaa);


//Calculate difference in current time step (Q)to add to base state vector
//Precondition: Xv is the 13x1 base state vector, u is a blank control vector (6x1), Pn is a diagonal covariance matrix (see predict_state_and_covariance)
//Postcondition: Q is correctly calculated
double* func_Q(double * xv,double *func_Q_zero_matrix, double*  pn,double delta_t, NSString* type, double* Q);


//Something about calculating the product of q and p (Shreenivaas, fill this out?)
//Precondition: qWr and v2qRES are 4-element arrays...qprodRES may not be necessary since that's also the return value...
//Postcondition: Returns qp
double* qprod(double* qWR, double* v2qRES, double* qprodRES);


//Create jacobian representing the change in state (minus features, it seems?)
//Precondition: dfv_by_dxv is allocated as a 13x13 array, xv is the 13x1 state matrix (allocated and set)
//              control vector (u) is allocated as a 6x1 array (all zeros)
//Postcondition: dfv_by_dxv set
double* dfv_by_dxv(double* dfv_by_dxv, double* xv, double* u, double dt, NSString* type);


//Return components of Jacobian
//Precondition: dqomagadt_by_domegaRES is allocated as a 4x3 (12 element) array
//Postcondition: dqomagadt_by_domegaRES elements are set properly
double* dqomegadt_by_domega(double* dqomegadt_by_domegaRES, double omega1, double omega2, double omega3, double delta_t);


//Method which comes under func_Q.
// Returns a 4 x 4 matrix as result
double* q2tr(double* qOld, double* q2trRES);


// Method which comes under func_Q
//Returns a 1 x 3 matrix
double* tr2pry(double* q2trRES, double* tr2rpyRES);


//Method which returns the Jacobian result taking Euler angles as input
//Returns the 4 x 3 matrix as result
double* dq_by_deuler(double* tr2rpyRES, double* dq_by_deulerRES);

#endif