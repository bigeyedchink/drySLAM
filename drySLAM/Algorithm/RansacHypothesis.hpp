//
//  RansacHypothesis.hpp
//  drySLAM
//
//  Created by Anthony Rodriguez on 4/17/15.
//  Copyright (c) 2015 ECE473 573 Common. All rights reserved.
//

#ifndef RANSACHYPOTHESIS_HPP
#define RANSACHYPOTHESIS_HPP

#include "SLAMModel.hpp"


//Main gateway
bool ransac_hypothesis(Filter* filter, NSMutableArray* features_info, Camera* cam);


//Generate state vector pattern
//Precondition: z_id, z_euc and state_vector_pattern are allocated but not initialized
//              z_id/z_euc end up being 2xn matrices
//              state_vector_pattern: (state_size x 4) matrix
//Postcondition: z_id, z_euc and state_vector_pattern are updated and set to correct values
//NOTE: z_id and z_euc are lists of x,y values stored as follows:   x1 x2 x3 x4
//                                                                  y1 y2 y3 y4
//       So, In order to store this we must create a list
//       of x values and a list of y values,
// then z is a 1d array stored in a way that looks like:            x1 x2 x3 x4 y1 y2 y3 y4
//***z_id and z_euc are stored as ints
bool generate_state_vector_pattern(NSMutableArray* features_info, int state_km1_size, double* state_vector_pattern, NSMutableArray* z_id, NSMutableArray* z_euc);


//Select potential match for high innovation inlier (I think)
//Precondition: zi, position and num_IC_matches are allocated.
//Postcondition: zi, position and num_IC_matches are allocated and set correctly
bool select_random_match(NSMutableArray* features_info,
                         int* zi, int* position, int* num_IC_matches);


//Predict measurements and count matches under a threshold
//Precondition: hypothesis_support, positions_li_inliers_sizes not initialized, positions_li_inliers is NULL
//              Positions_li_inliers arrays are boolean arrays (1s and 0s)
//Postcondition: hypothesis_support, positions_li_inliers_id, positions_li_inliers_id_size set correctly
//              If there are no id or euc inliers, keep that array set to NULL and set its size to 0
bool compute_hypothesis_support(double* xi, int state_size, Camera* cam, double* state_vector_pattern,
                                double* z_id, int z_id_length, double* z_euc, int threshold,
                                int* hypothesis_support,
                                NSMutableArray* positions_li_inliers_id,
                                NSMutableArray* positions_li_inliers_euc);


//Set appropriate features to be inliers
//Precondition: All variables are allocated and initialized
//Postcondition: Appropriate features_info elements have correct "low_innovation_inlier" values
bool set_as_most_supported_hypothesis(NSMutableArray* features_info,
                                      NSMutableArray* positions_li_inliers_id,
                                      NSMutableArray* positions_li_inliers_euc);

#endif
