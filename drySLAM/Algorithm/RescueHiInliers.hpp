//
//  RescueHiInliers.hpp
//  drySLAM
//
//  Created by Anthony Rodriguez on 4/17/15.
//  Copyright (c) 2015 ECE473 573 Common. All rights reserved.
//

#ifndef RESCUEHIINLIERS_HPP
#define RESCUEHIINLIERS_HPP

#include "SLAMModel.hpp"

#define CHI_2_INV_2_95 5.9915f

//Main gateway
//Dependencies: predict_camera_measurements and calculate_derivatives
bool rescue_hi_inliers(Filter* filter, NSMutableArray* features_info, Camera* cam);

//predict_camera_measurements() is in map_management.hpp
//calculate_derivatives is in searchICMatches.hpp

#endif
