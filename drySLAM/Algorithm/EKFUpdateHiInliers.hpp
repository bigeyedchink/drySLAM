//
//  EKFUpdateHiInliers.hpp
//  drySLAM
//
//  Created by Anthony Rodriguez on 4/17/15.
//  Copyright (c) 2015 ECE473 573 Common. All rights reserved.
//

#ifndef EKFUPDATEHIINLIERS_HPP
#define EKFUPDATEHIINLIERS_HPP

#include "SLAMModel.hpp"

//Main gateway:
bool ekf_update_hi_inliers(Filter* filter, NSMutableArray* features_info);

//update() is in EKFUpdateLiInliers.hpp

#endif