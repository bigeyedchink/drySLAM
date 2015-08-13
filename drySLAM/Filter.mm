//
//  Filter.m
//  drySLAM
//
//  Created by Anthony Rodriguez on 4/1/15.
//  Copyright (c) 2015 ECE473 573 Common. All rights reserved.
//

#import <Foundation/Foundation.h>
#include "Filter.hpp"
#include "SLAMModel.hpp"

@implementation Filter: NSObject{
    
    
}



-(id)init{
    Filter* filter = self;
    filter->x_k_k = NULL;
    filter->p_k_k = NULL;
    filter->state_size = 0;
    filter->std_a = std_a;
    filter->std_alpha = std_alpha;
    filter->std_z = 0;
    filter->x_k_km1 = NULL;
    filter->p_k_km1 = NULL;
    filter->state_km1_size = 0;
    filter->predicted_measurements = NULL;
    filter->predicted_measurements_size = 0;
    filter->H_predicted = NULL;
    filter->H_predicted_size = 0;
    filter->R_predicted = NULL;
    filter->R_predicted_size = 0;
    filter->S_predicted = NULL;
    filter->S_predicted_size = 0;
    filter->h = NULL;
    filter->h_size = 0;
    filter->type = [NSString stringWithFormat:@"constant_velocity"];
    return filter;
};

-(void)destroy{
    if (p_k_k!=NULL)
        free(p_k_k);
    if (p_k_km1!=NULL)
        free(p_k_km1);
    if (x_k_k!=NULL)
        free(x_k_k);
    if (x_k_km1!=NULL)
        free(x_k_km1);
    type = nil;
}


@end