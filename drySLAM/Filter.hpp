//
//  Filter.h
//  drySLAM
//
//  Created by Anthony Rodriguez on 4/1/15.
//  Copyright (c) 2015 ECE473 573 Common. All rights reserved.
//

#ifndef drySLAM_Filter_h
#define drySLAM_Filter_h
#import <Foundation/Foundation.h>

@interface Filter : NSObject{
    //State vector (n x 1 vector)
    @public double* x_k_k;
    //Covariance matrix (n x n square matrix)
    @public double* p_k_k;
    
    @public int state_size;
    
    @public double std_a;
    @public double std_alpha;
    @public double std_z;
    
    @public double* x_k_km1;
    @public double* p_k_km1;
    @public int state_km1_size;
    
    @public double* predicted_measurements;
    @public int predicted_measurements_size;
    
    @public double* H_predicted;
    @public int H_predicted_size;
    
    @public double* R_predicted;
    @public int R_predicted_size;
    
    @public double* S_predicted;
    @public int S_predicted_size;
    //S_Matching
    //z
    @public double* h;
    @public int h_size;

    @public NSString *type;
}

-(id)init;
-(void)destroy;



@end

#endif
