//
//  MapManagement.m
//  drySLAM
//
//  Created by Anthony Rodriguez on 4/17/15.
//  Copyright (c) 2015 ECE473 573 Common. All rights reserved.
//

#include "RansacHypothesis.hpp"
#include <stdlib.h>

bool ransac_hypothesis(Filter* filter, NSMutableArray* features_info, Camera* cam)
{
    int i;
    
    float p_at_least_one_spurious_free = 0.99;
    double threshold = filter->std_z;
    int n_hyp = 1000;
    int max_hypothesis_support = 0;
    
    int size = filter->state_km1_size;
    double *state_vector_pattern = (double *)malloc(size*4*sizeof(double));
    NSMutableArray *z_id = [[NSMutableArray alloc] init];
    NSMutableArray *z_euc = [[NSMutableArray alloc] init];
    generate_state_vector_pattern(features_info, size, state_vector_pattern, z_id, z_euc);
    
    for (i=0; i<n_hyp; i++)
    {
        double *zi;
        int position;
        int num_IC_matches;
        //select_random_match(features_info, zi, position, num_IC_matches);
        
        cv::Mat *x_k_km1_mat = new cv::Mat(1,size,CV_64FC1,filter->x_k_km1);
        cv::Mat *p_k_km1_mat = new cv::Mat(size,size,CV_64FC1,filter->p_k_km1);
        
        Feature *current = features_info[position-1];
        cv::Mat *hi = new cv::Mat();
        
        cv::Mat *Hi = new cv::Mat(size,size,CV_64FC1,current->H);
        cv::Mat *S = new cv::Mat(size,size,CV_64FC1,0.0);
        
        //compute_hypothesis_support(xi,
 
        double epsilon;
        int hypothesis_support;
        if (hypothesis_support>max_hypothesis_support)
        {
            max_hypothesis_support = hypothesis_support;
            
            NSMutableArray *position_li_inliers_id = [[NSMutableArray alloc] init];
            NSMutableArray *position_li_inliers_euc = [[NSMutableArray alloc] init];
            set_as_most_supported_hypothesis(features_info, position_li_inliers_id, position_li_inliers_euc);
            
            epsilon = 1-hypothesis_support/num_IC_matches;
            
            n_hyp = (int)(log(1-p_at_least_one_spurious_free))/log(1-(1-epsilon));
            if (n_hyp == 0)
                break;
        }
 
    }
    
    return false;
}

/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
/***************************************************************************************************************************************/

bool generate_state_vector_pattern(NSMutableArray* features_info, int state_size, double* state_vector_pattern, NSMutableArray* z_id, NSMutableArray* z_euc)
{
    if (features_info==NULL||state_size==-1||state_vector_pattern==NULL||z_id==NULL||z_euc==NULL)
        return false;
    
    int position=14;
    
    int i,j;
    //int J1=0;
    //int J2=0;
    
    for (i=0; i<state_size*4; ++i){
        state_vector_pattern[i] = 0.0;
    }
    
    NSMutableArray *z_id_x = [[NSMutableArray alloc] init];
    NSMutableArray *z_id_y = [[NSMutableArray alloc] init];
    NSMutableArray *z_euc_x = [[NSMutableArray alloc] init];
    NSMutableArray *z_euc_y = [[NSMutableArray alloc] init];
    
    for (i=0; i<[features_info count]; i++)
    {
        Feature *current=features_info[i];
        
        if ([current->type isEqualToString:@"inversedepth"])
        {
            if (current->z != NULL)
            {
                state_vector_pattern[(4*(position-1))+0] = 1;
                state_vector_pattern[(4*(position+0))+0] = 1;
                state_vector_pattern[(4*(position+1))+0] = 1;
                state_vector_pattern[(4*(position+2))+1] = 1;
                state_vector_pattern[(4*(position+3))+1] = 1;
                state_vector_pattern[(4*(position+4))+2] = 1;
                
                [z_id_x addObject:[NSNumber numberWithDouble:current->z[0]]];
                [z_id_y addObject:[NSNumber numberWithDouble:current->z[1]]];
                
                //J1++;
            }
            
            position += 6;
        }
        
        if ([current->type isEqualToString:@"cartesian"])
        {
            if (current->z != 0)
            {
                state_vector_pattern[(4*(position-1))+3] = 1;
                state_vector_pattern[(4*(position+0))+3] = 1;
                state_vector_pattern[(4*(position+1))+3] = 1;
                
                [z_euc_x addObject:[NSNumber numberWithDouble:current->z[0]]];
                [z_euc_y addObject:[NSNumber numberWithDouble:current->z[1]]];
                
                //J2++;
             }
            
             position += 3;
         }
        
        
    }
    
    for (i=0; i<[z_id_x count]; ++i){
        [z_id addObject: z_id_x[i]];
    }
    for (i=0; i<[z_id_y count]; ++i){
        [z_id addObject: z_id_y[i]];
    }
    for (i=0; i<[z_euc_x count]; ++i){
        [z_euc addObject: z_euc_x[i]];
    }
    for (i=0; i<[z_euc_y count]; ++i){
        [z_euc addObject: z_euc_y[i]];
    }
    /*
    if (J1 != 0)
    {
        for (j=0; j<J1; j++)
        {
            [z_id insertObject:[z_id_0 objectAtIndex:2*j+0]
                       atIndex:0*J1+j];
            [z_id insertObject:[z_id_0 objectAtIndex:2*j+1]
                       atIndex:1*J1+j];
        }
    }
    
    if (J2 != 0)
    {
        for (j=0; j<J2; j++)
        {
            [z_euc insertObject:[z_euc_0 objectAtIndex:2*j+0]
                        atIndex:0*J2+j];
            [z_euc insertObject:[z_euc_0 objectAtIndex:2*j+1]
                        atIndex:1*J2+j];
        }
    }*/
    
    return true;
}

/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
/***************************************************************************************************************************************/

bool select_random_match(NSMutableArray* features_info, int *zi, int* position, int* num_IC_matches)
{
    int map_size = (int)[features_info count];
    int individually_compatible[map_size];
    
    int i;
    int sum = 0;
    NSMutableArray *positions_individually_compatible = [[NSMutableArray alloc] init];
    for (i=0; i<map_size; i++)
    {
        Feature *current = features_info[i];
        
        if (current->individually_compatible)
        {
            individually_compatible[i] = 1;
            
            [positions_individually_compatible addObject: [NSNumber numberWithInt:i]];
        }
        else
            individually_compatible[i] = 0;
        
        sum += individually_compatible[i];
    }
    
    int random_match_position = (((double)rand()/((double)RAND_MAX))*sum);
    //position = positions_individually_compatible[random_match_position];
    
    Feature *current = (Feature*)features_info[random_match_position];
    zi[0] = current->z[0];
    zi[1] = current->z[1];
    
    *num_IC_matches = sum;
    *position = random_match_position;
    
    return true;
}

/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
/***************************************************************************************************************************************/

bool compute_hypothesis_support(double* xi, int state_size, Camera* cam, double* state_vector_pattern,
                                double* z_id, int z_id_length, double* z_euc, int threshold,
                                int* hypothesis_support,
                                NSMutableArray* positions_li_inliers_id,
                                NSMutableArray* positions_li_inliers_euc){
    *hypothesis_support = 0;
    bool success = true;
    if(z_id != NULL) {
        int n_id = z_id_length;
        
        NSMutableArray* ri = [[NSMutableArray alloc]init];
        NSMutableArray* anglesi = [[NSMutableArray alloc]init];
        NSMutableArray* rhoi = [[NSMutableArray alloc]init];
        for (int i=0; i<state_size; ++i){
            if (state_vector_pattern[4*i]==1)
                [ri addObject: [NSNumber numberWithDouble:xi[i]]];
            if (state_vector_pattern[4*i+1]==1)
                [anglesi addObject: [NSNumber numberWithDouble:xi[i]]];
            if (state_vector_pattern[4*i+2]==1)
                [rhoi addObject: [NSNumber numberWithDouble:xi[i]]];
        }
        
        //ri: count x 1 change to 3 x n_id
        double* ri_reshaped = (double*)malloc([ri count]*sizeof(double));
        int j = 0;
        for (int i=0; i<[ri count]; i = i+3){
            ri_reshaped[j] = (double)[(NSNumber*)ri[i] doubleValue];
            ri_reshaped[n_id+j] = (double)[(NSNumber*)ri[i+1] doubleValue];
            ri_reshaped[n_id*2+j] = (double)[(NSNumber*)ri[i+2] doubleValue];
            j++;
        }
        //reshape to 2 x n_id
        double* anglesi_reshaped = (double*)malloc([anglesi count]*sizeof(double));
        j=0;
        for (int i=0; i<[anglesi count]; i = i+2){
            anglesi_reshaped[j] = (double)[anglesi[i] doubleValue];
            anglesi_reshaped[n_id+j] = (double)[anglesi[i+1] doubleValue];
            ++j;
        }
        
        //3xn_id matrix
        double* mi = (double*)malloc(3*n_id*sizeof(double));
        //Compute azimuth-elevation angles
        for (int i=0; i<n_id; ++i){
            double theta = anglesi_reshaped[i];
            double phi = anglesi_reshaped[n_id+i];
            double cphi = cos(phi);
            mi[i] = cphi*sin(theta);        //first row
            mi[n_id+i] = -sin(phi);         //second row
            mi[n_id*2+i] = cphi*cos(theta); //third row
        }
        
        double* rwc = (double*)malloc(3*sizeof(double));
        rwc[0] = xi[0];
        rwc[1] = xi[1];
        rwc[2] = xi[2];
        
        double* xi_4To7 = (double*)malloc(4*sizeof(double));
        xi_4To7[0] = xi[3];
        xi_4To7[1] = xi[4];
        xi_4To7[2] = xi[5];
        xi_4To7[3] = xi[6];
        double* rotwc = (double*)malloc(3*3*sizeof(double));
        success = success&&q2r(xi_4To7, rotwc);
        
        double* ri_minus_rwc_x = (double*)malloc(n_id*sizeof(double));
        double* ri_minus_rwc_y = (double*)malloc(n_id*sizeof(double));
        double* ri_minus_rwc_z = (double*)malloc(n_id*sizeof(double));
        for (int i=0; i<n_id; ++i){
            ri_minus_rwc_x[i] = (ri_reshaped[i] - rwc[0])*[rhoi[i] doubleValue];
            ri_minus_rwc_y[i] = (ri_reshaped[n_id+i] - rwc[1])*[rhoi[i] doubleValue];
            ri_minus_rwc_z[i] = (ri_reshaped[n_id*2+i] - rwc[2])*[rhoi[i] doubleValue];
        }
        
        // 3 x n_id matrix
        double* ri_minus_rwc_by_rhoi = (double*)malloc(3*n_id*sizeof(double));
        for (int i=0; i<n_id; ++i){
            ri_minus_rwc_by_rhoi[i] = ri_minus_rwc_x[i] + mi[i];
            ri_minus_rwc_by_rhoi[n_id+i] = ri_minus_rwc_y[i] + mi[n_id+i];
            ri_minus_rwc_by_rhoi[n_id*2+i] = ri_minus_rwc_z[i] + mi[n_id*2+i];
        }
        
        cv::Mat rotcw_mat(3, 3, CV_64FC1, rotwc, cv::Mat::AUTO_STEP);
        //rotcw is the transpose of rotwc:
        cv::transpose(rotcw_mat, rotcw_mat);
        cv::Mat ri_minus_rwc_by_rhoi_plus_mi(3, n_id, CV_64FC1, ri_minus_rwc_by_rhoi, cv::Mat::AUTO_STEP);
        cv::Mat hc = rotcw_mat*ri_minus_rwc_by_rhoi_plus_mi;
        
        double u0 = cam->Cx;
        double v0 = cam->Cy;
        double f = cam->f;
        double ku = 1/cam->dx;
        
        double* h_image = (double*)malloc(2*n_id*sizeof(double));
        double* hcPtr = hc.ptr<double>(0);
        for (int i=0; i<n_id; ++i){
            double denominator = hcPtr[n_id*2+i];
            h_image[i] = (hcPtr[i]/denominator)*f*ku + u0;
            h_image[n_id+i] = (hcPtr[n_id+i]/denominator)*f*ku + v0;
        }
        
        double* h_distorted = (double*)malloc(2*n_id*sizeof(double));
        success = success&&distortFM(h_distorted, h_image, n_id, cam);
        
        double* nu = (double*)malloc(2*n_id*sizeof(double));
        for (int i=0; i<2*n_id; ++i){
            nu[i] = z_id[i] - h_distorted[i];
            NSLog(@"%lf ", nu[i]);
        }
        //double* residuals = (double*)malloc(n_id*sizeof(double));
        int inliersCounter = 0;
        for (int i=0; i<n_id; ++i){
            double residual = sqrt(nu[i]*nu[i] + nu[n_id+i]*nu[n_id+i]);
            if (residual<threshold){
                [positions_li_inliers_id addObject: [NSNumber numberWithBool:YES]];
                ++inliersCounter;
            }
            else
                [positions_li_inliers_id addObject: [NSNumber numberWithBool:NO]];
        }
        *hypothesis_support = inliersCounter;
        
        free(nu);
        free(h_distorted);
        free(h_image);
        hc.release();
        ri_minus_rwc_by_rhoi_plus_mi.release();
        rotcw_mat.release();
        free(ri_minus_rwc_by_rhoi);
        free(ri_minus_rwc_x);
        free(ri_minus_rwc_y);
        free(ri_minus_rwc_z);
        free(rotwc);
        free(rwc);
        free(mi);
        free(anglesi_reshaped);
        free(ri_reshaped);
        free(xi_4To7);
        rhoi = nil;
        ri = nil;
        anglesi = nil;
    }
    if (z_euc!=NULL){
        NSLog(@"Error: Cartesian coordinate systems are not supported yet!");
        return false;
    }
    
    return true;
}



/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
/***************************************************************************************************************************************/

bool set_as_most_supported_hypothesis(NSMutableArray* features_info,NSMutableArray* positions_li_inliers_id,NSMutableArray* positions_li_inliers_euc)
{
    if (features_info==NULL)
        return false;
    
    if (positions_li_inliers_id==NULL||positions_li_inliers_euc==NULL)
        return false;
    
    int j_id=0;
    int j_euc=0;
    int i;
    
    for (i=0;i<(int)[features_info count];i++)
    {
        Feature *current=features_info[i];
        
        if (current->z) //if ~isempty(features_info(i).z)
        {
            if ([current->type isEqualToString:@"cartesian"])
            {
                if (positions_li_inliers_euc[j_euc])
                    current->low_innovation_inlier=1;
                else
                    current->low_innovation_inlier=0;
                
                j_euc++;
            }
            
            if ([current->type isEqualToString:@"inversedepth"])
            {
                if (positions_li_inliers_euc[j_id])
                    current->low_innovation_inlier=1;
                else
                    current->low_innovation_inlier=0;
                
                j_id++;
            }
        }
    }
    
    return true;
}