//
//  SLAMModelController.m
//  drySLAM
//
//  Created by Anthony Rodriguez on 3/28/15.
//  Copyright (c) 2015 ECE473 573 Common. All rights reserved.
//

#import "SLAMModel.hpp"
#import "RANSAC.hpp"
#import "EKFilter.hpp"
#include "Algorithm/EKFPrediction.hpp"
#include "Algorithm/MapManagement.hpp"
#include "Algorithm/EKFUpdateHiInliers.hpp"
#include "Algorithm/EKFUpdateLiInliers.hpp"
#include "Algorithm/RansacHypothesis.hpp"
#include "Algorithm/RescueHiInliers.hpp"
#include "Algorithm/SearchICMatches.hpp"
#include <stdlib.h>
#include <math.h>
#include <vector>
using namespace cv;

@implementation SLAMModel

//Function definitions for handling the data structures
/*
+(bool) initializeCameraData:(struct cameraData*) cam{
    if (cam==NULL)
        return false;
    cam->k1 = 0.f;
    cam->k2 = 0.f;
    cam->nRows = 0;
    cam->nCols = 0;
    cam->Cx = 0.f;
    cam->Cy = 0.f;
    cam->f = 0.f;
    cam->dx = 0.f;
    cam->dy = 0.f;
    //cam->model = (char*)malloc(26*sizeof(char));
    //strcpy(cam->model, "two_distortion_parameters");
    return true;
}*/

/*+(bool) destroyCameraData:(struct cameraData*) cam{
    if (cam==NULL)
        return false;
    if (cam->model==NULL){
        NSLog(@"destroyCameraStruct error: cameraData model double-free");
        return false;
    }
    free(cam->model);
    return true;
}

+(bool) initializeFeatureData:(struct featureData**) feature{
    if (feature==NULL)
        return false;
    //(*feature)->patch_when_initialized = cv::Mat(41, 41, CV_8UC1, 255);
    //(*feature)->patch_when_matching = cv::Mat(13, 13, CV_8UC1, 255);
    (*feature)->r_wc_when_initialized = NULL;
    (*feature)->R_wc_when_initialized = NULL;
    (*feature)->uv_when_initialized = NULL;
    (*feature)->half_patch_size_when_initialized = 20;
    (*feature)->half_patch_size_when_matching = 6;
    (*feature)->times_predicted = 0;
    (*feature)->times_measured = 0;
    (*feature)->individually_compatible = 0;
    (*feature)->low_innovation_inlier = 0;
    (*feature)->high_innovation_inlier = 0;
    (*feature)->init_frame = 0;
    (*feature)->init_measurement[0] = 0;
    (*feature)->init_measurement[1] = 0;

    NSLog(@"*********Address %d************", &((*feature)->type));
    (*feature)->type = (char *)calloc(10, sizeof(char));
    strcpy((*feature)->type, "YAY");
    //(*feature)->type = (char*)malloc(20*sizeof(char));
    //strcpy((*feature)->type, "inversedepth");
    //feature->z = NULL;
    //feature->state_size = 0;
    //feature->measurement_size = 0;
    return true;
}

+(bool) destroyFeatureData:(struct featureData**) feature{
    if (feature==NULL)
        return false;
    //(*feature)->patch_when_initialized.release();
    //(*feature)->patch_when_matching.release();
    if ((*feature)->type==NULL){
        NSLog(@"Error destroyFeatureData: type double-free");
        return false;
    }
    free((*feature)->type);
    return true;
}

+(bool) initializeFilterData:(struct filterData*) filter{
    if (filter==NULL)
        return false;
    filter->x_k_k = NULL;
    filter->p_k_k = NULL;
    filter->state_size = 0;
    filter->std_a = 0;
    filter->std_alpha = 0;
    filter->std_z = 0;
    filter->x_k_km1 = NULL;
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
    return true;
}

+(bool) destroyFilterData:(struct filterData*) filter{
    if (filter==NULL)
        return false;
    return true;
}
*/
//This will probably be implemented in RANSAC
//finds features to delete and calls delete_a_feature
//-(void) deleteFeatures{
    //NSMutableArray* deletion_list = [[NSMutableArray alloc]init];
    //Find features to delete
    //iterate through deletion_list and delete them

        //int featureToDelete = 0;
        //double* x_k_k_new;
        //double* p_k_k_new;
        //RANSAC CALL:deleteOneFeature(features_info[i].type, featureToDelete, x_k_k, p_k_k, x_k_k_new, p_k_k_new);
        //free old xk, pk
        //set them to new xk,pk
    //a couple other things
//}

//Called by initialize_features. Done here because of computer vision tasks
//-(void)initialize_a_feature{
    
//}

//-(void)initialize_features{
    //while loop
        //CALL:initialize_a_feature()
//}

//Find correspondences in the search regions using normalized cross-correlation
//-(void)matching{
    
//}

//Updates features for low innovation inliers
//-(void)ransac_hypothesis{
    //RANSAC call generate_state_vector_pattern(features_info, features_info_size, x_k_km1, state_size, z_id, z_euc, state_vector_pattern);
    //for each hypothesis:
        //RANSAC call select_random_match(features_info, features_info_size, zi, position, IC_matches
        //RANSAC call compute_hypothesis_support
        //RANSAC call set_as_most_supported_hypothesis
    //end
//}


//Main gateway to run code (processes video at filePath)

-(void) processVideoFile:(NSString*)filePath
{
    int i,j;
    
    double v_0 = 0;
    double w_0 = 1e-15;
    double eps = 2e-52;
    double std_v_0 = 0.025;
    double std_w_0 = 0.025;
    
    double *x_k_k = (double *)malloc(13*1*sizeof(double));
    for (i=0; i<13; i++)
        x_k_k[i] = 0;
    x_k_k[3] = 1;
    x_k_k[7] = v_0;
    x_k_k[8] = v_0;
    x_k_k[9] = v_0;
    x_k_k[10] = w_0;
    x_k_k[11] = w_0;
    x_k_k[12] = w_0;
    
    double *p_k_k = (double *)malloc(13*13*sizeof(double));
    for (i=0; i<13; i++)
        for (j=0; j<13; j++)
            p_k_k[i*13+j] = 0;
    p_k_k[0] = eps;
    p_k_k[14] = eps;
    p_k_k[28] = eps;
    p_k_k[42] = eps;
    p_k_k[56] = eps;
    p_k_k[70] = eps;
    p_k_k[84] = eps;
    p_k_k[98] = pow(std_v_0, 2);
    p_k_k[112] = pow(std_v_0, 2);
    p_k_k[126] = pow(std_v_0, 2);
    p_k_k[140] = pow(std_w_0, 2);
    p_k_k[154] = pow(std_v_0, 2);
    p_k_k[168] = pow(std_v_0, 2);
    
    double *X = (double *)malloc(1*1000*sizeof(double));
    for (i=0; i<1000; i++)
        X[i] = ((double)rand()/(double)RAND_MAX)-0.5;
    
    double *Y = (double *)malloc(1*1000*sizeof(double));
    for (i=0; i<1000; i++)
        Y[i] = ((double)rand()/(double)RAND_MAX)-0.5;
    
    double *Z = (double *)malloc(1*1000*sizeof(double));
    for (i=0; i<1000; i++)
        Z[i] = ((double)rand()/(double)RAND_MAX)-0.5;
    
    double *theta = (double *)malloc(1*1000*sizeof(double));
    for (i=0; i<1000; i++)
        theta[i] = ((double)rand()/(double)RAND_MAX)-0.5;
    
    double *phi = (double *)malloc(1*1000*sizeof(double));
    for (i=0; i<1000; i++)
        phi[i] = ((double)rand()/(double)RAND_MAX)-0.5;
    
    double *lambda = (double *)malloc(1*1000*sizeof(double));
    for (i=0; i<1000; i++)
        lambda[i] = ((double)rand()/(double)RAND_MAX)-0.5;
    
    cv::Mat *a = new cv::Mat(1,6,CV_64FC1,0.0);
    double a_norm;
    for (i=0; i<1000; i++)
    {
        a->at<double>(0,0) = X[i];
        a->at<double>(0,1) = Y[i];
        a->at<double>(0,2) = Z[i];
        a->at<double>(0,3) = theta[i];
        a->at<double>(0,4) = phi[i];
        a->at<double>(0,5) = lambda[i];
        
        a_norm = cvNorm(a);
        for (j=0; j<6; j++)
            a[j] = a[j]/a_norm*sqrt(12.59158724374398);

        X[i] = a->at<double>(0,0);
        Y[i] = a->at<double>(0,1);
        Z[i] = a->at<double>(0,2);
        theta[i] = a->at<double>(0,3);
        phi[i] = a->at<double>(0,4);
        lambda[i] = a->at<double>(0,5);
    }
    
    for (j=0; j<1000; j++)
    {
        randSphere6D[0][j] = X[j];
        randSphere6D[1][j] = Y[j];
        randSphere6D[2][j] = Z[j];
        randSphere6D[3][j] = theta[j];
        randSphere6D[4][j] = phi[j];
        randSphere6D[5][j] = lambda[j];
    }
    
    int min_number_of_features_in_image = 25;
    
    NSMutableArray *measurements;
    NSMutableArray *predicted_measurements;
    
    //cvNamedWindow("Vedio",CV_WINDOW_AUTOSIZE);
    //CvCapture *capture = cvCreateFileCapture("dryslam.avi");
    //IplImage *frame;
    //frame = cvQueryFrame(capture);
    cv::Mat* img;
    int step;
    
    while (1)
    {
        if (!img)
            break;
        
        //step = cvGetCaptureProperty(capture, CV_CAP_PROP_POS_FRAMES);
        map_management(img, filter, cam, features_info, min_number_of_features_in_image, step-1);
        ekf_prediction(filter, features_info);
        //frame = cvQueryFrame(capture);
        search_IC_matches(filter, features_info, cam, img);
        ransac_hypothesis(filter, features_info, cam);
        ekf_update_li_inliers(filter, features_info);
        rescue_hi_inliers(filter, features_info, cam);
        ekf_update_hi_inliers(filter, features_info);
    }
    
    //cvReleaseCapture(&capture);
    //cvDestroyWindow("Vedio");
}

@end
