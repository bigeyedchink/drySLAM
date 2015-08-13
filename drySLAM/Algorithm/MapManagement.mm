//
//  MapManagement.mm
//  drySLAM
//
//  Created by Anthony Rodriguez on 4/17/15.
//  Copyright (c) 2015 ECE473 573 Common. All rights reserved.
//

#include "MapManagement.hpp"
#include "SearchICMatches.hpp"
#include <stdlib.h>
#include <math.h>
//#include <stdio>
#include <vector>

//Main function:
bool map_management(cv::Mat* im, Filter* filter, Camera* cam, NSMutableArray* features_info, int minNumberOfFeaturesInImage, int step)
{
    int i;
    
    deleteFeatures(features_info, filter);
    
    int measured = 0;
    for (i=0; i<(int)[features_info count]; i++)
    {
        Feature *current = features_info[i];
        
        if (current->low_innovation_inlier || current->high_innovation_inlier)
            measured++;
    }
    
    updateFeaturesInfo(features_info);
    
    inversedepth_2_cartesian(features_info, filter);
    
    if (measured == 0)
        initialize_features(step, cam, filter, features_info, minNumberOfFeaturesInImage, im);
    else
        if (measured < minNumberOfFeaturesInImage)
            initialize_features(step, cam, filter, features_info, minNumberOfFeaturesInImage-measured, im);
    
    return true;
}

/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
bool deleteOneFeature(NSMutableArray* features_info,Filter* filter,int featToDelete)
{
    Feature *last = features_info[featToDelete-1];
    int parToDelete;
    
    if ([last->type isEqualToString:@"cartesian"])
        parToDelete = 3;
    else
        parToDelete = 6;
    
    int indexFromWhichDelete = 14;
    int i, j;
    
    for (i=0; i<featToDelete-1; i++)
    {
        Feature *current = features_info[i];
        
        if ([current->type isEqualToString:@"inversedepth"])
            indexFromWhichDelete += 6;
        
        if ([current->type isEqualToString:@"cartesian"])
            indexFromWhichDelete += 3;
    }
    
    int size=filter->state_size;
    
    double *x_k_k_new = (double *)malloc(1*(size-parToDelete)*sizeof(double));
    for (j=0; j<indexFromWhichDelete-1; j++)
        x_k_k_new[j] = filter->x_k_k[j];
    for (j=indexFromWhichDelete+parToDelete-1; j<size; j++)
        x_k_k_new[j-parToDelete] = filter->x_k_k[j];
    
    double *p_k_k_new_0 = (double *)malloc(size*(size-parToDelete)*sizeof(double));
    for (i=0; i<size; i++)
    {
        for (j=0; j<indexFromWhichDelete-1; j++)
            p_k_k_new_0[i*(size-parToDelete)+j] = filter->p_k_k[i*size+j];
        for (j=indexFromWhichDelete+parToDelete-1; j<size; j++)
            p_k_k_new_0[i*(size-parToDelete)+j-parToDelete] = filter->p_k_k[i*size+j];
    }
    
    double *p_k_k_new = (double *)malloc((size-parToDelete)*(size-parToDelete)*sizeof(double));
    for (i=0; i<indexFromWhichDelete-1; i++)
    {
        for (j=0; j<size-parToDelete; j++)
            p_k_k_new[i*(size-parToDelete)+j] = p_k_k_new_0[i*(size-parToDelete)+j];
    }
    for (i=indexFromWhichDelete+parToDelete-1; i<size; i++)
    {
        for (j=0; j<size-parToDelete; j++)
            p_k_k_new[(i-parToDelete)*(size-parToDelete)+j] = p_k_k_new_0[i*(size-parToDelete)+j];
    }
    
    filter->state_size=size-parToDelete;
    
    for (i=0; i<size-parToDelete; i++)
    {
        filter->x_k_k[i] = x_k_k_new[i];
        
        for (j=0; j<size-parToDelete; j++)
            filter->p_k_k[i*(size-parToDelete)+j] = p_k_k_new[i*(size-parToDelete)+j];
    }
    
    free(x_k_k_new);
    free(p_k_k_new_0);
    free(p_k_k_new);
    
    return true;
}

bool deleteFeatures(NSMutableArray* features_info, Filter* filter)
{
    int i,j;
    
    if (features_info==NULL||filter==NULL)
    {
        return false;
    }
    
    NSMutableArray *deletionList=[[NSMutableArray alloc] init];
    
    for(i=0; i<(int)[features_info count]; i++)
    {
        Feature *current=features_info[i];
        
        if ((current->times_measured<0.5*(current->times_predicted))&&
            (current->times_predicted>5))
        {
            [deletionList addObject:[NSNumber numberWithInt:i]];
        }
    }
    
    if ([deletionList count])
    {
        for (j=(int)[deletionList count]-1;j>=0;j--)
        {
            cv::Mat x_k_k_new;
            cv::Mat p_k_k_new;
            
            deleteOneFeature(features_info,filter,[[deletionList objectAtIndex:j] intValue]);
            
            if ([[deletionList objectAtIndex:j] intValue]==0)
                [features_info removeObjectAtIndex:0];
            
            if ([[deletionList objectAtIndex:j] intValue]==[features_info count])
                [features_info removeObjectAtIndex:[features_info count]-1];
            
            if (([[deletionList objectAtIndex:j] intValue]!=[features_info count])&&
                ([[deletionList objectAtIndex:j] intValue]!=0))
                [features_info removeObjectAtIndex:[[deletionList objectAtIndex:j] intValue]];
            //filter->state_size -= 6;
        }
        
    }
    
    return true;
}

/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
/***************************************************************************************************************************************/


/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
/***************************************************************************************************************************************/

bool updateFeaturesInfo(NSMutableArray* features_info)
{
    int i;
    
    for (i=0; i<(int)[features_info count]; i++)
    {
        Feature *current = features_info[i];
        if (sizeof(current->h))     //~isempty(features_info(i).h)
            current->times_predicted = current->times_measured+1;
        if (current->low_innovation_inlier || current->high_innovation_inlier)
            current->times_measured++;
        
        current->individually_compatible = 0;
        current->low_innovation_inlier = 0;
        current->high_innovation_inlier = 0;
        
        current->h = NULL;
        current->z = NULL;
        current->H = NULL;
        current->S = NULL;
    }
    
    return true;
}

/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
void m(double theta, double phi, cv::Mat *m)
{
    m->at<double>(0,0) = cos(phi) * sin(theta);
    m->at<double>(1,0) = (-1)*sin(phi);
    m->at<double>(2,0) = cos(phi) * cos(theta);
}

bool inversedepth_2_cartesian(NSMutableArray* features_info, Filter* filter)
{
    float linearity_index_threshold = 0.1;
    
    int i,j;
    
    for (i=0; i<(int)[features_info count]; i++)
    {
        Feature *current = features_info[i];
        
        int initialPositionOfFeature = 14;
        
        if ([current->type isEqualToString:@"inversedepth"])
        {
            for (j=0; j<i-1; j++)
            {
                if ([current->type isEqualToString:@"cartesian"])
                    initialPositionOfFeature += 3;

                if ([current->type isEqualToString:@"inversedepth"])
                    initialPositionOfFeature += 6;
            }
        }
        
        double std_rho = sqrt(filter->p_k_k[(initialPositionOfFeature+4)*(filter->state_size)+initialPositionOfFeature+4]);
        double rho = filter->x_k_k[initialPositionOfFeature+4];
        double std_d = std_rho/pow(rho,2);
        
        double theta = filter->x_k_k[initialPositionOfFeature+2];
        double phi = filter->x_k_k[initialPositionOfFeature+3];
        cv::Mat *mi = new cv::Mat(3,1,CV_64FC1,0.0);
        m(theta, phi, mi);
        
        double *x_c1 = (double *)malloc(1*3*sizeof(double));
        for (i=0; i<3; i++)
            x_c1[i] = filter->x_k_k[initialPositionOfFeature+(i-1)];
        
        double *x_c2 = (double *)malloc(1*3*sizeof(double));
        for (i=0; i<3; i++)
            x_c2[i] = filter->x_k_k[i];
        
        double *p0 = (double *)malloc(1*6*sizeof(double));
        for (i=0; i<5; i++)
            p0[i] = filter->x_k_k[initialPositionOfFeature+(i-1)];

        double *ct = (double *)malloc(1*3*sizeof(double));
        inverseDepth2Cartesian(p0, ct);
        
        double sum1 = sqrt(pow(ct[0]-x_c1[0], 2) + pow(ct[1]-x_c1[1], 2) + pow(ct[2]-x_c1[2], 2));
        double sum2 = sqrt(pow(ct[0]-x_c2[0], 2) + pow(ct[1]-x_c2[1], 2) + pow(ct[2]-x_c2[2], 2));
        
        double d_c2p = sum2;
        double cos_alpha = ((ct[0]-x_c1[0])*(ct[0]-x_c2[0]) + (ct[1]-x_c1[1])*(ct[1]-x_c2[1]) +
                            (ct[2]-x_c1[2])*(ct[2]-x_c2[2]))/(sum1*sum2);
        
        double linearity_index = 4*std_d*cos_alpha/d_c2p;

        if (linearity_index < linearity_index_threshold)
        {
            int size = filter->state_size;
            
            double *X = (double *)malloc(1*(size-3)*sizeof(double));
            for (i=0; i<initialPositionOfFeature-1; i++)
                X[i] = filter->x_k_k[i];
            for (i=initialPositionOfFeature-1; i<initialPositionOfFeature+2; i++)
                X[i] = ct[i-(initialPositionOfFeature-1)];
            for (i=initialPositionOfFeature+2; i<size-3; i++)
                X[i] = filter->x_k_k[i+3];
            
            double *dm_dtheta = (double *)malloc(1*3*sizeof(double));
            dm_dtheta[0] = cos(phi)*cos(theta);
            dm_dtheta[1] = 0;
            dm_dtheta[2] = -cos(phi)*sin(theta);
            
            double *dm_dphi = (double *)malloc(1*3*sizeof(double));
            dm_dphi[0] = -sin(phi)*sin(theta);
            dm_dphi[1] = -cos(phi);
            dm_dphi[2] = -sin(phi)*cos(theta);
            
            double J[3][6] = {{1,0,0,1/rho*dm_dtheta[0],1/rho*dm_dphi[0],-mi->at<double>(0,0)/pow(rho, 2)},
                              {0,1,0,1/rho*dm_dtheta[1],1/rho*dm_dphi[1],-mi->at<double>(1,0)/pow(rho, 2)},
                              {0,0,1,1/rho*dm_dtheta[2],1/rho*dm_dphi[2],-mi->at<double>(2,0)/pow(rho, 2)}};
            
            double *J_all = (double *)malloc((size-3)*size*sizeof(double));
            for (i=0; i<initialPositionOfFeature-1; i++)
                for (j=0; j<initialPositionOfFeature-1; j++)
                {
                    if (i==j)
                        J_all[i*size+j] = 1;
                }
            for (i=initialPositionOfFeature-1; i<initialPositionOfFeature+2; i++)
                for (j=initialPositionOfFeature-1; j<initialPositionOfFeature+5; j++)
                    J_all[i*size+j] = J[i-(initialPositionOfFeature-1)][j-(initialPositionOfFeature-1)];
            for (i=initialPositionOfFeature+2; i<size-3; i++)
                for (j=initialPositionOfFeature+5; j<size; j++)
                {
                    if (i==j)
                        J_all[i*size+j] = 1;
                }
            
            cv::Mat *J_all_mat = new cv::Mat(size-3,size,CV_64FC1,J_all);
            cv::Mat *J_all_transpose = new cv::Mat(size,size-3,CV_64FC1,0.0);
            cv::transpose(*J_all_mat, *J_all_transpose);
            
            cv::Mat *P_0 = new cv::Mat(size,size,CV_64FC1,filter->p_k_k);
            cv::Mat *P = new cv::Mat(size-3,size-3,CV_64FC1,0.0);
            *P = (*J_all_mat) * (*P_0) * (*J_all_transpose);
            
            filter->state_size = size-3;
            
            for (i=0; i<filter->state_size; i++)
            {
                filter->x_k_k[i] = X[i];
                
                for (j=0; j<filter->state_size; j++)
                    filter->p_k_k[i*(filter->state_size)+j] = P->at<double>(i,j);
            }
            
            current->type = @"cartesian";
        }
    }
    
    return true;
}

/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
/***************************************************************************************************************************************/

bool initialize_features(int step, Camera* cam, Filter* filter, NSMutableArray* features_info, int numFeaturesToInitialize, cv::Mat* im)
{
    int max_attempts = 50;
    int attempts = 0;
    int initialized = 0;
    
    while ((initialized<numFeaturesToInitialize)&&
           (attempts<max_attempts))
    {
        attempts++;
        
        NSMutableArray* uv = [[NSMutableArray alloc]init];
        initialize_a_feature(filter, features_info, step, cam, im, uv);
        
        if ([uv count]!=0)
            initialized++;
        uv = nil;
    }
    
    return true;
}

/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
/***************************************************************************************************************************************/


bool initialize_a_feature(Filter* filter, NSMutableArray* features_info, int step, Camera* cam, cv::Mat* im_k, NSMutableArray* uv)
{
    
    int i,j;
    
    int half_patch_size_when_initialized = 20;
    int half_patch_size_when_matching = 6;
    int excluded_band = half_patch_size_when_initialized+1;
    int max_initialization_attempts = 1;
    int initializing_box_size[2] = {60,40};
    int initializing_box_semisize[2] = {30,20};
    int initial_rho = 1;
    int std_rho = 1;
    
    double std_pxl = filter->std_z;
    
    int rand_attempt = 1;
    int not_empty_box = 1;
    bool detected_new = false;
    
    predict_camera_measurements(features_info, filter, cam);
    
    double *uv_pred = (double *)malloc(1*(2*(int)[features_info count])*sizeof(double));
    for (i=0; i<(int)[features_info count]; i++)
    {
        Feature *current = (Feature *)features_info[i];
        if (current->h!=NULL){
            uv_pred[2*i+0] = current->h[0];
            uv_pred[2*i+1] = current->h[1];
        }
    }
    
    bool are_there_features = false;
    bool success = true;
    
    for (i=0; i<max_initialization_attempts; i++){
        if (detected_new)
            break;
    
        //double search_region_center[2] = {(double)rand()/((double)RAND_MAX), (double)rand()/((double)RAND_MAX)};
        double centerX = (double)rand()/((double)RAND_MAX);
        double centerY = (double)rand()/((double)RAND_MAX);
        
        centerX = round(centerX*(cam->nCols-2*excluded_band-2*initializing_box_semisize[0]))+excluded_band+initializing_box_semisize[0];
        centerY = round(centerY*(cam->nRows-2*excluded_band-2*initializing_box_semisize[1]))+excluded_band+initializing_box_semisize[1];
        
        int searchBoxMaxX = centerX + initializing_box_semisize[0];
        int searchBoxMinX = centerX - initializing_box_semisize[0];
        int searchBoxMaxY = centerY + initializing_box_semisize[1];
        int searchBoxMinY = centerY - initializing_box_semisize[1];
        int searchImRows = 2*initializing_box_semisize[1]+1;
        int searchImCols = 2*initializing_box_semisize[0]+1;
        
        //Build image search patch:
        uchar* searchIm_data = (uchar*)malloc(searchImRows*searchImCols*sizeof(uchar));
        uchar* imPtr = im_k->ptr<uchar>(0);
        for (int r=0; r<searchImRows; ++r){
            imPtr = im_k->ptr<uchar>(r+centerY) + (searchBoxMinX);
            for (int c=0; c<searchImCols; ++c){
                searchIm_data[r*searchImCols+c] = imPtr[c];
            }
        }
        cv::Mat searchIm(searchImRows, searchImCols, CV_8UC1, searchIm_data, cv::Mat::AUTO_STEP);
        std::vector<cv::KeyPoint> cs;
        cv::FAST(searchIm, cs, 100, true);
        
        bool are_there_corners = false;
        if (!cs.empty())
            are_there_corners = true;
        
        //Find centers of each
        for (int j=0; j<cs.size(); ++j){
            cs.at(j).pt.y = cs.at(j).pt.y + searchBoxMinY - 1;
            cs.at(j).pt.x = cs.at(j).pt.x + searchBoxMinX - 1;
        }
        
        if (cs.size()>0){
            int totalFeat = (int)[features_info count];
            std::vector<int> indicesInBox;
            for (j=0; j<totalFeat; ++j){
                int offset = 2*j;
                double row = uv_pred[offset];
                double col = uv_pred[offset+1];
                if (uv_pred[offset]>searchBoxMinY && uv_pred[offset]<searchBoxMaxY && uv_pred[offset+1]>searchBoxMinX && uv_pred[offset+1]<searchBoxMaxX){
                    indicesInBox.push_back(j);
                }
            }
            if (!indicesInBox.empty()){
                are_there_features = true;
            }
            else
                are_there_features = false;
        }
        
        double* uv_data = (double*)malloc(cs.size()*2*sizeof(double));
        int* uv_int_data = (int*)malloc(cs.size()*2*sizeof(double));
        if (are_there_corners && !are_there_features){
            detected_new = true;
            NSMutableArray* uv_row = [[NSMutableArray alloc]init];
            NSMutableArray* uv_col = [[NSMutableArray alloc]init];
            for (j=0; j<cs.size(); ++j){
                [uv_row addObject: [NSNumber numberWithFloat:cs.at(j).pt.y ]];
                [uv_col addObject: [NSNumber numberWithFloat:cs.at(j).pt.x ]];
                uv_data[j] = cs.at(j).pt.y;
                uv_data[(int)cs.size()+j] = cs.at(j).pt.x;
                uv_int_data[j] = (int)cs.at(j).pt.y;
                uv_int_data[(int)cs.size()+j] = (int)cs.at(j).pt.x;
            }
            for (j=0; j<[uv_row count]; ++j){
                [uv addObject:uv_row[j]];
            }
            for (j=0; j<[uv_col count]; ++j){
                [uv addObject:uv_col[j]];
            }
        }
        
        if ([uv count]!=0){
            double* newFeature = (double*)malloc(6*sizeof(double));
            success = success && add_features_inverse_depth(uv_data, (int)cs.size(), filter, cam, std_pxl, initial_rho, std_rho, newFeature);
            success = success && add_feature_to_info_vector(uv_int_data, im_k, filter->x_k_k, filter->state_size, features_info, step, newFeature);
            free(newFeature);
        }
        free(uv_data);
        free(uv_int_data);
        for (j=0; j<[features_info count]; ++j){
            Feature* feat = (Feature*)features_info[j];
            if (feat->h!=NULL){
                free(feat->h);
                feat->h = NULL;
            }
        }

    }
    return success;
}


/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
/***************************************************************************************************************************************/



/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
/***************************************************************************************************************************************/

bool addOneFeatureCovarianceInverseDepth(Filter *filter, double* uvd, double *Xv, double std_pxl, double std_rho, Camera* cam, double *P_RES)
{
    int i,j;
    
    double fku=cam->k->at<double>(0,0);
    double fkv=cam->k->at<double>(1,1);
    double U0=cam->k->at<double>(0,2);
    double V0=cam->k->at<double>(1,2);
    
    double *q_wc = (double *)malloc(1*4*sizeof(double));
    for (i=3; i<7; i++)
        q_wc[i-3] = *(Xv+i);
    
    double *R_wc = (double *)malloc(3*3*sizeof(double));
    q2r(q_wc, R_wc);
    
    double *uvu = (double *)malloc(1*2*sizeof(double));
    undistortFM(uvd, 1, cam, uvu);
    double uu = *(uvu+0);
    double vu = *(uvu+1);
    
    double x_c = -(U0-uu)/fku;
    double y_c = -(V0-vu)/fkv;
    double z_c = 1;
    
    cv::Mat *xyz_c = new cv::Mat(3,1,CV_64FC1,0.0);
    xyz_c->at<double>(0,0) = x_c;
    xyz_c->at<double>(1,0) = y_c;
    xyz_c->at<double>(2,0) = z_c;
    
    cv::Mat *R_wc_mat = new cv::Mat(3,3,CV_64FC1,R_wc, cv::Mat::AUTO_STEP);
    cv::Mat *xyz_w = new cv::Mat(3,1,CV_64FC1,0.0);
    (*xyz_w) = (*R_wc_mat) * (*xyz_c);
    double x_w = xyz_w->at<double>(0,0);
    double y_w = xyz_w->at<double>(1,0);
    double z_w = xyz_w->at<double>(2,0);
    
    cv::Mat *dtheta_dgw = new cv::Mat(1,3,CV_64FC1,0.0);
    dtheta_dgw->at<double>(0,0) = z_w/(pow(x_w,2)+pow(z_w,2));
    dtheta_dgw->at<double>(0,1) = 0;
    dtheta_dgw->at<double>(0,2) = -x_w/(pow(x_w,2)+pow(z_w,2));
    
    cv::Mat *dphi_dgw = new cv::Mat(1,3,CV_64FC1,0.0);
    dphi_dgw->at<double>(0,0) = (x_w*y_w)/ ((pow(x_w,2)+pow(y_w,2)+pow(z_w,2))*sqrt(pow(x_w,2)+pow(z_w,2)));
    dphi_dgw->at<double>(0,1) = -sqrt(pow(x_w,2)+pow(z_w,2))/ (pow(x_w,2)+pow(y_w,2)+pow(z_w,2));
    dphi_dgw->at<double>(0,2) = (z_w*y_w) / ((pow(x_w,2)+pow(y_w,2)+pow(z_w,2))*sqrt(pow(x_w,2)+pow(z_w,2)));
    
    cv::Mat *dgw_dqwr = new cv::Mat(3,4,CV_64FC1,0.0);
    dRqTimesABydq(q_wc, xyz_c, dgw_dqwr);
    
    cv::Mat *dtheta_dqwr = new cv::Mat(1,4,CV_64FC1,0.0);
    (*dtheta_dqwr) = (*dtheta_dgw) * (*dgw_dqwr);
    
    cv::Mat *dphi_dqwr = new cv::Mat(1,4,CV_64FC1,0.0);
    (*dphi_dqwr) = (*dphi_dgw) * (*dgw_dqwr);
    
    cv::Mat *dy_dqwr = new cv::Mat(6,4,CV_64FC1,0.0);
    for (j=0; j<4; j++)
    {
        dy_dqwr->at<double>(3,j) = dtheta_dqwr->at<double>(0,j);
        dy_dqwr->at<double>(4,j) = dphi_dqwr->at<double>(0,j);
    }
    
    cv::Mat *dy_drw = new cv::Mat(6,3,CV_64FC1,0.0);
    dy_drw->at<double>(0,0) = 1;
    dy_drw->at<double>(1,1) = 1;
    dy_drw->at<double>(2,2) = 1;
    
    cv::Mat *dy_dxv = new cv::Mat(6,13,CV_64FC1,0.0);
    for (i=0; i<6; i++)
    {
        for (j=0; j<3; j++)
            dy_dxv->at<double>(i,j) = dy_drw->at<double>(i,j);
        
        for (j=3; j<7; j++)
            dy_dxv->at<double>(i,j) = dy_dqwr->at<double>(i,j-3);
    }
    
    cv::Mat *dy_dxv_transpose = new cv::Mat(13,6,CV_64FC1,0.0);
    cv::transpose(*dy_dxv, *dy_dxv_transpose);
    
    cv::Mat *dyprima_dgw = new cv::Mat(5,3,CV_64FC1,0.0);
    for (j=0; j<3; j++)
    {
        dyprima_dgw->at<double>(3,j) = dtheta_dgw->at<double>(0,j);
        dyprima_dgw->at<double>(4,j) = dphi_dgw->at<double>(0,j);
    }
    
    cv::Mat *dgw_dgc = new cv::Mat(3,3,CV_64FC1,R_wc_mat);
    
    cv::Mat *dgc_dhu = new cv::Mat(2,3,CV_64FC1,0.0);
    dgc_dhu->at<double>(0,0) = 1/fku;
    dgc_dhu->at<double>(1,1) = 1/fkv;
    cv::transpose(*dgc_dhu, *dgc_dhu);
    
    cv::Mat *dhu_dhd = new cv::Mat(2,2,CV_64FC1,0.0);
    JacobUndistorFM(cam, uvd, dhu_dhd);
    
    cv::Mat *dyprima_dhd = new cv::Mat(5,2,CV_64FC1,0.0);
    (*dyprima_dhd) = (*dyprima_dgw) * (*dgw_dgc) * (*dgc_dhu) * (*dhu_dhd);
    
    cv::Mat *dy_dhd = new cv::Mat(6,3,CV_64FC1,0.0);
    for (i=0; i<5; i++)
        for (j=0; j<2; j++)
            dy_dhd->at<double>(i,j) = dyprima_dhd->at<double>(i,j);
    dy_dhd->at<double>(5,2) = 1;
    
    cv::Mat *dy_dhd_transpose = new cv::Mat(3,6,CV_64FC1,0.0);
    cv::transpose(*dy_dhd, *dy_dhd_transpose);
    
    cv::Mat *Ri = new cv::Mat(2,2,CV_64FC1,0.0);
    Ri->at<double>(0,0) = pow(std_pxl,2);
    Ri->at<double>(1,1) = pow(std_pxl,2);
    
    cv::Mat *Padd = new cv::Mat(3,3,CV_64FC1,0.0);
    Padd->at<double>(0,0) = Ri->at<double>(0,0);
    Padd->at<double>(1,1) = Ri->at<double>(1,1);
    Padd->at<double>(2,2) = pow(std_rho,2);
    
    int size = filter->state_size;
    
    cv::Mat *P_xv = new cv::Mat(13,13,CV_64FC1,0.0);
    for (i=0; i<13; i++)
        for (j=0; j<13; j++)
            P_xv->at<double>(i,j) = filter->p_k_k[i*size+j];
    
    cv::Mat *P_yxv = new cv::Mat(size-13,13,CV_64FC1,0.0);
    for (i=13; i<size; i++)
        for (j=0; j<13; j++)
            P_yxv->at<double>(i-13,j) = filter->p_k_k[i*size+j];
    
    cv::Mat *P_y = new cv::Mat(size-13,size-13,CV_64FC1,0.0);
    for (i=13; i<size; i++)
        for (j=13; j<size; j++)
            P_y->at<double>(i-13,j-13) = filter->p_k_k[i*size+j];
    
    cv::Mat *P_xvy = new cv::Mat(13,size-13,CV_64FC1,0.0);
    for (i=0; i<13; i++)
        for (j=13; j<size; j++)
            P_xvy->at<double>(i,j-13) = filter->p_k_k[i*size+j];
    
    cv::Mat *P_RES_1_3 = new cv::Mat(13,6,CV_64FC1,0.0);
    (*P_RES_1_3) = (*P_xv) * (*dy_dxv_transpose);
    
    cv::Mat *P_RES_2_3 = new cv::Mat(size-13,6,CV_64FC1,0.0);
    (*P_RES_2_3) = (*P_yxv) * (*dy_dxv_transpose);
    
    cv::Mat *P_RES_3_1 = new cv::Mat(6,13,CV_64FC1,0.0);
    (*P_RES_3_1) = (*dy_dxv) * (*P_xv);
    
    cv::Mat *P_RES_3_2 = new cv::Mat(6,size-13,CV_64FC1,0.0);
    (*P_RES_3_2) = (*dy_dxv) * (*P_xvy);
    
    cv::Mat *P_RES_3_3 = new cv::Mat(6,6,CV_64FC1,0.0);
    (*P_RES_3_3) = (*dy_dxv) * (*P_xv) * (*dy_dxv_transpose) + (*dy_dhd) * (*Padd) * (*dy_dhd_transpose);
    
    for (i=0; i<13; i++)
    {
        for (j=0; j<13; j++)
            P_RES[i*(size+6)+j] = P_xv->at<double>(i,j);
        
        for (j=13; j<size; j++)
            P_RES[i*(size+6)+j] = P_xvy->at<double>(i,j-13);
        
        for (j=size; j<size+6; j++)
            P_RES[i*(size+6)+j] = P_RES_1_3->at<double>(i,j-size);
    }
    
    for (i=13; i<size; i++)
    {
        for (j=0; j<13; j++)
            P_RES[i*(size+6)+j] = P_yxv->at<double>(i-13,j);
        
        for (j=13; j<size; j++)
            P_RES[i*(size+6)+j] = P_y->at<double>(i-13,j-13);
        
        for (j=size; j<size+6; j++)
            P_RES[i*(size+6)+j] = P_RES_2_3->at<double>(i-13,j-size);
    }
    
    for (i=size; i<size+6; i++)
    {
        for (j=0; j<13; j++)
            P_RES[i*(size+6)+j] = P_RES_3_1->at<double>(i-size,j);
        
        for (j=13; j<size; j++)
            P_RES[i*(size+6)+j] = P_RES_3_2->at<double>(i-size,j-13);
        
        for (j=size; j<size+6; j++)
            P_RES[i*(size+6)+j] = P_RES_3_3->at<double>(i-size,j-size);
    }
    
    free(q_wc);
    free(R_wc);
    free(uvu);
    
    xyz_c->release();
    xyz_w->release();
    R_wc_mat->release();
    dtheta_dgw->release();
    dphi_dgw->release();
    dgw_dqwr->release();
    dtheta_dqwr->release();
    dphi_dqwr->release();
    dy_dqwr->release();
    dy_drw->release();
    dy_dxv->release();
    dy_dxv_transpose->release();
    dyprima_dgw->release();
    dgw_dgc->release();
    dgc_dhu->release();
    dhu_dhd->release();
    dyprima_dhd->release();
    dy_dhd->release();
    dy_dhd_transpose->release();
    Ri->release();
    Padd->release();
    P_xv->release();
    P_yxv->release();
    P_y->release();
    P_xvy->release();
    P_RES_1_3->release();
    P_RES_2_3->release();
    P_RES_3_1->release();
    P_RES_3_2->release();
    
    return true;
}

bool add_features_inverse_depth(double* uvd, int nNewFeat, Filter *filter, Camera* cam, double std_pxl, double initial_rho, double std_rho, double* newFeature)
{
    int i,j;
    
    if (uvd == NULL)
        return true;
    
    else
    {
        double *Xv = (double *)malloc(1*13*sizeof(double));
        for (i=0; i<13; i++)
            Xv[i] = filter->x_k_k[i];
        
        int size = filter->state_size;
        
        hinv(uvd, Xv, cam, initial_rho, newFeature);
        double *X_RES = (double *)malloc(1*(size+6)*sizeof(double));
        for (i=0; i<size; i++)
            X_RES[i] = filter->x_k_k[i];
        for (i=size; i<size+6; i++)
            X_RES[i] = *(newFeature+(i-size));
        
        double *P_RES = (double *)malloc((size+6)*(size+6)*sizeof(double));
        addOneFeatureCovarianceInverseDepth(filter, uvd, Xv, std_pxl, std_rho, cam, P_RES);

        filter->state_size = size+6;
        
        free(filter->x_k_k);
        free(filter->p_k_k);

        filter->x_k_k = (double*)malloc(filter->state_size*sizeof(double));
        filter->p_k_k = (double*)malloc(filter->state_size*filter->state_size*sizeof(double));
                                
        for (i=0; i<size+6; i++)
        {
            filter->x_k_k[i] = X_RES[i];
            
            for (j=0; j<size+6; j++)
                filter->p_k_k[i*(size+6)+j] = P_RES[i*(size+6)+j];
        }

        free(Xv);
        free(X_RES);
        free(P_RES);
    }
    
    return true;
}


/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
/***************************************************************************************************************************************/

bool hinv(double* uvd, double* Xv, Camera* cam, double initialRho, double* newFeature)
{
    int i;
    
    double fku = cam->k->at<double>(0,0);
    double fkv = cam->k->at<double>(1,1);
    double U0 = cam->k->at<double>(0,2);
    double V0 = cam->k->at<double>(1,2);
    
    double *uvu = (double *)malloc(1*2*sizeof(double));
    undistortFM(uvd, 1, cam, uvu);
    double u = *(uvu+0);
    double v = *(uvu+1);
    
    double *r_w = (double *)malloc(1*3*sizeof(double));
    for (i=0; i<3; i++)
        r_w[i] = *(Xv+i);
    
    double *q_wr = (double *)malloc(1*4*sizeof(double));
    for (i=3; i<7; i++)
        q_wr[i-3] = *(Xv+i);
    
    double h_LR_x = -(U0-u)/fku;
    double h_LR_y = -(V0-v)/fkv;
    double h_LR_z = 1;
    
    cv::Mat *h_LR = new cv::Mat(3,1,CV_64FC1,0.0);
    h_LR->at<double>(0,0) = h_LR_x;
    h_LR->at<double>(1,0) = h_LR_y;
    h_LR->at<double>(2,0) = h_LR_z;
    
    double *q_wr_0 = (double *)malloc(3*3*sizeof(double));
    q2r(q_wr, q_wr_0);
    
    cv::Mat *q_wr_0_mat = new cv::Mat(3,3,CV_64FC1,q_wr_0);
    cv::Mat *n = new cv::Mat(3,1,CV_64FC1,0.0);
    (*n) = (*q_wr_0_mat)*(*h_LR);
    
    double nx = n->at<double>(0,0);
    double ny = n->at<double>(1,0);
    double nz = n->at<double>(2,0);
    
    for (i=0; i<3; i++)
        *(newFeature+i) = r_w[i];
    
    *(newFeature+3) = atan2(nx, nz);
    *(newFeature+4) = atan2(-ny, sqrt(pow(nx, 2)+pow(nz, 2)));
    *(newFeature+5) = initialRho;
    
    free(uvu);
    free(r_w);
    free(q_wr);
    free(q_wr_0);
    
    h_LR->release();
    q_wr_0_mat->release();
    
    return true;
 }


/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
/***************************************************************************************************************************************/

/***************************************************************************************************************************************/
/************************************************************************************************************Predict Camera Measurement*/
/***************************************************************************************************************************************/

bool q2r(double *q, double *R)
{
    if (q == NULL || R == NULL)
        return false;
    
    double x = *(q+1);
    double y = *(q+2);
    double z = *(q+3);
    double r = *(q+0);
    
    R[0] = r*r + x*x - y*y - z*z;
    R[1] = 2*(x*y-r*z);
    R[2] = 2*(z*x+r*y);
    R[3] = 2*(x*y+r*z);
    R[4] = r*r - x*x + y*y - z*z;
    R[5] = 2*(y*z-r*x);
    R[6] = 2*(z*x-r*y);
    R[7] = 2*(y*z+r*x);
    R[8] = r*r - x*x - y*y + z*z;
    
    return true;
}

/***************************************************************************************************************************************/

bool hu(double *uv_u, double *yi, int yi_numColums, Camera *cam)
{
    double u0 = cam->Cx;
    double v0 = cam->Cy;
    double f = cam->f;
    double ku = 1/cam->dx;
    double kv = 1/cam->dy;
    
    int i;
    
    for (i=0; i<yi_numColums; i++)
    {
        uv_u[0*yi_numColums+i] = u0 + (yi[0*yi_numColums+i]/yi[2*yi_numColums+i])*f*ku;
        uv_u[1*yi_numColums+i] = v0 + (yi[1*yi_numColums+i]/yi[2*yi_numColums+i])*f*kv;
    }
    
    return true;
}

/***************************************************************************************************************************************/

bool distortFM(double *uvd,double *uv,int uv_ncols,Camera* cam)
{
    double k1 = cam->k1;
    double k2 = cam->k2;
    double Cx = cam->Cx;
    double Cy = cam->Cy;
    double dx = cam->dx;
    double dy = cam->dy;
    
    double xu,yu,ru,rd,D,xd,yd;
    
    double f;
    double f_p;
    
    int i,j;
    
    for (i=0;i<uv_ncols;i++)
    {
        xu=(uv[i]-Cx)*dx;
        yu=(uv[i+uv_ncols]-Cy)*dy;
        
        ru=sqrt(pow(xu,2)+pow(yu,2));
        rd=ru/(1+k1*pow(ru,2)+k2*pow(ru,4));
        
        for (j=0;j<10;j++)
        {
            f=rd+k1*pow(rd,3)+k2*pow(rd,5)-ru;
            f_p=1+3*k1*pow(rd,2)+5*k2*pow(rd,4);
            rd=rd-f/f_p;
        }
        
        D=1+k1*pow(rd,2)+k2*pow(rd,4);
        xd=xu/D;
        yd=yu/D;
        
        uvd[i]=xd/dx+Cx;
        uvd[i+uv_ncols]=yd/dy+Cy;
    }
    
    return true;
}

/***************************************************************************************************************************************/

bool hiCartesian(Camera *cam, double* yi, double* t_wc, double* r_wc, double *zi)
{
    const double pi = 3.1415926;
    
    cv::Mat *r_wc_mat = new cv::Mat(3,3,CV_64FC1,r_wc);
    cv::Mat *r_cw = new cv::Mat(3,3,CV_64FC1,0.0);
    invert(*r_wc_mat, *r_cw);
    
    cv::Mat *yi_mat =  new cv::Mat(3,1,CV_64FC1,yi);
    cv::Mat *t_wc_mat = new cv::Mat(3,1,CV_64FC1,t_wc);
    cv::Mat *hrl_mat = new cv::Mat(3,1,CV_64FC1,0.0);
    *hrl_mat = (*r_cw) * ( *yi_mat - *t_wc_mat);
    
    double a0=atan2(hrl_mat->at<double>(0,0), hrl_mat->at<double>(2,0))*180/pi;
    double a1=atan2(hrl_mat->at<double>(1,0), hrl_mat->at<double>(2,0))*180/pi;
    
    if (a0<-60 || a0>+60 || a1<-60 || a1>+60)
    {
        zi = NULL;
        
        return true;
    }
    
    double *hrl = (double *)malloc(3*1*sizeof(double));
    hrl[0] = hrl_mat->at<double>(0,0);
    hrl[1] = hrl_mat->at<double>(1,0);
    hrl[2] = hrl_mat->at<double>(2,0);
    
    double *uv_u = (double *)malloc(2*1*sizeof(double));
    hu(uv_u, hrl, 1, cam);
    
    double *uv_d = (double *)malloc(2*1*sizeof(double));
    distortFM(uv_d, uv_u, 1, cam);
    
    if ((uv_d[0]>0) && (uv_d[0]<cam->nCols) &&
        (uv_d[1]>0) && (uv_d[1]<cam->nRows))
    {
        zi[0] = uv_d[0];
        zi[1] = uv_d[1];
        
        return true;
    }
    else
    {
        zi = NULL;
        
        return true;
    }
    
    r_wc_mat->release();
    r_cw->release();
    yi_mat->release();
    t_wc_mat->release();
    hrl_mat->release();
    
    free(hrl);
    free(uv_u);
    free(uv_d);
    
    return true;
}

/***************************************************************************************************************************************/

bool hiInverseDepth(Camera *cam, double* yinit, double* t_wc, double* r_wc, double* zi)
{
    const double pi = 3.1415926;
    
    cv::Mat *r_wc_mat = new cv::Mat(3,3,CV_64FC1,r_wc);
    cv::Mat *r_cw = new cv::Mat(3,3,CV_64FC1,0.0);
    transpose(*r_wc_mat, *r_cw);
    
    cv::Mat *yi = new cv::Mat(3,1,CV_64FC1,0.0);
    yi->at<double>(0,0) = *(yinit+0);
    yi->at<double>(1,0) = *(yinit+1);
    yi->at<double>(2,0) = *(yinit+2);
    
    double theta = *(yinit+3);
    double phi = *(yinit+4);
    double rho = *(yinit+5);
    
    cv::Mat *mi = new cv::Mat(3,1,CV_64FC1,0.0);
    m(theta, phi, mi);
    
    cv::Mat *hrl_mat = new cv::Mat(3,1,CV_64FC1,0.0);
    cv::Mat *t_wc_mat = new cv::Mat(3,1,CV_64FC1,t_wc);
    *hrl_mat = (*r_cw) * ((( *yi - *t_wc_mat )*rho)+ *mi);
    
    double a0=atan2(hrl_mat->at<double>(0,0), hrl_mat->at<double>(2,0))*180/pi;
    double a1=atan2(hrl_mat->at<double>(1,0), hrl_mat->at<double>(2,0))*180/pi;
    
    if (a0<-60 || a0>+60 || a1<-60 || a1>+60)
    {
        zi = NULL;
        
        return true;
    }
    
    double *hrl = (double *)malloc(3*1*sizeof(double));
    hrl[0] = hrl_mat->at<double>(0,0);
    hrl[1] = hrl_mat->at<double>(1,0);
    hrl[2] = hrl_mat->at<double>(2,0);
    
    double *uv_u = (double *)malloc(2*1*sizeof(double));
    hu(uv_u, hrl, 1, cam);
    
    double *uv_d = (double *)malloc(2*1*sizeof(double));
    distortFM(uv_d, uv_u, 1, cam);
    
    if ((uv_d[0]>0) && (uv_d[0]<cam->nCols) &&
        (uv_d[1]>0) && (uv_d[1]<cam->nRows))
    {
        zi[0] = uv_d[0];
        zi[1] = uv_d[1];
        
        return true;
    }
    else
    {
        zi = NULL;
        
        return true;
    }
    
    r_wc_mat->release();
    r_cw->release();
    yi->release();
    mi->release();
    hrl_mat->release();
    t_wc_mat->release();
    
    free(hrl);
    free(uv_u);
    free(uv_d);
    
    return true;
}

/***************************************************************************************************************************************/

bool predict_camera_measurements(NSMutableArray *features_info,Filter *filter,Camera *cam)
{
    int i,j;
    
    double *t_wc = (double *)malloc(1*3*sizeof(double));
    for (i=0; i<3; i++)
        t_wc[i] = filter->x_k_k[i];
    
    double *r_wc_0 = (double *)malloc(1*4*sizeof(double));
    for (i=3; i<7; i++)
        r_wc_0[i-3] = filter->x_k_k[i];
    
    double *r_wc = (double *)malloc(3*3*sizeof(double));
    q2r(r_wc_0, r_wc);
    
    double size = filter->state_size;
    double *features_0 = (double *)malloc(1*(size-13)*sizeof(double));
    int currFeaturesIdx = 0;
    
    for (i=13; i<size; i++)
        features_0[i-13] = filter->x_k_k[i];
    
    for (i=0; i<(int)[features_info count]; i++)
    {
        Feature *current=(Feature*)features_info[i];
        
        if ([current->type isEqualToString:@"cartesian"])
        {
            double *yi = (double *)malloc(1*3*sizeof(double));
            for (j=0; j<3; j++)
                yi[j] = features_0[j];
            
            double *features = (double *)malloc(1*(size-16)*sizeof(double));
            for (j=3;j<size-13;j++)
                features[j-3] = features_0[j];
            
            double *hi = (double *)malloc(1*2*sizeof(double));
            hiCartesian(cam, yi, t_wc, r_wc, hi);
            
            if (hi != NULL)
            {
                if (current->h == NULL)
                    current->h = (double*)malloc(1*2*sizeof(double));
                
                current->h[0] = hi[0];
                current->h[1] = hi[1];
            }
            
            free(features);
            free(hi);
        }
        
        if ([current->type isEqualToString:@"inversedepth"])
        {
            double *yi = (double *)malloc(1*6*sizeof(double));
            for (j=currFeaturesIdx; j<currFeaturesIdx+6; j++)
                yi[j-currFeaturesIdx] = features_0[j];
            currFeaturesIdx = currFeaturesIdx+6;
            
            
            double *features = (double *)malloc((size-currFeaturesIdx)*sizeof(double));
            
            for (j=currFeaturesIdx; j<size; j++)
                features[j-currFeaturesIdx] = features_0[j];
            
            double *hi = (double *)malloc(1*2*sizeof(double));
            hiInverseDepth(cam, yi, t_wc, r_wc, hi);
            
            if (hi != NULL)
            {
                if (current->h == NULL)
                    current->h = (double*)malloc(1*2*sizeof(double));
                
                current->h[0] = hi[0];
                current->h[1] = hi[1];
            }
            
            free(features);
            free(hi);
            free(yi);
        }
    }
    
    free(t_wc);
    free(r_wc_0);
    free(r_wc);
    free(features_0);
    
    return false;
}

/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
/***************************************************************************************************************************************/

 bool undistortFM(double *uvd, int uvd_nCols, Camera *cam, double *uvu)
{
    double Cx = cam->Cx;
    double Cy = cam->Cy;
    double k1 = cam->k1;
    double k2 = cam->k2;
    double dx = cam->dx;
    double dy = cam->dy;
    
   // double xd = ( uvd[0] - Cx )*dx;
    //double yd = ( uvd[1] - Cy )*dy;
    double* xd = (double*)malloc(2*uvd_nCols*sizeof(double));
    double* yd = (double*)malloc(2*uvd_nCols*sizeof(double));
    
    for (int a = 0; a < uvd_nCols; a++) {
        xd[a] = (uvd[(uvd_nCols*0)+a] - Cx)*dx;
        yd[a] = (uvd[(uvd_nCols*1)+a] - Cy)*dy;
    }
    
    
        
    //double rd = sqrt( pow(*xd, 2) + pow(*yd, 2) );
    
    double* rd = (double*)malloc(2*uvd_nCols*sizeof(double));
    
    for (int i = 0; i < uvd_nCols; i++) {
        rd[i] = sqrt((xd[i]*xd[i])+ yd[i]*yd[i]);
    }
    
    //double D = 1 + k1* pow(rd, 2) + k2* pow(rd, 4);
    double* D = (double*)malloc(2*uvd_nCols*sizeof(double));
    for (int j = 0; j < uvd_nCols; j++) {
        D[j] = 1 + k1*(rd[j]*rd[j]) + k2*(rd[j]*rd[j]*rd[j]*rd[j]);
    }
    
   // double xu = xd*D;
    //double yu = yd*D;
    
    double* xu = (double*)malloc(2*uvd_nCols*sizeof(double));
    double* yu = (double*)malloc(2*uvd_nCols*sizeof(double));
    
    for (int r = 0; r < uvd_nCols; r++)
    {
        xu[r] = xd[r]*D[r];
        yu[r] = yd[r]*D[r];
    }
    
    
   // uvu[0] = xu/dx + Cx;
    //uvu[1] = yu/dy + Cy;
    
    for (int d = 0; d < uvd_nCols; d++) {
        uvu[(uvd_nCols*0)+d] = (xu[d]/dx) + Cx;
        uvu[(uvd_nCols*1)+d] = (yu[d]/dy) + Cy;
    }
    
    free (xd);
    free(yd);
    free(rd);
    free(D);
    free(xu);
    free(yu);
    return true;
    
}

/***************************************************************************************************************************************/

void dR_by_dq0(double *q, double *dR_by_dq0RES)
{
    double q0 = *(q+0);
    double qx = *(q+1);
    double qy = *(q+2);
    double qz = *(q+3);
    
    dR_by_dq0RES[0] = +2*q0;
    dR_by_dq0RES[1] = -2*qz;
    dR_by_dq0RES[2] = +2*qy;
    dR_by_dq0RES[3] = +2*qz;
    dR_by_dq0RES[4] = +2*q0;
    dR_by_dq0RES[5] = -2*qx;
    dR_by_dq0RES[6] = -2*qy;
    dR_by_dq0RES[7] = +2*qx;
    dR_by_dq0RES[8] = +2*q0;
}

void dR_by_dqx(double *q, double *dR_by_dqxRES)
{
    double q0 = *(q+0);
    double qx = *(q+1);
    double qy = *(q+2);
    double qz = *(q+3);
    
    dR_by_dqxRES[0] = +2*qx;
    dR_by_dqxRES[1] = +2*qy;
    dR_by_dqxRES[2] = +2*qz;
    dR_by_dqxRES[3] = +2*qy;
    dR_by_dqxRES[4] = -2*qx;
    dR_by_dqxRES[5] = -2*q0;
    dR_by_dqxRES[6] = +2*qz;
    dR_by_dqxRES[7] = +2*q0;
    dR_by_dqxRES[8] = -2*qx;
    
}

void dR_by_dqy(double *q, double *dR_by_dqyRES)
{
    double q0 = *(q+0);
    double qx = *(q+1);
    double qy = *(q+2);
    double qz = *(q+3);
    
    dR_by_dqyRES[0] = -2*qy;
    dR_by_dqyRES[1] = +2*qx;
    dR_by_dqyRES[2] = +2*q0;
    dR_by_dqyRES[3] = +2*qx;
    dR_by_dqyRES[4] = +2*qy;
    dR_by_dqyRES[5] = +2*qz;
    dR_by_dqyRES[6] = -2*q0;
    dR_by_dqyRES[7] = +2*qz;
    dR_by_dqyRES[8] = -2*qy;
    
}

void dR_by_dqz(double *q, double *dR_by_dqzRES)
{
    double q0 = *(q+0);
    double qx = *(q+1);
    double qy = *(q+2);
    double qz = *(q+3);
    
    dR_by_dqzRES[0] = -2*qz;
    dR_by_dqzRES[1] = -2*q0;
    dR_by_dqzRES[2] = +2*qx;
    dR_by_dqzRES[3] = +2*q0;
    dR_by_dqzRES[4] = -2*qz;
    dR_by_dqzRES[5] = +2*qy;
    dR_by_dqzRES[6] = +2*qx;
    dR_by_dqzRES[7] = +2*qy;
    dR_by_dqzRES[8] = +2*qz;
    
}

bool dRqTimesABydq(double *q, cv::Mat* aMat, cv::Mat* dRqTimesABydqRES)

{
    if (aMat == NULL || dRqTimesABydqRES == NULL)
        return false;
    
    if (!aMat->data || !dRqTimesABydqRES->data)
        return false;
    
    int i;
    
    double *TempR1 = (double *)malloc(3*3*sizeof(double));
    cv::Mat Temp1 = cv::Mat(3,1,CV_64FC1,0.0);
    dR_by_dq0(q, TempR1);
    cv::Mat* TempR1_mat = new cv::Mat(3,3,CV_64FC1,TempR1);
    Temp1 = (*TempR1_mat) * (*aMat);
    
    double *TempR2 = (double *)malloc(3*3*sizeof(double));
    cv::Mat Temp2 = cv::Mat(3,1,CV_64FC1,0.0);
    dR_by_dqx(q, TempR2);
    cv::Mat* TempR2_mat = new cv::Mat(3,3,CV_64FC1,TempR2);
    Temp2 = (*TempR2_mat) * (*aMat);
    
    double *TempR3 = (double *)malloc(3*3*sizeof(double));
    cv::Mat Temp3 = cv::Mat(3,1,CV_64FC1,0.0);
    dR_by_dqy(q, TempR3);
    cv::Mat* TempR3_mat = new cv::Mat(3,3,CV_64FC1,TempR3);
    Temp3 = (*TempR3_mat) * (*aMat);
    
    double *TempR4 = (double *)malloc(3*3*sizeof(double));
    cv::Mat Temp4 = cv::Mat(3,1,CV_64FC1,0.0);
    dR_by_dqz(q, TempR4);
    cv::Mat* TempR4_mat = new cv::Mat(3,3,CV_64FC1,TempR4);
    Temp4 = (*TempR4_mat) * (*aMat);
    
    for (i=0;i<3;i++)
    {
        dRqTimesABydqRES->at<double>(i,0) = Temp1.at<double>(i,0);
        dRqTimesABydqRES->at<double>(i,1) = Temp2.at<double>(i,0);
        dRqTimesABydqRES->at<double>(i,2) = Temp3.at<double>(i,0);
        dRqTimesABydqRES->at<double>(i,3) = Temp4.at<double>(i,0);
    }
    
    free(TempR1);
    free(TempR2);
    free(TempR3);
    free(TempR4);
    
    Temp1.release();
    Temp2.release();
    Temp3.release();
    Temp4.release();
    
    TempR1_mat->release();
    TempR2_mat->release();
    TempR3_mat->release();
    TempR4_mat->release();
    
    return true;
}

/***************************************************************************************************************************************/

void JacobUndistorFM(Camera* cam, double* uvd, cv::Mat* J_unidistor)
{
    double k1 = cam->k1;
    double k2 = cam->k2;
    double Cx = cam->Cx;
    double Cy = cam->Cy;
    double dx = cam->dx;
    double dy = cam->dy;
    
    double ud = uvd[0];
    double vd = uvd[1];
    
    double xd = (ud-Cx)*dx;
    double yd = (vd-Cy)*dy;
    
    double rd2 = xd*xd+yd*yd;
    double rd4 = rd2*rd2;
    
    double uu_ud = (1 + k1*rd2 + k2*rd4)+(ud-Cx)*(k1 + 2*k2*rd2)*(2*(ud-Cx)*dx*dx);
    double vu_vd = (1 + k1*rd2 + k2*rd4)+(vd-Cy)*(k1 + 2*k2*rd2)*(2*(vd-Cy)*dy*dy);
    
    double uu_vd = (ud-Cx)*(k1 + 2*k2*rd2)*(2*(vd-Cy)*dy*dy);
    double vu_ud = (vd-Cy)*(k1 + 2*k2*rd2)*(2*(ud-Cx)*dx*dx);
    
    J_unidistor->at<double>(0,0)=uu_ud;
    J_unidistor->at<double>(0,1)=uu_vd;
    J_unidistor->at<double>(1,0)=vu_ud;
    J_unidistor->at<double>(1,1)=vu_vd;
}

/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
/***************************************************************************************************************************************/


/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
/***************************************************************************************************************************************/


/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
/***************************************************************************************************************************************/


/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
//TODO: Get rid of cv::Mats when they're not necessary


/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
/***************************************************************************************************************************************/

/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
/***************************************************************************************************************************************/







/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
/***************************************************************************************************************************************/

bool normJac(double* J, double* q){
    if (J==NULL || q==NULL)
        return false;
    
    double r = *q;
    double x = *(q+1);
    double y = *(q+2);
    double z = *(q+3);
    
    double xx = x*x;
    double yy = y*y;
    double zz = z*z;
    double rr = r*r;
    
    double J_coefficient = pow((rr + xx + yy + zz),(-3/2));
    
    J[0] = J_coefficient * (xx + yy + zz);
    J[1] = J_coefficient * (-r*x);
    J[2] = J_coefficient * (-r*y);
    J[3] = J_coefficient * (-r*z);
    J[4] = J_coefficient * (-x*r);
    J[5] = J_coefficient * (rr + yy + zz);
    J[6] = J_coefficient * (-x*y);
    J[7] = J_coefficient * (-x*z);
    J[8] = J_coefficient * (-y*r);
    J[9] = J_coefficient * (-y*x);
    J[10] = J_coefficient * (rr + xx+ zz);
    J[11] = J_coefficient * (-y*z);
    J[12] = J_coefficient * (-z*r);
    J[13] = J_coefficient * (-z*x);
    J[14] = J_coefficient * (-z*y);
    J[15] = J_coefficient * (rr + xx + yy);
    
    return true;
}

/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
/***************************************************************************************************************************************/


bool add_feature_to_info_vector(int* uv, cv::Mat* im_k, double* X_RES, int X_RES_size, NSMutableArray* features_info, int step, double *newFeature)
{
    int half_patch_size_when_initialized = 20;
    int half_patch_size_when_matching = 6;
    
    Feature *current = [[Feature alloc] init];
    
    int i,j;
    int ia = *(uv+1)-half_patch_size_when_initialized-1;
    int ib = *(uv+1)+half_patch_size_when_initialized;
    int ja = *(uv+0)-half_patch_size_when_initialized-1;
    int jb = *(uv+0)+half_patch_size_when_initialized;
    current->patch_when_initialized = new cv::Mat(2*half_patch_size_when_initialized+1,2*half_patch_size_when_initialized+1,CV_8UC1,0.0);
    for (i=ia; i<ib; i++)
    {
        for (j=ja; j<jb; j++)
        {
            current->patch_when_initialized->at<uchar>(i-ia,j-ja) = im_k->at<uchar>(i,j);
        }
    }
    
    current->patch_when_matching = new cv::Mat(2*half_patch_size_when_matching+1,2*half_patch_size_when_matching+1,CV_8UC1,0.0);
    
    current->r_wc_when_initialized = (double*)malloc(3*sizeof(double));
    current->r_wc_when_initialized[0] = *(X_RES+0);
    current->r_wc_when_initialized[1] = *(X_RES+1);
    current->r_wc_when_initialized[2] = *(X_RES+2);
    
    double *X_RES_0 = (double *)malloc(1*4*sizeof(double));
    X_RES_0[0] = X_RES[3];
    X_RES_0[1] = X_RES[4];
    X_RES_0[2] = X_RES[5];
    X_RES_0[3] = X_RES[6];
    if (current->R_wc_when_initialized==NULL){
        current->R_wc_when_initialized = (double*)malloc(16*sizeof(double));
    }
    bool success = q2r(X_RES_0, current->R_wc_when_initialized);
    
    current->uv_when_initialized[0] = *(uv+0);
    current->uv_when_initialized[1] = *(uv+1);
    
    current->half_patch_size_when_initialized = half_patch_size_when_initialized;
    current->half_patch_size_when_matching = half_patch_size_when_matching;
    
    current->times_predicted = 0;
    current->times_measured = 0;
    
    current->init_frame = step;
    current->init_measurement[0] = *(uv+0);
    current->init_measurement[1] = *(uv+1);
    
    current->type = [NSString stringWithFormat:@"inversedepth"];
    
    current->yi[0] = *(newFeature+0);
    current->yi[1] = *(newFeature+1);
    current->yi[2] = *(newFeature+2);
    current->yi[3] = *(newFeature+3);
    current->yi[4] = *(newFeature+4);
    current->yi[5] = *(newFeature+5);
    
    current->individually_compatible = 0;
    current->low_innovation_inlier = 0;
    current->high_innovation_inlier = 0;
    
    current->z = NULL;
    current->h = NULL;
    current->H = NULL;
    current->S = NULL;
    
    current->state_size = 6;
    current->measurement_size = 2;
    
    current->R[0][0] = 1;
    current->R[0][1] = 0;
    current->R[1][0] = 0;
    current->R[1][1] = 1;
    
    [features_info insertObject:current
                        atIndex:[features_info count]];
    
    return success;
}
