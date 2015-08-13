//
//  VisualGraphicsViewController.h
//  drySLAM
//
//  Created by Anthony Rodriguez on 3/26/15.
//  Copyright (c) 2015 ECE473 573 Common. All rights reserved.
//

#import <GLKit/GLKit.h>

@interface VisualGraphicsViewController : GLKViewController

//Draw the camera
//Quaternion variables (length 7)
-(void)drawCameraAtLocation:(float*)location locationArraySize:(int)location_size;

@end
