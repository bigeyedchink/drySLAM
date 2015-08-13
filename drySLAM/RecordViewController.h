//
//  RecordViewController.h
//  drySLAM
//
//  Created by Anthony Rodriguez on 3/8/15.
//  Copyright (c) 2015 ECE473 573 Common. All rights reserved.
//

//#import <opencv2/opencv.hpp>
#import <opencv2/highgui/highgui_c.h>
#import <opencv2/highgui/cap_ios.h>
#import <UIKit/UIKit.h>

@interface RecordViewController : UIViewController<CvVideoCameraDelegate, AVCaptureFileOutputRecordingDelegate, AVCaptureVideoDataOutputSampleBufferDelegate>{
    CvVideoCamera* videoCamera;
    BOOL cameraIsCapturing;
    AVCaptureSession* captureSession;
    AVCaptureMovieFileOutput* movieFileOutput;
    AVCaptureDeviceInput* videoInputDevice;
    //TODO: Implement this for real-time
    CvVideoWriter* videoWriter;
}

@property(nonatomic,retain) CvVideoCamera* videoCamera;
@property (retain) AVCaptureVideoPreviewLayer* previewLayer;

@property(weak,nonatomic) IBOutlet UIButton* backButton;
@property(weak,nonatomic) IBOutlet UIButton* startStopButton;
@property(weak,nonatomic) IBOutlet UIImageView* imageView;
@property(weak,nonatomic) IBOutlet UIView* AVView;
@property(weak,nonatomic) IBOutlet UILabel* feedbackLabel;



-(IBAction)backButtonPressed:(id)sender;
-(IBAction)startStopButtonPressed:(id)sender;

-(void) CameraSetOutputProperties;
-(AVCaptureDevice*) CameraWithPosition:(AVCaptureDevicePosition) Position;
-(IBAction)CameraToggleButtonPressed:(id)sender;


@end




