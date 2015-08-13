//
//  MainViewController.h (Repository version - TESTING)
//  drySLAM
//
//  Created by Anthony Rodriguez on 3/8/15.
//  Copyright (c) 2015 ECE473 573 Common. All rights reserved.
//

#import <UIKit/UIKit.h>

@interface MainViewController : UIViewController

@property(weak,nonatomic) IBOutlet UIButton* recordViewButton;
@property(weak,nonatomic) IBOutlet UIButton* playbackViewButton;
@property(weak,nonatomic) IBOutlet UIButton* processViewButton;

-(IBAction)recordViewButtonPressed:(id)sender;
-(IBAction)playbackViewButtonPressed:(id)sender;
-(IBAction)processViewButtonPressed:(id)sender;

@end
