//
//  ProcessViewController.h
//  drySLAM
//
//  Created by Anthony Rodriguez on 3/8/15.
//  Copyright (c) 2015 ECE473 573 Common. All rights reserved.
//

#import <UIKit/UIKit.h>

@interface ProcessViewController : UIViewController <UITableViewDelegate, UITableViewDataSource>

@property(weak,nonatomic) IBOutlet UIButton* backButton;
@property(weak,nonatomic) IBOutlet UIButton* processButton;
@property(weak,nonatomic) IBOutlet UILabel* feedbackLabel;
@property(weak,nonatomic) IBOutlet UITableView* tableView;

-(IBAction)backButtonPressed:(id)sender;
-(IBAction)processButtonPressed:(id)sender;

@end
