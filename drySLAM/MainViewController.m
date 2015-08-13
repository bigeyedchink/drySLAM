//
//  MainViewController.m
//  drySLAM
//
//  Created by Anthony Rodriguez on 3/8/15.
//  Copyright (c) 2015 ECE473 573 Common. All rights reserved.
//

#import "MainViewController.h"
#include "RecordViewController.h"
#include "ProcessViewController.h"
#include "PlaybackViewController.h"
#include "SLAMModelAdapter.h"

extern NSMutableArray* fileNamesArray;

@interface MainViewController ()

@end

@implementation MainViewController

RecordViewController* rvc;
ProcessViewController* pvc;
PlaybackViewController* pbvc;


- (void)viewDidLoad {
    [super viewDidLoad];
    [SLAMModelAdapter updateFileNamesArray];
    rvc = [[RecordViewController alloc]init];
    pvc = [[ProcessViewController alloc]init];
    pbvc = [[PlaybackViewController alloc]init];
}

-(IBAction)recordViewButtonPressed:(id)sender{
    [self presentViewController:rvc animated:YES completion:NULL];
}

-(IBAction)playbackViewButtonPressed:(id)sender{
    [self presentViewController:pbvc animated:YES completion:NULL];
}

-(IBAction)processViewButtonPressed:(id)sender{
    [self presentViewController:pvc animated:YES completion:NULL];
}

- (void)didReceiveMemoryWarning {
    [super didReceiveMemoryWarning];
    // Dispose of any resources that can be recreated.
}

/*
#pragma mark - Navigation

// In a storyboard-based application, you will often want to do a little preparation before navigation
- (void)prepareForSegue:(UIStoryboardSegue *)segue sender:(id)sender {
    // Get the new view controller using [segue destinationViewController].
    // Pass the selected object to the new view controller.
}
*/

@end
