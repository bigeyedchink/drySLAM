//
//  SLAMModelAdapter.h
//  drySLAM
//
//  Created by Anthony Rodriguez on 3/13/15.
//  Copyright (c) 2015 ECE473 573 Common. All rights reserved.
//

#ifndef drySLAM_SLAMModelAdapter_h
#define drySLAM_SLAMModelAdapter_h

NSMutableArray* fileNamesArray;
NSString* filePath;

@interface SLAMModelAdapter : NSObject

//Update array with the names of all files in directory
//TODO: Implement (right now we are using dummy values)
+(void)updateFileNamesArray;

//Call model from ViewControllers to process data
//TODO: Add arguments and return values
+(void)processDataWithRANSAC;


//Initialize filePath for program
+(void)setFilePath;



@end




#endif
