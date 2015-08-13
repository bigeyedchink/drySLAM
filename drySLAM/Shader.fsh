//
//  Shader.fsh
//  glTest
//
//  Created by Anthony Rodriguez on 3/26/15.
//  Copyright (c) 2015 ECE473 573 Common. All rights reserved.
//

varying lowp vec4 colorVarying;

void main()
{
    gl_FragColor = colorVarying;
}
