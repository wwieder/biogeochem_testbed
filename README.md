# biogeochem_testbed_1.1
Soil biogeochemical testbed 

Code Repository

Created by Will Wieder, Melannie Hartman, Ben Sulman, Emily Kyker-Snowman, Brooke Eastman: 

Updated 
Aug. 29, 2018
Feb 10, 2023

User’s manual and technical documentation for the biogeochemical testbed accepted in Global Change Biology, Oct 2017.
The biogeochemical testbed code base used in these simulations is included (commit 26b4630).

Updates to the code base address issues documented here and include modifications to CORPSE parameterization and the addition of a soil moisture scalar to MIMICS. Updates to the Example_Grid also simulate RCP4.5 and 8.5 through 2100 for each model.

Feb 2023 updates include: 
- Representation of coupled C-N biogeochemistry for MIMICS, 
- Representation of root exudation, currently set to zero for all simulation
- Switch to using input data from CLM5-SP with GSWP3 forcing
- This code base was used in single point simulations at the Fernow Experimental Forest by Eastman et al (2023) and global simulations by Wieder et al. (2023).

# Licence

The MIT License (MIT)

Copyright (c) 2017 Will Wieder, Melannie Hartman, & Ben Sulman:

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
