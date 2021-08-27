# biogeochem_testbed_1.1
Soil biogeochemical testbed 

Code Repository

Created by Will Wieder, Melannie Hartman, & Ben Sulman: 

Updated Aug. 29, 2018

User’s manual and technical documentation for the biogeochemical testbed accepted in Global Change Biology, Oct 2017.
The biogeochemical testbed code base used in these simulations is included (commit 26b4630).

Updates to the code base address issues documented here and include modifications to CORPSE parameterization and the addition of a soil moisture scalar to MIMICS. Updates to the Example_Grid also simulate RCP4.5 and 8.5 through 2100 for each model.

For questions, comments, or inquiries, please contact Will Wieder wwieder@ucar.edu

# Licence

The MIT License (MIT)

Copyright (c) 2017 Will Wieder, Melannie Hartman, & Ben Sulman:

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


## To get started using this repo:
Fork this repository to your personal repo (`fork` button in the upper right)

Then, from the local directory you want to work from you can use the following commands

#### This will clone your remote repository (now called `origin`) to your local directory
git clone https://github.com/`user`/biogeochem_testbed_1.1.git

cd biogeochem_testbed_1.1

#### Then you can link to the main, remote repo & call this `upstream`
git remote add upstream https://github.com/wwieder/biogeochem_testbed_1.1.git

git remote -v

#### Let's create a branch  
git checkout --no-track -b Testbed_CN origin/Testbed_CN

  - or maybe this would be better **not sure it works?**

git checkout -b --track Testbed_CN upstream/Testbed_CN

#### It's a good idea to periodically look for changes to the upstream master (or branch)
git status upstream/master

git pull upstream master

git merge

##### If there are conflicts you can needed you can rebase your repo 
git rebase upstream/master

#### Then you can bring in files you may want to use, create new directoryies, etc
mkdir POINT

mv `file_x` POINT/.

#### At some point you may want to add & commit new filese & directories to your local repo
git status   (are there any files to track or add?)

git add POINT 

git status    (now the directory `point` is staged for commit)

git commit    (this will open up a text file where you can add comments

#### You can get a nice visual on the status of your work
git log --decorate --oneline --graph

#### Create and check out new branch, best to start from the latest upstream/master
##### This example makes a branch from remote (not local) master _or_ branch
git checkout --no-track -b mybranch upstream/master

git fetch --all

#### Create a new branch has same name as remote branch & switch to it
git checkout --track upstream/Testbed_CN   
git checkout Testbed_CN

#### To update your local branch to the remote repo just try
git pull

#### to make a new branch within your own code...
##### This is a good idea if you're going to make code modifications!
git -b Testbed_CN_v2 Testbed_CN




