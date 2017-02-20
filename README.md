#Physics Research Code

This it my physics research code repo. It's mainly to help me organize and backup my code as well as have branches
that will run and compile on my OSX vs linux systems. The main code can be found in the `current` folder with other
folders for different parts of the code or useful scripts.  The main dependency is [cern root](https://root.cern.ch)
and I usually use the root [github repo](https://github.com/root-mirror/root) to get the latest version.


## Compile and run
All code is in the `current` folder. Inside current there is a cpp analysis program, a cpp skim program and an attempt at a python analysis program.
To make either the cpp or skim program go into the folder and run make.


### E1D analysis

### Skim

### Python

I'm using [cppyy](http://doc.pypy.org/en/latest/cppyy.html) to compile a c++ class which I then map out to as many processors as I want.  Still working on the reduce portion of the function. To run first `make` and then run `./main.py path/to/input/files path/to/output/files`.

There is also a `test.py` which runs all the code in pure python and as a cpp calculation but it seemed to very slow.
