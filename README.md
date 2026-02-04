# Frequent subtrajectory Mining under the Frechet Distance
It builds on these libraries:
- CGAL
- boost

To build the sampler run (provided that the above-mentioned libraries are available):

```
g++ -g ./chernoff_sampler.cpp -o ./chernoff_sampler -I./ext/ -I./src_untouched/ -std=c++20 -O3
```

Similarly, to build the matcher run:
```
g++ -g ./helloworld_pruning.cpp -o ./helloworld_pruning -I./ext/ -I./src_untouched/ -std=c++20 -O3
```
or to build the simplified version:
```
g++ -g ./helloworld_simplified.cpp -o ./helloworld_simplified -I./ext/ -I./src_untouched/ -std=c++20 -O3
```
To build subtrajectory clustering algo do the same:
```
g++ -g ./subtrajectory_clustering.cpp -o ./subtrajectory_clustering -I./ext/ -I./src_untouched/ -std=c++20 -O3
```

## Data format
The dataset with all trajectories must be saved in .txt files, where a row is ```x_coord y_coord trajectory_id```. 
```trajectory_id``` must be integer and trajectories' ids must be increasing and consecutives. 
Similarly, a file with trajectories upon which we build pathlets must be stored in the same format.

## Sampling 
Run 
```
./chernoff_sampler [OPTIONS] input_dataset output_dir
```
With these options:
-  ```-h,--help```                   Print this help message and exit
-  ```-e,--epsilon FLOAT```          The maximum allowed additive error on the frequency.
-  ```-d,--delta FLOAT```            Confidence parameter.
-  ```-r,--radius FLOAT```           The radius for the vc dimension.
-  ```-c,--cardinality INT```        Fixed Size for fixed size sample generation.
-  ```-s,--seed INT```               Seed for the random generator. Default is 0.
-  ```-m,--mode ENUM``` REQUIRED     Sampling mode: 0 for Fixed Size, 1 for Chernoff, 2 for VC, 3 for Rough VC Estimate.

The output sample file will be formatted according to the selected mode.
## Mining 
An example command is reported below
```
./helloworld_pruning sample_file pathlet_file output_file -r radius -f frequency_threshold
```
The script writes the identified pathlets to output_file. By default it will print all the frequent pathlets, to have the maximal ones only, post-process the output or change the corresponding line in the script. 

## Mining with simplification
An example command is reported below
```
./helloworld_simplified sample_file pathlet_file output_file -r radius -f frequency_threshold -s simplification_factor
```
Here ```simplification_factor``` is a floating point number between 0.0 and 1.0. 
The script writes the identified pathlets to output_file. By default it will print all the frequent pathlets, to have the maximal ones only, post-process the output or change the corresponding line in the script. 