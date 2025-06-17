# Frequent subtrajectory Mining under the Frechet Distance
It builds on these libraries:
- CGAL
- boost


To build the sampler run (provided that the above-mentioned libraries are available):

```
g++ -g ./chernoff_sampler.cpp -o ./chernoff_sampler -I./ext/ -I./src/ -std=c++20
```

Similarly, to build the matcher run:
```
g++ -g ./helloworld.cpp -o ./helloworld -I./ext/ -I./src/ -std=c++20
```