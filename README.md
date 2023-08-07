[![Visual demonstration](http://img.youtube.com/vi/xysaDVy2M3g/0.jpg)](http://www.youtube.com/watch?v=xysaDVy2M3g "Visual demonstration (YouTube)")

A "simple" NEE+MIS Kajiya-style pathtracer, where pseudorandom number generator, which is fed to the estimator, is temporally a low discrepancy sequence with blue noise spatial correlation. Scene is constructed with sphere primitives and procedural skybox, on which the sun disk is importance sampled. The code is by design a single file for simplicity. See commit history for an accidental semi step-by-step tutorial.
