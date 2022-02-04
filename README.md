DISCLAIMER: This is an archive of a group of diffusion MRI processing programs called mincdiffusion, posted here for historical reference. 

-This is not part of the minc project and does not meet its standards.

-This is not a mature project and not meant for external use. The current version is largely at the prototype/WIP stage. It is now not being maintained.

-The deconvolution and probabilistic deconvolution do not implement the now state-of-the-art constrained spherical deconvolution.

-The FACT integration used throughout all of the tractography is now not commonly used, and there is no interpolation in this implementation.

-The fiber ODF and maximum extraction is discretized at the acquisition angular resolution.

-The single fiber response function is hardcoded for an average healthy brain acquired with an outdated acquisition.

-This software works only for transverse images.
