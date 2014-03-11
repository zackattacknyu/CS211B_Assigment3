
About the "toytracer":

This directory contains the code for a very simple ray tracer called the "toytracer".  It is intended as a teaching tool, not as a serious piece of software.  Consequently, it was designed to be extremely small and simple, at the expense of efficiency.  The associated classes for vectors (Vec3) and colors (Color) are minimal, supporting nothing more than the basic operators, all of which are implemented as inline functions in a header file.

The toytracer can read in simple scene descriptions from a file and write out images in PPM format, which is a trivial image file format.  A PPM image consists of a simple header followed by a sequence of rgb triples, one byte per color component.

The toytracer provides a starting point for implementing various Monte Carlo techniques, such as stochastic ray tracing, bidirectional path tracing, and Metropolis path tracing.  It is not intended for ray tracing scenes consisting of more than a few primitive objects until it is augmented with a reasonable acceleration technique, such as spatial subdivision, or bounding volume hierarchies.

James Arvo
University of California, Irvine
April, 2008



