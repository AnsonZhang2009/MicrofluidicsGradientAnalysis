# MicrofluidicsGradientAnalysis
Hi👋 This is a R workflow designed to analyze the gradient of microfluidic devices through two methods: 
1. The Sobel-Kernel edge detection system to calculate gradient strength. 
2. Using the `quantize` function of the `magick` package to obtain the number of unique colors as an index of gradient strength.
## The Sobel-Kernel Edge Detection System
Current commonly used gradient detection systems for image processing include the Kobel-Kernel system. More information can be found [here](https://doi.org/10.3390/jimaging4060074).
## The `Quantize` Detection System
The quantize detection system was something the author of this Github came up originally, with help from the `magick` package. For more information, view the R script, which has extensive comments on every part of the process. 

