# RainGageInterpolation
Calculation and Interpolation of Bias correction to apply across daily NexRad data within a model domain.  

Historically NEXRAD statistics in South Florida show rain totals are overestimated for small rainfall events, and underestimated for large rainfall events relative to rain gauge network.  

Bias correction factors are calculated at rain gauges and krigged producing an annual bias rasters which are then multiplied by daily NexRad rasters created from pixel values.
