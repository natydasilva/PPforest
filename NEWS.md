# PPforest 0.2.0

* Changes in version 0.2.0
* This is the fifth release of the package, with minor changes.
* Added a predict() S3 method for PPforest objects.
* Included new data standardization methods: min-max, quantile, and scale. If the PPforest object was fitted using any of these standardizations, the predict() method will apply the same transformation to newdata.
* Argument class in PPforest version 0.1.3 was renamed as y 
* In version 0.1.3, the functions PPclasifly, baggtrees, and tres.pred were exported. In this version, they are now kept as internal functions.
* Fixed several bugs.
