# PPforest 0.2.1

* Fixed predictions that depended on the other rows predicted with them. Each tree's predicted groups were renumbered by the groups it predicted for those rows, so `predict()` on rows that do not include every group, for example a single group, returned the wrong labels. The same renumbering affected the training, test and out-of-bag predictions whenever a tree predicted no observation into some group, which becomes likely with many groups.
* Fixed the counting of tree votes when the votes skip a group: the votes for the groups after it were dropped and the winning label was wrong.
* `PPforest()` no longer fails with "subscript out of bounds" when no tree predicts one of the groups, and its out-of-bag and test errors compare each prediction with its own group.
* `predict()` works on a single row.

# PPforest 0.2.0

* Changes in version 0.2.0
* This is the fifth release of the package, with minor changes.
* Added a predict() S3 method for PPforest objects.
* Included new data standardization methods: min-max, quantile, and scale. If the PPforest object was fitted using any of these standardizations, the predict() method will apply the same transformation to newdata.
* Argument class in PPforest version 0.1.3 was renamed as y 
* In version 0.1.3, the functions PPclasifly, baggtrees, and tres.pred were exported. In this version, they are now kept as internal functions.
* Fixed several bugs.
