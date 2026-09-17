# ALON

Alon means "wave" in Tagalog, one of the over 170 languages spoken by Filipinos. Their ancestors were seafarers. But they did not tame the waves, for nature
cannot be subdued. They learned the language of the waves, and the waves learned theirs.

#

ALON (Abundance-based Latitudinal distributiON) is an R-based framework designed to assign microbial biogeography using global scale -omics data. 
We applied this framework to microbial hosts and viruses using metagenomes from the ocean. 

The pipeline consists of two stages: Per-latitude abundance prediction using generalized additive models (GAMs) and classification of latitudinal distributions.
The initial input requires a sample identifier, feature ID, latitude, and abundance. The final output is a tab-separated file with the feature ID, biogeographic 
classification, and the exclusivity status. The exclusivity status indicates whether the feature ID is found preferentially or exclusively in the assigned region.
