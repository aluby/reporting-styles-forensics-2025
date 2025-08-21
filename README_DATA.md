## This is the data dictionary for all the data sets used and mentioned in Bayesian Models for Analyzing the Impact of Examiner Variability and Response Category Options in Forensic Fingerprint Comparisons (2025).

Throughout the paper, we utilize data from the Hicklin et al. 2025 fingerprint black box study and from the Ulery et al. 2011 fingerprint black box study. Both data sets have been mutated to include the variables below, however the Ulery 2011 data set does not contain a variable for the 5 level conclusion scale as the original study only used a 3 level scale. Both data sets can be found in the data_changed folder, with the original, unmutated data living in the data_unchanged folder. The unchanged data sets include more variables than our mutated data sets, but information on the additional variables can be found in the appendix of the original papers.

The additional data sets that our model can be fit to include data from black box studies on shoe prints, handwriting, tire marks, another fingerprint data set from a study done by Busey et al, and bullet and cartridge case comparisons. The Busey et al. data set does not have a variable for the item difficulty, as examiners were not asked about this in the study. There are two unmutated bullet/cartridge case data sets, the KQ set representing known-questioned item pairs, where an examiner would compare an unknown bullet or cartridge case with a known item; this could be something like a bullet from a suspect's gun. The QQ set is more similar to what examiners are likely to see more often in their casework, as they may not have a known item to compare to, but may compare the given bullet/cartridge case to one from another crime scene. Both bullet data sets come from the same black box study, and after being mutated have been combined with a variable distinguishing between a KQ or QQ item pair. All 5 data sets contain similar variables to the Hicklin and Ulery fingerprint data sets, which are described below. Again, all supplementary data sets can be found in the data_unchanged folder, and more information on the variables that we did not include can be found in their respective supplementary materials sections, which have been linked below. All of the changed data sets that can be fit to our data are in the data_changed folder.

All information in the data sets is described below—this includes things like how many examiners were present in each study, how many item pairs were in each study, etc.

### Variables:

-   `exID` - examiner id

-   `pairID`- fingerprint image pair

-   `KQ/QQ` - only exists in the bullets data set; differentiates between a KQ or QQ item set, with KQ representing a known-questioned set and QQ representing a questioned-questioned set

-   `responseID` - each observation in the data set goes a unique response ID

-   `ground_truth` - the true mating of the item pair

-   `value_decision` - whether the image pair is of value or not; not in every mutated data set, but causes both response columns and reported difficulty column to be NA as th examiner would not have answered for that given item pair

-   `response3` - examiner's decision on a 3 level conclusion scale, with 1 representing an exclusion, 2 representing an inconclusive, and 3 representing an identification

-   `response5` - examiner's decision on a 5 level scale, with 1 representing exclusion, 2 inconclusive leaning exclusion, 3 inconclusive, 4 inconclusive leaning identification, and 5 identification

-   `response7` - only exists in the footwear and tire mark data sets; examiner's decision on a 7 level scale, with 1 representing exclusion, 2 almost exclusion, 3 inconclusive leaning exclusion, 4 inconclusive, 5 inconclusive leaning identification, 6 almost identification, and 7 identification

-   `reported_difficulty` - the perceived difficulty of each item pair according to the examiners, with 1 representing an easy item pair and 5 representing a difficult item pair

### Data set information:

| Data Set, with links to the original paper | Number of Examiners | Number of Item Pairs | Number of Responses | Additional info |
|---------------|---------------|---------------|---------------|---------------|
| [Hicklin](https://www.sciencedirect.com/science/article/pii/S0379073825000957#bib22) (2025) fingerprints | 156 | 300 | 14,224 | **base data set for results in our paper** |
| [Ulery](https://www.pnas.org/doi/10.1073/pnas.1018707108) (2011) fingerprints | 169 | 744 | 17,121 | **also used in our paper**; only 3 level conclusion scale; no value decision column |
| [Busey](https://www.sciencedirect.com/science/article/pii/S0379073822003735) (2022) fingerprints | 74 | 60 | 4,190 | no value decision column |
| [Hicklin](https://www.sciencedirect.com/science/article/pii/S0379073822002481#sec0085) (2022) footwear | 84 | 269 | 6,610 | has 7 level response scale |
| [Hicklin](https://www.pnas.org/doi/abs/10.1073/pnas.2119944119) (2022) hand writing | 86 | 200 | 7,196 | no value decision column |
| [Hicklin](https://www.sciencedirect.com/science/article/pii/S0379073824003694) (2024) bullets and cartridge cases | 49 | 620 | 1,581 | both KQ and QQ combined |
| [Richetelli](https://www.sciencedirect.com/science/article/pii/S0379073824000902?via%3Dihub) (2024) tire marks | 17 | 77 | 239 | has 7 level response scale |
