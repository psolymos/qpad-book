# Point count data analysis

> How to violate assumptions and get away with it

![](images/bsims.gif)

This book provides material for the workshop
*Analysis of point-count data in the presence of variable survey methodologies and detection error*
at the [AOS 2019 conference](https://amornithmeeting.org/)
by [Peter Solymos](http://peter.solymos.org).

The book and related materials in this repository is the basis of a
full day workshop (8 hours long with 3 breaks).

Prior exposure to [R](https://www.r-project.org/) language is necessary
(i.e. basic R object types and their manipulation, such as arrays, data frames, indexing)
because this is not covered as part of the course.
Check [this](_etc/R-basics.pdf) intro.

## Summary of course objectives

This course is aimed towards ornithologists analyzing field observations,
who are often faced by data heterogeneities due to
field sampling protocols changing from one project to another,
or through time over the lifespan of projects, or trying to combine
'legacy' data sets with new data collected by recording units.
Such heterogeneities can bias analyses when data sets are integrated
inadequately, or can lead to information loss when filtered and standardized to
common standards. Accounting for these issues is important for better
inference regarding status and trend of bird species and communities.

Analysts of such 'messy' data sets need to feel comfortable
with manipulating the data, need a full understanding the mechanics of the
models being used (i.e. critically interpreting the results and acknowledging
assumptions and limitations), and should be able to make informed choices when
faced with methodological challenges.

The course emphasizes critical thinking and active learning.
Participants will be asked to take part in the analysis:
first hand analytics experience from start to finish.
We will use publicly available data sets to demonstrate the data manipulation
and analysis. We will use freely available and open-source R packages.

The expected outcome of the course is a solid foundation for further
professional development via increased confidence in applying these methods
for field observations.


## Course outline

1. Introduction
2. Organizing and processing point count data
3. A primer in regression techniques

Short break

4. Behavioral complexities
5. The detection process

Luch break

6. Dealing with recordings
7. A closer look at assumptions

Short break

8. Understanding roadside surveys
9. Miscellaneous topics

## Todo

- Make sure that packages load (not just installed).
- Save `version=2` Rdata files as well for back compatibility.
- Add 'reset' button instead of `set.seed` slider.
- Increase D slider max to 20?
- Time/distance intervals to fix (provide only one, or change default)
- Repel distance: might be more activity around the observer.


