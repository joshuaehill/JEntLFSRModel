# Code Used Within the Paper *JEnt v2.2.0 LFSR Conditioning Analysis*

The LFSR used within JEnt v2.2.0 is linear, and thus can be modeled using matrix multiplication. This idea is explored in the paper *JEnt v2.2.0 LFSR Conditioning Analysis*, which was presented in the CMUF Entropy Working Group using [this presentation](https://www.untruth.org/~josh/sp80090b/20251209%20JEnt%20v2.2.0%20LFSR%20Conditioning%20Analysis%20Slides.pdf) (additional presentations are likely).

The code included in this repository can be used to generate the matrices $A$, $A_{\text{single}}$ , and $B$. Because this LFSR operation is the sum of two linear transforms, this LFSR operation is completely specified by its action on a basis.

For the purpose of demonstrating this behavior this code also performs stochastic testing to provide evidence that the described implementation of the LFSR processing as linear operations is equivalent to the version implemented within JEnt v2.2.0. Such testing is not formally required, but may help persuade readers who are not as comfortable with linear algebra.

This code also includes a program to model all the raw noise symbols in a range and then produce 1 million bytes of conditioned output per analyzed raw noise symbol by repeatedly LFSR processing the existing state. This is used to support the claim that the LFSR output is indistinguishable from pseudorandom values using the SP 800-90B testing tools.

This repository also includes the `filtermass-testing.tar.xz` file, which contains the JEnt v2.2.0 and v3.6.3 code that was used in filter mass testing.

The `hprime-results.tar.xz` file contains the raw results associated with the distribution (cited in Section 4.2.2 in the associated paper), but this file is too large for GitHub. This file can be found [here](https://untruth.org/~josh/sp80090b/jent-lfsr/hprime-results.tar.xz).
