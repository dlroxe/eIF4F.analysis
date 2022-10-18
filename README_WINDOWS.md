Please note the descriptions of Systems 1, 2, and 3 from [README.md](https://github.com/a3609640/eIF4F.analysis#readme).  Behavior
of the package on Windows was evaluated using these systems, at different points
during package development by David Roxe, as a favor to the authors.  The
following account includes some investigation regarding apparent differences
in current behavior on Windows, compared to past behavior.


On Nov 2. 2020, I submitted a descriptive change
[here](https://github.com/a3609640/EIF-analysis/commit/26a1b84082450ca219319ec134c83c8b53437b67)
to indicate that the package was observed to work on System 2.  In particular,
note that the then-current version of R was v4.0.3.

I learned of a report that the current version of this package could not
complete its execution (specifically its correlation analysis) within expected
time frames on Windows 10, even with powerful hardware.  The "expected time
frames" were the times consistently observed on Linux using System 1. My
observations on Windows 11 on System 3 matched the report on Windows 10.  I was
able to see every step complete fairly quickly, except for correlation analysis.

I used the "profiling" feature of RStudio to analyze the "lung"
portion of the correlation analysis (because this is much more tractable than
the "all" portion).  I found that two places accounted for nearly all CPU time:
ComplexHeatmap:Heatmap(), and .correlation_coefficient(), in COR.R.

Since correlation analysis was consistently observed to execute in timely
fashion on Linux, I attempted a controlled comparison.

First, to rule out any difficulty with RStudio on Windows, I evaluated direct
use of R-4.2.1 on Windows.  Although it was marginally faster for most tasks
performed by the package, correlation analysis still did not complete in a
reasonable time.  I concluded that RStudio was not the cause of the performance
difficulty.

An interesting feature of modern Windows is the "Windows Subsystem for Linux"
(WSL).  This feature permits the installation of various Linux distributions
within specialized environments managed by Windows.  Notably, Linux binaries
do not need to be recompiled in order to execute within WSL.  Windows is able
to execute them in the same way they execute natively on Linux.  Thus, WSL
provides an opportunity to rule out Windows as a culprit for the performance
difficulty, by comparing the behavior of Windows-native R and Linux-native R,
both executing under Windows.

It should be noted, however, that WSL is not perfect.  Although it can
decipher and execute the Linux binary format, it does not formally guarantee
equivalence to an actual Linux runtime environment.  Dynamically-linked Linux
libraries may be absent (or at least expected versions may be absent).  Thus,
certain runtime errors may arise during execution of Linux binaries under WSL
that may not arise during execution under Linux.

I installed the Ubuntu 22.04 WSL package, and configured it with several
several Ubuntu packages necessary to support the various dependencies of this R
package.  I then (with some trouble) completed an installation of this R
package in that WSL environment.  Note, WSL forced the use of R-4.1.2 instead of
R-4.2.1; this in turn forced use of BiocManager 3.1.4 instead of
BiocManager-3.1.5.  This in turn affected other package versions.

Under WSL, I was able to execute the package, but not reliably.  I sometimes
encountered runtime errors, particularly "recursive gc invocation".  I
attribute this to complications with WSL library compatibility.

Nevertheless, I was able to observe, repeatedly, execution of the entire package
(including correlation analysis) in approximately 20m.  Note that I made
[some adjustments to
Analysis.R](https://github.com/a3609640/eIF4F.analysis/commit/c6e25e5dd9f2cb4856429fb4186cff357a6251f5)
in order to facilitate this.  As an example of such an adjustment, it seems to
me that the Windows-native R can only fetch biological pathway analysis data
from a remote website using 'curl'; whereas the Linux-native R can only do so
using 'wget'.

To summarize:

  - The package worked on Windows 10 on System 2 in Nov 2020 (with
    Windows-native R 4.0.3), but correlation analysis seems to require
    inordinate execution time (if it completes execution at all) on Windows 10
    and Windows 11 in Oct 2022.
  - Because the performance problem has been seen recently on Windows 10 and 11,
    it does not seem as though the Windows version is the culprit.
  - Because the problem has been seen with and without RStudio, it
    does not seem as though RStudio is the culprit.
  - Profiling suggests that a great deal of CPU time is expended in
    .correlation_coefficient().  Although it's true that the package code has
    substantially changed ince Nov 2020, it should be noted that the body of
    .correlation_coefficient() [has not changed significantly](https://github.com/a3609640/EIF-analysis/blob/8ea7dada7e7d5ba732c6e10cb28215fdb6eb06e6/R/EIFanalysis.R#L5703).
    So, code changes in the package seem unlikely to be the culprit.
  - It is possible that the package works on Windows with versions R v4.0.3 and
    v4.1.2, but cannot work properly with Windows-native R at v4.2.1.  I did not
    evaluate behavior with a Windows-native R v4.1.2 binary.
  - R package versions may also vary between R v4.0.3 an R v.4.1.2.  It may be
    that these R package version differences are more blameworthy than R binary
    version differences.
  - It is possible that a difference between executable binaries for
    Windows-native R and Linux-native R is the culprit for poor performance
    on Windows-native R.

That is, I *do* suspect that a current pecularity of Windows, R, and this package
(and their manner of interaction) might explain the Windows performance
difficulty.  I *do not* suspect a code path within the R package itself.

Regarding the possibility of executing this package on Windows using WSL, as a
way to work-alike with Linux, that unfortunately does not seem to me like a
promising approach.

  - The installation procedure itself has many complexities beyond the existing
    procedure.  WSL installation itself is not yet entirely trivial.  After a
    Ubuntu-22.04 was installed, separate installs were necessary for several
    development packages.  The version of R available with default package
    repositories was too old; so system package repositories required
    adjustment.  R package installation required several retries, seemingly due
    to network interface difficulties.  Unprivileged installation of R packages
    did not work "out of the box".
  - As noted above, the WSL R installation is not stable; and despite attempts
    I am unable to provide details (much less a fix or workaround) for that at
    this time.
  - Even if the installation could be documented, and package execution could be
    made reliable, and execution time were therefore repeatably quick, resource
    consumption would still present a difficulty.  Specifically WSL is subject
    to a RAM limitation: It has access to at most 50% of system RAM.  Some web
    search results suggest "50% or 8GB, whichever is smaller", but at any rate I
    note that on System 3, which has 128G of RAM, WSL reports that it has 64G of
    RAM available.  Because 48G are recommended for exeuction of the package,
    under WSL the recommended hardware would include 96G of RAM (which would
    preclude the use of most PC laptops on the market today, because consumer
    PC laptops provide support for 64G or less.

My own instinct is that the correlation analysis portion of the package may not
be easily supportable on Linux at the present time.  I hasten to point out that
I have often have seen the package code other than correlation analysis complete
quickly.
