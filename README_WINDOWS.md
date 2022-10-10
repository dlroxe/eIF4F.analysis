Please note the descriptions of System 2 and System 3 from [README.md](https://github.com/a3609640/eIF4F.analysis#readme).  Behavior
of the package on Windows was evaluated using these systems, at different points
during package development.  The following is provided by David Roxe, who
undertook some evaluation as a favor to the authors.

On Nov 2. 2020, I submitted a descriptive change
[here](https://github.com/a3609640/EIF-analysis/commit/26a1b84082450ca219319ec134c83c8b53437b67)
to indicate that the package was observed to work on System 2.

Recently, I evaluated the package on Windows 11 on System 3.  I was made aware
of another evaluation on Windows 10, in which it was observed that correlation
analysis could not complete within a reasonable amount of time, even on a very
amply-resourced machine.  This observation was at odds with observed behavior
on Linux on System 1.  However, my observations on System 3 matched the report
of performance difficulties on Windows 10, so I investigated a little further.

I specifically observed a performance problem on Windows in correlation
analysis.  I used the "profiling" feature of RStudio to analyze the "lung"
portion of the correlation analysis (because this is much more tractable than
the "all" portion).  I found two places where substantial computing resources
were spent (to the point of dominating all other computational expenditures):

   - I found that roughly 50% of CPU time was spent in ComplexHeatmap:Heatmap().
   - I found that roughly 50% of CPU time was spent in
     .correlation_coefficient(), in COR.R.

For whatever reason, these functions are slow on Windows, yet perform adequately
on Linux.  I do not have an explanation for this.

I evaluated direct use of R-4.2.1 on Windows (without RStudio).  Although it was
marginally faster, the performance difficulties with correlation analysis
remained quite evident.

I installed Ubuntu 22.04 using the "Windows Subsystem for Linux" (WSL), plus
several Ubuntu packages necessary to support the various dependencies of this R
package.  I then (with some trouble) completed an installation of this R
package in that WSL environment.  Note, WSL forced the use of R-4.1.2 instead of
R-4.2.1; this in turn forced use of BiocManager 3.1.4 instead of
BiocManager-3.1.5.  This in turn affected other package versions.

Using only the bash command line provided by WSL, I was able to execute the
entire package in approximately 20m.  This required [some adjustments to
Analysis.R](https://github.com/a3609640/eIF4F.analysis/commit/c6e25e5dd9f2cb4856429fb4186cff357a6251f5).

Note that WSL executes a binary program exactly as it was compiled for Linux,
rather than a binary compiled specifically for Windows.

So:

  - The package worked on Windows 10 on System 2 in Nov 2020.
  - The package code and the versions of Windows, RStudio, R, and R packages
    have all changed since then.  It should be noted, however, that the
    body of .correlation_coefficient() [has not changed significantly](https://github.com/a3609640/EIF-analysis/blob/8ea7dada7e7d5ba732c6e10cb28215fdb6eb06e6/R/EIFanalysis.R#L5703).
  - Because the performance problem has been seen recently on Windows 10 and 11,
    it does not seem as though the Windows version is the culprit.
  - Likewise, because the problem has been seen with and without RStudio, it
    does not seem as though RStudio is the culprit.
  - Possible remaining culprits would seem to be the following:
    - The Windows-native binary may have a problem that the Linux-native binary
      does not have, even when both are executed on Windows.
    - The version differences (either in R itself or in R packages) between
      R-4.1.2 and R-4.2.1 may contain a performance regression.

I would like to be able to report a stable and reproducible procedure for
executing the package with WSL on Windows.  Among other things, this would
permit further investigation of the performance problem.  However, I have
encountered the following difficulties:

  - The installation procedure itself has many complexities beyond the existing
    procedure.  WSL installation itself is not yet entirely trivial.  After a
    Ubuntu-22.04 was installed, separate installs were necessary for several
    development packages.  The version of R available with default package
    repositories was too old; so system package repositories required
    adjustment.  R package installation required several retries, seemingly due
    to network interface difficulties.  Unprivileged installation of R packages
    did not work "out of the box".
  - After all that, the WSL R installation is not stable.  I have observed
    seemingly-transient errors with varying complaints.  I have also observed
    recurring complaints of "recursive gc invocation", which occur with low
    but annoying probability.  I have experienced multiple Blue Screens of
    Death for kernel security complaints (I think related to WSL's use of
    virtualization and/or the combination of applications I had open and/or
    my BIOS settings.)

Despite all that, I did observe a successful (and reasonably quick) package
execution on Windows.  This is all to say:

  - The expected evaluation times for the work seem accurate in principle.
  - The problems that cause long evaluation times seem specific to either
    a particular R version, or to R as it is compiled for Windows.
  - It is difficult to see how to recommend execution on Windows (at least for
    correlation analysis) for the time being, because the user seemingly has
    the choice between intolerably long execution time, or the unstable WSL
    environment.
