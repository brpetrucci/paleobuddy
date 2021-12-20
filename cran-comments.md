## Test environments
* local macOS Big Sur 11.6, R 4.1.1
* devtools::check_win (devel, release, oldrelease)
* rhub::check (all available platforms)
* GitHub actions (macOS latest, windows latest release; ubuntu latest release, devel and oldrel-1)

## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs

This is the first submission of this package, so there may be a first submission NOTE.

## Downstream dependencies
There are currently no downstream dependencies for this package

# Resubmission
This is a resubmission. In this version I have resolved comments by a CRAN member---great thanks on the feedback. The following have been changed:

* "Please always write package names, software names and API (application
programming interface) names in single quotes in title and description.
e.g: --> 'paleotree', 'APE', 'RPANDA', 'paleobuddy'"
 ** Changed all cases of package names in description to be wrapped in single quotes. No instances in title.
 
* "Please write references in the description of the DESCRIPTION file in
the form
authors (year) <doi:...>
authors (year) <arXiv:...>
authors (year, ISBN:...)
or if those are not available: authors (year) <https:...>
with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for
auto-linking.
(If you want to add a title as well please put it in quotes: "Title")"
 ** Changed reference of RPANDA to the first mentioned format. No other instances of references in DESCRIPTION

* "Please make sure that you do not change the user's options, par or
working directory. If you really have to do so within functions, please
ensure with an *immediate* call of on.exit() that the settings are reset
when the function is exited. e.g.:
...
oldpar <- par(no.readonly = TRUE)    # code line i
on.exit(par(oldpar))            # code line i + 1
...
par(mfrow=c(2,2))            # somewhere after
...
e.g.: R/draw.sim.R, R/sim.R
If you're not familiar with the function, please check ?on.exit. This
function makes it possible to restore options before exiting a function
even if the function breaks. Therefore it needs to be called immediately
after the option change within a function."
 ** Included first two lines mentioned before any change in par within functions, i.e. draw.sim.R (par(yaxt = "n")); sim.R (par(mfrow = c(1, 3))). No other instances of par or options changes within functions.


* "In examples and vignettes and demos: Please always make sure to reset to
user's options(), working directory or par() after you changed it:
e.g.:
oldpar <- par(mfrow = c(1,2))
...
par(oldpar)"
 ** Included the first line mentioned before any changes to par in examples, and the second before the end of any plotting, in examples for make.phylo.R and find.lineages.R


* "\dontrun{} should only be used if the example really cannot be executed
(e.g. because of missing additional software, missing API keys, ...) by
the user. That's why wrapping examples in \dontrun{} adds the comment
('# Not run:') as a warning for the user.
Does not seem necessary.

Please unwrap the examples if they are executable in < 5 sec, or replace
\dontrun{} with \donttest{}."
 ** Excluded \dontrun{} from any examples that took less than 5 sec in my local machine. Changed \dontrun{} to \donttest{} in other cases. Kept \dontrun{} in case of examples that will error on run, just used for illustration of wrong usage.
