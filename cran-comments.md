## Resubmission

This is a resubmission (version 0.1.1) addressing the UBSan report on the
M1-SAN check flagged by the CRAN team (email of 2026-07-16):

    dtweedie.cpp:28:13: runtime error: inf is outside the range of
    representable values of type 'int'
    dtweedie.cpp:36:13: runtime error: inf is outside the range of
    representable values of type 'int'

Cause: when `dtweedie()` was called with `x = +Inf`, the internal saddlepoint
index `j_max` became `+Inf`, and the subsequent conversion of `ceil()`/`floor()`
to `int` was undefined behaviour.

Fix: `x = +Inf` is now handled explicitly (density 0, log-density -Inf), and the
series expansion is only applied to finite positive `x`, so the value passed to
the `double`->`int` conversion is always finite. `NaN` inputs now propagate to
`NaN` as well. No user-facing behaviour changed for finite inputs.

## Test environments

* local macOS, R 4.5.2
* All 190 package tests pass (`[ FAIL 0 | WARN 0 | SKIP 0 | PASS 190 ]`).

## R CMD check results

0 errors | 0 warnings | 0 notes
