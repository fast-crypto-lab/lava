# lava
A Buchberger's algorithm

This programs finds groebner basis, and will produce a log very similar to magma's "F4" algorithm.
To use this program, simply make, then ./f4_magma-test <InputFileName>. Input format is demostrated in sample files, a particular note is that
"," after each input polynomial is strictly required. Also this code so far only works on F_2, so it is required no duplicated variable
exist in the same monomial. Any incorrect input may result in unexpected behavior.
The log is always printed, if you don't want to see it, just go in the code and comment some cout. Adding a enable_log switch is not on 
our priority list, sorry.
To make, m4ri library is required, please moldify the makefile if you installed it other than the default path.

Mutant strategy branch:
This branch uses a mutant-first strategy, and results in a speed faster than magma.
