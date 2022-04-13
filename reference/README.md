A light translation of the original odin code provided by Joe, compared manually against the original code with a fixed seed to confirm that we do get the correct answer.

Notable departures from the original:

* moved `beta` inside the variables, as with the mode version
* no longer run the model from the beginning but restart from the last state
* no longer initialise a new odin object at each time step but reuse a single one

This combination makes it much faster, taking 1.7s vs 27s originally, but has negligible effect on the solution.

Running as-is gives a `log_likelihood` of -1714.20949841969
