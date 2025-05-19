# After running the above simulations, I realized we don't need separate
# optimization runs for each daylong recording, since at this point the
# optimization is pure random search. So we can merge the three sets of
# simulation runs and then re-do the search for best fits.

recordings = c("0344_000913","0833_010606","0054_000603")
simTypes = c("nonInteractive","a2interactive","bidirectional")

