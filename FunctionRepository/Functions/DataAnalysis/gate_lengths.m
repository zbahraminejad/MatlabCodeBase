function shorttraces=gate_lengths(tracestats,motherstats,minlengthtrace,minlengthmother)
shorttraces=tracestats(:,3)<=minlengthtrace | motherstats(:,3)<=minlengthmother;