#!/usr/local/bin/perl

my $KTUP = 12;

require "mkpl.pl";

my @fnames;
my @flength;
my @ftypes;

sub blast_init {
  my ($indexfile) = @_;
  open F0, $indexfile;
  while (<F0>) {
    chomp;
    if ($_ =~ /^(\d+) (\d+) (\S+)/) {
      my $idx = $1;
      $fnames[$idx] = $3;
      $flength[$idx] = $2;
      my $type = "Feature";
      if ($fnames[$idx] =~ /^(\w):([\w\W]+),antisense$/) {
        if ($1 eq "G") { $type = "Gene"; }
        elsif ($1 eq "P") { $type = "Promoter"; }
        elsif ($1 eq "O") { $type = "Origin"; }
        elsif ($1 eq "R") { $type = "Regulatory"; }
        elsif ($1 eq "T") { $type = "Terminator"; }
        elsif ($1 eq "E") { $type = "Cutter"; }
        elsif ($1 eq "f") { $type = "ExactFeature"; }
      }
      elsif ($fnames[$idx] =~ /^(\w):([\w\W]+)$/) {
        if ($1 eq "G") { $type = "Gene"; }
        elsif ($1 eq "P") { $type = "Promoter"; }
        elsif ($1 eq "O") { $type = "Origin"; }
        elsif ($1 eq "R") { $type = "Regulatory"; }
        elsif ($1 eq "T") { $type = "Terminator"; }
        elsif ($1 eq "E") { $type = "Cutter"; }
        elsif ($1 eq "f") { $type = "ExactFeature"; }
      }
      $ftypes[$idx] = $type;
    }
  }
  close F0;
}

# a fragment train is a set of consecutive fragments with perfect
# percent identity, more than 20% of the total size of the gene

sub is_fragment_train {
  my ($f,$h,$m,$i,$d) = @_;
  my $l = $flength[$f];
  if ($i == 0 && $d == 0 && $h >= 0.2*$l) {
    return 1;
  }
  else { return 0; }
}

# scoring function:
#   +1   for each exact match
#   +0   for missing nucleotides at the start or end of sequence
#   +0.3 for each nucleotide in a mutated fragment
#        (i.e. assume roughly 70% of the fragment are mutated)
#   -0.1 for each deleted nucleotide
#   +0   for each inserted nucleotide in a gene and if total inserts
#        is less than a threshold that's related to gene size
#        (i.e. we don't penalize small inserts for genes)
#   -0.1 for each inserted nucleotide not in a gene
#
# $i: inserts
# $d: deletes
# $m: mutations (actually, # nucleotides from missing fragments)

$E_THRESH  = 0.25;
my $WT_THRESH = 0.05;

sub percent_identity {
  my ($seqsize,$f,$h,$m,$i,$d) = @_;
  my $l = $flength[$f];
  if (0 && ($f == 425 || $f == 426)) { print "  pi $l,$h,$m,$i,$d\n"; }
  my $s = $h;
  if ($ftypes[$f] eq "ExactFeature" || $ftypes[$f] eq "Cutter") {
    if ($s == $l && $i == 0 && $d == 0 && $m == 0) { return 0; }
    else { return 1; }
  }
  # missing nucleotides
  $s += ($l-($s+$m))*0;
  # mutations
  $s += $m*0.3;
  # insertions
  $s += $i*0;
  if ($ftypes[$f] ne "Gene") {
    $s += $i*(-0.1);
  }
  elsif ($i > $l*3/4) {
    $s += $i*(-0.1);
  }
  # deletions
  $s += $d*(-0.1);
  my $e = 1-$s/$l;
  if (0 && ($f == 425 || $f == 426)) {
    print "  $f: $l,$h,$m,$i,$d --> $s, $l = $e\n";
  }
  return $e;
}

my @features_to_insert;

sub insertFeatures {
  my ($plasmidid, $seqsize) = @_;

  my @results;
  @features_to_insert =
    sort {($a->{"startp"}+(-1)*(1-$a->{"clockwise"})*$a->{"length"}) <=>
          ($b->{"startp"}+(-1)*(1-$b->{"clockwise"})*$b->{"length"})
	 } @features_to_insert;

  # go through all the genes and features, for engulfed genes or
  # features with similar names, only keep those with the best scores
  my @x;
  for (my $i=0; $i<=$#features_to_insert; $i++) {
    my $r0 = $features_to_insert[$i];
    my %f0 = %$r0;
    my $s0 = $f0{"startp"}+(-1)*(1-$f0{"clockwise"})*$f0{"length"};
    if ($s0 > $seqsize) { next; }
    if ($f0{"type"} ne "Cutter") {
      my $n0 = $f0{"name"};
      my $e0 = $f0{"score"};
      my $found = 0;
      for (my $j=0; $j<=$#features_to_insert; $j++) {
        my $r1 = $features_to_insert[$j];
        my %f1 = %$r1;
        my $n1 = $f1{"name"};
        my $s1 = $f1{"startp"}+(-1)*(1-$f1{"clockwise"})*$f1{"length"};
	if ($s1 > $s0+$f0{"length"}) { last; }
	# check for interleaving gene
	if ($f1{"type"} eq "Gene" &&
	    (($s1 >= $s0 && $s1 <= $s0+$f0{"length"}) ||
	     ($s1 <= $s0 && $s1+$f1{"length"} >= $s0) ||
	     ($s1 >= $s0 && $s1+$f1{"length"} <= $s0+$f0{"length"}) ||
	     ($s1 <= $s0 && $s1+$f1{"length"} >= $s0+$f0{"length"}))) {
	  my $e1 = $f1{"score"};
	  if ($e0 > $e1) {
	    $found = 1;
	    last;
	  }
	}
	# check for similar feature overlap: first two cases are for
	# interleaved features, the last two are for features
	# engulfing one another
	if ($f1{"type"} eq $f0{"type"} &&
	    ($f1{"name"} =~ /$f0{"name"}/ ||
	     $f0{"name"} =~ /$f1{"name"}/) &&
	    (($s1 >= $s0 && $s1 <= $s0+$f0{"length"}) ||
	     ($s1 <= $s0 && $s1+$f1{"length"} >= $s0) ||
	     ($s1 >= $s0 && $s1+$f1{"length"} <= $s0+$f0{"length"}) ||
	     ($s1 <= $s0 && $s1+$f1{"length"} >= $s0+$f0{"length"}))) {
	  my $e1 = $f1{"score"};
	  if ($e0 > $e1) {
	    $found = 1;
	    last;
	  }
	}
      }
      if ($found) { next; }
    }
    push @x, \%f0;
  }

  @features_to_insert = @x;
  for (my $i=0; $i<=$#features_to_insert; $i++) {
    my $r = $features_to_insert[$i];
    my %f = %$r;
    if ($f{"plasmidid"} == $plasmidid) {
      if ($INSERTFEATURE) {
        if ($VECTORMODE) {
          insertVectorFeatures
            ($f{"plasmidid"},
             $f{"type"},
             $f{"name"},
             $f{"startp"},
             $f{"length"},
             $f{"clockwise"});
	}
	else {
          insertPlasmidFeatures
            ($f{"plasmidid"},
             $f{"type"},
             $f{"name"},
             $f{"startp"},
             $f{"length"},
             $f{"clockwise"});
	}
      }
      else {
        my $ft = $f{"type"};
        my $fn = $f{"name"};
        my $fs = $f{"startp"};
        my $fl = $f{"length"};
        my $fc = $f{"clockwise"};
        my $fx = $f{"score"};
        push @results, "$ft,$fn,$fs,$fl,$fc,$fx";
        if ($VERBOSE) {
          print "$ft,$fn,$fs,$fl,$fc,$fx\n";
	}
      }
    }
  }
  undef @features_to_insert;
  return @results;
}

sub check_feature {
  my ($plasmidid, $seqsize, $fidx,
      $hits, $mutations, $inserts, $deletes,
      $startf, $startp, $stopp) = @_;
  my $escore = percent_identity
    ($seqsize, $fidx, $hits, $mutations, $inserts, $deletes);
  my $ftrain = is_fragment_train
    ($fidx, $hits, $mutations, $inserts, $deletes);
  if ($escore < $E_THRESH || $ftrain) {
    my $name;
    my $type = $ftypes[$fidx];
    my $length = $stopp-$startp+1;
    my $clockwise = 1;
    if ($fnames[$fidx] =~ /^(\w):([\w\W]+),antisense$/) {
      $name = $2;
      $clockwise = 0;
    }
    elsif ($fnames[$fidx] =~ /^(\w):([\w\W]+)$/) {
      $name = $2;
    }
    else {
      $name = $fnames[$fidx];
    }
    if ($type eq "ExactFeature" || $type eq "Cutter") {
      # correct size
      $stopp = $startp+$flength[$fidx]-1;
    }
    if ($type eq "ExactFeature") {
      $type = "Feature";
    }
    if ($type eq "Cutter" && $name =~ /^([\w\W]+),(\d+)\/(\d+)/) {
      if ($clockwise) { $startp += ($2-1); }
      else { $stopp -= ($2-1); }
      $length = 1;
      $name = $1;
    }
    if ($type eq "Gene") {
      if ($escore >= $E_THRESH && $ftrain) {
        my $s0 = $startf*$KTUP;
        my $s1 = $s0 + $length;
        if ($clockwise) {
          $name = $name."(".$s0."-".$s1.")";
        }
        else {
          $name = $name."(".$s1."-".$s0.")";
        }
      }
      elsif ($inserts > 2*$KTUP) {
        $name = $name."(w/ gaps)";
      }
      elsif ($escore >= $WT_THRESH || $deletes > $KTUP) {
        $name = $name."(variant)";
      }
    }
    if ($escore < $E_THRESH || ($ftrain && $type eq "Gene")) {
      my %feat;
      if ($clockwise) {
        $feat{"plasmidid"} = $plasmidid;
        $feat{"type"} = $type;
        $feat{"name"} = $name;
        $feat{"startp"} = $startp;
        $feat{"length"} = $length;
        $feat{"clockwise"} = $clockwise;
        $feat{"score"} = $escore;
      }
      else {
        $feat{"plasmidid"} = $plasmidid;
        $feat{"type"} = $type;
        $feat{"name"} = $name;
        $feat{"startp"} = $stopp;
        $feat{"length"} = $length;
        $feat{"clockwise"} = $clockwise;
        $feat{"score"} = $escore;
      }
      push @features_to_insert, \%feat;
      return 1;
    }
  }
  return 0;
}

# take all the fragments of a feature, build potential trains, then
# report fragment trains with high percent identity

sub process_frags {
  my ($plasmidid, $seqsize, @results) = @_;
  
  my $feature;
  my @hits;
  my @inserts;
  my @deletes;
  my @mutations;
  my @trains;
  my @short;

  my $FRAGMENT = 0;
  my $POSITION = 1;

  for (my $i=0; $i<=$#results; $i++) {
    if ($results[$i] =~ /^(\d+) (\d+) (\d+) (\d+)/) {
      if (0 && ($1 == 425 || $1 == 426)) {
        print "$results[$i] ($fnames[$1]) $3\n";
      }
      # 1: feature index
      # 2: fragment index
      # 3: position
      # 4: shift
      if (!defined $feature) { $feature = $1; }
      elsif ($feature != $1) { die "mismatching features\n"; }
      my $fragment = $2;
      my $position = $3;
      if ($4 > 0) { $position += $4; }

      my $extended = 0;
      my $createnew = 1;
            
      my $fl = $flength[$feature];
      my $hitscore = $fl-$fragment*$KTUP;
      if ($hitscore > $KTUP) { $hitscore = $KTUP; }

      # first pass, try extend trains with consecutive fragments
      for (my $j=0; $j<=$#trains; $j++) {
        my $r = $trains[$j];
	my @train = @$r;
	my $lastr = $train[$#train];
	if ($lastr->[$FRAGMENT] < $fragment &&
	    $lastr->[$POSITION]+$KTUP <= $position) {
          # can extend this train
	  my $fdiff = $fragment-$lastr->[$FRAGMENT];
	  my $pdiff = $position-($lastr->[$POSITION]+$KTUP);

	  # extend train with consecutive fragment
	  if ($lastr->[$FRAGMENT] == $fragment-1 && $pdiff == 0) {
            if (0 && ($feature == 425 || $feature == 426)) {
	      print "ext consecutive ($j) $feature $fragment $position\n";
	    }
            my @v;
	    $v[$FRAGMENT] = $fragment;
	    $v[$POSITION] = $position;
	    push @train, \@v;
	    $trains[$j] = \@train;

	    $hits[$j] += $hitscore;
            $extended = 1;
	    $createnew = 0;
	  }
	}
      }

      # second pass, try extend a train after just mutations or with
      # insertions
      if ($extended == 0) {
        for (my $j=0; $j<=$#trains; $j++) {
          my $r = $trains[$j];
	  my @train = @$r;
	  my $lastr = $train[$#train];
	  if ($short[$j] == 0 &&
	      $lastr->[$FRAGMENT] < $fragment &&
	      $lastr->[$POSITION]+$KTUP <= $position) {
            # can extend this train
	    my $fdiff = $fragment-$lastr->[$FRAGMENT];
	    my $pdiff = $position-($lastr->[$POSITION]+$KTUP);

            # no inserts or deletes, just mutations
	    if ($pdiff == $KTUP*($fdiff-1)) {
              if (0 && ($feature == 425 || $feature == 426)) {
	        print "ext inframe ($j) $feature $fragment $position\n";
              }
              my @v;
	      $v[$FRAGMENT] = $fragment;
	      $v[$POSITION] = $position;
	      push @train, \@v;
	      $trains[$j] = \@train;
	      
	      $hits[$j] += $hitscore;
	      $mutations[$j] += ($pdiff);
              $extended = 1;
	      $createnew = 0;
	    }

            # there are inserts, and insert size is smaller than X
	    # percent of the feature length (i.e. if insert size is
	    # too big, just consider the feature to be broken up)
	    elsif ($KTUP*($fdiff-1)<= $pdiff &&
	           $pdiff < int($flength[$feature]*0.75)) {

              # there are inserts, duplicate the earlier portion of
	      # the sequence (before insertion) to two copies and
	      # extend one copy
            
	      my $escore = percent_identity
	        ($seqsize, $feature, $hits[$j], $mutations[$j],
	         $inserts[$j], $deletes[$j]);
              if ($escore < $E_THRESH) {
	        my @train2 = @train;
	        push @trains, \@train2;
	        $hits[$#trains] = $hits[$j];
	        $mutations[$#trains] = $mutations[$j];
	        $inserts[$#trains] = $inserts[$j];
	        $deletes[$#trains] = $deletes[$j];
	        $short[$#trains] = 1;
                if (0 && ($feature == 425 || $feature == 426)) {
	          print "dupinsert ($j->$#trains) $escore\n";
		}
              }

              if (0 && ($feature == 425 || $feature == 426)) {
	        print "extend inserts ($j) $feature $fragment $position\n";
	      }
              my @v;
	      $v[$FRAGMENT] = $fragment;
	      $v[$POSITION] = $position;
	      push @train, \@v;
	      $trains[$j] = \@train;

	      $hits[$j] += $hitscore;

              # an insert should really only cause one fragment to be
	      # missing, but there may be some mutations around the
	      # splicing sites so we conservatively assume two
	      # fragments may be missing due to the insert. all other
	      # missing fragments must have mutations.

	      if ($fdiff > 3) {
	        # ($fdiff-1)-2 is the number of fragments that may
		# have mutations, assuming we allow two fragments to
		# be mutated due to the insert
	        $mutations[$j] += (($fdiff-1)-2)*$KTUP;
	        $inserts[$j] += ($pdiff-$KTUP*($fdiff-1));
	      }
	      else {
	        $inserts[$j] += ($pdiff-$KTUP*($fdiff-1));
	      }
              $extended = 1;
	      # still want to create a new train from scratch
            }
	  }
	}
      }
     
      # third pass, try extend a train with deletions
      for (my $j=0; $j<=$#trains && $extended == 0; $j++) {
        my $r = $trains[$j];
	my @train = @$r;
	my $lastr = $train[$#train];
	if ($lastr->[$FRAGMENT] < $fragment &&
	    $lastr->[$POSITION]+$KTUP <= $position) {
          # can extend this train
	  my $fdiff = $fragment-$lastr->[$FRAGMENT];
	  my $pdiff = $position-($lastr->[$POSITION]+$KTUP);

	  if ($short[$j] == 0) {
	    my $new_hits = $hits[$j]+$hitscore;
	    my $g = ($pdiff-$KTUP*($fdiff-1));
	    my $new_mutations = $mutations[$j];
	    my $new_inserts = $inserts[$j];
	    my $new_deletes = $deletes[$j]+abs($g);

	    # hypothetical hits is current hits plus all remaining
	    # fragments
	    my $hypo_hits =
	      $KTUP*($#train+1)+$flength[$feature]-$fragment*$KTUP;
            my $escore = percent_identity
	      ($seqsize, $feature, $hypo_hits, $new_mutations,
	       $new_inserts, $new_deletes);
            if ($escore < $E_THRESH) {
	      my @train2 = @train;
              my @v;
	      $v[$FRAGMENT] = $fragment;
	      $v[$POSITION] = $position;
	      push @train2, \@v;
	      push @trains, \@train2;
	      my $lastf = $lastr->[$FRAGMENT];

	      $hits[$#trains] = $new_hits;
	      $mutations[$#trains] = $new_mutations;
	      $inserts[$#trains] = $new_inserts;
	      $deletes[$#trains] = $new_deletes;
	      $short[$#trains] = 0;

	      # print "duplicate ($j($lastf)->$#trains) $feature $fragment $position: $new_hits $new_mutations $new_inserts $new_deletes ($escore)\n";
	    }
	  }
	}
      }

      if ($createnew) {
        my @train3;
        my @v;
        $v[$FRAGMENT] = $fragment;
        $v[$POSITION] = $position;
        push @train3, \@v;
        push @trains, \@train3;

        $hits[$#trains] = $hitscore;
        $mutations[$#trains] = 0;
        $inserts[$#trains] = 0;
        $deletes[$#trains] = 0;
	if ($extended) {
	  # already extended a train, make new train a short train
          $short[$#trains] = 1;
	}
	else { $short[$#trains] = 0; }
	
        if (0 && ($feature == 425 || $feature == 426)) {
	  if ($extended) {
            print "short train ($#trains) $feature $fragment $position\n";
	  } 
	  else {
            print "new train ($#trains) $feature $fragment $position\n";
	  }
	}
      }
    }
  }

  my $lastp0;
  my $lastp1;

  for (my $i=0; $i<=$#trains; $i++) {
    my $r = $trains[$i];
    my @train = @$r;
    my $firstr = $train[0];
    my $lastr = $train[$#train];
    if (defined $lastp0) {
      if ($lastp0 == $firstr->[$POSITION] &&
          $lastp1 == $lastr->[$POSITION]+$KTUP-1) {
        next;
      }
      # insider knowledge: with the way we iterate the @trains array
      # in this loop, small fragments that are part of a longer train
      # will come after the long train
      if ($lastp0 <= $firstr->[$POSITION] &&
          $lastp1 >= $lastr->[$POSITION]+$KTUP-1) {
        next;
      }
    }
    if (0 && ($feature == 425 || $feature == 426)) {
      print "scoring $fnames[$feature] ($feature) $i\n";
    }
    my $fl = $flength[$feature];
    my $stopp = $lastr->[$POSITION];
    if ($lastr->[$FRAGMENT] == int((($fl-1)/$KTUP))) { # last fragment
      if (($fl%$KTUP) == 0) { $stopp += ($KTUP-1); }
      else { $stopp += (($fl%$KTUP)-1); }
    }
    else { $stopp += ($KTUP-1); }

    my $ir = check_feature
      ($plasmidid, $seqsize, $feature,
       $hits[$i], $mutations[$i], $inserts[$i], $deletes[$i],
       $firstr->[$FRAGMENT], $firstr->[$POSITION], $stopp);
    if ($ir) {
      $lastp0 = $firstr->[$POSITION];
      $lastp1 = $stopp;
    }
  }
}

true;

