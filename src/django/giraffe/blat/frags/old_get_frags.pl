#!/usr/local/bin/perl

use IO::Handle;
use Time::HiRes qw (gettimeofday);
STDERR->autoflush(1);
STDOUT->autoflush(1);

require "blastpl.pl";

$FINDORF = 1;
$VERBOSE = 1;
$BM      = 0;

sub get_features {
    my ( $r_files, $r_plids, $blprefix ) = @_;

    my @files = @$r_files;
    my @plids = @$r_plids;
    my @results;

    my ( $bm_s, $bm_ms ) = gettimeofday;
    if ($BM) { print "$bm_s,$bm_ms - pre-orf\n"; }

    if ($FINDORF) {
        for ( my $fi = 0 ; $fi <= $#files ; $fi++ ) {
            my $plasmidid = $plids[$fi];
            push @results, find_orfs( $files[$fi], $plasmidid );
        }
    }

    my ( $bm_s, $bm_ms ) = gettimeofday;
    if ($BM) { print "$bm_s,$bm_ms - orf\n"; }

    push @results, blast_sequences( \@files, $r_plids, $blprefix );
    return @results;
}

sub find_orfs {
    my ( $file, $plasmidid ) = @_;
    my @r = `/usr/addgene/events/findorfs < $file`;
    my @results;

    for ( my $i = 0 ; $i <= $#r ; $i++ ) {
        if ($VERBOSE) { print "$r[$i]"; }
        chomp $r[$i];
        my @v         = split " ", $r[$i];
        my $type      = "Orf";
        my $name      = $v[0];
        my $position  = $v[1];
        my $length    = $v[2];
        my $clockwise = 1;

        if ( $name =~ /^Orf\d$/ ) {
            $type = "Orf";
        }
        elsif ( $name =~ /^Orf(\d+),antisense$/ ) {
            $type      = "Orf";
            $position  = $position + $length - 1;
            $clockwise = 0;
            $name      = "Orf$1";
        }
        else {
            die "bad features list for $file\n";
        }
        push @results, "$type,$name,$position,$length,$clockwise,1";
    }
    return @results;
}

sub blast_sequences {
    my ( $r_files, $r_plids, $blprefix ) = @_;

    my @files = @$r_files;
    my @plids = @$r_plids;

    my $bldata = $blprefix . "data";
    my $args = join " ", @files;

    my ( $bm_s, $bm_ms ) = gettimeofday;
    if ($BM) { print "$bm_s,$bm_ms - blastinit\n"; }

    blast_init( $blprefix . "index" );

    my ( $bm_s, $bm_ms ) = gettimeofday;
    if ($BM) { print "$bm_s,$bm_ms - findfrags\n"; }

    my @r = `/usr/addgene/events/findfrags $bldata $args`;

    my ( $bm_s, $bm_ms ) = gettimeofday;
    if ($BM) { print "$bm_s,$bm_ms - gotfrags\n"; }

    my @results;
    my $seqtotal;
    my $fn;
    my @frags;
    my @empty;

    my $pli = 0;

    for ( my $ri = 0 ; $ri <= $#r ; $ri++ ) {
        my $rl = $r[$ri];
        chomp $rl;
        if ( $rl =~ /^======(\d+) ([\w\W]+)/ ) {
            if ( defined $fn ) {
                while ( $files[$pli] ne $fn && $pli <= $#files ) { $pli++; }
                if ( $pli > $#files ) { die "result does not match input\n"; }
                my $plasmidid = $plids[$pli];
                if ($VERBOSE) { print "$fn $files[$pli]\n"; }
                $pli++;

                my %fresults;
                my @ff;
                for ( $fi = 0 ; $fi <= $#frags ; $fi++ ) {
                    if ( $frags[$fi] =~ /^(\d+) (\d+) (\d+)/ ) {
                        my $f = $fresults{$1};
                        my @v;
                        if ( defined $f ) { @v = @$f; }
                        else {
                            push @ff, $1;
                        }
                        push @v, $frags[$fi];
                        $fresults{$1} = \@v;
                    }
                }
                my @fv = split /[\s.(\+]/, $fn;
                my @pl = split /\//,       $fv[0];

                for ( $fi = 0 ; $fi <= $#ff ; $fi++ ) {
                    my $f = $fresults{ $ff[$fi] };
                    my @v = @$f;
                    process_frags( $plasmidid, $seqtotal, @v );
                }
                push @results, insertFeatures( $plasmidid, $seqtotal );
            }
            $seqtotal = $1;
            $fn       = $2;
            @frags    = @empty;
        }
        else {
            push @frags, $rl;
        }
    }

    my ( $bm_s, $bm_ms ) = gettimeofday;
    if ($BM) { print "$bm_s,$bm_ms - processed frags\n"; }

    if ( defined $fn ) {
        while ( $files[$pli] ne $fn && $pli <= $#files ) { $pli++; }
        if ( $pli > $#files ) { die "result does not match input\n"; }
        my $plasmidid = $plids[$pli];
        if ($VERBOSE) { print "$fn $files[$pli]\n"; }
        $pli++;

        my %fresults;
        my @ff;
        for ( $fi = 0 ; $fi <= $#frags ; $fi++ ) {
            if ( $frags[$fi] =~ /^(\d+) (\d+) (\d+)/ ) {
                my $f = $fresults{$1};
                my @v;
                if ( defined $f ) { @v = @$f; }
                else {
                    push @ff, $1;
                }
                push @v, $frags[$fi];
                $fresults{$1} = \@v;
            }
        }
        my @fv = split /[\s.(\+]/, $fn;
        my @pl = split /\//,       $fv[0];

        for ( $fi = 0 ; $fi <= $#ff ; $fi++ ) {
            my $f = $fresults{ $ff[$fi] };
            my @v = @$f;
            process_frags( $plasmidid, $seqtotal, @v );
        }
        push @results, insertFeatures( $plasmidid, $seqtotal );
    }

    my ( $bm_s, $bm_ms ) = gettimeofday;
    if ($BM) { print "$bm_s,$bm_ms - done\n"; }

    return @results;
}

true;
