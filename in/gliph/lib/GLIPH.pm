package GLIPH;
use strict;

# Author:  Jacob Glanville 
# Contact: jakeg@stanford.edu

# Constructor ----------------------------------------------
sub new {
  my ($class) = @_;
  my $self = {};
  $self->{tcrfile}               = ""; # tcr file
  $self->{corename}              = "cag"; 
  $self->{unique_h3s}            = {}; 
  $self->{redundant_h3s}         = [];

  $self->{refdb_hash_h3s}        = {};
  $self->{refdb_unique_h3s}      = [];
  $self->{refdb_nseqs}           =  0;
  $self->{refdb_kmers}           = {};
  $self->{refdb_redundant_nseqs} = [];

  $self->{nseqs}                 =  0;
  $self->{nseqs_redundant}       =  0;
  $self->{sample_kmers}          = {};
  $self->{sample_highcount_kmer_array} = [];

  $self->{kmer_sim_log}          = "";
  $self->{clone_network}         = "";
  $self->{convergence_table}     = "";
  $self->{global_similarity_table}= "";
  $self->{localConvergenceMotifList} ="";
  
  $self->{globalConverganceCutoff} = "";

  bless $self,'GLIPH';
  return $self;
}

# Methods --------------------------------------------------

sub loadTCRs {
  my($self,$textfile)=@_;

  $self->{corename}=$textfile;
  $self->{corename}=~s/.txt$//;

  open(FILE,$textfile);
  my @lines=<FILE>;
  chomp(@lines);
  close(FILE);

  for(my $x=0;$x<scalar(@lines);$x++){
    if($lines[$x]=~m/^C[AC-WY][AC-WY][AC-WY][AC-WY]*F/){
      my @fields=split(/\t/,$lines[$x]);
      push @{$self->{redundant_h3s}},$fields[0];
      if(defined($self->{unique_h3s}{$fields[0]})){
        $self->{unique_h3s}{$fields[0]}++;
      }else{
        $self->{unique_h3s}{$fields[0]}=1;
      }
    }
  }
  $self->{nseqs}=scalar(keys %{$self->{unique_h3s}});
  $self->{nseqs_redundant}=scalar(@{$self->{redundant_h3s}});

}

sub loadRefDB {
  my($self,$refdb,$discontinuous)=@_;

  if(-f $refdb){
    open(FILE,$refdb);
    my @lines=<FILE>;
    chomp(@lines);
    close(FILE);
    my $count=0;
    #my %refdb_hash_h3s=();
    for(my $x=0;$x<scalar(@lines);$x++){
      if($lines[$x]=~m/>/){
        $count++;
        my @fields=split(/;/,$lines[$x]);
        my $h3=$fields[4];
        # $self->{refdb_redundant_nseqs} = [];
        if(defined($self->{refdb_hash_h3s}{$fields[0]})){
          $self->{refdb_hash_h3s}{$fields[0]}++;
        }else{
          $self->{refdb_hash_h3s}{$fields[0]}=1;
        }
      }
    }
  }
  @{$self->{refdb_unique_h3s}} = keys %{$self->{refdb_hash_h3s}};
  $self->{refdb_nseqs} = scalar(@{$self->{refdb_unique_h3s}});
 
  %{$self->{refdb_kmers}}=$self->getAllKmerLengths(\@{$self->{refdb_unique_h3s}},$discontinuous);

}

sub stochasticMotifSample {
  my($self,$samplingDepth,$discontinuous)=@_;

  print "Subsampling (depth " . $self->{nseqs} . "): ";
  my $iter=0;
  for(my $sim=0;$sim<$samplingDepth;$sim++){
    my @random_subsample=$self->randomSubsample(\@{$self->{refdb_unique_h3s}},$self->{nseqs_redundant});
    print "\tMONKEY obtained " . scalar(@random_subsample) . " subsample $sim from reference database containing " . scalar(@{$self->{refdb_unique_h3s}}) . " reads\n";
    my %refdb_subsample_kmers=$self->getAllKmerLengths(\@random_subsample,$discontinuous);
    print "\tMONKEY count " . $refdb_subsample_kmers{"AGD"} . "\n";
    $self->appendToKmerLog(\@{$self->{sample_highcount_kmer_array}},\%refdb_subsample_kmers,"sim-$sim",$self->{kmer_sim_log});

    if($iter>49){
      print "\n";
      print "Subsampling (depth " . $self->{nseqs} . "): ";
      $iter=0;
    }
    print "#";
    $| = 1;
    $iter++;
  }
}

sub analyzeKmerLog {
  my($self,$simdepth,$minfoldchange,$minp)=@_;
  #my($kmer_sim_log,$localConvergenceMotifList,$samplingDepth,
  #   $nseqs,$motif_min_foldchange,$motif_minp,$refdb_kmers,$refdb_nseqs)=@_;
  #my($logfile,$localConvergenceMotifList,$simdepth,$seqspersim,$minfoldchange,$minp,$refdb_kmers,$refdb_nseqs)=@_;

  # $self->{kmer_sim_log};
  # $self->{localConvergenceMotifList};
  # $samplingDepth
  my $seqspersim = $self->{nseqs};
  # $motif_min_foldchange
  # $motif_minp
  # $self->{refdb_kmers}
  # $self->{refdb_nseqs}

  my $localConvergenceMotifList = $self->{localConvergenceMotifList};
  open(LOG,">$localConvergenceMotifList");

  open(FILE,$self->{kmer_sim_log});
  my @lines=<FILE>;
  chomp(@lines);
  close(FILE);

  my @motifs=split(/ /,$lines[0]);
  my @discovery_sample_counts=split(/ /,$lines[1]);

  print LOG "Motif\tCounts\tavgRef\ttopRef\tOvE\tp-value\n";
  print "Motif\tCounts\tavgRef\ttopRef\tOvE\tp-value\n";

  for(my $m=1;$m<scalar(@motifs);$m++){

    # get number of simulations at the level observed in discovery sample
    # get highest
    # get average (median?)
    my $highest=0;
    my $average=0;
    my $odds_as_enriched_as_discovery=0;
    for(my $sim=2;$sim<scalar(@lines);$sim++){
      my @fields=split(/ /,$lines[$sim]);
      if($fields[$m]>=$discovery_sample_counts[$m]){
        $odds_as_enriched_as_discovery += 1/$simdepth;
      }
      if($fields[$m]>$highest){
        $highest=$fields[$m];
      }
      $average += ($fields[$m]/(scalar(@lines)-2));
    }

    # get observed vs expected
    my $ove = 0;
    if($average>0){
      $ove=$discovery_sample_counts[$m]/$average;
    }else{
      # if found 0, assume you just missed it. pseudocount of 1
      $ove = 1 / ($simdepth * $seqspersim);
      #$ove=">" . $discovery_sample_counts[$m];
    }
    $ove=( int($ove * 1000) / 1000);
    $average=( int($average * 100) / 100);
    # calculate fisher's exact by generating a confusion matrix
    # of (counts in $seqspersim) vs ( $$refdb_kmers{kmer} in $refdb_nseqs) 
    #               motif    !motif
    #  discovery    n11      n12   | n1p
    #  refdb        n21      n22   | n2p
    #              --------------
    #               np1      np2   npp

    # 3 100 3 100
    #my $n11 = $discovery_sample_counts[$m];
    #my $n1p = $seqspersim;
    #my $np1 = $discovery_sample_counts[$m];
    #if(defined($$refdb_kmers{$motifs[$m]})){
    #  $np1 += $$refdb_kmers{$motifs[$m]};
    #}
    #my $npp = $seqspersim + $refdb_nseqs;
    #print "$n11 $n1p $np1 $npp | $refdb_nseqs\n";
    #my $left_value = calculateStatistic( n11=>$n11,
    #                                     n1p=>$n1p,
    #                                     np1=>$np1,
    #                                     npp=>$npp);

    if($odds_as_enriched_as_discovery<$minp){
      if($odds_as_enriched_as_discovery == 0){
        $odds_as_enriched_as_discovery=(1/$simdepth);
      }
      my $this_minfoldchange=$self->motifMinFoldPerCounts($discovery_sample_counts[$m]);
      # if($ove>=$minfoldchange){
      if($ove>=$this_minfoldchange){ #$minfoldchange){
        print LOG        $motifs[$m]
          . "\t" . $discovery_sample_counts[$m]
          . "\t" . $average
          . "\t" . $highest
          . "\t" . $ove
          . "\t" . $odds_as_enriched_as_discovery
          . "\n"; #. "\tfisher=$left_value\n";
        print        $motifs[$m]
          . "\t" . $discovery_sample_counts[$m]
          . "\t" . $average
          . "\t" . $highest
          . "\t" . $ove
          . "\t" . $odds_as_enriched_as_discovery
          . "\n";
      }
    }
  }
  close(LOG);
}

sub motifMinFoldPerCounts {
  my($self,$motif_count)=@_;
  # minimum fold-change for the motif to count
  # there is a motif count to fold-change relationship
  if($motif_count == 2){
    return 1000;
  }elsif($motif_count == 3){
    return 100;
  }elsif($motif_count == 4){
    return 10;#50;
  }elsif($motif_count > 4){
    return 10;#25;
  }
}

sub obtainLocalMotifs {
  my($self,$kmer_mindepth,$discontinuous)=@_;

  my @unique_cdrs = keys %{$self->{unique_h3s}};

  # already have my @unique_cdrs=keys %unique_h3s; # monkey
  my %sample_kmers=$self->getAllKmerLengths(\@unique_cdrs,$discontinuous);
  my @sample_kmer_array=keys %sample_kmers;
  print "  kmers obtained:             " . scalar(@sample_kmer_array) . "\n";
  my @sample_highcount_kmer_array = $self->getMinDepthKmerArray(\%sample_kmers,$kmer_mindepth);
  print "  mindepth>=$kmer_mindepth kmers obtained: " . scalar(@sample_highcount_kmer_array) . "\n";

  # store the kmers above the mindepth
  @{$self->{sample_highcount_kmer_array}}=@sample_highcount_kmer_array;

  # print kmer scores to the kmer subsampling log
  $self->appendToKmerLog(\@sample_highcount_kmer_array,\%sample_kmers,"Discovery",$self->{kmer_sim_log}); 
}

sub appendToKmerLog {
  # $self->appendToKmerLog(\@{$self->{sample_highcount_kmer_array}},\%refdb_subsample_kmers,"sim-$sim",$self->{kmer_sim_log});
  my($self,$sample_highcount_kmer_array,$sample_kmers,$label,$kmer_sim_log)=@_;
  if($label eq "Discovery"){
    open(LOG,">$kmer_sim_log");
    print LOG "Sample";
    for(my $x=0;$x<scalar(@$sample_highcount_kmer_array);$x++){
      print LOG " " . $$sample_highcount_kmer_array[$x];
    }
    print LOG "\n";
    close LOG;
  }

  open(LOG,">>$kmer_sim_log");

  print LOG $label;
  print "\tMONKEY\tlabel=($label)\n";
  for(my $x=0;$x<scalar(@$sample_highcount_kmer_array);$x++){
    if($$sample_highcount_kmer_array[$x] eq "AGD"){
      print "\tMONKEY\ttest=" . $$sample_highcount_kmer_array[$x] . "\t" . $$sample_kmers{$$sample_highcount_kmer_array[$x]} . "\n";
    }
    if(defined($$sample_kmers{$$sample_highcount_kmer_array[$x]})){
      print LOG " " . $$sample_kmers{$$sample_highcount_kmer_array[$x]};
    }else{
      print LOG " 0";
    }
  }
  print LOG "\n";
  close(LOG);
}



sub getMinDepthKmerArray {
  my($self,$kmer_hash,$mindepth)=@_;
  my @all_kmers=keys %$kmer_hash;
  my @selected_kmers=();
  for(my $x=0;$x<scalar(@all_kmers);$x++){
    if($$kmer_hash{$all_kmers[$x]}>=$mindepth){
      push @selected_kmers,$all_kmers[$x];
    }
  }
  return @selected_kmers;
}

sub getAllKmerLengths {
  my($self,$unique_cdrs,$discontinuous)=@_;
  my %sample_kmers=();

  my %sample_2mers=$self->getKmers(\@$unique_cdrs,2);
  my %sample_3mers=$self->getKmers(\@$unique_cdrs,3);
  my %sample_4mers=$self->getKmers(\@$unique_cdrs,4);

  if($discontinuous){
    my %sample_xxox4mers=$self->getKmers(\@$unique_cdrs,"xx.x");
    my %sample_xoxx4mers=$self->getKmers(\@$unique_cdrs,"x.xx");
    %sample_kmers=(%sample_2mers,%sample_3mers,%sample_4mers,%sample_xxox4mers,%sample_xoxx4mers);
  }else{
    %sample_kmers=(%sample_2mers,%sample_3mers,%sample_4mers); 
  }

  return %sample_kmers;
}

sub getKmers {
  my($self,$h3_array,$kmer_size)=@_;
  my %kmer_hash=();
  my $nseqs=scalar(@$h3_array);
  my $mask_type=$kmer_size;
  for(my $s=0;$s<$nseqs;$s++){
    my $seq=$$h3_array[$s];
       $seq=~s/^...//;
       $seq=~s/...$//;

    my $length = length($seq);

    # dealing with x.x, x..x, or x...x
    my $mask=0;
    if($mask_type =~ m/x/){
      $mask=1;
      $kmer_size=length($mask_type);
    }

    # run search
    for(my $p=0;($p+$kmer_size)<=$length;$p++){
      my $kmer = substr($seq,$p,$kmer_size);
      # dealing with x.x, x..x, or x...x
      if($mask){
        my @chars=split(/ */,$kmer);
        if($mask_type eq "xx.x"){
          $chars[2]=".";
        }elsif($mask_type eq "x.xx"){
          $chars[1]=".";
        }
        $kmer=join("",@chars);
      }

      if(defined($kmer_hash{$kmer})){
        $kmer_hash{$kmer}++;
      }else{
        $kmer_hash{$kmer}=1;
      }
    }
  }
  return %kmer_hash;
}





sub setGlobalDistCutoff {
  my($self,$globalConverganceCutoff)=@_;

  $self->{globalConverganceCutoff}=$globalConverganceCutoff;

#depth dist1 dist2 dist3
#25   0.0016  0.00239 0.012
#50   0.0004  0.0032  0.0274
#75   0.0008  0.00626 0.03133
#100  0.00059 0.0076  0.04969
#125  0.00176 0.01024 0.06328 (0.01 cutoff for dist 2)
#150  0.00119 0.0122  0.06913
#200  0.00225 0.01725 0.0896
#300  0.00333 0.02553 0.11913
#400  0.00457 0.03215 0.14322
#500  0.00449 0.03854 0.16704
#600  0.00695 0.0478  0.18053 (0.05 cutoff for dist 2 - 650)
#700  0.00671 0.05374 0.20252
#800  0.00737 0.05976 0.21313
#900  0.00873 0.06694 0.22701
#1000 0.00996 0.07113 0.23726 (0.01 cutoff for dist 1)
#2000 0.0193  0.12085 0.29814
#3000 0.02916 0.15676 0.3221
#4000 0.03842 0.18685 0.32965
#5000 0.04548 0.2138 0.32751   (0.05 cutoff for dist 1  - 5500)
#6000 0.05694 0.23055 0.32172

  # nseqs<125 = 0.01, nseqs<650 = 0.05, nseqs<1000 = 0.01, nseqs<5500 = 0.05
  if( $self->{globalConverganceCutoff} eq ""){
    if($self->{nseqs}<125){
      $self->{globalConverganceCutoff}=2;
    }elsif($self->{nseqs}<5500){
      $self->{globalConverganceCutoff}=1;
    }else{
      $self->{globalConverganceCutoff}=0; # technically it gets riskier at this point
    }
  }
}

sub getGlobalDistCutoff {
  my($self)=@_;
  return $self->{globalConverganceCutoff};
}

sub makeDepthFigTable { 
  my($self,$sampling_depth)=@_;

  print "  simulated stochastic resampling of depth $self->{nseqs} ";
  print "unique seqs out of $self->{nseqs_redundant} seqs from $self->{refdb}\n";
  $self->getGlobalConvergenceDistribution(\@{$self->{redundant_h3s}},"Sample",$self->{global_similarity_table});
 
  for(my $s=0;$s<$sampling_depth;$s++){
    my @random_subsample=$self->randomSubsample(\@{$self->{refdb_unique_h3s}},$self->{nseqs_redundant});
    $self->getGlobalConvergenceDistribution(\@random_subsample,"Sim$s",$self->{global_similarity_table});
  }
  exit;
}

sub randomSubsample {
  my($self,$array,$depth)=@_;
  my @id_array=();
  my @random_subsample=();

  unless(defined($depth)){
    $depth=scalar(@$array);
  }
  if($depth>scalar(@$array)){
    $depth=scalar(@$array);
  }

  $self->fisher_yates_shuffle(\@$array);

  for(my $s=0;$s<$depth;$s++){
    push @random_subsample,$$array[$s];
  }
  return @random_subsample;
}

sub fisher_yates_shuffle {
  my ($self,$array) = @_;
  my $i = @$array;
  while ( --$i ) {
    my $j = int rand( $i+1 );
    @$array[$i,$j] = @$array[$j,$i];
  }
}





sub getGlobalConvergenceDistribution {
  my($self,$sequences,$samplename,$global_similarity_table)=@_;

  open(FILE,">>$global_similarity_table");

  my @counts=(0,0,0,0,0, 0,0,0,0,0,
              0,0,0,0,0, 0,0,0,0,0,
              0,0,0,0,0, 0,0,0,0,0,
              0,0,0,0,0, 0,0,0,0,0);
  for(my $x=0;$x<scalar(@$sequences);$x++){
    my $distance=length($$sequences[$x]);
    for(my $y=0;$y<scalar(@$sequences);$y++){
      if($x != $y){
        if(length($$sequences[$x]) == length($$sequences[$y])){
          my $this_dist=$self->getHammingDist($$sequences[$x],$$sequences[$y]);
          if($this_dist<$distance){
            $distance=$this_dist;
          }
        }
      }
    }
    $counts[$distance]++;
  }

  # print the result
  print $samplename;
  print FILE $samplename;
  for(my $x=0;$x<13;$x++){
    print "\t" . $counts[$x];
    print FILE "\t" . $counts[$x];
  }
  print "\n";
  print FILE "\n";
  return @counts;

  close(FILE);
}

sub getHammingDist {
  my($self,$seq1,$seq2)=@_;
  my @chars1=split(/ */,$seq1);
  my @chars2=split(/ */,$seq2);

  my $mismatch_columns=0;

  for(my $c=0;$c<scalar(@chars1);$c++){
    if($chars1[$c] ne $chars2[$c]){
      $mismatch_columns++;
    }
  }
  return $mismatch_columns;
}



sub setOutputFiles {
  my($self,$samplingDepth,$motif_minp,$motif_min_foldchange)=@_;

  $self->{kmer_sim_log}            = $self->{corename} . "-kmer_resample_" . $samplingDepth . "_log.txt";
  $self->{clone_network}           = $self->{corename} . "-clone-network.txt";
  $self->{convergence_table}       = $self->{corename} . "-convergence-groups.txt";
  $self->{global_similarity_table} = $self->{corename} . "-global-similarity.txt";

  $self->{localConvergenceMotifList} = $self->{corename} . "-kmer_resample_" . $samplingDepth . "_minp"
                             . $motif_minp . "_ove" . $motif_min_foldchange . ".txt";

}

# if(defined(${$self->{seqnames}}{$seqname})){
#    if(defined(${$self->{msa}}[${$self->{seqnames}}{$seqname}][$position])){
#      return ${$self->{msa}}[${$self->{seqnames}}{$seqname}][$position];

1;
