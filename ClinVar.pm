=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

 Stephen Kazakoff <sh.kazakoff@gmail.com>
    
=cut

=head1 NAME

 ClinVar

=head1 SYNOPSIS

 mv ClinVar.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin ClinVar,/path/to/output/b37/clinvar_alleles.single.b37.vcf.gz,/path/to/output/b37/clinvar_alleles.multi.b37.vcf.gz
 ./vep -i variations.vcf --plugin ClinVar,/path/to/output/b38/clinvar_alleles.single.b38.vcf.gz,/path/to/output/b38/clinvar_alleles.multi.b38.vcf.gz

=head1 DESCRIPTION

 A VEP plugin that retrieves ClinVar annotation from the VCF output
 generated using the MacArthur Lab ClinVar XML processing pipeline,
 available here:

 https://github.com/macarthur-lab/clinvar

 If you use this plugin or the MacArthur Lab ClinVar XML processing
 pipeline, please cite:

 Zhang X, Minikel EV, O'Donnell-Luria AH et al. ClinVar data parsing
 [version 1; referees: 2 approved]. Wellcome Open Res 2017, 2:33 (doi:
 10.12688/wellcomeopenres.11640.1)

 This plugin simply retrieves matches using reference coordinates and
 an alternate allele sequence. Simple variants are ordered first, and
 any additional variation types are added to a column called
 'ADDITIONAL_VARIATION_TYPES'.

 A column called 'TRANSCRIPT_MATCH' is also added to indicate whether
 the current HGVS transcript ID matches to any of the ClinVar HGVS_C
 (or HGVS_N, if present) annotation.

 The tabix utility must be installed in your path to use this plugin.


=cut

package ClinVar;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::VEP qw(parse_line);

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  # test tabix
  `tabix --version` || die("ERROR: tabix does not seem to be in your path\n");

  $self->expand_left(0);
  $self->expand_right(0);

  $self->get_user_params();

  $self->{_prefix} = 'ClinVar';

  my ($single, $multi) = @{ $self->files() };

  my $single_headers = $self->parse_info_fields($single);

  unless (scalar(keys %$single_headers) > 1) {
    die("ERROR: Could not read column headers from $single\n");
  }

  if ($multi) {
    my $multi_headers = $self->parse_info_fields($multi);

    if (keys %$multi_headers != keys %$single_headers) {
      die("ERROR: The input files provided do not contain identical INFO fields\n");
    }
  }

  my %headers = %$single_headers;

  $headers{$_} = $_ for qw(TRANSCRIPT_MATCH ADDITIONAL_VARIATION_TYPES);

  $self->{headers} = \%headers;

  return $self;
}

sub parse_info_fields {
  my ($self, $file) = @_;

  my %headers;

  open(my $fh, '-|', "tabix -H $file") or die $!;

  while (<$fh>) {
    chomp;

    next unless /
      \#\#INFO=
      <
        ID=(?<id>.*?),
        Number=(?<number>.*?),
        Type=(?<type>.*?),
        Description="(?<description>.*?)"
      >
    /x;

    $headers{$+{id}} = $+{description};
  }

  close $fh;

  return \%headers;
}

sub feature_types {
  return ['Feature', 'Intergenic'];
}

sub get_header_info {
  my $self = shift;

  my %headers = %{ $self->{headers} };

  $headers{ join('_', $self->{_prefix}, $_) } = delete $headers{$_} for keys %headers;

  return \%headers;
}

sub run {
  my ($self, $vfoa) = @_;

  my $vf = $vfoa->variation_feature;

  my @results = map {
    $_->{result}
  }
  grep {
    $_->{start} eq $vf->{start} &&
    $_->{end} eq $vf->{end} &&
    $_->{alt} eq $vfoa->variation_feature_seq
  }
  @{ $self->get_data($vf->{chr}, $vf->{start} - 1, $vf->{end}) };

  return {} unless @results;

  # create a copy of the variants, ensuring any simple variants are listed first

  my @simple_variants;
  my @complex_variants;

  for (@results) {
    if ($_->{VARIATION_TYPE} eq 'Variant') {
      push(@simple_variants, $_);
    }
    else {
      push(@complex_variants, $_);
    }
  }

  my @variants = (@simple_variants, @complex_variants);

  # use the first variant in the list to report results

  my $result = shift @variants;

  # work on a copy of the results with the 'ClinVar_' prefix

  my %h = map { join('_', $self->{_prefix}, $_) => $result->{$_} } keys %$result;

  # add any additional variation types found

  if (@variants) {

    my %add_vars = map { $_->{VARIATION_TYPE} => undef } @variants;

    $h{ClinVar_ADDITIONAL_VARIATION_TYPES} = join(';', sort keys %add_vars);
  }

  # check if the HGVS transcript exists in a set of ClinVar HGVS values

  if ($vfoa->isa('Bio::EnsEMBL::Variation::TranscriptVariationAllele')) {

    if ($vfoa->hgvs_transcript) {

      my @hgvs_c = split(/\|/, $result->{'HGVS_C'}) if $result->{'HGVS_C'};
      my @hgvs_n = split(/\|/, $result->{'HGVS_N'}) if $result->{'HGVS_N'};

      my %hgvs = map { $_ => undef } (@hgvs_c, @hgvs_n);

      $h{ClinVar_TRANSCRIPT_MATCH} = 'YES' if exists $hgvs{$vfoa->hgvs_transcript};
    }
  }

  return \%h;
}

sub parse_data {
  my ($self, $line) = @_;

  my ($vf) = @{ parse_line({format => 'vcf', minimal => 1}, $line) };

  my ($ref, $alt) = split(/\//, $vf->allele_string);

  my %result = map { split /=/ } grep { /=/ } split(/;/, (split(/\t/, $line))[-1]);

  return {
    chr => $vf->{chr},
    start => $vf->{start},
    end => $vf->{end},
    ref => $ref,
    alt => $alt,
    result => \%result,
  };
}

sub get_start {
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}

1;

