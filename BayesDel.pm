=head1 LICENSE

Copyright 2020 QIMR Berghofer Medical Research Institute

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

Stephen Kazakoff <Stephen.Kazakoff@qimrberghofer.edu.au>

=cut

=head1 NAME

gnomAD

=head1 SYNOPSIS

mv BayesDel.pm ~/.vep/Plugins
./vep -i variations.vcf --plugin BayesDel,/path/to/BayesDel_nsfp33a_noAF

=head1 DESCRIPTION

A VEP plugin that retrieves BayesDel annotation

curl -sL http://www.bjfenglab.org/download/get_GRCh37_inst.sh | bash
curl -sL http://www.bjfenglab.org/download/get_GRCh38_inst.sh | bash

=cut

package BayesDel;

use strict;
use warnings;

use Data::Dumper;

use FileHandle;

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub new {
  my $class = shift;
  my $self = $class->SUPER::new(@_);

  my $dir = shift @{ $self->params };
  die("ERROR: BayesDel directory not specified\n") unless $dir;
  die("ERROR: BayesDel directory not found\n") unless -d $dir;

  my $params = $self->params_to_hash();
  $self->{params} = $params;

  $params->{per_bp} //= 3;
  $params->{min} //= -1.5;
  $params->{step} //= 0.01;

  my %offsets;

  if ($params->{per_bp} == 3) {
    my @alleles = qw/A C G T/;

    for my $ref_allele (@alleles) {
      my $offset = 0;

      for my $alt_allele (@alleles) {
        next if $ref_allele eq $alt_allele;

        $offsets{$ref_allele}{$alt_allele} = $offset++;
      }
    }
  }

  unless (%offsets) {
    die("ERROR: Could not retrieve offset values with per_bp=$params->{per_bp}\n");
  }

  $self->{offsets} = \%offsets;

  my %handles;
  for my $chrom (1..22, qw/X Y M/) {

    my $file = "$dir/$chrom";

    next unless -e $file;

    my $fh = FileHandle->new($file) or die("ERROR: Could not open file $file\n");

    binmode($fh) or die("ERROR: Could not binmode file $file\n");

    $handles{$chrom} = $fh;
  }

  $self->{file_handles} = \%handles;

  return $self;
}

sub feature_types {
  return ['Feature', 'Intergenic'];
}

sub get_header_info {
  return { BayesDel => 'BayesDel score' };
}


sub run {
  my ($self, $tva) = @_;

  my $vf = $tva->variation_feature;

  return {} unless $vf->{start} eq $vf->{end};

  # get allele, reverse comp if needed
  my $allele = $tva->variation_feature_seq;
  reverse_comp(\$allele) if $vf->{strand} < 0;
  
  return {} unless $allele =~ /^[ACGT]$/;

  my $chr = ($vf->{chr} =~ /MT/) ? 'M' : $vf->{chr};
  my $fh = $self->{file_handles}{$chr};

  return {} unless $fh;

  my $bp = $vf->{start};
  my $offset = $self->{offsets}{$vf->ref_allele_string}{$allele};

  my $params = $self->params();

  my $ret = seek($fh, (($bp - 1) * $params->{per_bp}) + $offset, 0);
  return {} unless $ret == 1;

  my $chars = read($fh, my $buf, 1);
  return {} unless $chars == 1;

  my $score;

  unless (ord($buf) == 255) {
    $score = ord($buf) * $params->{step} + $params->{min};
  }

  return { BayesDel => $score };
}

1;
