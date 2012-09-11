#!/usr/bin/env perl

use strict;
use warnings;
use File::Spec;
use FindBin;
use Getopt::Long;

sub print_usage {
  die "Usage: $0 [--configuration=Debug|Release] target0 [target1 ...]\n";
}

print_usage unless @ARGV;

my $config = 'Debug';
my $help;

GetOptions(
  'configuration=s' => \$config,
  'help' => \$help,
);

print_usage if $help;

my %targets;
for my $target (@ARGV) {
  my ($dir, $target_name) = split /:/, $target;
  $dir =~ s!/$!!;
  if ($target_name and (ref $targets{$dir} or not defined $targets{$dir})) {
    push @{ $targets{$dir} }, $target_name;
  } else {
    $targets{$dir} = 'ALL';
  }
}

my $common_file = File::Spec->catfile($FindBin::Dir, 'common.gypi');
for my $target_dir (keys %targets) {
  my $gyp_file = File::Spec->catfile(
    $FindBin::Dir, $target_dir, "$target_dir.gyp");

  my @targets = ref $targets{$target_dir} ? @{ $targets{$target_dir} } : ();

  my $gyp_status = system 'gyp',
    '--depth' => $FindBin::Dir,
    '--format' => 'make',
    '--include' => $common_file,
    $gyp_file;
  die 'gyp terminated with an error status' if $gyp_status != 0;
  my $make_status = system 'make', "BUILDTYPE=$config", @targets;
  die 'make terminated with an error status' if $make_status != 0;
}
