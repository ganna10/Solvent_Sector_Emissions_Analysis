#! /usr/bin/env perl
# Get differences in max O3 on 2nd day in each mechanism and each model run, output to tex file
# Version 0: Jane Coates 23/2/2015
# Version 1: Jane Coates 12/3/2015 adding differences between mechanisms of each speciation

use strict;
use diagnostics;
use MECCA;
use PDL;
use PDL::NiceSlice;

#my $base = "/work/users/jco/Solvent_Emissions";
my $base = "/local/home/coates/Solvent_Emissions";
my @mechanisms = qw(MCM MOZART RADM2);
my @speciations = qw( DE94 EMEP GR05 GR95 IPCC TNO UK08 UK98 );
my @runs = qw( Solvents_Only mean_NO_source_Solvents_Only );
my (%data, %speciation_data);

foreach my $mechanism (@mechanisms) {
    foreach my $speciation (@speciations) {
        foreach my $run (@runs) {
            my $dir = "$base/$mechanism/${speciation}_$run";
            my $mecca = MECCA->new("$dir/boxmodel");
            my $O3 = $mecca->tracer("O3") * 1e9;
            $O3 = $O3(72:144); #get 2nd day mixing ratios
            $data{$mechanism}{$run}{$speciation} = $O3->max;
            $speciation_data{$speciation}{$run}{$mechanism} = $O3->max;
        }
    }
}

my (%min_max_mech, %min_max_spec);
foreach my $mechanism (sort keys %data) {
    foreach my $run (sort keys %{$data{$mechanism}}) {
        (my $print_run = $run) =~ s/_/ /g;
        $print_run =~ s/(\w+)/\u$1/g;
        my @sorted = reverse sort { $data{$mechanism}{$run}{$a} <=> $data{$mechanism}{$run}{$b} } keys %{$data{$mechanism}{$run}};
        $min_max_mech{$mechanism}{$print_run}{"Maximum"}{$sorted[0]} = $data{$mechanism}{$run}{$sorted[0]}; 
        $min_max_mech{$mechanism}{$print_run}{"Minimum"}{$sorted[-1]} = $data{$mechanism}{$run}{$sorted[-1]}; 
        $min_max_mech{$mechanism}{$print_run}{"Difference"}{"difference"} = sprintf "%.2f", $data{$mechanism}{$run}{$sorted[0]} - $data{$mechanism}{$run}{$sorted[-1]};
    }
};

foreach my $speciation (sort keys %speciation_data) {
    foreach my $run (sort keys %{$speciation_data{$speciation}}) {
        (my $print_run = $run) =~ s/_/ /g;
        $print_run =~ s/(\w+)/\u$1/g;
        my @sorted = reverse sort { $speciation_data{$speciation}{$run}{$a} <=> $speciation_data{$speciation}{$run}{$b} } keys %{$speciation_data{$speciation}{$run}};
        $min_max_spec{$speciation}{$print_run}{"Maximum"}{$sorted[0]} = $speciation_data{$speciation}{$run}{$sorted[0]}; 
        $min_max_spec{$speciation}{$print_run}{"Minimum"}{$sorted[-1]} = $speciation_data{$speciation}{$run}{$sorted[-1]}; 
        $min_max_spec{$speciation}{$print_run}{"Difference"}{"difference"} = sprintf "%.2f", $speciation_data{$speciation}{$run}{$sorted[0]} - $speciation_data{$speciation}{$run}{$sorted[-1]};
    }
}

my $output = "2nd_day_differences.tex";
open my $out, '>:encoding(utf-8)', $output or die "Can't open $output : $!";
print $out "\\documentclass{report}\n";
print $out "\\usepackage{multirow}\n";
print $out "\\usepackage{booktabs}\n";
print $out "\\begin{document}\n";
print $out "\n\\begin{table}\n";
print $out "\t\\begin{center}\n";
print $out "\t\t\\makebox[\\textwidth][c]{\n";
print $out "\t\t\\begin{tabular}{lllll}\n";
print $out "\t\t\t\\toprule\n";
print $out "\t\t\t\\multirow{2}{*}{\\textbf{Mechanism}} & \\multirow{2}{*}{\\textbf{Run}} & \\textbf{Difference in O3} & \\textbf{Maximum O3} & \\textbf{Least O3} \\\\ \n";
print $out "\t\t\t & & \\textbf{Mixing Ratios} & \\textbf{Speciation} & \\textbf{Speciation} \\\\ \\toprule\n";

foreach my $mechanism (sort keys %min_max_mech) {
    foreach my $run (sort keys %{$min_max_mech{$mechanism}}) {
        print $out "\t\t\t$mechanism & $run & $min_max_mech{$mechanism}{$run}{'Difference'}{'difference'} & ", keys %{$min_max_mech{$mechanism}{$run}{"Maximum"}}, " & ", keys %{$min_max_mech{$mechanism}{$run}{"Minimum"}}, "\\\\ \n";
    }
    print $out "\t\t\t\\midrule\n" unless ($mechanism eq "RADM2");
}
print $out "\t\t\t\\bottomrule\n";
print $out "\t\t\\end{tabular}}\n";
print $out "\n\t\t\\caption{\\textbf{Second Day Differences in O3 mixing ratios between speciations}}\n";
print $out "\t\t\\end{center}\n";
print $out "\\end{table}\n";
print $out "\\newpage\n";
print $out "\n\\begin{table}\n";
print $out "\t\\begin{center}\n";
print $out "\t\t\\makebox[\\textwidth][c]{\n";
print $out "\t\t\\begin{tabular}{lllll}\n";
print $out "\t\t\t\\toprule\n";
print $out "\t\t\t\\multirow{2}{*}{\\textbf{Speciation}} & \\multirow{2}{*}{\\textbf{Run}} & \\textbf{Difference in O3} & \\textbf{Maximum O3} & \\textbf{Least O3} \\\\ \n";
print $out "\t\t\t & & \\textbf{Mixing Ratios} & \\textbf{Mechanism} & \\textbf{Mechanism} \\\\ \\toprule\n";

foreach my $speciation (sort keys %min_max_spec) {
    foreach my $run (sort keys %{$min_max_spec{$speciation}}) {
        print $out "\t\t\t$speciation & $run & $min_max_spec{$speciation}{$run}{'Difference'}{'difference'} & ", keys %{$min_max_spec{$speciation}{$run}{"Maximum"}}, " & ", keys %{$min_max_spec{$speciation}{$run}{"Minimum"}}, "\\\\ \n";
    }
    print $out "\t\t\t\\midrule\n" unless ($speciation eq "UK98");
}
print $out "\t\t\t\\bottomrule\n";
print $out "\t\t\\end{tabular}}\n";
print $out "\n\t\t\\caption{\\textbf{Second Day Differences in O3 mixing ratios between mechanisms}}\n";
print $out "\t\t\\end{center}\n";
print $out "\\end{table}\n";
print $out "\n\\end{document}\n";
close $out;
