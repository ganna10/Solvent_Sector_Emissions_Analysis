#! /usr/bin/env perl

use strict;
use diagnostics;
use MECCA;
use KPP;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/Solvent_Emissions/MCM";
my (%families, %weights, %data);
$families{"HO2x"} = [ qw( HO2 HO2NO2 ) ];

opendir DIR, $base or die "Can't open $base : $!";
my @dirs = grep { $_ =~ /EMEP_tagged_|EMEP_Solvents_Only/ } readdir DIR;
close DIR;

foreach my $dir (@dirs) {
    print $dir, "\n";
    my $boxmodel = "$base/$dir/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $eqn = "$base/$dir/gas.eqn";
    my $kpp = KPP->new($eqn);
    my $ro2_file = "$base/$dir/RO2_species.txt";
    my @no2_reservoirs = get_no2_reservoirs($kpp, $ro2_file);
    $families{$dir} = [ qw( O3 O O1D NO2 NO3 N2O5 HO2NO2 ), @no2_reservoirs ]; ##full Ox family for dir
    $weights{$dir} = { NO3 => 2, N2O5 => 3 }; 
    my $label;
    if ($dir =~ /Solvents_Only/) {
        $label = "No_tagging";
    } else {
        ($label = $dir) =~ s/EMEP_(.*?)/$1/;
    }
    $data{$label} = get_data($mecca, $kpp, $dir);
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
);

$R->run(q` data = data.frame() `);
foreach my $run (sort keys %data) {
    print "$run: $data{$run} \n";
    $R->set('run', $run);
    $R->set('prod', $data{$run});
    $R->run(q` pre = data.frame(Dummy = c(1)) `,
            q` pre[run] = prod `,
            q` pre = gather(pre, Run, Prod, -Dummy) `,
            q` data = rbind(data, pre) `);
    #$R->run(q` pre = data.frame(Run = run) `);
#    foreach my $category (sort keys %{$data{$run}}) {
#        print "\t$category: $data{$run}{$category}\n";
#        $R->set('category', $category);
#        $R->set('prod', $data{$run}{$category});
#        $R->run(q` pre[category] = prod `);
#    }
#    $R->run(q` pre = gather(pre, Category, Prod, -Run) `,
#            q` data = rbind(data, pre) `,
#    );
}
$R->run(q` data = data[, -1] `);
my $p = $R->run(q` print(data) `);
print $p, "\n";

$R->run(q` p = ggplot(data, aes(x = Run, y = Prod)) `,
        q` p = p + geom_bar(stat = "identity") `,
);

$R->run(q` CairoPDF(file = "MCM_Ox_tagged_non_tagged.pdf") `,
        q` print(p) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($mecca, $kpp, $dir) = @_;
    
    my @loop = ("HO2x");
    my (%production_rates, %consumption_rates, $prod, $cons);
    foreach my $species (@loop) {
        $kpp->family({
                name    => $species,
                members => $families{$species},
                weights => $weights{$species},
        });
        my $producers = $kpp->producing($species);
        my $producer_yields = $kpp->effect_on($species, $producers);
        my $consumers = $kpp->consuming($species);
        my $consumer_yields = $kpp->effect_on($species, $consumers);
        print "No consumers for $species\n" if (@$consumers == 0);
        print "No producers for $species\n" if (@$producers == 0);

        for (0..$#$producers) {
            my $reaction = $producers->[$_];
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
            next if ($rate->sum == 0);
            $prod += $rate->sum;
            #my ($number, $parent) = split /_/, $reaction;
            #if (defined $parent) {
            #$production_rates{$parent} += $rate->sum;
            #} else {
            #my ($reactants, $products) = split / = /, $kpp->reaction_string($reaction);
            #$production_rates{$reactants} += $rate->sum;
            #}
       }
 
       for (0..$#$consumers) {
           my $reaction = $consumers->[$_];
           my $reaction_number = $kpp->reaction_number($reaction);
           my $rate = $mecca->rate($reaction_number) * $consumer_yields->[$_];
           next if ($rate->sum == 0);
            $cons += $rate->sum;
           #my ($number, $parent) = split /_/, $reaction;
           #if (defined $parent) {
           #$consumption_rates{$parent} += $rate->sum;
           #} else {
           #my ($reactants, $products) = split / = /, $kpp->reaction_string($reaction);
           #$consumption_rates{$reactants} += $rate->sum;
           #}
       }
    }
    $prod += $cons;
    #remove_common_processes(\%production_rates, \%consumption_rates);
    #my %common_processes;
    #$common_processes{$_} = 1 foreach ( grep { defined $production_rates{$_} } keys %consumption_rates );
    #$production_rates{$_} += $consumption_rates{$_} foreach (keys %common_processes);
    #my %final;
    #foreach my $item (keys %production_rates) {
        #print "$item: ", $production_rates{$item}->sum, "\n";
        #my $category = get_category($item);
        #print "$item: $category\n";
        #$final{$category} = $production_rates{$item};
        #}
    #print "$_: $final{$_}\n" foreach (keys %final);
    #return \%final;
    return $prod;
}

sub remove_common_processes {
    my ($production, $consumption) = @_;
    my %common_processes;
    $common_processes{$_} = 1 for (grep { defined $production->{$_} } keys %$consumption) ;

    foreach my $process (keys %common_processes) {
        my $net_effect = $production->{$process} + $consumption->{$process};
        #print $process, $net_effect->nelem, "\n";
        if ($net_effect->sum > 0) { #if net production remove consumption processes, retain net production
            if (which($net_effect < 0)->nelem > 0) { #which gets indices of non-0 values, nelem get nr of elements
                #print "which if $process $net_effect\n";
                $production->{$process} .= $net_effect;
                $production->{$process}->where($net_effect < 0) .= 0; #where returns output pdl corresponding to input pdl and meeting the condition
                $consumption->{$process} .= $net_effect;
                $consumption->{$process}->where($net_effect > 0) .= 0; #where returns output pdl corresponding to input pdl and meeting the condition
                next;
            } 
            $production->{$process} = $net_effect;
            delete $consumption->{$process};
        } else { #net consumption
            if (which($net_effect > 0)->nelem > 0) {
                $production->{$process} .= $net_effect;
                $production->{$process}->where($net_effect < 0) .= 0;
                $consumption->{$process} .= $net_effect;
                $consumption->{$process}->where($net_effect > 0) .= 0;
                next;
            }
            $consumption->{$process} = $net_effect;
            delete $production->{$process};
        }
    }
} 

sub get_category {
    my ($process, $mechanism, $speciation) = @_;
    my $category;
    if ($process =~ /CO \+ OH/) {
        $category = "CO";
    } elsif ($process =~ /\+/) {
        $category = "Inorganic";
    } elsif ($process eq "CH4") {
        $category = "Methane";
    } elsif ($process =~ /^TOL/) {
        $category = "Toluene";
    } elsif ($process =~ /XYL/) {
        $category = "Xylenes";
    } elsif ($process =~ /^TM/) {
        $category = "Trimethylbenzenes";
    } elsif ($process eq "BENZENE") {
        $category = "Benzene";
    } elsif ($process =~ /BENZ|DIME|STYRENE|TOL$/) {
        $category = "Other aromatics";
    } elsif ($process eq "ETH" or $process eq "C2H6") {
        $category = "Ethane";
    } elsif ($process eq "OL2" or $process eq "C2H4") {
        $category = "Ethene";
    } elsif ($process =~ /NC4|IC4/) {
        $category = "Butanes";
    } elsif ($process eq "HC8" or $process =~ /NC6|NC7|NC8|NC9|NC1\d|M\dPE|M\dHEX|CHEX/) {
        $category = "Higher alkanes";
    } elsif ($process =~ /CH3OH|C2H5OH|IPROPOL|NBUTOL|NPROPOL|BUT2OL|IBUTOL|MIBKAOH|C6H5CH2OH|ETHGLY|PROPGLY/) {
        $category = "Alcohols";
    } elsif ($process =~ /CH3OCH3|BUOX2ETOH|PR2OHMOX|EOX2EOL|MO2EOL/) {
        $category = "Ethers";
    } elsif ($process =~ /ACET$/) {
        $category = "Esters";
    } elsif ($process =~ /TRICLETH|CH3CCL3|CH2CL2|TCE/) {
        $category = "Chlorinated";
    } elsif ($process =~ /HCOOH|CH3COOH|PROPACID|CH3CO2H/) {
        $category = "Acids";
    } elsif ($process eq "KET" or $process =~ /MEK|CH3COCH3|MIBK|CYHEXONE/) {
        $category = "Ketones";
    } elsif ($process eq "ALD" or $process =~ /MACR|CH3CHO|C2H5CHO|C3H7CHO|IPRCHO|C4H9CHO|ACR|C4ALDB|HCHO|^CH2O/) {
        $category = "Aldehydes";
    } elsif ($process =~ /C5H8|ISO|BIGENE|BUT1ENE|C2H2/) {
        $category = "Other alkenes, alkynes, dienes";
    } elsif ($process eq "HC5" or $process eq "BIGALK" or $process =~ /NC5|IC5|NEOP/) {
        $category = "Pentanes";
    } elsif ($process eq "HC3" or $process eq "C3H8") {
        $category = "Propane";
    } elsif ($process eq "OLT" or $process eq "C3H6") {
        $category = "Propene";
    } elsif ($process =~ /C10H16|LIMONENE|PINENE|OLI/) {
        $category = "Terpenes";
    } elsif ($process =~ /Inorganic|Emissions/) {
        $category = $process;
    } else {
        $category = "Inorganic";
        #print "No category found for $process in $mechanism, $speciation\n";
    }
    return $category;
}

sub get_no2_reservoirs { #get species that are produced when radical species react with NO2
    my ($kpp, $file) = @_; 
    open my $in, '<:encoding(utf-8)', $file or die $!; 
    my @ro2;
    for (<$in>) {
        push @ro2, split /\s+/, $_; 
    }
    close $in;
    my @no2_reservoirs;
    foreach my $ro2 (@ro2) {
        my ($reactions) = $kpp->reacting_with($ro2, 'NO2');
        foreach my $reaction (@$reactions) {
            my ($products) = $kpp->products($reaction);
            if (@$products == 1) {
                push @no2_reservoirs, $products->[0];
            }   
        }   
    }   
    return @no2_reservoirs;
} 
