#! /usr/bin/env perl
# analysis of net rate of reactive carbon loss during 2-methoxyethanol (MO2EOL) degradation in each mechanism just TNO speciation
# Version 0: Jane Coates 4/3/2015

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/Solvent_Emissions";
my $mecca = MECCA->new("$base/RADM2/TNO_tagged_solvents_only/boxmodel");
my $times = $mecca->time;
my $NTIME = $mecca->time->nelem;
my $dt = $mecca->dt->at(0);
my $N_PER_DAY = 43200 / $dt;
my $N_DAYS = int $NTIME / $N_PER_DAY;

my @mechanisms = qw( MCM MOZART RADM2 );
my (%n_carbon, %data);

foreach my $mechanism (@mechanisms) {
    my ($boxmodel, $mecca, $eqnfile, $kpp, $carbon_file, $parent);
    if ($mechanism =~ /MOZ|RAD/) {
        $boxmodel = "$base/$mechanism/TNO_tagged_solvents_only/boxmodel";
        $mecca = MECCA->new($boxmodel); 
        $eqnfile = "$base/$mechanism/TNO_tagged_solvents_only/gas.eqn";
        $kpp = KPP->new($eqnfile);
        $carbon_file = "$base/$mechanism/TNO_tagged_solvents_only/carbons.txt";
        $n_carbon{$mechanism} = get_carbons($mechanism, $carbon_file);
    } else { 
        $boxmodel = "$base/$mechanism/TNO_tagged_solvents_only_esters/boxmodel";
        $mecca = MECCA->new($boxmodel); 
        $eqnfile = "$base/$mechanism/TNO_tagged_solvents_only_esters/gas.eqn";
        $kpp = KPP->new($eqnfile);
        $carbon_file = "$base/$mechanism/TNO_tagged_solvents_only_esters/carbons.txt";
        $n_carbon{$mechanism} = get_carbons($mechanism, $carbon_file);
        $parent = "MO2EOL";
    }
    if ($mechanism =~ /RA/) {
        $parent = "HC8";
    } elsif ($mechanism =~ /MOZ/) {
        $parent = "BIGALK";
    }
    $data{$mechanism} = get_data($kpp, $mecca, $n_carbon{$mechanism}, $parent);
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
);

$R->set('Time', [("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7")]);
$R->run(q` plot.data = data.frame() `);
foreach my $run (sort keys %data) {
    $R->set('mechanism', $run);
    $R->set('rate', [map { $_ } $data{$run}->dog]);
    $R->run(q` data = data.frame(Time) `);
    $R->run(q` data[mechanism] = rate `);
    $R->run(q` data = gather(data, Mechanism, Rate, -Time) `,
            q` plot.data = rbind(plot.data, data) `,
    );
}
#my $p = $R->run(q` print(data) `);
#print "$p\n";
$R->run(q` my.colours = c(  "MCM" = "#dc3522", "MOZART" = "#cc9900", "RADM2" = "#035c28") `);

$R->run(q` plot = ggplot(data = plot.data, aes(x = Time, y = Rate, colour = Mechanism, group = Mechanism)) `,
        q` plot = plot + geom_line() `,
        q` plot = plot + geom_point() `,
        q` plot = plot + theme_bw() `,
        q` plot = plot + ylab("Net Carbon Loss Rate (molecules cm-3 s-1)") `,
        q` plot = plot + theme(panel.grid = element_blank()) `,
        q` plot = plot + theme(legend.key = element_blank()) `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + theme(axis.title.y = element_text(face = "bold")) `,
        q` plot = plot + theme(axis.title.x = element_blank()) `,
        q` plot = plot + theme(axis.text.x = element_text(angle = 45, face = "bold", vjust = 0.5)) `,
        q` plot = plot + theme(legend.justification = c(1, 0)) `,
        q` plot = plot + theme(panel.border = element_rect(colour = "black")) `,
        q` plot = plot + theme(legend.position = c(1, 0)) `,
        q` plot = plot + theme(strip.background = element_blank()) `,
        q` plot = plot + theme(strip.text = element_text(face = "bold")) `,
        q` plot = plot + scale_colour_manual(values = my.colours) `,
);

$R->run(q` CairoPDF(file = "MO2EOL_net_reactive_carbon_loss.pdf", width = 7.5, height = 5.3) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($kpp, $mecca, $carbons, $VOC) = @_;
    my %carbons = %$carbons;
    my %carbon_loss_rate;
        
    ##get all reactions and loop over them
    my $all_reactions = $kpp->all_reactions();

    for (0..$#$all_reactions) { #get rates for all producing reactions
        my $reaction = $all_reactions->[$_];
        my ($r_number, $parent) = split /_/, $reaction; #remove tag from reaction number
        next unless (defined $parent and $parent eq $VOC);
        my $reaction_string = $kpp->reaction_string($reaction);
        $reaction_string =~ s/_(.*?)\b//g;
        next if ($reaction_string eq "CO + OH = HO2"); 
        next if (exists $carbon_loss_rate{$reaction_string});
        my ($net_carbon) = get_total_C($reaction_string, $carbons, $kpp);
        next if ($net_carbon == 0);
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $net_carbon * $mecca->rate($reaction_number); 
        next if ($rate->sum == 0); # do not include reactions that do not occur 
        $carbon_loss_rate{$reaction_string} += $rate(1:$NTIME-2);
    }

    foreach my $reaction (keys %carbon_loss_rate) {
        if ($VOC =~ /MO2EOL/) {
            $carbon_loss_rate{$reaction} = $carbon_loss_rate{$reaction} * 0.5; ###MO2EOL represents 2 VOC with 50:50 split
        } elsif ($VOC =~ /HC8/) {
            $carbon_loss_rate{$reaction} = $carbon_loss_rate{$reaction} * 2.8e-3; ###ratio of MO2EOL emissions to total HC8 emissions / 2
        } elsif ($VOC =~ /BIGALK/) {
            $carbon_loss_rate{$reaction} = $carbon_loss_rate{$reaction} * 1.7e-3; ###ratio of MO2EOL emissions to total BIGALK emissions / 2
        }
    }

    my $overall_carbon_loss_rate = 0;
    $overall_carbon_loss_rate += $carbon_loss_rate{$_} foreach (keys %carbon_loss_rate);
    $overall_carbon_loss_rate = $overall_carbon_loss_rate->reshape($N_PER_DAY, $N_DAYS);
    $overall_carbon_loss_rate = $overall_carbon_loss_rate->sumover;
    $overall_carbon_loss_rate = $overall_carbon_loss_rate(0:13:2);
    return $overall_carbon_loss_rate;
}

sub get_total_C {
    my ($reaction_string, $carbons, $kpp) = @_;
    my ($reactant_c, $product_c, @reactants, @products);

    my @inorganic = qw( hv OH HO2 O3 NO NO2 NO3 H2O HNO3 H2 PAROP O CO2 XO2 XO2N OHOP UNITY CL SO3 SO2 H2O2 O2 HONO NULL );
    my ($reactants, $products) = split / = /, $reaction_string;
    push @reactants, split / \+ /, $reactants;
    push @products, split / \+ /, $products;
    
    foreach my $reactant (@reactants) {
        next if ($reactant ~~ @inorganic);
        $reactant_c += get_species_carbon($reactant, $carbons);
    }
    
    return 0 unless (defined $reactant_c);
    foreach my $product (@products) {
        my ($yield, $item);
        if ($product =~ /^[0-9]|^\.[0-9]/) {
            ($yield, $item) = split / /, $product;
            next if ($item ~~ @inorganic);
            $product_c += $yield * get_species_carbon($item, $carbons);
        } else {
            next if ($product ~~ @inorganic);
            $product_c += get_species_carbon($product, $carbons);
        } 
    }
    $product_c = 0 unless (defined $product_c);
    return $product_c - $reactant_c;
}

sub get_species_carbon {
    my ($species, $carbons) = @_;
    my %carbons = %$carbons;
    my $carbon;
    if (exists $carbons{$species}) {
        $carbon = $carbons{$species};
    } else {
        print "No C found for species: $species\n";
    }
    return $carbon;
}

sub get_carbons {
    my ($run, $file) = @_;
    my $carbons;
    if ($run =~ /MCM/) {
        $carbons = mcm_n_carbon($file);
    } elsif ($run =~ /MOZART/) {
        $carbons = mozart_n_carbon($file);
    } elsif ($run =~ /CRI|RADM2|RACM|CB/) {
        $carbons = carbons_others($file);
    } else {
        print "$run doesn't match\n";
    }
    return $carbons;
}

sub mcm_n_carbon {
    my ($file) = @_;

    open my $in, '<:encoding(utf-8)', $file or die $!;
    my @lines = (<$in>);
    close $in;
    my %carbons;
    foreach my $line (@lines) {
        chomp $line;
        my ($species, $smile) = split ' ', $line;
        my $C_number = 0;
        if ($smile =~ /\./) {
            $C_number = 8;
        } else {
            $C_number++ while ($smile =~ m/C/gi);
        }
        $carbons{$species} = $C_number;
    }
    return \%carbons;
}

sub mozart_n_carbon {
    my ($file) = @_;

    open my $in, '<:encoding(utf-8)', $file or die $!;
    my $words = join ',', (<$in>);
    close $in;
    my ($string) = $words =~ /Solution(.*?)End\sSolution/s;
    $string =~ s/^\s+,//;
    $string =~ s/,\s+$//;
    $string =~ s/\s+/ /g;
    $string =~ s/RO2(.*?)->//;
    $string =~ s/ROOH(.*?)->//;
    my @species = split ',', $string;
    my %carbons;
    foreach my $species (@species) {
        $species =~ s/^\s+|\s+$//g;
        my $C_number = 0;
        if ($species !~ /->/ and $species !~ /(C[0-9])/) {
            $C_number ++ while ($species =~ m/C/g);
            $carbons{$species} = $C_number;
        } elsif ($species !~ /->/ and $species =~ /(C[0-9])/) { 
            my ($c_nr) = $species =~ /(C[0-9]+)/s;
            $c_nr =~ s/C//; 
            $C_number = $c_nr;
            $carbons{$species} = $C_number;
        } else {
            my ($mech, $molecule) = split ' -> ', $species;
            $mech =~ s/^\s+|\s+$//g;
            if ($molecule =~ /(C[0-9]+)/) { 
                my ($c_nr) = $molecule =~ /(C[0-9]+)/s;
                $c_nr =~ s/C//; 
                $C_number = $c_nr;
                $carbons{$mech} = $C_number;
            } else {
                $C_number ++ while ($molecule =~ m/C/g);
                $carbons{$mech} = $C_number;
            }
        }
    } 
    return \%carbons;
}

sub carbons_others { #get C-number from file names that have species and C# separated by space
    my ($file) = @_;

    open my $in, '<:encoding(utf-8)', $file or die "Cannot open file $file: $!";
    my (@lines) = (<$in>);
    close $in;
    my %carbons;
    foreach my $line (@lines) {
        my ($species, $C_number) = split '\s', $line;
        $carbons{$species} = $C_number;
    }
    return \%carbons;
}
