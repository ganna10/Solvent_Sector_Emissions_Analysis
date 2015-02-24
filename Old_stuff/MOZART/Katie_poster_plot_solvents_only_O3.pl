#!/usr/bin/perl
# plot O3 concentrations for Katie's poster
# Version 0: Jane Coates 7/11/2014

use strict;
use diagnostics;
use MECCA;
use PDL;
use PDL::NiceSlice;
use PDL::NetCDF;
use Statistics::R;

#Create x-axis for plot in hours
my $run = "/local/home/coates/Solvent_Emissions/MOZART/TNO_Solvents_Only/boxmodel";
my $mecca = MECCA->new($run); 
my $ntime = $mecca->time->nelem; #number of time points
my $times = $mecca->time;
$times -= $times->at(0);
$times = $times(1:$ntime-2);
$times /= 86400;
my $times_list = join ":", $times->dog;
my @time_axis = split /:/, $times_list;

my %concentration;
my $species = $ARGV[0];

#TNO data
my $tno_conc_file = "/local/home/coates/Solvent_Emissions/MOZART/TNO_Solvents_Only/boxmodel/mecca1_tracer.nc";
my $tno_nc = PDL::NetCDF->new($tno_conc_file);
print "TNO: ";
($concentration{"TNO"}) = get_concentration($species, $tno_nc, $ntime); 

#IPCC data
my $ipcc_conc_file = "/local/home/coates/Solvent_Emissions/MOZART/IPCC_Solvents_Only/boxmodel/mecca1_tracer.nc";
my $ipcc_nc = PDL::NetCDF->new($ipcc_conc_file);
print "IPCC: ";
($concentration{"IPCC"}) = get_concentration($species, $ipcc_nc, $ntime); 

#EMEP data
my $emep_conc_file = "/local/home/coates/Solvent_Emissions/MOZART/EMEP_Solvents_Only/boxmodel/mecca1_tracer.nc";
my $emep_nc = PDL::NetCDF->new($emep_conc_file);
print "EMEP: ";
($concentration{"EMEP"}) = get_concentration($species, $emep_nc, $ntime); 

#DE94 data
#my $de94_conc_file = "/local/home/coates/Solvent_Emissions/MOZART/DE94_Solvents_Only/boxmodel/mecca1_tracer.nc";
#my $de94_nc = PDL::NetCDF->new($de94_conc_file);
#($concentration{"DE94"}) = get_concentration($species, $de94_nc, $ntime); 

#GR95 data
#my $gr95_conc_file = "/local/home/coates/Solvent_Emissions/MOZART/GR95_Solvents_Only/boxmodel/mecca1_tracer.nc";
#my $gr95_nc = PDL::NetCDF->new($gr95_conc_file);
#($concentration{"GR95"}) = get_concentration($species, $gr95_nc, $ntime); 

#GR05 data
#my $gr05_conc_file = "/local/home/coates/Solvent_Emissions/MOZART/GR05_Solvents_Only/boxmodel/mecca1_tracer.nc";
#my $gr05_nc = PDL::NetCDF->new($gr05_conc_file);
#($concentration{"GR05"}) = get_concentration($species, $gr05_nc, $ntime); 

#UK98 data
#my $uk98_conc_file = "/local/home/coates/Solvent_Emissions/MOZART/UK98_Solvents_Only/boxmodel/mecca1_tracer.nc";
#my $uk98_nc = PDL::NetCDF->new($uk98_conc_file);
#($concentration{"UK98"}) = get_concentration($species, $uk98_nc, $ntime); 

#UK08 data
my $uk08_conc_file = "/local/home/coates/Solvent_Emissions/MOZART/UK08_Solvents_Only/boxmodel/mecca1_tracer.nc";
my $uk08_nc = PDL::NetCDF->new($uk08_conc_file);
print "UK08: ";
($concentration{"UK08"}) = get_concentration($species, $uk08_nc, $ntime); 


concentration_plot(\@time_axis, \%concentration, $ntime);

sub get_concentration {
    my ($species, $nc) = @_;
    
    my $concs += $nc->get($species)->squeeze;
    $concs = $concs(1:$ntime-2) * 1e9 * 2; #convert ppb to ug/m3
    print $concs->max, "\n";
    my $concs_list = join ":", $concs->dog;
    my @conc_array = split /:/, $concs_list;
    return \@conc_array;
}

sub concentration_plot {
    my ($time, $data) = @_;
    my %data = %$data;

    my $R = Statistics::R->new();
    $R->run(q` library(ggplot2) `);
    $R->run(q` library(scales) `);
    $R->run(q` library(grid) `);
    $R->run(q` library(gridExtra) `);
    $R->run(q` library(gtable) `);
    $R->run(q` library(Cairo) `);
    $R->set('time', [@$time]);   
    $R->set('main.title', "$species - Solvents Only Speciation");
    $R->set('file.name', "${species}_Solvents_Only_Speciation.png");
    $R->set('rep.number', 4); #number of model setups

    $R->run(q` times = rep(time, rep.number) `,
            q` Solvent_Speciation = {} `,
            q` Mixing.ratio = {} `,
    );

    foreach my $item (sort keys %data) {
        my $R_name = $R->set('name', $item);
        my $R_data = $R->set('mixing.ratio', [@{$data{$item}}]);
        $R->run(q` Solvent_Speciation = cbind(Solvent_Speciation, rep(name, length(time))) `,
                q` Mixing.ratio = cbind(Mixing.ratio, mixing.ratio) `,
        );
    }

    if ($species eq "O3") { #set y-axis & breaks
        $R->set('y.max', 200);
        $R->set('y.breaks', 25);
    } elsif ($species eq "CO") {
        $R->set('y.max', 250);
        $R->set('y.breaks', 50);
    } elsif ($species eq "OH") {
        $R->set('y.max', 7e-4) ;
        $R->set('y.breaks', 1e-4);
    }

    #create dataframe after converting the matrices above to vectors
    $R->run(q` Solvent_Speciation = c(Solvent_Speciation) `,
            q` Mixing.ratio = c(Mixing.ratio) `,
            q` data = data.frame(times, Solvent_Speciation, Mixing.ratio) `,
            q` data$Solvent_Speciation = factor(data$Solvent_Speciation, levels = c("IPCC", "EMEP", "TNO", "UK08")) `,
    ); 
    $R->run(q` my.colours = c("DE94" = "#010101", "EMEP" = "#f37d22", "GR05" = "#008c47", "GR98" = "#1859a9", "IPCC" = "#ed2d2e", "TNO" = "#662c91", "UK08" = "#12b2b2", "UK98" = "#b33893", "#a11d20") `);
    $R->run(q` plot.lines = function () { list( ylab(expression(bold(paste("Concentration (", mu,g, "/", m^3, ")")))),
                                                xlab("Time (days)"), 
                                                ggtitle(main.title), 
                                                geom_line(size = 3), 
                                                scale_x_continuous(limits = c(0, 7), breaks = seq(0, 7, 1)), 
                                                theme_bw(), 
                                                theme(axis.title.x = element_text(size = 50, face = "bold")), 
                                                theme(axis.title.y = element_text(size = 50)), 
                                                theme(panel.grid.major = element_blank()) ,
                                                theme(panel.grid.minor = element_blank()) ,
                                                scale_y_continuous(limits = c(75, y.max), breaks = seq(75, y.max, y.breaks)), 
                                                theme(axis.text.x = element_text(size = 35), axis.text.y = element_text(size = 35), legend.title = element_text(size = 45, face = "bold"), legend.key.size = unit(5, "cm"), legend.text = element_text(size = 40), legend.key = element_blank(), plot.title = element_text(size = 60, face = "bold")), 
                                                scale_colour_manual(values = my.colours) ) } `); 

    $R->run(q` plot1 = ggplot(data = data, aes(x = times, y = Mixing.ratio, colour = Solvent_Speciation)) `,
            q` plot1 = plot1 + plot.lines() `, 

            q` CairoPNG(file = file.name, width = 2000, height = 1500) `,
            q` print(plot1) `,
            q` dev.off() `,
    );

    $R->stop();
} 
