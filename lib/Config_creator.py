import os
import shutil
from lib.Logger import log_error


def create_config(config):
    circos_path = shutil.which('circos')
    if circos_path == '':
        log_error("No circos executable found")

    circos_config_path = circos_path.replace('bin/circos', 'etc')

    template = f"""
    <colors>
        <<include {circos_config_path}/colors.conf>>
    </colors>
    
    <fonts>
        <<include {circos_config_path}/fonts.conf>>
    </fonts>
    
    show_ticks			= yes
    show_tick_labels	= yes
    
    
    <ticks>
        radius		= dims(ideogram,radius_outer)
        multiplier	= 0.001
        color		= black
        thickness	= 2p
        size		= 15p
        
        <tick>
            spacing			= 10u
            show_label		= yes
            label_size		= 20p
            label_offset	= 10p
            format			= %d
            suffix = " kb"
        </tick>
        
        <tick>
            spacing	= 1u
            color	= grey
            size	= 10p
        </tick>
    </ticks>
    
    karyotype = circos_input_files/karyotype.txt
    
    <ideogram>
        <spacing>
            default = 0.005r
            break = 0u
        </spacing>
        fill			= yes
        radius			= 0.60r
        thickness		= 30p
        stroke_thickness	= 2
        stroke_color	= black
        
        show_label		= no
        label_radius	= dims(ideogram,radius) + 0.075r
        label_size		= 24
        label_parallel	= yes
    </ideogram>
    
    chromosomes_units = 1000
    chromosomes_display_default = yes
    
    <image>
        <<include {circos_config_path}/image.conf>>
    </image>
    
    <<include {circos_config_path}/colors_fonts_patterns.conf>>
    
    <<include {circos_config_path}/housekeeping.conf>>
    
    <links>
        <link>
            ribbon			= yes
            color			= green_a4
            stroke			= yes
            stroke_color	= green_a4
            stroke_thickness	= 0.5
            file			= circos_input_files/rRNA_links
            radius			= 0.9r
            bezier_radius	= 0.3r
            thickness		= 1
        </link>
        <link>
            ribbon			= yes
            color			= red_a4
            stroke			= yes
            stroke_color	= red_a4
            stroke_thickness	= 0.5
            file			= circos_input_files/tRNA_links
            radius			= 0.9r
            bezier_radius	= 0.3r
            thickness		= 1
        </link>
        <link>
            ribbon			= yes
            color			= blue_a4
            stroke			= yes
            stroke_color	= blue_a4
            stroke_thickness	= 0.5
            file			= circos_input_files/prot_links
            radius			= 0.9r
            bezier_radius	= 0.3r
            thickness		= 1
        </link>
        <link>
            ribbon			= yes
            color			= black_a5
            stroke			= yes
            stroke_color	= black_a5
            stroke_thickness	= 0.5
            file			= circos_input_files/duplication_links
            radius			= 0.9r
            bezier_radius	= 0.3r
            thickness		= 1
        </link>
    </links>
    
    <highlights>
        <highlight>
            file = circos_input_files/plus_prot_names_highlights
            r0 = 1.0r-30p
            r1 = 1.0r-100p
            fill_color = blue
            stroke_color = blue
            stroke_thickness = 2
        </highlight>
        <highlight>
            file = circos_input_files/plus_tRNA_names_highlights
            r0 = 1.0r-30p
            r1 = 1.0r-100p
            fill_color = red
            stroke_color = red
            stroke_thickness = 2
        </highlight>
        <highlight>
            file = circos_input_files/plus_rRNA_names_highlights
            r0 = 1.0r-30p
            r1 = 1.0r-100p
            fill_color = green
            stroke_color = green
            stroke_thickness = 2
        </highlight>
        <highlight>
            file = circos_input_files/minus_prot_names_highlights
            r0 = 1.0r
            r1 = 1.0r+70p
            fill_color = blue
            stroke_color = blue
            stroke_thickness = 2
        </highlight>
        <highlight>
            file = circos_input_files/minus_tRNA_names_highlights
            r0 = 1.0r
            r1 = 1.0r+70p
            fill_color = red
            stroke_color = red
            stroke_thickness = 2
        </highlight>
        <highlight>
            file = circos_input_files/minus_rRNA_names_highlights
            r0 = 1.0r
            r1 = 1.0r+70p
            fill_color = green
            stroke_color = green
            stroke_thickness = 2
        </highlight>
    </highlights>
    
    <plots>
        <plot>
            type = text
            file = circos_input_files/minus_names
            color = black
            r0 = 1.0r+70p
            r1 = 1.0r+400p
            
            show_links = yes
            link_dims = 10p,50p,200p,50p,10p
            link_thickness = 2p
            link_color = black
            
            label_snuggle = yes
            max_snuggle_distance = 1r
            
            label_size = 20p
            label_font = condensed
            padding = 0p
            rpadding = 0p
        </plot>
        <plot>
            type = text
            file = circos_input_files/plus_names
            color = black
            r0 = 1.0r+0p
            r1 = 1.0r+550p
            
            show_links = yes
            link_dims = 10p,50p,50p,50p,10p
            link_thickness = 2p
            link_color = black
            
            label_snuggle = yes
            max_snuggle_distance = 1r
            
            label_size = 20p
            label_font = condensed
            padding = 0p
            rpadding = 0p
        </plot>
    </plots>
    """

    with open(os.path.join(config['output_dir'], 'circos.config'), 'w') as file:
        file.write(template)

