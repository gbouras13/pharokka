from pycirclize import Circos
from pycirclize.parser import Gff
from pycirclize.parser import Genbank
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import os
import numpy as np




# Load GFF file

def create_plot(out_dir, prefix,  interval, annotations, title_size, plot_title, truncate, outfile, dpi, label_size):
    # read in gff file
    gff_file =  os.path.join(out_dir, prefix + ".gff")
    gff = Gff(gff_file)

    # Load Genbank file
    gbk_file = os.path.join(out_dir, prefix + ".gbk")
    # get only to range of gff - as by default gbk takes all contigs, gff only the first
    gbk = Genbank(gbk_file, max_range = gff.range_size)

    # instantiate circos
    circos = Circos(sectors={gbk.name: gbk.range_size})

    # title if not blank

    circos.text(plot_title, size=int(title_size), r=190)

    sector = circos.get_sector(gbk.name)
    cds_track = sector.add_track((70, 80))
    cds_track.axis(fc="#EEEEEE", ec="none")

    # Plot forward CDS
    # cds_track.genomic_features(
    #     gff.extract_features("CDS", target_strand=1),
    #     plotstyle="arrow",
    #     r_lim=(75, 80),
    #     fc="salmon",
    # )


#### plot each PHROG fwd and reverse

## colours 

#4deeea              
#74ee15
#ffe700
#f000ff
#001eff
#8900ff
#ff008d
#1c1c1c
#ffffff    


# unknown 

    unk_col = '#4deeea'

    unk_fwd_list = []
    for f in gff.extract_features("CDS", target_strand=1):
        if f.qualifiers.get("function")[0] == 'unknown function':
            unk_fwd_list.append(f)

    cds_track.genomic_features(
        unk_fwd_list,
        plotstyle="arrow",
        r_lim=(75, 80),
        fc=unk_col,
    )

    unk_rev_list = []
    for f in gff.extract_features("CDS", target_strand=-1):
        if f.qualifiers.get("function")[0] == 'unknown function':
            unk_rev_list.append(f)

    cds_track.genomic_features(
        unk_rev_list,
        plotstyle="arrow",
        r_lim=(70, 75),
        fc=unk_col,
    )

# other

    other_col = '#4deeea'

    fwd_list = []
    for f in gff.extract_features("CDS", target_strand=1):
        if f.qualifiers.get("function")[0] == 'other':
            fwd_list.append(f)

    cds_track.genomic_features(
        fwd_list,
        plotstyle="arrow",
        r_lim=(75, 80),
        fc=other_col,
    )

    rev_list = []
    for f in gff.extract_features("CDS", target_strand=-1):
        if f.qualifiers.get("function")[0] == 'other':
            rev_list.append(f)

    cds_track.genomic_features(
        rev_list,
        plotstyle="arrow",
        r_lim=(70, 75),
        fc=other_col,
    )

# tail 

    tail_col = '#74ee15'

    fwd_list = []
    for f in gff.extract_features("CDS", target_strand=1):
        if f.qualifiers.get("function")[0] == 'tail':
            fwd_list.append(f)

    cds_track.genomic_features(
        fwd_list,
        plotstyle="arrow",
        r_lim=(75, 80),
        fc=tail_col,
    )

    rev_list = []
    for f in gff.extract_features("CDS", target_strand=-1):
        if f.qualifiers.get("function")[0] == 'tail':
            rev_list.append(f)

    cds_track.genomic_features(
        rev_list,
        plotstyle="arrow",
        r_lim=(70, 75),
        fc=tail_col,
    )

# transcription regulation

    transcription_col = '#ffe700'

    fwd_list = []
    for f in gff.extract_features("CDS", target_strand=1):
        if f.qualifiers.get("function")[0] == 'transcription regulation':
            fwd_list.append(f)

    cds_track.genomic_features(
        fwd_list,
        plotstyle="arrow",
        r_lim=(75, 80),
        fc=transcription_col,
    )

    rev_list = []
    for f in gff.extract_features("CDS", target_strand=-1):
        if f.qualifiers.get("function")[0] == 'transcription regulation':
            rev_list.append(f)

    cds_track.genomic_features(
        rev_list,
        plotstyle="arrow",
        r_lim=(70, 75),
        fc=transcription_col,
    )

# DNA, RNA and nucleotide metabolism

    dna_col = '#f000ff'

    fwd_list = []
    for f in gff.extract_features("CDS", target_strand=1):
        if f.qualifiers.get("function")[0] == 'DNA':
            fwd_list.append(f)

    cds_track.genomic_features(
        fwd_list,
        plotstyle="arrow",
        r_lim=(75, 80),
        fc=dna_col,
    )

    rev_list = []
    for f in gff.extract_features("CDS", target_strand=-1):
        if f.qualifiers.get("function")[0] == 'DNA':
            rev_list.append(f)

    cds_track.genomic_features(
        rev_list,
        plotstyle="arrow",
        r_lim=(70, 75),
        fc=dna_col,
    )

#lysis

    lysis_col = '#001eff'

    fwd_list = []
    for f in gff.extract_features("CDS", target_strand=1):
        if f.qualifiers.get("function")[0] == 'lysis':
            fwd_list.append(f)

    cds_track.genomic_features(
        fwd_list,
        plotstyle="arrow",
        r_lim=(75, 80),
        fc=lysis_col,
    )

    rev_list = []
    for f in gff.extract_features("CDS", target_strand=-1):
        if f.qualifiers.get("function")[0] == 'lysis':
            rev_list.append(f)

    cds_track.genomic_features(
        rev_list,
        plotstyle="arrow",
        r_lim=(70, 75),
        fc=lysis_col,
    )

# moron, auxiliary metabolic gene and host takeover

    moron_col = '#8900ff'

    fwd_list = []
    for f in gff.extract_features("CDS", target_strand=1):
        if f.qualifiers.get("function")[0] == 'moron':
            fwd_list.append(f)

    cds_track.genomic_features(
        fwd_list,
        plotstyle="arrow",
        r_lim=(75, 80),
        fc=moron_col,
    )

    rev_list = []
    for f in gff.extract_features("CDS", target_strand=-1):
        if f.qualifiers.get("function")[0] == 'moron':
            rev_list.append(f)

    cds_track.genomic_features(
        rev_list,
        plotstyle="arrow",
        r_lim=(70, 75),
        fc=moron_col,
    )




# integration and excision

    int_col = '#1c1c1c'

    fwd_list = []
    for f in gff.extract_features("CDS", target_strand=1):
        if f.qualifiers.get("function")[0] == 'integration and excision':
            fwd_list.append(f)

    cds_track.genomic_features(
        fwd_list,
        plotstyle="arrow",
        r_lim=(75, 80),
        fc=other_col,
    )

    rev_list = []
    for f in gff.extract_features("CDS", target_strand=-1):
        if f.qualifiers.get("function")[0] == 'integration and excision':
            rev_list.append(f)

    cds_track.genomic_features(
        rev_list,
        plotstyle="arrow",
        r_lim=(70, 75),
        fc=int_col,
    )


# head and packaging

    head_col = '#ff008d'

    fwd_list = []
    for f in gff.extract_features("CDS", target_strand=1):
        if f.qualifiers.get("function")[0] == 'head and packaging':
            fwd_list.append(f)

    cds_track.genomic_features(
        fwd_list,
        plotstyle="arrow",
        r_lim=(75, 80),
        fc=head_col,
    )

    rev_list = []
    for f in gff.extract_features("CDS", target_strand=-1):
        if f.qualifiers.get("function")[0] == 'head and packaging':
            rev_list.append(f)

    cds_track.genomic_features(
        rev_list,
        plotstyle="arrow",
        r_lim=(70, 75),
        fc=head_col,
    )

# connector

    con_col = '#808080'

    fwd_list = []
    for f in gff.extract_features("CDS", target_strand=1):
        if f.qualifiers.get("function")[0] == 'connector':
            fwd_list.append(f)

    cds_track.genomic_features(
        fwd_list,
        plotstyle="arrow",
        r_lim=(75, 80),
        fc=con_col,
    )

    rev_list = []
    for f in gff.extract_features("CDS", target_strand=-1):
        if f.qualifiers.get("function")[0] == 'connector':
            rev_list.append(f)

    cds_track.genomic_features(
        rev_list,
        plotstyle="arrow",
        r_lim=(70, 75),
        fc=con_col,
    )




    # trunction 
    truncate = int(truncate)


    # Extract CDS product labels
    pos_list, labels, length_list = [], [], []
    for f in gff.extract_features("CDS"):
        start, end = int(str(f.location.end)), int(str(f.location.start))
        pos = (start + end) / 2.
        length = end - start
        label = f.qualifiers.get("product", [""])[0]
        if label == "" or label.startswith("hypothetical") or label.startswith("unknown") :
            continue
        if len(label) > truncate:
            label = label[:truncate] + "..."
        pos_list.append(pos)
        labels.append(label)
        length_list.append(length)




    # Plot forward CDS

    # Plot reverse CDS
    # cds_track.genomic_features(
    #     gff.extract_features("CDS", target_strand=-1),
    #     plotstyle="arrow",
    #     r_lim=(70, 75),
    #     fc="skyblue",
    # )


    # cds_track.genomic_features(
    #     gff.extract_features("CDS", target_strand=1),
    #     plotstyle="arrow",
    #     r_lim=(75, 80),
    #     fc="salmon",
    # )


    #### thin out annotations
    annotations = float(annotations)


    # if sparse > 1:
    #     pos_list = pos_list[::sparse]
    #     labels = labels[::sparse]

    if annotations == 0:
        print("by inputting --annotations 0 you have chosen to plot no annotations. Continuing.")
    elif annotations == 0:
        print("by inputting --annotations 1 you have chosen to plot all annotations. Continuing.")
    elif annotations > 1:
        print("You have input a --annotations value greater than 1. Setting to 1 (will plot all annotations). Continuing.")
        annotations = 1
    elif annotations < 0: # sparse < 0
        print("You have input a --annotations value less than 1. Setting to 0 (will plot no annotations). Continuing.")
        annotations = 0

    ####### running the sparsity

    median_length = np.quantile(length_list, annotations)
    # Create an empty list to store the filtered indices
    filtered_indices = []

    # Loop through the indices of the length_list
    for i in range(len(length_list)):
        # If the length at this index is greater than or equal to the median, add the index to filtered_indices
        if length_list[i] < median_length:
            filtered_indices.append(i)

    # Use the filtered indices to create new lists for pos_list, labels, and length_list
    pos_list = [pos_list[i] for i in filtered_indices]
    labels = [labels[i] for i in filtered_indices]
    length_list = [length_list[i] for i in filtered_indices]



    # Plot CDS product labels on outer position
    cds_track.xticks(
        pos_list,
        labels,
        label_orientation="vertical",
        show_bottom_line=True,
        label_size=label_size,
        line_kws=dict(ec="grey"),
    )

        # Plot GC content
    gc_content_track = sector.add_track((42.5, 60))

    pos_list, gc_contents = gbk.calc_gc_content()
    gc_contents = gc_contents - gbk.calc_genome_gc_content()
    positive_gc_contents = np.where(gc_contents > 0, gc_contents, 0)
    negative_gc_contents = np.where(gc_contents < 0, gc_contents, 0)
    abs_max_gc_content = np.max(np.abs(gc_contents))
    vmin, vmax = -abs_max_gc_content, abs_max_gc_content
    gc_content_track.fill_between(
        pos_list, positive_gc_contents, 0, vmin=vmin, vmax=vmax, color="black"
    )
    gc_content_track.fill_between(
        pos_list, negative_gc_contents, 0, vmin=vmin, vmax=vmax, color="grey"
    )

    # Plot GC skew
    gc_skew_track = sector.add_track((25, 42.5))

    pos_list, gc_skews = gbk.calc_gc_skew()
    positive_gc_skews = np.where(gc_skews > 0, gc_skews, 0)
    negative_gc_skews = np.where(gc_skews < 0, gc_skews, 0)
    abs_max_gc_skew = np.max(np.abs(gc_skews))
    vmin, vmax = -abs_max_gc_skew, abs_max_gc_skew
    gc_skew_track.fill_between(
        pos_list, positive_gc_skews, 0, vmin=vmin, vmax=vmax, color="green"
    )
    gc_skew_track.fill_between(
        pos_list, negative_gc_skews, 0, vmin=vmin, vmax=vmax, color="purple"
    )

    label_size = int(label_size)

    # Plot xticks & intervals on inner position
    cds_track.xticks_by_interval(
        interval=int(interval),
        outer=False,
        show_bottom_line=False,
        label_formatter=lambda v: f"{v/ 1000:.0f} Kb", # no decimal place
        label_orientation="vertical",
        line_kws=dict(ec="grey"),
        label_size=8
    )

    # # Add legend
    handle_phrogs = [
        Patch(color=unk_col, label="Unknown/Other Function"),
        Patch(color=tail_col, label="Tail"),
        Patch(color=transcription_col, label="Transcription Regulation"),
        Patch(color=dna_col, label="DNA/RNA & nucleotide \n metabolism"),
        Patch(color=lysis_col, label="Lysis"),
        Patch(color=moron_col, label="Moron, auxiliary metabolic \n gene & host takeover"),
        Patch(color=int_col, label="Integration & excision"),
        Patch(color=head_col, label="Head & packaging"),
        Patch(color=con_col, label="Connector")
    ]




    fig = circos.plotfig()

      
    line_legend = circos.ax.legend(
    handles=handle_phrogs,
    bbox_to_anchor=(0.92, 1.2),
    fontsize=9.5,
    loc="center",
    title="PHROG CDS",
    handlelength=2,
)
    circos.ax.add_artist(line_legend)

    handle_gc_content = [
        Line2D([], [], color="black", label="Positive GC Content", marker="^", ms=6, ls="None"),
        Line2D([], [], color="grey", label="Negative GC Content", marker="v", ms=6, ls="None"),
    ]

    handle_gc_skew = [
        Line2D([], [], color="green", label="Positive GC Skew", marker="^", ms=6, ls="None"),
        Line2D([], [], color="purple", label="Negative GC Skew", marker="v", ms=6, ls="None")
    ]

    # shrink plot a bit (0.8)
    box = circos.ax.get_position()
    circos.ax.set_position([box.x0, box.y0, box.width * 0.65, box.height*0.9])

    gc_legend_cont = circos.ax.legend(
    handles=handle_gc_content,
    bbox_to_anchor=(0.08, 1.25),
    loc="center",
    fontsize=9.5,
    title="GC Content (Outer Ring)",
    handlelength=2,
    )

    circos.ax.add_artist(gc_legend_cont)

    gc_legend_skew = circos.ax.legend(
    handles=handle_gc_skew,
    bbox_to_anchor=(0.08, 1.15),
    loc="center",
    fontsize=9.5,
    title="GC Skew (Inner Ring)",
    handlelength=2,
    )
    
    circos.ax.add_artist(gc_legend_skew)

    dpi = int(dpi)
    
    fig.savefig(outfile, dpi=dpi)