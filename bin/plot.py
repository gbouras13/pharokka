from pycirclize import Circos
from pycirclize.parser import Gff
from pycirclize.parser import Genbank
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import os
import numpy as np



# Load GFF file

def create_plot(out_dir, prefix, plot_name, interval, sparse, title_size):
    # read in gff file
    gff_file =  os.path.join(out_dir, prefix + ".gff")
    gff = Gff(gff_file)

    # Load Genbank file
    gbk_file = os.path.join(out_dir, prefix + ".gbk")
    # get only to range of gff - as by default gbk takes all contigs, gff only the first
    gbk = Genbank(gbk_file, max_range = gff.range_size)

    # instantiate circos
    circos = Circos(sectors={gbk.name: gbk.range_size})
    # remove text name as cant fix overlap for now
    circos.text(plot_name, size=int(title_size), r=190)

    sector = circos.get_sector(gbk.name)
    cds_track = sector.add_track((90, 100))
    cds_track.axis(fc="#EEEEEE", ec="none")
    # Plot forward CDS
    cds_track.genomic_features(
        gff.extract_features("CDS", target_strand=1),
        plotstyle="arrow",
        r_lim=(95, 100),
        fc="salmon",
    )
    # Plot reverse CDS
    cds_track.genomic_features(
        gff.extract_features("CDS", target_strand=-1),
        plotstyle="arrow",
        r_lim=(90, 95),
        fc="skyblue",
    )
    # Extract CDS product labels
    pos_list, labels = [], []
    for f in gff.extract_features("CDS"):
        start, end = int(str(f.location.end)), int(str(f.location.start))
        pos = (start + end) / 2
        label = f.qualifiers.get("product", [""])[0]
        if label == "" or label.startswith("hypothetical") or label.startswith("unknown") :
            continue
        if len(label) > 25:
            label = label[:25] + "..."
        pos_list.append(pos)
        labels.append(label)

    #### thin out annotations

    sparse = int(sparse)

    if sparse > 1:
        pos_list = pos_list[::sparse]
        labels = labels[::sparse]

    # Plot CDS product labels on outer position
    cds_track.xticks(
        pos_list,
        labels,
        label_orientation="vertical",
        show_bottom_line=True,
        label_size=7,
        line_kws=dict(ec="grey"),
    )

        # Plot GC content
    gc_content_track = sector.add_track((50, 65))

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
    gc_skew_track = sector.add_track((35, 50))

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

    # Plot xticks & intervals on inner position
    cds_track.xticks_by_interval(
        interval=int(interval),
        outer=False,
        show_bottom_line=False,
        label_formatter=lambda v: f"{v/ 1000:.1f} Kb",
        label_orientation="vertical",
        line_kws=dict(ec="grey"),
    )

    # # Add legend
    handle = [
        Patch(color="salmon", label="Forward CDS"),
        Patch(color="skyblue", label="Reverse CDS"),
        Line2D([], [], color="black", label="Positive GC Content", marker="^", ms=6, ls="None"),
        Line2D([], [], color="grey", label="Negative GC Content", marker="v", ms=6, ls="None"),
        Line2D([], [], color="green", label="Positive GC Skew", marker="^", ms=6, ls="None"),
        Line2D([], [], color="purple", label="Negative GC Skew", marker="v", ms=6, ls="None")
    ]

    fig = circos.plotfig()

    _ = circos.ax.legend(handles=handle, 
                         bbox_to_anchor=(1.3, 1.3),  
                         loc="center",
                         fontsize=9)
    
    fig.savefig(os.path.join(out_dir,  str(prefix) + "_plot.png"), dpi=600)