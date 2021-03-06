#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse

from ChromoDoter.Parsers.PSL import CollectionPSL
from ChromoDoter.Parsers.LAST import CollectionLast
from ChromoDoter.Collections.General import IdList, SynDict
from ChromoDoter.Routines import CircosRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_synteny_psl", action="store", dest="input_synteny_psl", required=True,
                    help="Input psl file created by halSyntheny")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file with links for circos")

parser.add_argument("-w", "--white_target_id_file", action="store", dest="white_target_id_file",
                    help="File with target scaffold ids from white list or corresponding comma-separated list."
                         "NOTE: filtering is done BEFORE renaming by synonyms!")
parser.add_argument("-b", "--black_target_id_file", action="store", dest="black_target_id_file",
                    help="File with target scaffold ids from black list or corresponding comma-separated list."
                         "NOTE: filtering is done BEFORE renaming by synonyms!")

parser.add_argument("-x", "--white_query_id_file", action="store", dest="white_query_id_file",
                    help="File with query scaffold ids from white list or corresponding comma-separated list."
                         "NOTE: filtering is done BEFORE renaming by synonyms!")
parser.add_argument("-c", "--black_query_id_file", action="store", dest="black_query_id_file",
                    help="File with query scaffold ids from black list or corresponding comma-separated list."
                         "NOTE: filtering is done BEFORE renaming by synonyms!")

parser.add_argument("-s", "--query_syn_file", action="store", dest="query_syn_file",
                    help="File with query scaffold id synonyms")
parser.add_argument("--query_syn_file_key_column", action="store", dest="query_syn_file_key_column",
                    default=0, type=int,
                    help="Column(0-based) with key(current id) for scaffolds in query synonym file")
parser.add_argument("--query_syn_file_value_column", action="store", dest="query_syn_file_value_column",
                    default=1, type=int,
                    help="Column(0-based) with value(synonym id) for scaffolds in query synonym file synonym")

parser.add_argument("-y", "--target_syn_file", action="store", dest="target_syn_file",
                    help="File with target scaffold id synonyms")
parser.add_argument("--target_syn_file_key_column", action="store", dest="target_syn_file_key_column",
                    default=0, type=int,
                    help="Column(0-based) with key(current id) for scaffolds in target synonym file")
parser.add_argument("--target_syn_file_value_column", action="store", dest="target_syn_file_value_column",
                    default=1, type=int,
                    help="Column(0-based) with value(synonym id) for scaffolds in target synonym file synonym")
"""
parser.add_argument("-u", "--target_order_file", action="store", dest="target_order_file",
                    help="File with order of target scaffolds on the plot or corresponding comma-separated list."
                         "NOTE: ordering is done AFTER renaming by synonyms!")
parser.add_argument("-z", "--query_order_file", action="store", dest="query_order_file",
                    help="File with order of query scaffolds on the plot or corresponding comma-separated list."
                         "NOTE: ordering is done AFTER renaming by synonyms!")

parser.add_argument("-l", "--target_label", action="store", dest="target_label",
                    help="Label for target genome(X axis)")
parser.add_argument("-r", "--query_label", action="store", dest="query_label",
                    help="Label for query genome(Y axis)")
parser.add_argument("-t", "--title", action="store", dest="title",
                    help="Title of dot plot")
parser.add_argument("-e", "--extensions", action="store", dest="extensions", type=lambda x: x.split(","),
                    default=["png", ],
                    help="Comma-separated list of extensions for histogram files. Default: png")

parser.add_argument("-d", "--dpi", action="store", dest="dpi", type=int, default=400,
                    help="DPI of figure. Default: 400")
parser.add_argument("-f", "--figsize", action="store", dest="figsize", type=lambda s: list(map(int, s.split(","))),
                    default=(12, 12),
                    help="Size of figure in inches(two comma-separated ints). "
                         "Ignored if --auto_scale_figure is set. Default: 12,12")
parser.add_argument("--auto_scale_figure", action="store_true", dest="auto_scale_figure", default=False,
                    help="Autoscale figure size based on total length of scaffolds (independently for query and target)"
                         ". Depends on --mbp_per_inch/--figure_height/--figure_width options. Default: False")
parser.add_argument("--mbp_per_inch", action="store", dest="mbp_per_inch", type=float, default=None,
                    help="Megabasepairs per inch. This option is used for autoscaling of figure size. Default: not set")
parser.add_argument("--figure_width", action="store", dest="figure_width", type=int, default=None,
                    help="Figure width. This option is used for autoscaling of figure size."
                         " Ignored if --mbp_per_inch is set. Default: not set")
parser.add_argument("--figure_height", action="store", dest="figure_height", type=int, default=None,
                    help="Figure height. This option is used for autoscaling of figure size."
                         " Ignored if --mbp_per_inch or --figure_width is set.  Default: not set")
parser.add_argument("-a", "--antialiasing", action="store_true", dest="antialiasing", default=False,
                    help="Enable antialiasing. Use this option only for small sequences, i.e segments of chromosomes ")

parser.add_argument("--linewidth", action="store", dest="linewidth", type=float,
                    default=0.01,
                    help="Width of alignment lines. Default: 0.01")
parser.add_argument("--gridwidth", action="store", dest="gridwidth", type=float,
                    default=1,
                    help="Width of grid lines. Default: 1")
parser.add_argument("--scaffold_label_fontsize", action="store", dest="scaffold_label_fontsize", type=float,
                    default=13,
                    help="Fontsize for scaffold labels. Default: 13")
parser.add_argument("--grid_color", action="store", dest="grid_color", default='black',
                    help="Color of grid lines. Default: 'black'")
parser.add_argument("--hide_grid", action="store_true", dest="hide_grid", default=False,
                    help="Hide grid. Default: False")
parser.add_argument("--bar_color", action="store", dest="bar_color", default='grey',
                    help="Color of bars. Default: 'grey'")
parser.add_argument("--same_strand_color", action="store", dest="same_strand_color", default='blue',
                    help="Color of alignment line if query and target are in same strand. Default: 'blue'")
parser.add_argument("--diff_strand_color", action="store", dest="diff_strand_color", default='red',
                    help="Color of alignment line if query and target are in different strand. Default: 'red'")

parser.add_argument("--target_scaffold_labels_angle", action="store", dest="target_scaffold_labels_angle",
                    default=45, type=int,
                    help="Angle for labels of target scaffolds. Default: 45")
parser.add_argument("--query_scaffold_labels_angle", action="store", dest="query_scaffold_labels_angle",
                    default=0, type=int,
                    help="Angle for labels of query scaffolds. Default: 0")
parser.add_argument("--hide_target_labels", action="store_true", dest="hide_target_labels", default=False,
                    help="Hide labels of target scaffolds. Default: False")
parser.add_argument("--hide_query_labels", action="store_true", dest="hide_query_labels", default=False,
                    help="Hide labels of query scaffolds. Default: False")
parser.add_argument("--bottom_offset", action="store", dest="bottom_offset",
                    default=0.1, type=float,
                    help="Bottom offset for subplot. Default: 0.1")
parser.add_argument("--top_offset", action="store", dest="top_offset",
                    default=0.9, type=float,
                    help="Top offset for subplot. Default: 0.9")
parser.add_argument("--left_offset", action="store", dest="left_offset",
                    default=0.1, type=float,
                    help="Left offset for subplot. Default: 0.1")
parser.add_argument("--right_offset", action="store", dest="right_offset",
                    default=0.9, type=float,
                    help="Right offset for subplot. Default: 0.9")
parser.add_argument("--x_axis_visible", action="store_true", dest="x_axis_visible",
                    default=False,
                    help="Make X axis visible. Default: False")
parser.add_argument("--y_axis_visible", action="store_true", dest="y_axis_visible",
                    default=False,
                    help="Make Y axis visible. Default: False")
parser.add_argument("--axes_label_distance", action="store", dest="axes_label_distance",
                    default=12, type=float,
                    help="Distance between axes and its labels. Default: 12")
parser.add_argument("--remove_scaffolds_absent_in_target_ordered_list", action="store_true",
                    dest="remove_scaffolds_absent_in_target_ordered_list",
                    default=False,
                    help="Remove scaffolds that are absent in target orderlist. Default: False")
parser.add_argument("--remove_scaffolds_absent_in_query_ordered_list", action="store_true",
                    dest="remove_scaffolds_absent_in_query_ordered_list",
                    default=False,
                    help="Remove scaffolds that are absent in scaffold orderlist. Default: False")
parser.add_argument("--show_length_ticks", action="store_true",
                    dest="show_length_ticks",
                    default=False,
                    help="Show scaffold length ticks. Default: False")
parser.add_argument("--show_tick_grid", action="store_true",
                    dest="show_tick_grid",
                    default=False,
                    help="Show tick grid. Default: False")
parser.add_argument("--tick_step", action="store", dest="tick_step",
                    default=10000000, type=int,
                    help="Scaffold length tick step. Default: 10 000 000, i.e 10 Mbp")
parser.add_argument("--tick_unit", action="store", dest="tick_unit",
                    default=1000000, type=int,
                    help="Scaffold length tick unit. Default: 1 000 000, i.e. Mbp")
"""
args = parser.parse_args()

if args.white_target_id_file:
    target_white_list = IdList(filename=args.white_target_id_file) if os.path.isfile(args.white_target_id_file) else IdList(args.white_target_id_file.split(","))
else:
    target_white_list = IdList()

if args.black_target_id_file:
    target_black_list = IdList(filename=args.black_target_id_file) if os.path.isfile(args.black_target_id_file) else IdList(args.black_target_id_file.split(","))
else:
    target_black_list = IdList()

if args.white_query_id_file:
    query_white_list = IdList(filename=args.white_query_id_file) if os.path.isfile(args.white_query_id_file) else IdList(args.white_query_id_file.split(","))
    #print query_white_list
else:
    query_white_list = IdList()

if args.black_query_id_file:
    query_black_list = IdList(filename=args.black_query_id_file) if os.path.isfile(args.black_query_id_file) else IdList(args.black_query_id_file.split(","))
else:
    query_black_list = IdList()

query_syn_dict = SynDict(filename=args.query_syn_file,
                         key_index=args.query_syn_file_key_column,
                         value_index=args.query_syn_file_value_column)
target_syn_dict = SynDict(filename=args.target_syn_file,
                          key_index=args.target_syn_file_key_column,
                          value_index=args.target_syn_file_value_column)

psl_collection = CollectionPSL(args.input_synteny_psl,
                               target_white_list=target_white_list,
                               target_black_list=target_black_list,
                               query_white_list=query_white_list,
                               query_black_list=query_black_list,
                               query_syn_dict=query_syn_dict,
                               target_syn_dict=target_syn_dict,
                               format="psl")

circos_links_df = CircosRoutines.convert_synteny_psl_to_circos_links(psl_collection)
circos_links_df.to_csv(args.output, sep=" ", index=False, header=False)
