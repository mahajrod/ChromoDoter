__author__ = 'mahajrod'

import datetime
import numpy as np
#import matplotlib

import pandas as pd

#matplotlib.use('Agg')

import matplotlib.pyplot as plt
plt.ioff()
from matplotlib.patches import Rectangle, Circle
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection


from ChromoDoter.Parsers.LAST import CollectionLast
from ChromoDoter.Parsers.PSL import CollectionPSL


class DrawingRoutines:

    @staticmethod
    def millions(x, pos):
        return '%1.1fMbp' % (x*1e-6)

    @staticmethod
    def billions(x, pos):
        return '%1.1fGbp' % (x*1e-9)

    @staticmethod
    def get_filtered_scaffold_list(raw_scaffold_list,
                                   scaffold_black_list=[],
                                   sort_scaffolds=False,
                                   scaffold_ordered_list=None,
                                   scaffold_white_list=[],
                                   scaffold_length_dict=None,
                                   length_column_index=0,
                                   min_length=None,
                                   remove_scaffolds_absent_in_ordered_list=False):

        white_set = set(scaffold_white_list)
        black_set = set(scaffold_black_list)
        length_set = None

        scaffold_set = set(raw_scaffold_list)

        if scaffold_length_dict and min_length:
            if isinstance(scaffold_length_dict, pd.DataFrame):
                length_df = scaffold_length_dict
            else:
                length_df = pd.DataFrame.from_dict(scaffold_length_dict, orient="index")
            length_set = set(length_df[length_df.iloc[:, length_column_index] >= min_length].index)

        if white_set:
            scaffold_set = scaffold_set & white_set

        if black_set:
            scaffold_set = scaffold_set - black_set

        if length_set is not None:
            scaffold_set = scaffold_set & length_set

        scaffold_list = list(scaffold_set)

        if sort_scaffolds:
            scaffold_list.sort()

        final_scaffold_list = []

        if scaffold_ordered_list:
            for entry in scaffold_ordered_list:
                if entry in scaffold_list:
                    final_scaffold_list.append(entry)
                    scaffold_list.remove(entry)
                else:
                    print("WARNING!!!Entry(%s) from order list is absent in list of scaffolds!" % entry)
            if not remove_scaffolds_absent_in_ordered_list:
                final_scaffold_list = final_scaffold_list + scaffold_list
        else:
            final_scaffold_list = scaffold_list

        return final_scaffold_list

    def draw_dot_plot_per_scaffold_from_last_self_alignment(self, last_collection,
                                                            output_prefix=None,
                                                            extension_list=("png", ),
                                                            scaffold_black_list=(), scaffold_white_list=(),
                                                            scaffold_reverse_list=(),
                                                            figsize=(16, 16), dpi=300,
                                                            grid_color='black',
                                                            bar_color='grey',
                                                            same_strand_color='red',
                                                            diff_strand_color='blue',
                                                            title=None,
                                                            label=None,
                                                            gridwidth=1,
                                                            show_grid=True,
                                                            linewidth=0.01,
                                                            scaffold_label_fontsize=13,
                                                            axes_label_fontstyle="italic",
                                                            axes_label_weight="normal",
                                                            axes_label_fontsize=15,
                                                            axes_label_distance=6,
                                                            antialiased_lines=None,
                                                            target_scaffold_labels_angle=45,
                                                            query_scaffold_labels_angle=0,
                                                            show_labels=True,
                                                            bottom_offset=0.1,
                                                            top_offset=0.9,
                                                            left_offset=0.1,
                                                            right_offset=0.9,
                                                            x_axis_visible=False,
                                                            y_axis_visible=False,
                                                            x_labelpad=None,
                                                            y_labelpad=None
                                                            ):

        for scaffold in scaffold_white_list:
            self.draw_dot_plot_from_last_alignment(last_collection,
                                                  output_prefix="%s.%s" % (output_prefix, scaffold),
                                                  extension_list=extension_list,
                                                  target_black_list=scaffold_black_list, target_white_list=[scaffold,],
                                                  target_ordered_list=[scaffold,], target_reverse_list=scaffold_reverse_list,
                                                  query_black_list=scaffold_black_list, query_white_list=[scaffold,],
                                                  query_ordered_list=[scaffold,], query_reverse_list=scaffold_reverse_list,
                                                  figsize=figsize, dpi=dpi,
                                                  grid_color=grid_color,
                                                  bar_color=bar_color,
                                                  same_strand_color=same_strand_color,
                                                  diff_strand_color=diff_strand_color,
                                                  title=title,
                                                  target_label=label,
                                                  query_label=label,
                                                  gridwidth=gridwidth,
                                                  show_grid=show_grid,
                                                  linewidth=linewidth,
                                                  scaffold_label_fontsize=scaffold_label_fontsize,
                                                  axes_label_fontstyle=axes_label_fontstyle,
                                                  axes_label_weight=axes_label_weight,
                                                  axes_label_fontsize=axes_label_fontsize,
                                                  axes_label_distance=axes_label_distance,
                                                  antialiased_lines=antialiased_lines,
                                                  target_scaffold_labels_angle=target_scaffold_labels_angle,
                                                  query_scaffold_labels_angle=query_scaffold_labels_angle,
                                                  show_target_labels=show_labels,
                                                  show_query_labels=show_labels,
                                                  bottom_offset=bottom_offset,
                                                  top_offset=top_offset,
                                                  left_offset=left_offset,
                                                  right_offset=right_offset,
                                                  x_axis_visible=x_axis_visible,
                                                  y_axis_visible=y_axis_visible,
                                                  x_labelpad=x_labelpad,
                                                  y_labelpad=y_labelpad
                                                   )

    def draw_dot_plot_per_scaffold_vs_all_from_last_self_alignment(self, last_collection,
                                                            output_prefix=None,
                                                            extension_list=("png", ),
                                                            scaffold_black_list=(), scaffold_white_list=(),
                                                            scaffold_reverse_list=(),
                                                            ordered_list=(),
                                                            figsize=(16, 16), dpi=300,
                                                            grid_color='black',
                                                            bar_color='grey',
                                                            same_strand_color='red',
                                                            diff_strand_color='blue',
                                                            title=None,
                                                            label=None,
                                                            gridwidth=1,
                                                            show_grid=True,
                                                            linewidth=0.01,
                                                            scaffold_label_fontsize=13,
                                                            axes_label_fontstyle="italic",
                                                            axes_label_weight="normal",
                                                            axes_label_fontsize=15,
                                                            axes_label_distance=6,
                                                            antialiased_lines=None,
                                                            target_scaffold_labels_angle=45,
                                                            query_scaffold_labels_angle=0,
                                                            show_labels=True,
                                                            bottom_offset=0.1,
                                                            top_offset=0.9,
                                                            left_offset=0.1,
                                                            right_offset=0.9,
                                                            x_axis_visible=False,
                                                            y_axis_visible=False,
                                                            x_labelpad=None,
                                                            y_labelpad=None
                                                           ):

        for scaffold in scaffold_white_list:
            self.draw_dot_plot_from_last_alignment(last_collection,
                                                   output_prefix="%s.%s" % (output_prefix, scaffold),
                                                   extension_list=extension_list,
                                                   target_black_list=scaffold_black_list, target_white_list=scaffold_white_list,
                                                   target_ordered_list=ordered_list, target_reverse_list=scaffold_reverse_list,
                                                   query_black_list=scaffold_black_list, query_white_list=[scaffold,],
                                                   query_ordered_list=[scaffold,], query_reverse_list=scaffold_reverse_list,
                                                   figsize=figsize, dpi=dpi,
                                                   grid_color=grid_color,
                                                   bar_color=bar_color,
                                                   same_strand_color=same_strand_color,
                                                   diff_strand_color=diff_strand_color,
                                                   title=title,
                                                   target_label=label,
                                                   query_label=label,
                                                   gridwidth=gridwidth,
                                                   show_grid=show_grid,
                                                   linewidth=linewidth,
                                                   scaffold_label_fontsize=scaffold_label_fontsize,
                                                   axes_label_fontstyle=axes_label_fontstyle,
                                                   axes_label_weight=axes_label_weight,
                                                   axes_label_fontsize=axes_label_fontsize,
                                                   axes_label_distance=axes_label_distance,
                                                   antialiased_lines=antialiased_lines,
                                                   target_scaffold_labels_angle=target_scaffold_labels_angle,
                                                   query_scaffold_labels_angle=query_scaffold_labels_angle,
                                                   show_target_labels=show_labels,
                                                   show_query_labels=show_labels,
                                                   bottom_offset=bottom_offset,
                                                   top_offset=top_offset,
                                                   left_offset=left_offset,
                                                   right_offset=right_offset,
                                                   x_axis_visible=x_axis_visible,
                                                   y_axis_visible=y_axis_visible,
                                                   x_labelpad=x_labelpad,
                                                   y_labelpad=y_labelpad
                                                   )

    def draw_dot_plot_from_last_alignment(self, last_collection,
                                          output_prefix=None,
                                          extension_list=("png", ),
                                          target_black_list=(), target_white_list=(),
                                          target_ordered_list=(), target_reverse_list=(),
                                          query_black_list=(), query_white_list=(),
                                          query_ordered_list=(), query_reverse_list=(),
                                          remove_scaffolds_absent_in_target_ordered_list=False,
                                          remove_scaffolds_absent_in_query_ordered_list=False,
                                          figsize=(16, 16), dpi=300,
                                          auto_scale_figure=False,
                                          mbp_per_inch=None,
                                          figure_height=None,
                                          figure_width=None,
                                          grid_color='black',
                                          bar_color='grey',
                                          same_strand_color='red',
                                          diff_strand_color='blue',
                                          title=None,
                                          target_label=None,
                                          query_label=None,
                                          gridwidth=1,
                                          show_grid=True,
                                          linewidth=0.01,
                                          scaffold_label_fontsize=13,
                                          axes_label_fontstyle="italic",
                                          axes_label_weight="normal",
                                          axes_label_fontsize=15,
                                          axes_label_distance=6,
                                          antialiased_lines=None,
                                          target_scaffold_labels_angle=45,
                                          query_scaffold_labels_angle=0,
                                          show_target_labels=True,
                                          show_query_labels=True,
                                          bottom_offset=0.1,
                                          top_offset=0.9,
                                          left_offset=0.1,
                                          right_offset=0.9,
                                          x_axis_visible=False,
                                          y_axis_visible=False,
                                          x_labelpad=None,
                                          y_labelpad=None,
                                          show_length_ticks=False,
                                          show_tick_grid=False,
                                          tick_step=10000000,
                                          tick_unit=1000000):

        target_scaffold_list = self.get_filtered_scaffold_list(last_collection.target_scaffold_list,
                                                               scaffold_black_list=target_black_list,
                                                               sort_scaffolds=False,
                                                               scaffold_ordered_list=target_ordered_list,
                                                               scaffold_white_list=target_white_list,
                                                               remove_scaffolds_absent_in_ordered_list=remove_scaffolds_absent_in_target_ordered_list)

        query_scaffold_list = self.get_filtered_scaffold_list(last_collection.query_scaffold_list,
                                                              scaffold_black_list=query_black_list,
                                                              sort_scaffolds=False,
                                                              scaffold_ordered_list=query_ordered_list,
                                                              scaffold_white_list=query_white_list,
                                                              remove_scaffolds_absent_in_ordered_list=remove_scaffolds_absent_in_query_ordered_list)

        target_length_df = last_collection.target_scaffold_lengths.loc[target_scaffold_list]
        target_length_df["cum_end"] = target_length_df["length"].cumsum()
        target_length_df["cum_start"] = target_length_df["cum_end"] - target_length_df["length"]

        total_target_len = target_length_df["length"].sum()

        query_length_df = last_collection.query_scaffold_lengths.loc[query_scaffold_list]
        query_length_df["cum_end"] = query_length_df["length"].cumsum()
        query_length_df["cum_start"] = query_length_df["cum_end"] - query_length_df["length"]

        total_query_len = query_length_df["length"].sum()

        bar_width_fraction = 0.04
        bar_width = int(max(total_query_len, total_target_len) * bar_width_fraction)

        print("%s\tDrawing..." % str(datetime.datetime.now()))
        print("%s\t\tInitializing figure..." % str(datetime.datetime.now()))

        query_to_target_ratio = float(total_query_len) / float(total_target_len)

        if auto_scale_figure:
            if mbp_per_inch is not None:
                figsi = (int(total_target_len / mbp_per_inch / 1000000), int(total_query_len / mbp_per_inch/1000000))
            elif figure_width is not None:
                figsi = (figure_width, int(query_to_target_ratio * figure_width))
            elif figure_height is not None:
                figsi = (int(figure_height / query_to_target_ratio), figure_height)
        else:
            figsi = figsize

        figure = plt.figure(figsize=figsi, dpi=dpi)
        ax = plt.subplot(1, 1, 1)

        print("%s\t\tInitializing figure finished..." % str(datetime.datetime.now()))
        print("%s\t\tDrawing grid..." % str(datetime.datetime.now()))

        ax.add_patch(Rectangle((0, total_query_len), total_target_len, bar_width, color=bar_color))  # top bar
        ax.add_patch(Rectangle((0, -bar_width), total_target_len, bar_width, color=bar_color))       # bottom bar
        ax.add_patch(Rectangle((-bar_width, 0), bar_width, total_query_len, color=bar_color))        # left bar
        ax.add_patch(Rectangle((total_target_len, 0), bar_width, total_query_len, color=bar_color))  # right bar

        """
        @staticmethod
        def create_tick_formatter_function(max_value, tick_type="nucleotide"):
            max_val = max_value * 1.1
            if tick_type == "nucleotide":
                if max_val // (10 ** 9) > 2:
                    def tick_formater(x, pos):
                        return '%1.1f Gbp' % (x * 1e-9)
                elif max_val // (10 ** 6) > 200:
                    def tick_formater(x, pos):
                        return '%.0f Mbp' % (x * 1e-6)
                elif max_val // (10 ** 6) > 2:
                    def tick_formater(x, pos):
                        return '%.1f Mbp' % (x * 1e-6)
                elif max_val // (10 ** 3) > 2:
                    def tick_formater(x, pos):
                        return '%.1f kbp' % (x * 1e-3)
                else:
                    def tick_formater(x, pos):
                        return '%i bp' % (int(x))
    
                return FuncFormatter(tick_formater)
    
            else:
                raise ValueError("ERROR!!! Tick formter for %s is not implemented yet!" % tick_type)
        """

        if show_grid:
            for query_cum_start in query_length_df["cum_start"]:
                ax.add_line(Line2D((-bar_width, total_target_len+bar_width), (query_cum_start, query_cum_start),
                                   color=grid_color, linewidth=gridwidth))

            ax.add_line(Line2D((-bar_width, total_target_len+bar_width), (total_query_len, total_query_len),
                               color=grid_color, linewidth=gridwidth))

            for target_cum_start in target_length_df["cum_start"]:
                ax.add_line(Line2D((target_cum_start, target_cum_start), (-bar_width, total_query_len + bar_width),
                                   color=grid_color, linewidth=gridwidth))
            ax.add_line(Line2D((total_target_len, total_target_len), (-bar_width, total_query_len + bar_width),
                               color=grid_color, linewidth=gridwidth))
        else:
            for query_cum_start in query_length_df["cum_start"]:
                ax.add_line(Line2D((-bar_width, 0), (query_cum_start, query_cum_start),
                                   color=grid_color, linewidth=gridwidth))
                ax.add_line(Line2D((total_target_len, total_target_len + bar_width), (query_cum_start, query_cum_start),
                                   color=grid_color, linewidth=gridwidth))

            ax.add_line(Line2D((-bar_width, 0), (total_query_len, total_query_len),
                               color=grid_color, linewidth=gridwidth))
            ax.add_line(Line2D((total_target_len, total_target_len + bar_width), (total_query_len, total_query_len),
                               color=grid_color, linewidth=gridwidth))

            for target_cum_start in target_length_df["cum_start"]:
                ax.add_line(Line2D((target_cum_start, target_cum_start), (-bar_width, 0),
                                   color=grid_color, linewidth=gridwidth))
                ax.add_line(Line2D((target_cum_start, target_cum_start), (total_query_len, total_query_len + bar_width),
                                   color=grid_color, linewidth=gridwidth))

            ax.add_line(Line2D((total_target_len, total_target_len), (-bar_width, 0),
                               color=grid_color, linewidth=gridwidth))
            ax.add_line(Line2D((total_target_len, total_target_len), (total_query_len, total_query_len + bar_width),
                               color=grid_color, linewidth=gridwidth))

        if show_length_ticks:
            for target_scaffold in target_length_df.index:
                tick_labels = np.arange(0, target_length_df.loc[target_scaffold, "length"], tick_step)
                tick_list = list(tick_labels + target_length_df.loc[target_scaffold, "cum_start"])
                tick_labels = list(map(str, tick_labels // tick_unit))
                tick_number = len (tick_list)
                if len(tick_list) > 1:
                    for tick in tick_list[1:]:
                        ax.add_line(Line2D((tick, tick), (-bar_width/4, 0),
                                           color="black", linewidth=gridwidth / 4))
                        ax.add_line(Line2D((tick, tick), (total_query_len + bar_width / 4, total_query_len),
                                           color="black", linewidth=gridwidth / 4))
                        if show_tick_grid:
                            ax.add_line(Line2D((tick, tick), (0, total_query_len),
                                               color="black", linewidth=gridwidth / 8))
                    if len(tick_list) >= 5:
                        counter = 0
                        for tick, tick_label in zip(tick_list[::5][1:], tick_labels[::5][1:]):
                            counter += 1

                            ax.add_line(Line2D((tick, tick), (-bar_width * 0.4, 0),
                                               color="darkred", linewidth=gridwidth / 2))
                            ax.add_line(Line2D((tick, tick), (total_query_len + bar_width * 0.4, total_query_len),
                                               color="darkred", linewidth=gridwidth / 2))
                            if show_tick_grid:
                                ax.add_line(Line2D((tick, tick), (0, total_query_len),
                                                   color="darkred", linewidth=gridwidth / 8))
                            if tick_number - counter * 5 >= 2:
                                ax.text(tick, -bar_width * 0.5, tick_label,
                                        fontsize=5, #scaffold_label_fontsize/3,
                                        horizontalalignment='center',
                                        verticalalignment='top', )
                                ax.text(tick, total_query_len + bar_width * 0.5, tick_label,
                                        fontsize=5,  # scaffold_label_fontsize/3,
                                        horizontalalignment='center',
                                        verticalalignment='bottom', )

            for query_scaffold in query_length_df.index:
                tick_labels = np.arange(0, query_length_df.loc[query_scaffold, "length"], tick_step)
                tick_list = list(tick_labels + query_length_df.loc[query_scaffold, "cum_start"])
                tick_labels = list(map(str, tick_labels // tick_unit))

                if len(tick_list) > 1:
                    for tick in tick_list[1:]:
                        ax.add_line(Line2D((-bar_width/4, 0), (tick, tick),
                                           color="black", linewidth=gridwidth / 4))
                        ax.add_line(Line2D((total_target_len + bar_width / 4, total_target_len), (tick, tick),
                                           color="black", linewidth=gridwidth / 4))
                        if show_tick_grid:
                            ax.add_line(Line2D((0, total_target_len), (tick, tick),
                                               color="black", linewidth=gridwidth / 8))
                    if len(tick_list) >= 5:
                        for tick, tick_label in zip(tick_list[::5][1:], tick_labels[::5][1:]):
                            ax.add_line(Line2D((-bar_width * 0.4, 0), (tick, tick),
                                               color="darkred", linewidth=gridwidth / 2))
                            ax.add_line(Line2D((total_target_len + bar_width * 0.4, total_target_len), (tick, tick),
                                               color="darkred", linewidth=gridwidth / 2))
                            if show_tick_grid:
                                ax.add_line(Line2D((0, total_target_len), (tick, tick),
                                                   color="darkred", linewidth=gridwidth / 8))
                            ax.text(-bar_width * 0.5, tick, tick_label,
                                    fontsize=5,  # scaffold_label_fontsize/3,
                                    horizontalalignment='right',
                                    verticalalignment='center', )
                            ax.text(total_target_len + bar_width * 0.5, tick, tick_label,
                                    fontsize=5,  # scaffold_label_fontsize/3,
                                    horizontalalignment='left',
                                    verticalalignment='center', )

        print("%s\t\tDrawing grid finished..." % str(datetime.datetime.now()))
        print("%s\t\tAdding labels..." % str(datetime.datetime.now()))

        if show_target_labels:
            for target_scaffold_id in target_scaffold_list:
                ax.text((target_length_df.loc[target_scaffold_id]["cum_start"] + target_length_df.loc[target_scaffold_id]["cum_end"])/2,
                        total_query_len + 1.5 * bar_width, target_scaffold_id, fontsize=scaffold_label_fontsize,
                        rotation=target_scaffold_labels_angle,
                        horizontalalignment='left',
                        verticalalignment='bottom',)
                ax.text((target_length_df.loc[target_scaffold_id]["cum_start"] + target_length_df.loc[target_scaffold_id]["cum_end"])/2,
                        -1.5 * bar_width, target_scaffold_id, fontsize=scaffold_label_fontsize,
                        rotation=target_scaffold_labels_angle,
                        horizontalalignment='right',
                        verticalalignment='top',)
        if show_query_labels:
            for query_scaffold_id in query_scaffold_list:
                ax.text(total_target_len + 1.5 * bar_width,
                        (query_length_df.loc[query_scaffold_id]["cum_start"] + query_length_df.loc[query_scaffold_id]["cum_end"])/2,
                         query_scaffold_id, fontsize=scaffold_label_fontsize,
                        rotation=query_scaffold_labels_angle,
                        horizontalalignment='left',
                        verticalalignment='center' if query_scaffold_labels_angle == 0 else 'bottom')
                ax.text(-1.5 * bar_width,
                        (query_length_df.loc[query_scaffold_id]["cum_start"] + query_length_df.loc[query_scaffold_id]["cum_end"])/2,
                         query_scaffold_id, fontsize=scaffold_label_fontsize,
                        rotation=query_scaffold_labels_angle,
                        horizontalalignment='right',
                        verticalalignment='center' if query_scaffold_labels_angle == 0 else 'top')

        if title:
            plt.title(title)
        if target_label:
            ax.text(total_target_len/2,
                    -bar_width*axes_label_distance,
                    target_label,
                    fontsize=axes_label_fontsize,
                    fontstyle=axes_label_fontstyle,
                    fontweight=axes_label_weight,
                    horizontalalignment='center',
                    verticalalignment='center')
        if query_label:
            ax.text(-bar_width*axes_label_distance,
                    total_query_len/2,
                    query_label,
                    fontsize=axes_label_fontsize,
                    fontstyle=axes_label_fontstyle,
                    fontweight=axes_label_weight,
                    rotation=90,
                    horizontalalignment='center',
                    verticalalignment='center')

        plt.xlim(xmin=-bar_width * 2, xmax=total_target_len + 2 * bar_width)
        plt.ylim(ymin=-bar_width * 2, ymax=total_query_len + 2 * bar_width)

        if not x_axis_visible:
            ax.spines['bottom'].set_color('none')
            ax.spines['top'].set_color('none')

        if not y_axis_visible:
            ax.spines['right'].set_color('none')
            ax.spines['left'].set_color('none')

        ax.get_yaxis().set_visible(y_axis_visible)
        ax.get_xaxis().set_visible(x_axis_visible)

        print("%s\t\tAdding labels finished..." % str(datetime.datetime.now()))
        print("%s\t\tDrawing alignments..." % str(datetime.datetime.now()))

        if isinstance(last_collection, CollectionLast):
            def get_strand_specific_records(collection, target_scaffold_id, query_scaffold_id):
                return (collection.records[collection.records[collection.target_id_syn].isin([target_scaffold_id])
                                            & collection.records[collection.query_id_syn].isin([query_scaffold_id])
                                            & (collection.records[collection.query_strand_syn] == collection.records[collection.target_strand_syn])],
                        collection.records[collection.records[collection.target_id_syn].isin([target_scaffold_id])
                                            & collection.records[collection.query_id_syn].isin([query_scaffold_id])
                                            & (collection.records[collection.query_strand_syn] != collection.records[collection.target_strand_syn])])

        elif isinstance(last_collection, CollectionPSL):
            def get_strand_specific_records(collection, target_scaffold_id, query_scaffold_id):
                return (collection.records[collection.records[collection.target_id_syn].isin([target_scaffold_id])
                                           & collection.records[collection.query_id_syn].isin([query_scaffold_id])
                                           & (collection.records[collection.query_strand_syn] == "++")],
                        collection.records[collection.records[collection.target_id_syn].isin([target_scaffold_id])
                                           & collection.records[collection.query_id_syn].isin([query_scaffold_id])
                                           & (collection.records[collection.query_strand_syn] == "+-")])
        else:
            raise ValueError("ERROR!!! Unknown collection type (neither CollectionPSL nor CollectionLast)")

        def line_segments_generator(dataframe):
            for row_tuple in dataframe.itertuples(index=False):
                yield (row_tuple[:2], row_tuple[2:])

        for query_scaffold_id in query_scaffold_list:
            for target_scaffold_id in target_scaffold_list:

                same_strand_records_df, diff_strand_records_df = get_strand_specific_records(last_collection, target_scaffold_id, query_scaffold_id)
                """
                same_strand_records_df = \
                    last_collection.records[last_collection.records["target_id"].isin([target_scaffold_id])
                                            & last_collection.records["query_id"].isin([query_scaffold_id])
                                            & (last_collection.records["query_strand"] == last_collection.records["target_strand"])]
    
    
                
                diff_strand_records_df = \
                    last_collection.records[last_collection.records["target_id"].isin([target_scaffold_id])
                                            & last_collection.records["query_id"].isin([query_scaffold_id])
                                            & (last_collection.records["query_strand"] != last_collection.records["target_strand"])]
                """
                if not same_strand_records_df.empty:
                    data = pd.DataFrame()
                    """
                    data["x1"] = same_strand_records["target_start"] + target_length_df.loc[target_scaffold_id]["cum_start"]
                    data["y1"] = same_strand_records["query_start"] + query_length_df.loc[query_scaffold_id]["cum_start"]
    
                    data["x2"] = data["x1"] + same_strand_records["target_hit_len"] - 1
                    data["y2"] = data["y1"] + same_strand_records["query_hit_len"] - 1
                    """
                    data["x1"] = same_strand_records_df[last_collection.target_start_syn] + target_length_df.loc[target_scaffold_id]["cum_start"]
                    data["y1"] = same_strand_records_df[last_collection.query_start_syn] + query_length_df.loc[query_scaffold_id]["cum_start"]
                    if isinstance(last_collection, CollectionLast):
                        data["x2"] = data["x1"] + same_strand_records_df[last_collection.target_hit_len_syn] - 1
                        data["y2"] = data["y1"] + same_strand_records_df[last_collection.query_hit_len_syn] - 1
                    elif isinstance(last_collection, CollectionPSL):
                        data["x2"] = same_strand_records_df[last_collection.target_end_syn] + target_length_df.loc[target_scaffold_id]["cum_start"] - 1
                        data["y2"] = same_strand_records_df[last_collection.query_end_syn] + query_length_df.loc[query_scaffold_id]["cum_start"] - 1

                    lines = LineCollection(line_segments_generator(data), colors=same_strand_color, linestyle='solid',
                                           linewidths=linewidth, antialiased=antialiased_lines)

                    ax.add_collection(lines)

                if not diff_strand_records_df.empty:
                    data = pd.DataFrame()
                    data["x1"] = diff_strand_records_df[last_collection.target_start_syn] + target_length_df.loc[target_scaffold_id]["cum_start"]
                    if isinstance(last_collection, CollectionLast):
                        # in last tab format coordinates for minus strand start from the end of sequence
                        data["y1"] = query_length_df.loc[query_scaffold_id]["length"] - diff_strand_records_df[last_collection.query_start_syn] + query_length_df.loc[query_scaffold_id]["cum_start"]
                        data["x2"] = data["x1"] + diff_strand_records_df[last_collection.target_hit_len_syn] - 1
                        data["y2"] = data["y1"] - diff_strand_records_df[last_collection.query_hit_len_syn] + 1

                    elif isinstance(last_collection, CollectionPSL):
                        # in PSL format coordinates for minus strand start from the start of sequence
                        data["y1"] = diff_strand_records_df[last_collection.query_end_syn] + query_length_df.loc[query_scaffold_id]["cum_start"] - 1
                        data["x2"] = diff_strand_records_df[last_collection.target_end_syn] + target_length_df.loc[target_scaffold_id]["cum_start"] - 1
                        data["y2"] = diff_strand_records_df[last_collection.query_start_syn] + query_length_df.loc[query_scaffold_id]["cum_start"] - 1
                    #print(data)
                    lines = LineCollection(line_segments_generator(data), colors=diff_strand_color, linestyle='solid',
                                           linewidths=linewidth, antialiased=antialiased_lines)
                    #print(data)
                    #print(diff_strand_records_df[["target_start", "target_hit_len", "query_start", "query_hit_len"]])
                    ax.add_collection(lines)

        print("%s\t\tDrawng alignments finished..." % str(datetime.datetime.now()))
        print("%s\tDrawing finished..." % str(datetime.datetime.now()))

        plt.subplots_adjust(left=left_offset, bottom=bottom_offset, right=right_offset, top=top_offset)

        if output_prefix:
            print("%s\tWriting to file..." % str(datetime.datetime.now()))
            for extension in extension_list:
                plt.savefig("%s.%s" % (output_prefix, extension))

            print("%s\tWriting to file finished..." % str(datetime.datetime.now()))
