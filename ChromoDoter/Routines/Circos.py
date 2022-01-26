__author__ = 'mahajrod'

#import datetime
#import numpy as np


#import pandas as pd

#from ChromoDoter.Parsers.LAST import CollectionLast
#from ChromoDoter.Parsers.PSL import CollectionPSL


class CircosRoutines:

    @staticmethod
    def convert_synteny_psl_to_circos_links(collection_psl):

        links_df = collection_psl.records[[collection_psl.query_id_syn, collection_psl.query_start_syn, collection_psl.query_end_syn,
                                          collection_psl.target_id_syn, collection_psl.target_start_syn, collection_psl.target_end_syn,
                                          collection_psl.query_strand_syn]]
        links_df.loc[links_df[collection_psl.query_strand_syn] == "+-", collection_psl.target_start_syn], \
        links_df.loc[links_df[collection_psl.query_strand_syn] == "+-", collection_psl.target_end_syn] = \
                    links_df.loc[links_df[collection_psl.query_strand_syn] == "+-", collection_psl.target_end_syn], \
                    links_df.loc[links_df[collection_psl.query_strand_syn] == "+-", collection_psl.target_start_syn]

        links_df[collection_psl.query_strand_syn] = "#" + links_df[collection_psl.query_strand_syn]

        return links_df[[collection_psl.query_id_syn, collection_psl.query_start_syn, collection_psl.query_end_syn,
                         collection_psl.target_id_syn, collection_psl.target_start_syn, collection_psl.target_end_syn,
                         collection_psl.query_strand_syn]]


