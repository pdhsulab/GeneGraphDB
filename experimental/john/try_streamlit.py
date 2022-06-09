import pandas as pd
import streamlit as st
from streamlit_agraph import Config, Edge, Node, agraph

from db_inference.simple_neo4j import SimpleNeo4j
from db_inference.batch_compute_colocalization_neo4j import get_colocalization_scores

import pandas as pd
import numpy as np
import json
import plotly.express as px
from pyvis.network import Network
import networkx as nx
from db_inference import simple_sql_db
import prot_db.bigtable_constants as btc

from google.cloud.bigtable.row_set import RowSet
from google.cloud.bigtable import row_filters


# HARDCODED_BAIT_P30 = "f3aeb9d18ae62eb49a"  # very slow
# HARDCODED_BAIT_P30 = "6bc2b0f7f986cd84bf"  # semi slow
# HARDCODED_BAIT_P30 = "3d83934184a2717769" # fast
HARDCODED_BAIT_P30 = "4046582228d4fb75f3" # slow


def main():
    st.title('Protein Graph Demo')

    colocalization_scores = get_colocalization_scores([HARDCODED_BAIT_P30], 'test_bait')
    colocalization_df = pd.DataFrame.from_dict(colocalization_scores, orient='index')
    st.dataframe(colocalization_df)

    agraph_ret = build_graph(colocalization_df)
    st.markdown(agraph_ret)

    if (agraph_ret is not None) and ('action' in agraph_ret) and (agraph_ret['action'] == 'onClickNode'):
        query_p30 = agraph_ret['node']
        annotations, sources = get_info(query_p30)

        # display source data
        source_studies = []
        for source in sources:
            source_studies.extend(source.keys())
        source_markdown = '## MGNify study sources:\n'
        for study in sorted(set(source_studies)):
            source_markdown += f"* {study}\n"
        st.markdown(source_markdown)

        annotation_df = pd.DataFrame(annotations)
        st.dataframe(annotation_df)


def get_bt_row_keys(p30_id):
    seq_db = simple_sql_db.SequenceSqlDb()
    sql_db = simple_sql_db.SimpleSqlDb()

    row_keys = []
    tgt_p100s = [row['p100'] for row in sql_db.get_p30_cluster_members(p30_id)]
    for pid in tgt_p100s:
        aa_seq = seq_db.get_sequence(pid)
        if aa_seq.endswith("*"):
            aa_seq = aa_seq[:-1]
        row_key = btc.row_key(aa_seq)
        row_keys.append(row_key)
    return row_keys


def get_bt_rows(row_keys):
    table = btc.get_table(cloud=True)
    row_filter = row_filters.CellsColumnLimitFilter(1)
    row_set = RowSet()
    for row_key in row_keys:
        row_set.add_row_key(row_key)
    rows = table.read_rows(row_set=row_set, filter_=row_filter)
    return rows

def parse_annotation(anno):
    parsed = {}
    for line in anno.split(";"):
        k, v = line.split("=")
        parsed[k] = v
    return parsed

def get_info(query_p30):
    row_keys = get_bt_row_keys(query_p30)
    rows = get_bt_rows(row_keys)

    sources = []
    annotations = []
    for row in rows:
        sources.append(btc.get_mgnify_study_to_analysis(row))
        annotations.append(btc.get_annotations_from_row(row))

    annotations = [parse_annotation(anno) for anno in annotations if anno is not None]
    return annotations, sources


def build_graph(colocalization_df):
    nodes = []
    edges = []

    bait_id = colocalization_df.iloc[0]['bait_p30']
    bait_num_p100s = int(colocalization_df.iloc[0]['num_bait_p100s'])
    nodes.append( Node(id=bait_id,
                       label="bait",
                       size=bait_num_p100s))

    for _, row in colocalization_df.iterrows():
        nodes.append(
            Node(id=row['tgt_p30'],
                 label=row['tgt_p30'],
                 size=int(row['num_tgt_p100s']))
        )

        edges.append(
            Edge(source=bait_id,
                 target=row['tgt_p30'],
                 label=str(row['bait_colocalization']),
                 strokeWidth=10.0 * row['bait_colocalization'],
            )
        )

        # edges.append(
        #     Edge(source=row['tgt_p30'],
        #          target=bait_id,
        #          label=str(row['tgt_colocalization']),
        #          strokeWidth=10.0 * row['tgt_colocalization'],
        #          type="CURVE_SMOOTH"
        #     )
        # )

    config = Config(width=500,
                    height=500,
                    directed=True,
                    nodeHighlightBehavior=True,
                    highlightColor="#F7A7A6", # or "blue"
                    collapsible=True,
                    node={'labelProperty':'label', 'size': 1000},
                    link={'labelProperty': 'label', 'renderLabel': False},
                    # **kwargs e.g. node_size=1000 or node_color="blue"
                    )

    return_value = agraph(nodes=nodes,
                          edges=edges,
                          config=config)
    return return_value


if __name__ == "__main__":
    main()
