from genegraphdb import *
import csv
import os

def get_runtime_summarystats(comment, prot_time = 0, infile_name = "ggdb_load_stats.csv", outfile_name = "ggdb_summary_stats.csv"):
    outfile = open(outfile_name, "a")
    with open(infile_name, newline='') as f:
        reader = csv.reader(f)
        next(reader)
        for line in f:
            line = line.strip('\n').split(',')
            cur_load_time = float(line[1])
            p2p_edge_time = line[2]
            if p2p_edge_time != "null":
                cur_p2p_edge_time = float(p2p_edge_time)
            else:
                cur_p2p_edge_time = prot_time
            break
        try:
            next(reader)
            for line in f:
                line = line.strip('\n').split(',')
                cur_load_time += float(line[1])
                if p2p_edge_time != "null":
                    cur_p2p_edge_time += float(line[2])
        except StopIteration:
            print("end of csv reached")
        if p2p_edge_time == "null":
            cur_load_time += cur_p2p_edge_time
        print(str(cur_load_time) + "," + str(cur_p2p_edge_time) + "," + comment, file=outfile)
