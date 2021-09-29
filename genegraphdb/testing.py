from genegraphdb import *
import csv
import os
# To do: ordering of arguments is super unscalable: fix
def get_runtime_summarystats(comment, prot_time = 0, samples_path = "", infile_name = "ggdb_load_stats.csv", outfile_name = "ggdb_summary_stats.csv"):
    if samples_path != "":
        infile_name = samples_path + infile_name
        outfile_name = samples_path + outfile_name
    print(infile_name) #REMOVE THIS
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
