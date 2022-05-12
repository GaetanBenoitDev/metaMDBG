#include <stdio.h>
#include <stdlib.h>
#include "abpoa.h"
#include "abpoa_graph.h"
#include "utils.h"

/* example of dot file for graphviz */
/*
   digraph test1 {
       a -> b -> c;
       a -> {x y};
       b [shape=box];
       c [label="hello\nworld",color=blue,fontsize=24,
         fontname="Palatino-Italic",fontcolor=red,style=filled];
       a -> z [label="hi", weight=100];
       x -> z [label="multi-line\nlabel"];
       edge [style=dashed,color=red];
       b -> x;
       {rank=same; b x}
   }
   graph test2 {
       a -- b -- c [style=dashed];
       a -- {x y};
       x -- c [w=10.0];
       x -- y [w=5.0,len=3];
   }
*/


extern char ab_nt256_table[256];
// base (index, rank, node_id)
// A (1, 1, 2) A: base 1: index 1: rank 2: node_id
void abpoa_dump_pog(abpoa_t *ab, abpoa_para_t *abpt) {
    char PROG[20] = "abpoa"; int font_size=24;

    abpoa_graph_t *abg = ab->abg;
    if (abg->is_topological_sorted == 0) abpoa_topological_sort(abg, abpt);

    // all settings
    // char node_color[5][10] = {"purple3", "red3", "seagreen4", "gold2", "gray"}; // ACGTN
    char node_color[5][10] = {"pink1", "red1", "gold2", "seagreen4", "gray"}; // ACGTN
    // float dpi_size = 3000, graph_width = 100, graph_height = 6; 
    float node_width=1;
    char rankdir[5] = "LR", node_style[10]="filled", node_fixedsize[10]="true", node_shape[10]="circle";
    int show_aligned_mismatch = 1;

    int i, j, id, index, out_id; char base;
    char **node_label = (char**)_err_malloc(abg->node_n * sizeof(char*));
    for (i = 0; i < abg->node_n; ++i) node_label[i] = (char*)_err_malloc(sizeof(char) * 128);
 
    char *dot_fn = (char*)malloc(strlen(abpt->out_pog) + 10);
    strcpy(dot_fn, abpt->out_pog);
    FILE *fp = xopen(strcat(dot_fn, ".dot"), "w");
    fprintf(fp, "// %s graph dot file.\n// %d nodes.\n", PROG, abg->node_n);
    // fprintf(fp, "digraph ABPOA_graph {\n\tgraph [dpi=%f]; size=\"%f,%f\";\n\trankdir=\"%s\";\n\tnode [width=%f, style=%s, fixedsize=%s, shape=%s];\n", dpi_size, graph_width, graph_height, rankdir, node_width, node_style, node_fixedsize, node_shape);
    fprintf(fp, "digraph ABPOA_graph {\n\tgraph [rankdir=\"%s\"];\n\tnode [width=%f, style=%s, fixedsize=%s, shape=%s];\n", rankdir, node_width, node_style, node_fixedsize, node_shape);

    for (i = 0; i < abg->node_n; ++i) {
        id = abpoa_graph_index_to_node_id(abg, i);
        index = i;
        if (id == ABPOA_SRC_NODE_ID) {
            base = 'S';
            //sprintf(node_label[id], "\"%c\n(%d,%d,%d)\"", base, index, rank, id);
            // only show seq
            sprintf(node_label[id], "\"%c\n%d\"", base,index);
            fprintf(fp, "%s [color=%s, fontsize=%d]\n", node_label[id], node_color[4], font_size);
        } else if (id == ABPOA_SINK_NODE_ID) {
            base = 'E';
            //sprintf(node_label[id], "\"%c\n(%d,%d,%d)\"", base, index, rank, id);
            // only show seq
            sprintf(node_label[id], "\"%c\n%d\"", base,index);
            fprintf(fp, "%s [color=%s, fontsize=%d]\n", node_label[id], node_color[4], font_size);
        } else {
            base = ab_nt256_table[abg->node[id].base];
            //sprintf(node_label[id], "\"%c\n(%d,%d,%d)\"", base, index, rank, id);
            // only show seq
            sprintf(node_label[id], "\"%c\n%d\"", base,index);
            fprintf(fp, "%s [color=%s, fontsize=%d]\n", node_label[id], node_color[abg->node[id].base], font_size);
        }
    }
    int x_index = -1;
    for (i = 0; i < abg->node_n; ++i) {
        id = abpoa_graph_index_to_node_id(abg, i);
        // out_edge
        for (j = 0; j < abg->node[id].out_edge_n; ++j) {
            out_id = abg->node[id].out_id[j];
            fprintf(fp, "\t%s -> %s [label=\"%d\", penwidth=%d]\n", node_label[id], node_label[out_id], abg->node[id].out_weight[j], abg->node[id].out_weight[j]+1);
        }
        if (abg->node[id].aligned_node_n > 0) {
            fprintf(fp, "\t{rank=same; %s ", node_label[id]);
            for (j = 0; j < abg->node[id].aligned_node_n; ++j)
                fprintf(fp, "%s ", node_label[abg->node[id].aligned_node_id[j]]);
            fprintf(fp, "};\n");
            if (show_aligned_mismatch) {
                if (i > x_index) {
                    x_index = i;
                    // mismatch dashed line
                    fprintf(fp, "\t{ edge [style=dashed, arrowhead=none]; %s ", node_label[id]);
                    for (j = 0; j < abg->node[id].aligned_node_n; ++j) {
                        fprintf(fp, "-> %s ", node_label[abg->node[id].aligned_node_id[j]]);
                        index = abpoa_graph_node_id_to_index(abg, abg->node[id].aligned_node_id[j]);
                        x_index = index > x_index ? index : x_index;
                    }
                    fprintf(fp, "}\n");
                }
            }
        }
    }
    fprintf(fp, "}\n");

    for (i = 0; i < abg->node_n; ++i) free(node_label[i]); free(node_label);
    err_fclose(fp);

    char cmd[1024];
    char *type = strrchr(abpt->out_pog, '.');
    if (strcmp(type+1, "pdf") != 0 && strcmp(type+1, "png") != 0)
        err_fatal_simple("POG can only be dump to .pdf/.png file");
    sprintf(cmd, "dot %s -T%s > %s", dot_fn, type+1, abpt->out_pog);
    free(dot_fn);
    if (system(cmd) != 0) err_fatal(__func__, "Fail to plot %s DAG.", PROG);
}
