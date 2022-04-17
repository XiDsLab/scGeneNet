using NetworkInference
using Test
using DelimitedFiles

# These tests use a dataset generated from the 10-node Yeast1 network from http://gnw.sourceforge.net/
# GeneNetWeaver: In silico benchmark generation and performance profiling of network inference methods.
# Schaffter T, Marbach D, and Floreano D. Bioinformatics, 27(16):2263-70, 2011.

println("Getting nodes...")
data_path = joinpath("D:\\Project\\Pro_VMPLN\\code\\simulation_1.27.1\\data_PIDC", "nor")
data_path_save = joinpath("D:\\Project\\Pro_VMPLN\\code\\simulation_1.27.1\\result_PIDC", "nor")
##
p_vec = [1]
m_vec = [1,2]
net_vec = [1,2,3,4]
#net_vec = [1,2]
ARI_vec = [1,2,3]
##
for p_index = p_vec
    for m_index = m_vec
        for net_index = net_vec
            for ARI_index = ARI_vec
                for rep = 1:20
                    for g = 1:3
                        if net_index == 1
                            data_file_path = joinpath(data_path, string("norobs_mat_group",g,"_sam1000_dim300_drop",m_index,"_nettypeER-graph_ARI%3A",ARI_index,"_clustermethod1_rep",rep,".txt"))
                        end
                        if net_index == 2
                            data_file_path = joinpath(data_path, string("norobs_mat_group",g,"_sam1000_dim300_drop",m_index,"_nettypeHUB-graph_ARI%3A",ARI_index,"_clustermethod1_rep",rep,".txt"))
                        end
                        if net_index == 3
                            data_file_path = joinpath(data_path, string("norobs_mat_group",g,"_sam1000_dim300_drop",m_index,"_nettypeAFF-graph_ARI%3A",ARI_index,"_clustermethod1_rep",rep,".txt"))
                        end
                        if net_index == 4
                            data_file_path = joinpath(data_path, string("norobs_mat_group",g,"_sam1000_dim300_drop",m_index,"_nettypePA-graph_ARI%3A",ARI_index,"_clustermethod1_rep",rep,".txt"))
                        end
                        println(string("p_index",p_index,"_m_index",m_index,"_net_index",net_index,"_ARI_index",ARI_index,"_rep",rep,"_group",g))
                        ##
                        nodes = get_nodes(data_file_path)
                        ##
                        pidc_network = InferredNetwork(PIDCNetworkInference(), nodes)
                        ##
                        if net_index == 1
                            data_file_path_save = joinpath(data_path_save, string("PIDCnor_res_group",g,"_sam1000_dim300_drop",m_index,"_nettypeER-graph_ARI%3A",ARI_index,"_clustermethod1_rep",rep,".txt"))
                        end
                        if net_index == 2
                            data_file_path_save = joinpath(data_path_save, string("PIDCnor_res_group",g,"_sam1000_dim300_drop",m_index,"_nettypeHUB-graph_ARI%3A",ARI_index,"_clustermethod1_rep",rep,".txt"))
                        end
                        if net_index == 3
                            data_file_path_save = joinpath(data_path_save, string("PIDCnor_res_group",g,"_sam1000_dim300_drop",m_index,"_nettypeAFF-graph_ARI%3A",ARI_index,"_clustermethod1_rep",rep,".txt"))
                        end
                        if net_index == 4
                            data_file_path_save = joinpath(data_path_save, string("PIDCnor_res_group",g,"_sam1000_dim300_drop",m_index,"_nettypePA-graph_ARI%3A",ARI_index,"_clustermethod1_rep",rep,".txt"))
                        end
                        ##
                        write_network_file(data_file_path_save, pidc_network)
                    end
                end
            end
        end
    end
end
