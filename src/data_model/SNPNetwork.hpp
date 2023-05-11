//
// Created by juli on 25.05.22.
//

#ifndef GENEPISEEKER_SNPNETWORK_HPP
#define GENEPISEEKER_SNPNETWORK_HPP

/**
 * @brief contains the SNP-SNP-interaction network and provides access to the graph in different ways
 *
 * this implementation uses an internal igraph representation of the network
 */

#include "SNP_t.hpp"
#include <vector>
#include <igraph.h>
#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <type_traits>

namespace epi {

    class SNPNetwork_base {

    protected:
        // boost graph stuff
        struct node_data {
            SNP snp;
        };
        struct edge_data {
            // std::set<uint16_t> labels{};
            std::uint64_t labels = 0;
        };

        typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS, node_data, edge_data> network_repr;
        typedef boost::graph_traits<network_repr>::vertex_descriptor vertex_descriptor;
        typedef boost::graph_traits<network_repr>::vertex_iterator vertex_iter;
        typedef boost::property_map<network_repr, SNP node_data::*>::type vertex_snp_map_type;

        network_repr network_boost;
        std::unordered_map<SNP, vertex_descriptor> node_boost_map;

        std::vector<SNP_t> get_adjacent_snps(vertex_descriptor v_desc);
        void request_igraph_repr();

        std::unordered_set<uint_fast64_t> edge_collision_list;

        std::vector<std::string> edge_labels;
        std::unordered_map<std::string, uint16_t> edge_label_map;
        uint16_t next_edge_label_id = 0;


        // igraph stuff (needed for some graph algorithms)
        igraph_t network_igraph;
        bool igraph_valid = false;
        bool igraph_created = false;
        // std::unordered_map<SNP, igraph_integer_t> node_snp_map;
        std::vector<vertex_descriptor> node_igraph_map;
        std::unordered_map<vertex_descriptor, igraph_integer_t> node_igraph_rev_map;

        // for non-boot graphs
        bool boost_mode = false;
        std::unordered_map<SNP, std::vector<SNP>> needl_adjacency_list;
        std::unordered_map<uint_fast64_t, edge_data> needl_edge_list;
        size_t needl_num_nodes = 0;
        size_t needl_num_edges = 0;
        std::vector<SNP> needl_igraph_map;
        std::unordered_map<SNP, igraph_integer_t> needl_igraph_rev_map;

    public:

        SNPNetwork_base(bool boost_mode);
        virtual ~SNPNetwork_base();

        void clear();

        template<class Iterator>
        void set_nodes(const Iterator _begin, const Iterator _end);


        template<class Iterator>
        void add_nodes(const Iterator _begin, const Iterator _end);
        void add_node(const SNP_t &snp);

        template<class Iterator>
        void add_edges(const Iterator _begin, const Iterator _end, std::string edge_tag = "");
        void add_edge(const SNPEdge &edge, std::string edge_tag);
        void add_edge(const SNPEdge &edge, uint16_t edge_tag = 0);

        uint16_t get_edge_tag_id(const std::string& edge_tag);
        std::vector<std::string> get_all_edge_tags();


        template<class Iterator>
        void remove_edges(const Iterator _begin, const Iterator _end);
        void remove_edge(const SNPEdge & edge);


        template<class Iterator>
        void remove_nodes(const Iterator _begin, const Iterator _end);
        void remove_node(const SNP_t &node);

        void replace_nodes(const std::vector<std::pair<SNP_t, SNP_t>> &pairs);

        size_t num_nodes();
        size_t num_edges();

        std::vector<SNP_t> get_network_snps();
        std::vector<SNP_t> get_adjacent_snps(const SNP_t& snp);
        std::set<std::string> get_edge_labels(const SNPEdge &edge);
        std::set<uint16_t> get_edge_label_ids(const SNPEdge &edge);

        template<class Iterator>
        std::vector<SNP_t> get_adjacent_snps(const Iterator _begin, const Iterator _end);
        std::vector<SNP_t> get_adjacent_snps(const std::vector<SNP_t> &snps);
        std::vector<SNP_t> get_adjacent_snps(const std::unordered_set<SNP_t, SNP_t::SNPHash>& snps);


        std::vector<std::pair<SNPEdge, std::set<uint16_t>>> get_adjacent_edges(const SNP_t& snp);
        std::vector<std::pair<SNPEdge, std::set<uint16_t>>> get_adjacent_edges(vertex_descriptor v_desc);


        std::vector<std::pair<SNP_t, std::vector<SNP_t>>> get_adjacency_list();
        long get_network_diameter();
        bool is_connected();
        size_t get_degree(const SNP_t& snp);
        bool edge_exists(const SNPEdge & edge);
        void clear_edges();

        std::vector<std::vector<SNP_t>> cluster_leiden(double resolution, double beta, unsigned long max_leiden_steps);
        std::vector<SNP_t> get_articulation_points();

        // SNP_t get_node(igraph_integer_t i);
        // igraph_integer_t get_node_index(const SNP_t &snp);
        bool contains_node(const SNP_t &snp);


        class all_iterator : public std::iterator<std::forward_iterator_tag, SNP_t, long, const SNP_t *, SNP_t> {
            vertex_iter curr;
            vertex_snp_map_type map;
            std::unordered_map<SNP, std::vector<SNP>>::iterator needl_curr;
            bool boost_mode;
        public:
            all_iterator(vertex_iter _curr, vertex_snp_map_type _map);
            all_iterator(std::unordered_map<SNP, std::vector<SNP>>::iterator _needl_curr);
            all_iterator& operator++();
            all_iterator operator++(int);

            bool operator==(all_iterator other) const;
            bool operator!=(all_iterator other) const;

            SNP_t operator*() const;
        };

        class all_iterable {
        public:
            all_iterable(vertex_iter __begin, vertex_iter __end, vertex_snp_map_type __map);
            all_iterable(std::unordered_map<SNP, std::vector<SNP>>::iterator __begin, std::unordered_map<SNP, std::vector<SNP>>::iterator __end);
            all_iterator begin();
            all_iterator end();
        private:
            vertex_iter _begin, _end;
            vertex_snp_map_type _map;
            std::unordered_map<SNP, std::vector<SNP>>::iterator _needl_begin, _needl_end;
            bool _boost_mode;
        };

        all_iterable all();
    };

    class SNPNetwork : public SNPNetwork_base {
    public:
        SNPNetwork(bool boost_mode = false);
        SNPNetwork(const SNPNetwork & other);
        SNPNetwork& operator=(const SNPNetwork&) = default;
    };


} // epi

#include "SNPNetwork.ipp"

#ifdef HEADER_ONLY
#include "SNPNetwork.cpp"
#endif

#endif //GENEPISEEKER_SNPNETWORK_HPP
