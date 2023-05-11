//
// Created by juli on 25.05.22.
//

#include "SNPNetwork.hpp"

#include <utility>
#include "../util/types.hpp"
#include "../util/TimeLogger.hpp"
#include "SNPStorage.hpp"

#if (IGRAPH_THREAD_SAFE != 1)
#error Thread-safe version if igraph library needed
#endif

namespace epi {
    SNPNetwork_base::SNPNetwork_base(bool boost_mode) {
        get_edge_tag_id("");
        this->boost_mode = boost_mode;
    }

    SNPNetwork_base::~SNPNetwork_base() {
        if (igraph_created) igraph_destroy(&network_igraph);
    }

    void SNPNetwork_base::clear() {
        if (igraph_created) igraph_destroy(&network_igraph);
        igraph_valid = false;
        igraph_created = false;

        network_boost.clear();
        node_boost_map.clear();
        edge_collision_list.clear();

        needl_adjacency_list.clear();
        needl_edge_list.clear();
        needl_num_edges = 0;
        needl_num_nodes = 0;
    }

    void SNPNetwork_base::add_node(const SNP_t &snp) {
        if (boost_mode) {
            if (node_boost_map.find(snp.value) == node_boost_map.end()) {
                auto descriptor = boost::add_vertex({snp.value}, network_boost);
                node_boost_map.insert({snp.value, descriptor});
                igraph_valid = false;
            }
        } else {
            if (needl_adjacency_list.find(snp.value) == needl_adjacency_list.end()) {
                needl_adjacency_list.insert({ snp.value, {}});
                igraph_valid = false;
                needl_num_nodes++;
            }
        }
    }

    void SNPNetwork_base::add_edge(const SNPEdge &edge, std::string edge_tag) {
        // find edge tag
        auto tid = get_edge_tag_id(edge_tag);
        add_edge(edge, tid);
    }

    void SNPNetwork_base::add_edge(const SNPEdge &edge, uint16_t edge_tag) {
        if (edge.bits.snp1 == edge.bits.snp2) return;

        if (boost_mode) {
            if (edge_collision_list.find(edge.val) != edge_collision_list.end()) {
                auto snp1_it = node_boost_map.find(edge.bits.snp1);
                auto snp2_it = node_boost_map.find(edge.bits.snp2);

                auto edge_val_map = boost::get(&edge_data::labels, network_boost);
                auto boost_edge = boost::edge(snp1_it->second, snp2_it->second, network_boost).first;

                edge_val_map[boost_edge] |= 1 << edge_tag;
                return;
            }

            auto desc1 = node_boost_map.find(edge.bits.snp1);
            auto desc2 = node_boost_map.find(edge.bits.snp2);
            if (desc1 == node_boost_map.end() ||
                desc2 == node_boost_map.end()) {
                throw epi::Error("Tried to insert edge to SNP that is not present in the graph.");
            }

            // if (!boost::in_edge_set(network_boost, desc1->second, desc2->second)) {
            uint64_t labels = 0;
            labels |= 1 << edge_tag;
            boost::add_edge(desc1->second, desc2->second, {labels}, network_boost);
            edge_collision_list.insert(edge.val);
            igraph_valid = false;
            // }
        } else {
            auto edge_repr = needl_edge_list.find(edge.val);
            if (edge_repr != needl_edge_list.end()) {
                edge_repr->second.labels |= 1 << edge_tag;
                return;
            }

            auto desc1 = needl_adjacency_list.find(edge.bits.snp1);
            auto desc2 = needl_adjacency_list.find(edge.bits.snp2);
            if (desc1 == needl_adjacency_list.end() ||
                desc2 == needl_adjacency_list.end()) {
                throw epi::Error("Tried to insert edge to SNP that is not present in the graph.");
            }

            uint64_t labels = 0;
            labels |= 1 << edge_tag;
            needl_edge_list.insert({ edge.val, { labels }});
            desc1->second.push_back(edge.bits.snp2);
            desc2->second.push_back(edge.bits.snp1);
            ++needl_num_edges;
        }
    }


    size_t SNPNetwork_base::num_nodes() {
        if (boost_mode) {
            return boost::num_vertices(network_boost);
        } else return needl_num_nodes;
    }

    size_t SNPNetwork_base::num_edges() {
        if (boost_mode) {
            return boost::num_edges(network_boost);
        } else return needl_num_edges;
    }

    std::vector<SNP_t> SNPNetwork_base::get_network_snps() {
        std::vector<SNP_t> snps;
        if (boost_mode) {
            for (auto &node_item: node_boost_map) {
                snps.push_back(SNP_t(node_item.first));
            }
        } else {
            for (auto &node_item : needl_adjacency_list) {
                snps.push_back(SNP_t(node_item.first));
            }
        }

        return snps;
    }

    std::vector<SNP_t> SNPNetwork_base::get_adjacent_snps(const SNP_t& snp) {
        if(boost_mode) {
            auto item = node_boost_map.find(snp.value);
            if (item == node_boost_map.end()) {
                throw epi::Error("Requested SNP is not in the network.");
            }
            return get_adjacent_snps(item->second);
        } else {
            auto item = needl_adjacency_list.find(snp.value);
            if (item == needl_adjacency_list.end()) {
                throw epi::Error("Requested SNP is not in the network.");
            }
            std::vector<SNP_t> adjacent_snps;
            for (auto adjacent : item->second) {
                adjacent_snps.push_back(SNP_t(adjacent));
            }
            return adjacent_snps;
        }
    }

    std::vector<SNP_t> SNPNetwork_base::get_adjacent_snps(vertex_descriptor v_desc) {
        std::vector<SNP_t> adjacent_snps;
        auto snp_val_map = boost::get(&node_data::snp, network_boost);
        auto it = boost::adjacent_vertices(v_desc, network_boost);
        for (; it.first != it.second; it.first++) {
            adjacent_snps.push_back(SNP_t(snp_val_map[*it.first]));
        }

        return adjacent_snps;
    }

    std::vector<std::pair<SNPEdge, std::set<uint16_t>>> SNPNetwork_base::get_adjacent_edges(const SNP_t& snp) {
        if (boost_mode) {
            auto item = node_boost_map.find(snp.value);
            if (item == node_boost_map.end()) {
                throw epi::Error("Requested SNP is not in the network.");
            }
            return get_adjacent_edges(item->second);
        } else {
            auto adjacent_snps = get_adjacent_snps(snp);
            std::vector<std::pair<SNPEdge, std::set<uint16_t>>> adjacent_edges;
            for (auto & adjacent : adjacent_snps) {
                SNPEdge edge (snp, adjacent);
                auto edge_repr = needl_edge_list.find(edge.val);

                auto labels_ids = edge_repr->second.labels;
                std::set<uint16_t> labels;
                for (uint64_t i = 0; i < edge_labels.size(); i++) {
                    uint64_t bitmask = 1 << i;
                    if (bitmask & labels_ids) labels.insert(i);
                }

                adjacent_edges.emplace_back(edge, labels);
            }

            return adjacent_edges;
        }
    }

    std::vector<std::pair<SNPEdge, std::set<uint16_t>>> SNPNetwork_base::get_adjacent_edges(vertex_descriptor v_desc) {
        std::vector<std::pair<SNPEdge, std::set<uint16_t>>> adjacent_edges;

        auto snp_val_map = boost::get(&node_data::snp, network_boost);
        auto edge_val_map = boost::get(&edge_data::labels, network_boost);
        auto it = boost::out_edges(v_desc, network_boost);
        for (; it.first != it.second; it.first++) {
            vertex_descriptor adjacent_snp = it.first->m_source == v_desc ? it.first->m_target : it.first->m_source;
            auto labels_ids = edge_val_map[*it.first];
            std::set<uint16_t> labels;
            for (uint64_t i = 0; i < edge_labels.size(); i++) {
                uint64_t bitmask = 1 << i;
                if (bitmask & labels_ids) labels.insert(i);
            }

            adjacent_edges.emplace_back(SNPEdge(SNP_t(snp_val_map[v_desc]), SNP_t(snp_val_map[adjacent_snp])), labels);
        }

        return adjacent_edges;
    }

    std::vector<std::pair<SNP_t, std::vector<SNP_t>>> SNPNetwork_base::get_adjacency_list() {
        std::vector<std::pair<SNP_t, std::vector<SNP_t>>> adjacency_list;
        adjacency_list.reserve(num_nodes());

        if (boost_mode) {
            auto node_map = boost::get(&node_data::snp, network_boost);
            auto it = boost::vertices(network_boost);
            for (; it.first != it.second; it.first++) {
                auto neighbours = get_adjacent_snps(*it.first);
                adjacency_list.push_back({SNP_t(node_map[*it.first]), neighbours});
            }
        } else {
            for (auto & snp1 : needl_adjacency_list) {
                std::vector<SNP_t> adjacent;
                adjacent.reserve(snp1.second.size());
                for (auto snp2 : snp1.second) {
                    adjacent.push_back(SNP_t(snp2));
                }

                adjacency_list.emplace_back(SNP_t(snp1.first), adjacent);
            }
        }

        return adjacency_list;
    }

    std::vector<std::vector<SNP_t>> SNPNetwork_base::cluster_leiden(double resolution, double beta, unsigned long max_leiden_steps)
    {
        request_igraph_repr();

        igraph_vector_t membership;
        igraph_vector_init_seq(&membership, 0, num_nodes() - 1);
        igraph_integer_t num_clusters = num_nodes();

        igraph_real_t prev_score, score = -std::numeric_limits<double>::max();
        size_t num_leiden_steps = 0;
        do {
            prev_score = score;
            igraph_community_leiden(&network_igraph, nullptr, nullptr, resolution, beta, num_leiden_steps == 0, &membership, &num_clusters,
                                    &score);
            num_leiden_steps ++;
        } while(score > prev_score && num_leiden_steps <= max_leiden_steps);

        std::vector<std::vector<SNP_t>> clusters(num_clusters);
        if (boost_mode) {
            auto node_map = boost::get(&node_data::snp, network_boost);
            for (size_t i = 0; i < num_nodes(); i++) {
                size_t group_id = VECTOR(membership)[i];
                auto node_id = node_igraph_map[i];
                auto snp_id = SNP_t(node_map[node_id]);
                clusters[group_id].push_back(snp_id);
                // clusters[group_id].push_back(SNP_t(node_map[node_igraph_map[i]]));
            }
        } else {
            for (size_t i = 0; i < num_nodes(); i++) {
                size_t group_id = VECTOR(membership)[i];
                auto node_id = needl_igraph_map[i];
                auto snp_id = SNP_t(node_id);
                clusters[group_id].push_back(snp_id);
                // clusters[group_id].push_back(SNP_t(node_map[node_igraph_map[i]]));
            }
        }

        igraph_vector_destroy(&membership);

        return clusters;
    }

    bool SNPNetwork_base::contains_node(const SNP_t& snp) {
        if (boost_mode) {
            auto item = node_boost_map.find(snp.value);
            return item != node_boost_map.end();
        } else {
            auto item = needl_adjacency_list.find(snp.value);
            return item != needl_adjacency_list.end();
        }
    }

    long SNPNetwork_base::get_network_diameter() {
        request_igraph_repr();
        igraph_real_t diameter;
        igraph_diameter(&network_igraph, &diameter, nullptr, nullptr, nullptr, false, true);

        long diameter_num = long(diameter);
        if (diameter == IGRAPH_NAN || diameter == IGRAPH_INFINITY) diameter_num = -1;

        return diameter_num;
    }

    bool SNPNetwork_base::is_connected() {
        request_igraph_repr();

        igraph_bool_t is_connected;
        igraph_is_connected(&network_igraph, &is_connected, IGRAPH_WEAK);

        return is_connected;
    }

    size_t SNPNetwork_base::get_degree(const SNP_t& snp) {
        if (boost_mode) {
            auto item = node_boost_map.find(snp.value);
            if (item == node_boost_map.end()) {
                throw epi::Error("Requested SNP is not in the network.");
            }

            return boost::out_degree(item->second, network_boost);
        } else {
            auto item = needl_adjacency_list.find(snp.value);
            if (item == needl_adjacency_list.end()) {
                throw epi::Error("Requested SNP is not in the network.");
            }

            return item->second.size();
        }
    }

    void SNPNetwork_base::remove_edge(const SNPEdge &edge) {
        if (boost_mode) {
            auto snp1_it = node_boost_map.find(edge.bits.snp1);
            auto snp2_it = node_boost_map.find(edge.bits.snp2);
            if (snp1_it == node_boost_map.end() || snp2_it == node_boost_map.end()) {
                throw epi::Error("Requested SNP is not in the network.");
            }

            boost::remove_edge(snp1_it->second, snp2_it->second, network_boost);
            auto coll_it = edge_collision_list.find(edge.val);
            if (coll_it != edge_collision_list.end()) edge_collision_list.erase(coll_it);
        } else {
            auto edge_repr = needl_edge_list.find(edge.val);
            if (edge_repr == needl_edge_list.end()){
                throw epi::Error("Requested edge is not in the network.");
            }
            needl_edge_list.erase(edge_repr);

            auto snp1 = needl_adjacency_list.find(edge.bits.snp1);
            auto snp2 = needl_adjacency_list.find(edge.bits.snp2);

            auto snp2_pos = std::find(snp1->second.begin(), snp1->second.end(), edge.bits.snp2);
            auto snp1_pos = std::find(snp2->second.begin(), snp2->second.end(), edge.bits.snp1);
            snp1->second.erase(snp2_pos);
            snp2->second.erase(snp1_pos);

            --needl_num_edges;
        }

        igraph_valid = false;
    }

    void SNPNetwork_base::remove_node(const SNP_t &snp) {
        if (boost_mode) {
            auto item = node_boost_map.find(snp.value);
            if (item == node_boost_map.end()) {
                throw epi::Error("Requested SNP is not in the network.");
            }

            auto snp_val_map = boost::get(&node_data::snp, network_boost);
            auto it = boost::adjacent_vertices(item->second, network_boost);
            for (; it.first != it.second; it.first++) {
                SNPEdge out_edge(snp, SNP_t(snp_val_map[*it.first]));
                auto it_out = edge_collision_list.find(out_edge.val);
                if (it_out != edge_collision_list.end()) {
                    edge_collision_list.erase(it_out);
                }
            }

            boost::clear_vertex(item->second, network_boost);
            boost::remove_vertex(item->second, network_boost);

            // node boost map was invalidated because of vertex deletion --> rebuild
            node_boost_map.clear();
            auto snp_val_map_new = boost::get(&node_data::snp, network_boost);
            auto all_snps = boost::vertices(network_boost);
            for (; all_snps.first != all_snps.second; all_snps.first++) {
                node_boost_map.insert({snp_val_map_new[*all_snps.first], *all_snps.first});
            }
        } else {
            auto item = needl_adjacency_list.find(snp.value);
            if (item == needl_adjacency_list.end()) {
                throw epi::Error("Requested SNP is not in the network.");
            }

            for (auto adjacent : item->second) {
                auto edge_repr = needl_edge_list.find(SNPEdge(snp, SNP_t(adjacent)).val);
                needl_edge_list.erase(edge_repr);

                auto adjacent_snp = needl_adjacency_list.find(adjacent);
                auto adjacent_pos = std::find(adjacent_snp->second.begin(), adjacent_snp->second.end(), snp.value);
                adjacent_snp->second.erase(adjacent_pos);
                --needl_num_edges;
            }
            needl_adjacency_list.erase(item);
            --needl_num_nodes;
        }

        igraph_valid = false;
    }

    std::vector<SNP_t> SNPNetwork_base::get_adjacent_snps(const std::vector<SNP_t> &snps) {
        std::unordered_set<SNP_t, SNP_t::SNPHash> all_adjacent;
        for (auto & snp : snps) {
            auto adjacent = get_adjacent_snps(snp);
            all_adjacent.insert(adjacent.begin(), adjacent.end());
        }
        std::vector<SNP_t> all_adjacent_vec;
        for (auto & snp : all_adjacent) {
            auto it = std::find(snps.begin(), snps.end(), snp);
            if (it == snps.end()) all_adjacent_vec.push_back(snp);
        }
        return  all_adjacent_vec;
    }

    std::vector<SNP_t> SNPNetwork_base::get_adjacent_snps(const std::unordered_set<SNP_t, SNP_t::SNPHash>& snps) {
        std::unordered_set<SNP_t, SNP_t::SNPHash> all_adjacent;
        for (auto & snp : snps) {
            auto adjacent = get_adjacent_snps(snp);
            all_adjacent.insert(adjacent.begin(), adjacent.end());
        }
        std::vector<SNP_t> all_adjacent_vec;
        for (auto & snp : all_adjacent) {
            auto it = snps.find(snp);
            if (it == snps.end()) all_adjacent_vec.push_back(snp);
        }
        return  all_adjacent_vec;
    }

    std::vector<SNP_t> SNPNetwork_base::get_articulation_points() {
        if (boost_mode) {
            std::vector<vertex_descriptor> articulation_vertices;
            boost::articulation_points(network_boost, std::back_inserter(articulation_vertices));
            auto snp_val_map = boost::get(&node_data::snp, network_boost);

            std::vector<SNP_t> articulation_snps;
            for (auto &v: articulation_vertices) {
                articulation_snps.push_back(SNP_t(snp_val_map[v]));
            }

            return articulation_snps;
        } else {
            throw epi::Error("Articulation points can only be calculated in boost_mode");
            return{};
        }
    }


    void SNPNetwork_base::request_igraph_repr() {
        if (!igraph_valid) {
            if (igraph_created) {
                // need to destroy previous igraph first
                igraph_destroy(&network_igraph);
            }

            // create clone of boost graph
            igraph_empty(&network_igraph, num_nodes(), false);
            if (boost_mode) {
                node_igraph_rev_map.clear();
                node_igraph_map.clear();
                auto v_it = boost::vertices(network_boost);
                igraph_integer_t i = 0;
                for (; v_it.first != v_it.second; v_it.first++) {
                    node_igraph_rev_map.insert({*v_it.first, i++});
                    node_igraph_map.push_back(*v_it.first);
                }
            } else {
                needl_igraph_map.clear();
                needl_igraph_rev_map.clear();

                igraph_integer_t i = 0;
                for (auto & snp : needl_adjacency_list) {
                    needl_igraph_map.push_back(snp.first);
                    needl_igraph_rev_map.insert({ snp.first, i++ });
                }
            }

            igraph_vector_t edge_list;
            igraph_vector_init(&edge_list, num_edges() * 2);
            size_t vector_pos = 0;
            if (boost_mode) {
                auto it = boost::edges(network_boost);
                for (; it.first != it.second; it.first++) {
                    VECTOR(edge_list)[vector_pos++] = node_igraph_rev_map[it.first->m_source];
                    VECTOR(edge_list)[vector_pos++] = node_igraph_rev_map[it.first->m_target];
                }
            } else {
                for (auto & edge : needl_edge_list) {
                    SNPEdge edge_t;
                    edge_t.val = edge.first;
                    VECTOR(edge_list)[vector_pos++] = needl_igraph_rev_map[edge_t.bits.snp1];
                    VECTOR(edge_list)[vector_pos++] = needl_igraph_rev_map[edge_t.bits.snp2];
                }
            }

            igraph_add_edges(&network_igraph, &edge_list, nullptr);
            igraph_vector_destroy(&edge_list);

            igraph_valid = true;
            igraph_created = true;
        }
    }


    SNPNetwork_base::all_iterator::all_iterator(vertex_iter _curr, vertex_snp_map_type _map) {
        curr = _curr;
        map = _map;
        boost_mode = true;
    }

    SNPNetwork_base::all_iterator::all_iterator(std::unordered_map<SNP, std::vector<SNP>>::iterator _needl_curr) {
        boost_mode = false;
        needl_curr = _needl_curr;
    }

    SNPNetwork_base::all_iterator& SNPNetwork_base::all_iterator::operator++() {
        if (boost_mode) curr++;
        else needl_curr ++;
        return *this;
    }

    SNPNetwork_base::all_iterator SNPNetwork_base::all_iterator::operator++(int) {
        all_iterator retval = *this;
        ++(*this);
        return retval;
    }


    /*
    SNPNetwork_base::all_iterator& SNPNetwork_base::all_iterator::operator--() {
        if (boost_mode) curr--;
        else needl_curr--;
        return *this;
    }

    SNPNetwork_base::all_iterator SNPNetwork_base::all_iterator::operator--(int) {
        all_iterator retval = *this;
        --(*this);
        return retval;
    }
     */


    bool SNPNetwork_base::all_iterator::operator!=(SNPNetwork_base::all_iterator other) const {
        if (boost_mode) return curr != other.curr;
        else return needl_curr != other.needl_curr;
    }

    bool SNPNetwork_base::all_iterator::operator==(SNPNetwork_base::all_iterator other) const  {
        if (boost_mode) return curr == other.curr;
        else return needl_curr == other.needl_curr;
    }


    SNP_t SNPNetwork_base::all_iterator::operator*() const {
        if (boost_mode) return SNP_t(map[*curr]);
        else return SNP_t(needl_curr->first);
    }

    SNPNetwork_base::all_iterable::all_iterable(vertex_iter __begin, vertex_iter __end,
                                           vertex_snp_map_type __map) {
        _begin = __begin;
        _end = __end;
        _map = __map;
        _boost_mode = true;
    }

    SNPNetwork_base::all_iterable::all_iterable(std::unordered_map<SNP, std::vector<SNP>>::iterator __begin, std::unordered_map<SNP, std::vector<SNP>>::iterator __end) {
        _needl_begin = __begin;
        _needl_end = __end;
        _boost_mode = false;
    }

    SNPNetwork_base::all_iterator SNPNetwork_base::all_iterable::begin() {
        if (_boost_mode)
            return all_iterator(_begin, _map);
        else
            return all_iterator(_needl_begin);
    }

    SNPNetwork_base::all_iterator SNPNetwork_base::all_iterable::end() {
        if (_boost_mode)
            return all_iterator(_end, _map);
        else
            return all_iterator(_needl_end);
    }

    SNPNetwork_base::all_iterable SNPNetwork_base::all() {
        if (boost_mode) {
            auto vertices = boost::vertices(network_boost);
            auto snp_val_map = boost::get(&node_data::snp, network_boost);
            return {vertices.first, vertices.second, snp_val_map};
        } else {
            return { needl_adjacency_list.begin(), needl_adjacency_list.end() };
        }
    }

    void SNPNetwork_base::replace_nodes(const std::vector<std::pair<SNP_t, SNP_t>> &pairs) {
        igraph_valid = false;

        if (boost_mode) {
            auto node_map = boost::get(&node_data::snp, network_boost);
            for (auto &pair: pairs) {
                // it is not necessary that both of the snps already exist in the network
                auto item1 = node_boost_map.find(pair.first.value);
                auto item2 = node_boost_map.find(pair.second.value);

                bool item1_in_net = item1 != node_boost_map.end();
                bool item2_in_net = item2 != node_boost_map.end();

                vertex_descriptor snp1_desc, snp2_desc;

                if (item1_in_net) {
                    // change snp information of that specific node
                    node_map[item1->second] = pair.second.value;
                    snp1_desc = item1->second;
                }

                if (item2_in_net) {
                    // change snp information of that specific node
                    node_map[item2->second] = pair.first.value;
                    snp2_desc = item2->second;
                }

                // replace entries in node_boost_map afterwards to not invalidate iterators
                if (item1_in_net) node_boost_map.erase(pair.first.value);
                if (item2_in_net) node_boost_map.erase(pair.second.value);

                if (item1_in_net) node_boost_map.insert({pair.second.value, snp1_desc});
                if (item2_in_net) node_boost_map.insert({pair.first.value, snp2_desc});
            }

            // regenerate the edge collision list
            edge_collision_list.clear();

            auto edges = boost::edges(network_boost);
            for (; edges.first != edges.second; edges.first++) {
                SNP_t snp1(node_map[edges.first->m_source]);
                SNP_t snp2(node_map[edges.first->m_target]);
                SNPEdge edge(snp1, snp2);
                edge_collision_list.insert(edge.val);
            }
        } else {
            for (auto &pair: pairs) {
                auto num_edges_before = num_edges();
                // it is not necessary that both of the snps already exist in the network
                auto item1 = needl_adjacency_list.find(pair.first.value);
                auto item2 = needl_adjacency_list.find(pair.second.value);

                bool item1_in_net = item1 != needl_adjacency_list.end();
                bool item2_in_net = item2 != needl_adjacency_list.end();

                if (!item1_in_net) {
                    add_node(pair.first);
                    item1 = needl_adjacency_list.find(pair.first.value);
                }
                if (!item2_in_net) {
                    add_node(pair.second);
                    item2 = needl_adjacency_list.find(pair.second.value);
                }

                // get previous adjacency data
                std::vector<SNP> adjacent1, adjacent2;

                if (item1_in_net) {
                    adjacent1 = { item1->second.begin(), item1->second.end() };
                }

                if (item2_in_net) {
                    adjacent2 = { item2->second.begin(), item2->second.end() };
                }


                // remove both nodes
                remove_node(pair.first);
                remove_node(pair.second);

                // add back nodes
                if (item2_in_net) {
                    add_node(pair.first);
                }
                if (item1_in_net) {
                    add_node(pair.second);
                }

                // insert edges
                if (item1_in_net) {
                    for (auto adj : adjacent1) {
                        if (adj == pair.second.value) adj = pair.first.value;
                        SNPEdge new_edge (SNP_t(pair.second.value), SNP_t(adj));
                        add_edge(new_edge);
                    }
                }

                if (item2_in_net) {
                    for (auto adj : adjacent2) {
                        if (adj == pair.first.value) adj = pair.second.value;
                        SNPEdge new_edge (SNP_t(pair.first.value), SNP_t(adj));
                        add_edge(new_edge);
                    }

                }

                if (num_edges() != num_edges_before) {
                    std::cout << "Error!!!"<< item1_in_net << ", " << item2_in_net << std::endl;
                }
            }
        }
    }

    bool SNPNetwork_base::edge_exists(const SNPEdge &edge) {
        if (boost_mode) {
            return edge_collision_list.find(edge.val) != edge_collision_list.end();
        } else {
            return needl_edge_list.find(edge.val) != needl_edge_list.end();
        }
    }

    std::set<std::string> SNPNetwork_base::get_edge_labels(const SNPEdge &edge) {
        uint64_t labels_ids;
        if (boost_mode) {
            auto snp1_it = node_boost_map.find(edge.bits.snp1);
            auto snp2_it = node_boost_map.find(edge.bits.snp2);
            if (snp1_it == node_boost_map.end() || snp2_it == node_boost_map.end()) {
                throw epi::Error("Requested SNP is not in the network.");
            }

            auto edge_val_map = boost::get(&edge_data::labels, network_boost);
            auto boost_edge = boost::edge(snp1_it->second, snp2_it->second, network_boost).first;

            labels_ids = edge_val_map[boost_edge];
        } else {
           auto edge_repr = needl_edge_list.find(edge.val);
           if (edge_repr == needl_edge_list.end()) {
               throw epi::Error("Requested edge is not in the network.");
           }
           labels_ids = edge_repr->second.labels;
        }

        std::set<std::string> labels;
        for (uint64_t i = 0; i < edge_labels.size(); i++) {
            uint64_t bitmask = 1 << i;
            if (bitmask & labels_ids) labels.insert(edge_labels[i]);
        }

        return labels;
    }

    std::set<uint16_t> SNPNetwork_base::get_edge_label_ids(const SNPEdge &edge) {
        uint64_t labels_ids;
        if (boost_mode) {
            auto snp1_it = node_boost_map.find(edge.bits.snp1);
            auto snp2_it = node_boost_map.find(edge.bits.snp2);
            if (snp1_it == node_boost_map.end() || snp2_it == node_boost_map.end()) {
                throw epi::Error("Requested SNP is not in the network.");
            }

            auto edge_val_map = boost::get(&edge_data::labels, network_boost);
            auto boost_edge = boost::edge(snp1_it->second, snp2_it->second, network_boost).first;

            labels_ids = edge_val_map[boost_edge];
        } else {
            auto edge_repr = needl_edge_list.find(edge.val);
            if (edge_repr == needl_edge_list.end()) {
                throw epi::Error("Requested edge is not in the network.");
            }
            labels_ids = edge_repr->second.labels;
        }

        std::set<uint16_t> labels;
        for (uint64_t i = 0; i < edge_labels.size(); i++) {
            uint64_t bitmask = 1 << i;
            if (bitmask & labels_ids) labels.insert(i);
        }

        return labels;
    }

    uint16_t SNPNetwork_base::get_edge_tag_id(const std::string& edge_tag) {
        auto item = edge_label_map.find(edge_tag);
        if (item != edge_label_map.end()) return item->second;

        // check that we did not reach max number of labels
        edge_label_map.insert({ edge_tag, next_edge_label_id });
        edge_labels.push_back(edge_tag);

        if (edge_labels.size() >= 64) {
            throw epi::Error("Tried to insert more edge labels than allowed. Currently, only a total of 64 edge labels is allowed.");
        }
        return next_edge_label_id ++;
    }

    std::vector<std::string> SNPNetwork_base::get_all_edge_tags() {
        return { edge_labels.begin(), edge_labels.end() };
    }

    void SNPNetwork_base::clear_edges() {
        if (boost_mode) {
            edge_collision_list.clear();
            auto ptrs = boost::vertices(network_boost);
            for (; ptrs.first != ptrs.second; ++ptrs.first) {
                boost::clear_vertex(*ptrs.first, network_boost);
            }
        } else {
            needl_edge_list.clear();
            needl_num_edges = 0;
            for (auto &snp : needl_adjacency_list) {
                snp.second.clear();
            }
        }
        igraph_valid = false;
    }

    SNPNetwork::SNPNetwork(const SNPNetwork &other) : SNPNetwork_base(other) {
        igraph_valid = false;
        igraph_created = false;
    }

    SNPNetwork::SNPNetwork(bool boost_mode) : SNPNetwork_base(boost_mode) {}
} // epi