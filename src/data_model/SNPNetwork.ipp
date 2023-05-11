namespace epi {
    template<class Iterator>
    void SNPNetwork_base::set_nodes(const Iterator _begin, const Iterator _end) {
        clear();

        add_nodes(_begin, _end);
    }

    template<class Iterator>
    void SNPNetwork_base::add_nodes(const Iterator _begin, const Iterator _end) {
        for (auto it = _begin; it != _end; it++) {
            add_node(*it);
        }
    }

    template<class Iterator>
    void SNPNetwork_base::remove_nodes(const Iterator _begin, const Iterator _end) {
        for (auto it = _begin; it != _end; it++) {
            remove_node(*it);
        }
    }

    template<class Iterator>
    void SNPNetwork_base::add_edges(const Iterator _begin, const Iterator _end, std::string edge_tag) {
        auto tid = get_edge_tag_id(edge_tag);
        for (auto it = _begin; it != _end; it++) {
                add_edge(*it, tid);
        }
    }

    template<class Iterator>
    void SNPNetwork_base::remove_edges(const Iterator _begin, const Iterator _end) {
        for (auto it = _begin; it != _end; it++) {
            remove_edge(*it);
        }
    }


    template<class Iterator>
    std::vector<SNP_t> SNPNetwork_base::get_adjacent_snps(const Iterator _begin, const Iterator _end) {
        std::vector<SNP_t> snps;
        std::unordered_set<SNP_t, SNP_t::SNPHash> all_adjacent;
        for (auto it = _begin; it != _end; it++) {
            snps.push_back(*it);
            auto adjacent = get_adjacent_snps(*it);
            all_adjacent.insert(adjacent.begin(), adjacent.end());
        }

        std::vector<SNP_t> all_adjacent_vec;
        for (auto & snp : all_adjacent) {
            auto it = std::find(snps.begin(), snps.end(), snp);
            if (it == snps.end()) all_adjacent_vec.push_back(snp);
        }
        return  all_adjacent_vec;
    }
}