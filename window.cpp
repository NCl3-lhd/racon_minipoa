/*!
 * @file window.cpp
 *
 * @brief Window class source file
 */

#include <algorithm>
#include <iostream>

#include "window.hpp"

namespace racon {

std::shared_ptr<Window> createWindow(uint64_t id, uint32_t rank, WindowType type,
    const char* backbone, uint32_t backbone_length, const char* quality,
    uint32_t quality_length) {

    if (backbone_length == 0 || backbone_length != quality_length) {
        fprintf(stderr, "[racon::createWindow] error: "
            "empty backbone sequence/unequal quality length!\n");
        exit(1);
    }

    return std::shared_ptr<Window>(new Window(id, rank, type, backbone,
        backbone_length, quality, quality_length));
}

Window::Window(uint64_t id, uint32_t rank, WindowType type, const char* backbone,
    uint32_t backbone_length, const char* quality, uint32_t quality_length)
        : id_(id), rank_(rank), type_(type), consensus_(), sequences_(),
        qualities_(), positions_() {

    sequences_.emplace_back(backbone, backbone_length);
    qualities_.emplace_back(quality, quality_length);
    positions_.emplace_back(0, 0);
}

Window::~Window() {
}

void Window::add_layer(const char* sequence, uint32_t sequence_length,
    const char* quality, uint32_t quality_length, uint32_t begin, uint32_t end) {

    if (sequence_length == 0 || begin == end) {
        return;
    }

    if (quality != nullptr && sequence_length != quality_length) {
        fprintf(stderr, "[racon::Window::add_layer] error: "
            "unequal quality size!\n");
        exit(1);
    }
    if (begin >= end || begin > sequences_.front().second || end > sequences_.front().second) {
        fprintf(stderr, "[racon::Window::add_layer] error: "
            "layer begin and end positions are invalid!\n");
        exit(1);
    }

    sequences_.emplace_back(sequence, sequence_length);
    qualities_.emplace_back(quality, quality_length);
    positions_.emplace_back(begin, end);
}

bool Window::generate_consensus(para_t* para, bool trim) {

    if (sequences_.size() < 3) {
        consensus_ = std::string(sequences_.front().first, sequences_.front().second);
        return false;
    }

    graph* DAG = new graph();
    DAG->init(para);
    aligned_buff_t* mpool = new aligned_buff_t;
    
    std::vector<res_t> res = alignment(para, DAG, nullptr, 0, sequences_.front().first, sequences_.front().second, mpool);
    DAG->add_path(para->m, 0, res, 1);
    DAG->topsort(para, 0);

    std::vector<uint32_t> rank;
    rank.reserve(sequences_.size());
    for (uint32_t i = 0; i < sequences_.size(); ++i) {
        rank.emplace_back(i);
    }

    std::sort(rank.begin() + 1, rank.end(), [&](uint32_t lhs, uint32_t rhs) {
        return positions_[lhs].first < positions_[rhs].first; });
    // perform minipoa
    uint32_t offset = 0.01 * sequences_.front().second;
    for (uint32_t j = 1; j < sequences_.size(); ++j) {
        uint32_t i = rank[j];

        // spoa::Alignment alignment;
        std::vector<res_t> res;
        int sink_id = -1;
        if (positions_[i].first < offset && positions_[i].second >
            sequences_.front().second - offset) {
            res = poa(para, DAG, 0, 1, i, sequences_[i].first, sequences_[i].second, mpool, para->ab_band);
            // alignment = alignment_engine->align(sequences_[i].first,
            //     sequences_[i].second, graph);
            sink_id = 1;
        }
        else {

            int beg_id = positions_[i].first + 2;
            int end_id = positions_[i].second + 1 == sequences_.front().second ? 1 : positions_[i].second + 1 + 2;
            res = poa(para, DAG, beg_id, end_id, i, sequences_[i].first + 1, sequences_[i].second - 1, mpool, para->ab_band);
            if(end_id == 1) sink_id = 1;
            // std::vector<int32_t> mapping;
            // auto subgraph = graph->subgraph(positions_[i].first,
            //     positions_[i].second, mapping);
            // alignment = alignment_engine->align(sequences_[i].first,
            //     sequences_[i].second, subgraph);
            // subgraph->update_alignment(alignment, mapping);
        }
        
        // add res to graph
        DAG->add_path(para->m, i, res, sink_id);
        DAG->topsort(para, 0);
        // if (qualities_[i].first == nullptr) {
        //     graph->add_alignment(alignment, sequences_[i].first,
        //         sequences_[i].second);
        // }
        // else {
        //     graph->add_alignment(alignment, sequences_[i].first,
        //         sequences_[i].second, qualities_[i].first,
        //         qualities_[i].second);
        // }
    }

    DAG->build_consensus(true);
    consensus_ = DAG->cons;
    std::vector<int> coverages = DAG->coverages;

    if (type_ == WindowType::kTGS && trim) {
        uint32_t average_coverage = (sequences_.size() - 1) / 2;

        int32_t begin = 0, end = consensus_.size() - 1;
        for (; begin < static_cast<int32_t>(consensus_.size()); ++begin) {
            if (coverages[begin] >= average_coverage) {
                break;
            }
        }
        for (; end >= 0; --end) {
            if (coverages[end] >= average_coverage) {
                break;
            }
        }

        if (begin >= end) {
            fprintf(stderr, "[racon::Window::generate_consensus] warning: "
                "contig %lu might be chimeric in window %u!\n", id_, rank_);
        }
        else {
            consensus_ = consensus_.substr(begin, end - begin + 1);
        }
    }
    delete mpool;
    mpool = nullptr;  // 防止后续误用
    delete DAG;
    DAG = nullptr;  // 防止后续误用
    return true;
}

}
