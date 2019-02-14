#ifndef FSA_CONTIG_GRAPH_HPP
#define FSA_CONTIG_GRAPH_HPP

#include <array>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <list>
#include <deque>

#include "utility.hpp"

#include "contig_link_store.hpp"
#include "sequence.hpp"

class ContigEdge;
class ContigNode {
public:
    using ID = Seq::EndId;

public:
    ContigNode(ID id) : id_(id) {}

    ID Id() const { return id_; }
    static ID ReverseId(ID id)  { return Seq::ReverseEndId(id); }
    size_t InDegree() const { return in_edges_.size(); }
    size_t OutDegree() const { return out_edges_.size(); }

    ID id_;
    std::vector<ContigEdge*> in_edges_;
    std::vector<ContigEdge*> out_edges_;

    ContigEdge* best_in_edge_{ nullptr };
    ContigEdge* best_out_edge_{ nullptr };
};

class ContigEdge {
public:
    using Hash = ArrayHash<int,2>;
    using Compare = ArrayCompare<int,2>;
    using ID = std::array<ContigNode::ID, 2>;   

public:
    ContigEdge(ContigNode* in, ContigNode* out) : in_node_(in), out_node_(out) {}

    ID Id() const { return std::array<int, 2>{in_node_->id_, out_node_->id_}; }
    static ID ReverseId(ID id)  {
        return std::array<int, 2>{ContigNode::ReverseId(id[1]), ContigNode::ReverseId(id[0]) };
    }

    std::vector<Seq::Area> GetSeqArea();
    size_t LinkLength() { return link_->TotalLength(); }

    size_t Support() { return link_->groups.size() >= 1 ? link_->groups[0].size() : 0; }

    ContigNode* in_node_;
    ContigNode* out_node_;

    ContigLink* link_;
};

class ContigGraph {
public:
    ContigGraph(ContigLinkStore &links);
    ~ContigGraph();

    void Create();

    void AddLink(ContigLink& link);
    void AddEdge(Seq::EndId s, Seq::EndId t, ContigLink& link);

    void CheckRepeat();

    void IdentifyPaths(const std::string &method);
    const std::unordered_set<Seq::Id>& GetContained() const { return contained_; }

    std::deque<ContigNode*> ExtendPath(ContigNode* e, std::unordered_set<ContigNode*> &visitied, const std::string &method);

    template<typename I, typename O>
    std::deque<ContigNode*> ExtendPathWithMethod(ContigNode* e, std::unordered_set<ContigNode*> &visited, I get_in_node, O get_out_node);

    ContigEdge* GetEdge(ContigNode* in, ContigNode* out) {
        for (ContigEdge* e : in->out_edges_) {
            if (e->out_node_ == out) return e;
        }
        return nullptr;
    }

    ContigEdge* ReverseEdge(ContigEdge* e) { return edges_[ContigEdge::ReverseId(e->Id())]; }
    ContigNode* ReverseNode(ContigNode* n) { return nodes_[ContigNode::ReverseId(n->Id())]; }

    void CalucateBest(const std::string &method);
    void Output(const std::string &fname);

    const std::list<std::deque<ContigNode*>>& GetPaths() { return paths_; }
    std::string ConstructContig(const std::list<ContigEdge*> &path);

protected:
    ContigLinkStore &links_;
    std::unordered_map<int, ContigNode*> nodes_;
    std::unordered_map<ContigEdge::ID, ContigEdge*, ContigEdge::Hash, ContigEdge::Compare> edges_;
    
    std::list<std::deque<ContigNode*>> paths_;
    std::unordered_set<Seq::Id> contained_;
};

#endif // FSA_CONTIG_GRAPH_HPP  
