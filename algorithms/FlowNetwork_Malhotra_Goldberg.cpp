#include <iostream>
#include <stdexcept>
#include <cassert>
#include <limits.h>
#include <optional>
#include <type_traits>
#include <vector>
#include <queue>


//Flow Network - addEdge(from, to, cap)
//Malhotra/Goldberg(network), getNetwork()

namespace NFlow{

template<typename TFlow>
class TNetwork {
private:
    struct TEdge_;

public:
    typedef unsigned int TVertex;
    typedef unsigned int TVertexNumber;
    typedef unsigned int TEdgeNum;

    class TEdgeIterator {
    friend class TNetwork;

    public:
        TFlow getFlow() const {
            return getEdge_().flow;
        }

        TFlow getCapacity() const {
            return getEdge_().capacity;
        }

        TFlow getResudialCapacity() const {
            return getCapacity() - getFlow();
        }

        TVertex getFinish() const {
            return getEdge_().finish;
        }

        void pushFlow(TFlow flow_value) {
            const auto edge_num = network_->graph_[vertex_][edge_num_];
            auto& edges_ = network_->edges_;
            if (edges_[edge_num].flow + flow_value > edges_[edge_num].capacity) {
                throw std::logic_error("Edge's flow is bigger than capacity");
            }
            edges_[edge_num].flow     += flow_value;
            edges_[edge_num ^ 1].flow -= flow_value;
        }

        TEdgeIterator& operator++() {
            if (edge_num_ < network_->graph_[vertex_].size()) {
                ++edge_num_;
            }
            return *this;
        }

        bool isEnd() const {
            return edge_num_ == network_->graph_[vertex_].size();
        }

    private:
        typedef unsigned int TEdgeNum_;

        TNetwork* network_;
        TVertex   vertex_;
        TEdgeNum_ edge_num_;

        TEdgeIterator(TNetwork* network, TVertex vertex) :
            network_(network),
            vertex_(vertex),
            edge_num_(0)
        {}

        const TEdge_& getEdge_() const {
            if (isEnd()) {
                throw std::out_of_range("Iterator out of range");
            }
            const auto edge_num = network_->graph_[vertex_][edge_num_];
            return network_->edges_[edge_num];
        }
    };

    TNetwork(TVertexNumber vertex_number, TVertex source, TVertex sink) :
        vertex_number_(vertex_number),
        source_(source),
        sink_(sink)
    {
        if (source >= vertex_number || sink   >= vertex_number) {
            throw std::out_of_range("Source or sink index is too large");
        }
        if (source == sink) {
            throw std::logic_error("Source and sink are the same");
        }
        graph_.resize(vertex_number_);
    }

    void addEdge(TVertex start, TVertex finish, TFlow capacity) {
        // add forward edge
        graph_[start].push_back(edges_.size());
        edges_.emplace_back(finish, /* flow = */ 0, capacity);
        // add backward edge
        graph_[finish].push_back(edges_.size());
        edges_.emplace_back(start,  /* flow = */ 0, /* capacity = */ 0);
    }

    TEdgeIterator getEdgeIterator(TVertex vertex) {
        return TEdgeIterator(this, vertex);
    }

    TVertexNumber getVertexNumber() const {
        return vertex_number_;
    }

    TVertex getSource() const {
        return source_;
    }

    TVertex getSink() const {
        return sink_;
    }

    TFlow getFlowValue() const {
        TFlow flow = 0;

        for (auto edge_num : graph_[source_]) {
            const auto& edge = edges_[edge_num];
            flow += edge.flow;
        }

        return flow;
    }

private:
    struct TEdge_ {
        TVertex finish;
        TFlow   flow;
        TFlow   capacity;

        TEdge_(TVertex finish, TFlow flow, TFlow capacity) :
            finish(finish),
            flow(flow),
            capacity(capacity)
        {}
    };

    std::vector< std::vector<TEdgeNum> > graph_;
    std::vector<TEdge_> edges_;
    TVertex vertex_number_;
    TVertex source_;
    TVertex sink_;
};

} // end of namespace NFlow


namespace NMalhotra {

template<typename TFlow>
class TMalhotra {
public:
    typedef NFlow::TNetwork<TFlow> TNetwork;

    TMalhotra(const TNetwork& network) :
        network_(network)
    {
        const auto vertex_number = network.getVertexNumber();
        incoming_potential_.resize(vertex_number);
        outcoming_potential_.resize(vertex_number);
        is_available_.resize(vertex_number);
        graph_.resize(vertex_number);
        reversed_graph_.resize(vertex_number);

        findMaxFlow_();
    }

    const TNetwork& getNetwork() const {
        return network_;
    }

private:
    typedef typename TNetwork::TVertex       TVertex_;
    typedef typename TNetwork::TVertexNumber TVertexNumber_;
    typedef typename TNetwork::TEdgeNum      TEdgeNum_;
    typedef typename TNetwork::TEdgeIterator TEdgeIterator_;
    typedef unsigned int                TDist_;
    typedef std::make_unsigned_t<TFlow> TPotential_;

    struct Edge_ {
        TVertex_       finish;
        TEdgeIterator_ network_edge;
        Edge_(TVertex_ finish, TEdgeIterator_ network_edge) :
            finish(finish),
            network_edge(network_edge)
        {}
    };

    TNetwork network_;
    std::vector<TPotential_>          incoming_potential_;
    std::vector<TPotential_>          outcoming_potential_;
    std::vector<bool>                 is_available_;
    std::vector< std::vector<Edge_> > graph_;
    std::vector< std::vector<Edge_> > reversed_graph_;

    TPotential_ getPotential_(TVertex_ vertex) {
        if (vertex == network_.getSource()) {
            return outcoming_potential_[vertex];
        }
        if (vertex == network_.getSink()) {
            return incoming_potential_[vertex];
        }
        return std::min(incoming_potential_[vertex], outcoming_potential_[vertex]);
    }

    TVertex_ getMinPotentialVertex_() {
        TVertex_ min_potential_vertex = network_.getSource();
        for (TVertexNumber_ vertex = 0; vertex < network_.getVertexNumber(); ++vertex) {
            if (is_available_[vertex] && getPotential_(vertex) < getPotential_(min_potential_vertex)) {
                min_potential_vertex = vertex;
            }
        }
        return min_potential_vertex;
    }

    void removeZeroPotentialVertex_(TVertex_ vertex) {
        is_available_[vertex] = false;
        for (const auto edge : graph_[vertex]) {
            incoming_potential_[edge.finish] -= edge.network_edge.getResudialCapacity();
        }
        for (const auto edge : reversed_graph_[vertex]) {
            outcoming_potential_[edge.finish] -= edge.network_edge.getResudialCapacity();
        }
    }

    void findMaxFlow_() {
        while(build_graph_()) {
            removeIncorrectEdges_();
            calcPotential_();
            const auto source = network_.getSource();
            const auto sink   = network_.getSink();
            while(std::min(getPotential_(source), getPotential_(sink)) > 0) {
                const auto min_potential_vertex = getMinPotentialVertex_();
                if (getPotential_(min_potential_vertex) == 0) {
                    removeZeroPotentialVertex_(min_potential_vertex);
                } else {
                    pushFlow_(min_potential_vertex);
                }
            }
        }
    }

    bool build_graph_() {
        const auto INF           = std::numeric_limits<TDist_>::max();
        const auto vertex_number = network_.getVertexNumber();
        const auto source        = network_.getSource();
        const auto sink          = network_.getSink();
        for (TVertexNumber_ vertex = 0; vertex < vertex_number; ++vertex) {
            is_available_[vertex] = false;
            graph_[vertex].clear();
            reversed_graph_[vertex].clear();
        }
        std::vector<TDist_> dist(vertex_number, INF);
        dist[source] = 0;

        std::queue<TVertex_> queue;
        queue.push(source);

        while(!queue.empty()) {
            const auto cur_vertex = queue.front();
            queue.pop();
            for (auto it = network_.getEdgeIterator(cur_vertex); !it.isEnd(); ++it) {
                const auto cur_finish = it.getFinish();
                if (it.getResudialCapacity() > 0){
                    if (dist[cur_finish] == INF) {
                        dist[cur_finish] = dist[cur_vertex] + 1;
                        queue.push(cur_finish);
                    }

                    if (dist[cur_finish] == dist[cur_vertex] + 1) {
                        graph_[cur_vertex].emplace_back(cur_finish, it);
                        reversed_graph_[cur_finish].emplace_back(cur_vertex, it);
                    }
                }
            }
        }

        if (dist[sink] == INF) {
            return false;
        }

        dist.assign(vertex_number, INF);
        dist[sink] = 0;
        queue.push(sink);

        while(!queue.empty()) {
            const auto cur_vertex = queue.front();
            queue.pop();
            for (const auto& edge : reversed_graph_[cur_vertex]) {
                if (dist[edge.finish] == INF) {
                    dist[edge.finish] = dist[cur_vertex] + 1;
                    queue.push(edge.finish);
                }
            }
        }

        for (TVertexNumber_ vertex = 0; vertex < vertex_number; ++vertex) {
            is_available_[vertex] = dist[vertex] != INF;
        }

        return true;
    }

    void removeIncorrectEdges_() {
        for (TVertexNumber_ vertex = 0; vertex < network_.getVertexNumber(); ++vertex) {
            if (!is_available_[vertex]) {
                graph_[vertex].clear();
                reversed_graph_[vertex].clear();
            } else {
                TEdgeNum_ edge_num = 0;
                auto& graph = graph_[vertex];
                while(edge_num < graph.size()) {
                    if (!is_available_[graph[edge_num].finish]) {
                        std::swap(graph[edge_num], graph.back());
                        graph.pop_back();
                    } else {
                        ++edge_num;
                    }
                }

                auto& reversed_graph = reversed_graph_[vertex];
                while(edge_num < reversed_graph.size()) {
                    if (!is_available_[reversed_graph[edge_num].finish]) {
                        std::swap(reversed_graph[edge_num], reversed_graph.back());
                        reversed_graph.pop_back();
                    } else {
                        ++edge_num;
                    }
                }
            }
        }
    }

    void calcPotential_() {
        for (TVertexNumber_ vertex = 0; vertex < network_.getVertexNumber(); ++vertex) {
            incoming_potential_[vertex]  = 0;
            outcoming_potential_[vertex] = 0;
            for (const auto& edge : reversed_graph_[vertex]) {
                incoming_potential_[vertex] += edge.network_edge.getResudialCapacity();
            }
            for (const auto& edge : graph_[vertex]) {
                outcoming_potential_[vertex] += edge.network_edge.getResudialCapacity();
            }
        }
    }

    void pushFlow_(TVertex_ min_potential_vertex) {
        TFlow flow_value = getPotential_(min_potential_vertex);
        std::queue< std::pair<TVertex_, TFlow> > queue;
        queue.push({min_potential_vertex, flow_value});

        while(!queue.empty()) {
            auto [cur_vertex, flow] = queue.front();
            queue.pop();
            if (cur_vertex == network_.getSink()) {
                continue;
            }
            auto& graph = graph_[cur_vertex];
            while(flow) {
                auto& cur_edge = graph.back();
                if (!is_available_[cur_edge.finish] || cur_edge.network_edge.getResudialCapacity() == 0) {
                    graph.pop_back();
                } else {
                    TFlow cur_flow = std::min(flow, cur_edge.network_edge.getResudialCapacity());
                    cur_edge.network_edge.pushFlow(cur_flow);
                    outcoming_potential_[cur_vertex]     -= cur_flow;
                    incoming_potential_[cur_edge.finish] -= cur_flow;
                    flow                                 -= cur_flow;
                    queue.push({cur_edge.finish, cur_flow});
                }
            }
        }

        queue.push({min_potential_vertex, flow_value});

        while(!queue.empty()) {
            auto [cur_vertex, flow] = queue.front();
            queue.pop();
            if (cur_vertex == network_.getSource()) {
                continue;
            }
            auto& graph = reversed_graph_[cur_vertex];
            while(flow) {
                auto& cur_edge = graph.back();
                if (!is_available_[cur_edge.finish] || cur_edge.network_edge.getResudialCapacity() == 0) {
                    graph.pop_back();
                } else {
                    TFlow cur_flow = std::min(flow, cur_edge.network_edge.getResudialCapacity());
                    cur_edge.network_edge.pushFlow(cur_flow);
                    incoming_potential_[cur_vertex]       -= cur_flow;
                    outcoming_potential_[cur_edge.finish] -= cur_flow;
                    flow                                  -= cur_flow;
                    queue.push({cur_edge.finish, cur_flow});
                }
            }
        }
    }
};

} // end of namespace NMalhotra


namespace NGoldberg {

template<typename TFlow>
class TGoldberg {
public:
    typedef NFlow::TNetwork<TFlow> TNetwork;

    TGoldberg(const TNetwork& network) :
        network_(network)
    {
        const auto vertex_number = network.getVertexNumber();
        height_.resize(vertex_number);
        potential_.resize(vertex_number);
        for (TVertexNumber_ vertex = 0; vertex < vertex_number; ++vertex) {
            edge_iterator_.push_back(network_.getEdgeIterator(vertex));
        }
        findMaxFlow_();
    }

    const TNetwork& getNetwork() const {
        return network_;
    }

private:
    typedef typename TNetwork::TVertex       TVertex_;
    typedef typename TNetwork::TVertexNumber TVertexNumber_;
    typedef typename TNetwork::TEdgeIterator TEdgeIterator_;
    typedef unsigned int THeight_;
    typedef std::make_unsigned_t<TFlow> TPotential_;

    TNetwork network_;
    std::vector<THeight_>       height_;
    std::vector<TPotential_>    potential_;
    std::vector<TEdgeIterator_> edge_iterator_;
    std::queue<TVertex_>        overflowed_vertexes_;

    void pushFlow(TVertex_ vertex, TEdgeIterator_ edge) {
        const TFlow flow = std::min(potential_[vertex], (TPotential_)edge.getResudialCapacity());
        const auto source = network_.getSource();
        const auto sink   = network_.getSink();
        if (vertex != source && vertex != sink) {
            potential_[vertex] -= flow;
        }
        const auto finish = edge.getFinish();
        if (finish != source && finish != sink) {
            potential_[finish] += flow;
        }
        edge.pushFlow(flow);
    }

    void relabel(TVertex_ vertex) {
        THeight_ new_height = std::numeric_limits<THeight_>::max();
        for (auto it = network_.getEdgeIterator(vertex); !it.isEnd(); ++it) {
            if (it.getResudialCapacity() > 0) {
                new_height = std::min(new_height, height_[it.getFinish()] + 1);
            }
        }
        height_[vertex] = new_height;
    }

    void discharge(TVertex_ vertex) {
        auto& edge = edge_iterator_[vertex];
        while(potential_[vertex] > 0) {
            if (edge.isEnd()) {
                relabel(vertex);
                edge = network_.getEdgeIterator(vertex);
            } else {
                const TVertex_ finish = edge.getFinish();
                if (edge.getResudialCapacity() > 0 && height_[vertex] == height_[finish] + 1) {
                    const bool was_overflowed = potential_[finish] > 0;
                    pushFlow(vertex, edge);
                    if (!was_overflowed && potential_[finish] > 0) {
                        overflowed_vertexes_.push(finish);
                    }
                } else {
                    ++edge;
                }
            }
        }
    }

    void findMaxFlow_() {
        for (TVertexNumber_ vertex = 0; vertex < network_.getVertexNumber(); ++vertex) {
            potential_[vertex] = 0;
            height_[vertex]    = 0;
        }

        const auto source = network_.getSource();
        height_[source] = network_.getVertexNumber();

        for (auto it = network_.getEdgeIterator(source); !it.isEnd(); ++it) {
            const auto cur_finish = it.getFinish();
            const auto sink       = network_.getSink();
            const auto flow       = it.getResudialCapacity();
            it.pushFlow(flow);
            if (cur_finish != sink) {
                potential_[cur_finish] += flow;
            }

            if (potential_[cur_finish] > 0) {
                overflowed_vertexes_.push(cur_finish);
            }
        }

        while(!overflowed_vertexes_.empty()) {
            const auto cur_vertex = overflowed_vertexes_.front();
            overflowed_vertexes_.pop();
            discharge(cur_vertex);
        }
    }
};

} // end of namespace NGoldberg


int main(){
    int n;
    std::cin >> n;
    std::vector<int> cost(n);
    for (int i = 0; i < n; i++) {
        std::cin >> cost[i];
    }

    const int INF = std::numeric_limits<int>::max();
    const unsigned int source = n;
    const unsigned int sink = n + 1;
    NFlow::TNetwork<int> network(n + 2, source, sink);

    for (int vertex = 0; vertex < n; vertex++) {
        int cnt;
        std::cin >> cnt;
        while(cnt--) {
            int parent;
            std::cin >> parent;
            parent--;
            network.addEdge(vertex, parent, INF);
        }
    }

    int result = 0;

    for (int vertex = 0; vertex < n; vertex++) {
        if (cost[vertex] > 0) {
            result += cost[vertex];
            network.addEdge(source, vertex, cost[vertex]);
        } else {
            network.addEdge(vertex, sink, -cost[vertex]);
        }
    }

    NMalhotra::TMalhotra malhotra(network);
    const auto malhotra_result_network = malhotra.getNetwork();

    NGoldberg::TGoldberg goldberg(network);
    const auto goldberg_result_network = goldberg.getNetwork();

    assert(malhotra_result_network.getFlowValue() == goldberg_result_network.getFlowValue());

    std::cout << result - malhotra_result_network.getFlowValue();
}
