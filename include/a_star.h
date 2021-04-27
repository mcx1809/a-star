#ifndef A_STAR_H_
#define A_STAR_H_

#include <boost/intrusive/set.hpp>
#include <boost/pool/object_pool.hpp>
#include <cstdint>
#include <functional>
#include <iterator>
#include <memory>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

template <typename Node>
class DynamicNodePool {
 public:
  DynamicNodePool() { Reset(); }

  void SetMaxSize(std::size_t size) { max_nodes_size_ = size; }

  Node* Get() {
    if (nodes_size_ < max_nodes_size_) {
      if (nodes_size_ == nodes_size_threshold_) {
        nodes_->set_next_size(2 * nodes_->get_next_size());
        nodes_size_threshold_ += nodes_->get_next_size();
      }
      nodes_size_++;

      auto* n = nodes_->construct();
      n->id = nodes_size_;

      return n;
    } else {
      return nullptr;
    }
  }

  void Reset() {
    nodes_.reset(new boost::object_pool<Node>());

    auto n = 4 * 1024 / sizeof(Node);
    n = n ? n : 32;
    nodes_->set_next_size(n);
    nodes_size_ = 0;
    nodes_size_threshold_ = nodes_->get_next_size();
  }

 private:
  std::unique_ptr<boost::object_pool<Node>> nodes_;
  std::size_t nodes_size_;
  std::size_t nodes_size_threshold_;
  std::size_t max_nodes_size_ = 0;
};

template <typename Node>
class PreAllocNodePool {
 public:
  PreAllocNodePool() { Reset(); }

  void SetMaxSize(std::size_t size) { nodes_.resize(size); }

  Node* Get() {
    if (nodes_size_ < nodes_.size()) {
      auto* n = &nodes_[nodes_size_++];
      n->id = nodes_size_;

      return n;
    } else {
      return nullptr;
    }
  }

  void Reset() { nodes_size_ = 0; }

 private:
  std::vector<Node> nodes_;
  std::size_t nodes_size_;
};

template <typename Vertex, typename CostScalar, bool DYNAMIC_NODE_POOL = true>
class AStar {
  struct NodeMapTag;

  using NodeMapBaseHook = boost::intrusive::set_base_hook<
      boost::intrusive::tag<NodeMapTag>,
      boost::intrusive::link_mode<boost::intrusive::normal_link>>;

  struct CostMapTag;

  using CostMapBaseHook = boost::intrusive::set_base_hook<
      boost::intrusive::tag<CostMapTag>,
      boost::intrusive::link_mode<boost::intrusive::normal_link>>;

  struct Node;

  struct NodeMapKeyOfValue {
    using type = Vertex;

    const type& operator()(const Node& n) const { return n.v; }
  };

  using NodeMap =
      boost::intrusive::set<Node, boost::intrusive::base_hook<NodeMapBaseHook>,
                            boost::intrusive::key_of_value<NodeMapKeyOfValue>>;

  using CostMap =
      boost::intrusive::set<Node, boost::intrusive::base_hook<CostMapBaseHook>>;

  struct Node : public NodeMapBaseHook, public CostMapBaseHook {
    bool operator<(const Node& o) const {
      auto f = g + h;
      auto of = o.g + o.h;
      return f < of ? true : (of < f ? false : id < o.id);
    }

    std::size_t id;

    Node* parent;

    Vertex v;
    CostScalar g;
    CostScalar h;
  };

 public:
  class Iterator
      : public std::iterator<std::input_iterator_tag,
                             std::tuple<Vertex, CostScalar, CostScalar>> {
    friend class AStar;

    Iterator(const Node* n = nullptr) { n_ = n; }

   public:
    Iterator& operator=(const Iterator&) = default;

    bool operator!=(const Iterator& o) const { return n_ != o.n_; }

    bool operator==(const Iterator& o) const { return n_ == o.n_; }

    Iterator& operator++() {
      n_ = n_->parent;
      return *this;
    }

    Iterator& operator++(int) {
      n_ = n_->parent;
      return *this;
    }

    std::tuple<const Vertex*, CostScalar, CostScalar> operator*() const {
      return std::make_tuple(&n_->v, n_->g, n_->h);
    }

   private:
    const Node* n_;
  };

 public:
  AStar() {
    on_each_n_ = [this](const Vertex& v, bool recalc_h) {
      if (close_list_.find(v) == close_list_.end()) {
        auto it = open_list_.find(v);
        if (it != open_list_.end()) {
          auto n = &*it;
          auto g = cur_->g + calc_d_(cur_->v, v);

          if (calc_h_) {
            auto h = recalc_h ? calc_h_(v) : n->h;
            if (g + h < n->g + n->h) {
              cost_list_.erase(cost_list_.iterator_to(*n));
              n->parent = cur_;
              n->v = v;
              n->g = g;
              n->h = h;
              cost_list_.insert(*n);
            }
          } else {
            if (g < n->g) {
              cost_list_.erase(cost_list_.iterator_to(*n));
              n->parent = cur_;
              n->v = v;
              n->g = g;
              cost_list_.insert(*n);
            }
          }
        } else {
          AddOpenList(cur_, v, cur_->g + calc_d_(cur_->v, v));
        }
      }
    };
  }

  void SetHeuristicFn(const std::function<CostScalar(const Vertex&)>& fn) {
    calc_h_ = fn;
  }

  void SetHeuristicFn(std::function<CostScalar(const Vertex&)>&& fn) {
    calc_h_ = std::move(fn);
  }

  void SetDistanceFn(
      const std::function<CostScalar(const Vertex&, const Vertex&)>& fn) {
    calc_d_ = fn;
  }

  void SetDistanceFn(
      std::function<CostScalar(const Vertex&, const Vertex&)>&& fn) {
    calc_d_ = std::move(fn);
  }

  void SetNeighborsFn(
      const std::function<
          void(const Vertex&, const std::function<void(const Vertex&, bool)>&)>&
          fn) {
    get_n_ = fn;
  }

  void SetNeighborsFn(
      std::function<void(const Vertex&,
                         const std::function<void(const Vertex&, bool)>&)>&&
          fn) {
    get_n_ = std::move(fn);
  }

  void SetMaxNodes(std::size_t nodes) {
    Clear();
    node_pool_.SetMaxSize(nodes);
  }

  template <typename GoalFunc>
  CostScalar CalcPathCostWithFn(const Vertex& start,
                                const GoalFunc& is_goal_fn) {
    Clear();

    AddOpenList(nullptr, start, CostScalar(0));

    while (true) {
      cur_ = OpenListFront();
      if (cur_) {
        if (!is_goal_fn(cur_->v)) {
          RemoveOpenList(cur_);
          AddCloseList(cur_);
          GetNeighbors();
        } else {
          break;
        }
      } else {
        return CostScalar(0);
      }
    }

    return cur_->g;
  }

  CostScalar CalcPathCost(const Vertex& start, const Vertex& goal) {
    return CalcPathCostWithFn(
        start, [&](const Vertex& v) { return !(v < goal || goal < v); });
  }

  CostScalar ReCalcPathCost(const Vertex& goal) {
    if (reversed_) {
      Reverse();
    }

    auto it = close_list_.find(goal);
    if (it != close_list_.end()) {
      cur_ = &*it;
    } else {
      if (calc_h_) {
        RefreshOpenList();
      }

      while (true) {
        cur_ = OpenListFront();
        if (cur_) {
          if (cur_->v < goal || goal < cur_->v) {
            RemoveOpenList(cur_);
            AddCloseList(cur_);
            GetNeighbors();
          } else {
            break;
          }
        } else {
          return CostScalar(0);
        }
      }
    }

    return cur_->g;
  }

  template <typename GoalFunc>
  std::tuple<std::size_t, CostScalar> SearchWithFn(const Vertex& start,
                                                   const GoalFunc& is_goal_fn) {
    CalcPathCostWithFn(start, is_goal_fn);
    return cur_ ? Reverse() : std::make_tuple(std::size_t(0), CostScalar(0));
  }

  std::tuple<std::size_t, CostScalar> Search(const Vertex& start,
                                             const Vertex& goal) {
    return SearchWithFn(
        start, [&](const Vertex& v) { return !(v < goal || goal < v); });
  }

  Iterator begin() { return Iterator(cur_); }

  Iterator end() { return Iterator(); }

 private:
  void Clear() {
    open_list_.clear();
    cost_list_.clear();
    close_list_.clear();

    node_pool_.Reset();

    cur_ = nullptr;
  }

  void AddOpenList(Node* parent, const Vertex& v, CostScalar g) {
    auto* n = node_pool_.Get();

    if (n) {
      n->parent = parent;
      n->v = v;
      n->g = g;
      if (calc_h_) {
        n->h = calc_h_(v);
      }

      AddOpenList(n);
    }
  }

  void AddOpenList(Node* n) {
    open_list_.insert(*n);
    cost_list_.insert(*n);
  }

  Node* OpenListFront() {
    auto it = cost_list_.begin();
    if (it != cost_list_.end()) {
      return &*it;
    } else {
      return nullptr;
    }
  }

  void RemoveOpenList(Node* n) {
    cost_list_.erase(cost_list_.iterator_to(*n));
    open_list_.erase(open_list_.iterator_to(*n));
  }

  void RefreshOpenList() {
    cost_list_.clear();

    for (auto& n : open_list_) {
      n.h = calc_h_(n.v);
      cost_list_.insert(n);
    }
  }

  void AddCloseList(Node* n) { close_list_.insert(*n); }

  void GetNeighbors() {
    if (get_n_) {
      get_n_(cur_->v, on_each_n_);
    }
  }

  std::tuple<std::size_t, CostScalar> Reverse() {
    std::size_t r = 1;
    auto g = cur_->g;

    for (Node *prev = nullptr, *cur = cur_;; ++r) {
      auto* next = cur->parent;

      cur->parent = prev;
      prev = cur;
      if (next) {
        cur = next;
      } else {
        cur_ = cur;
        break;
      }
    }

    reversed_ = !reversed_;

    return std::make_tuple(r, g);
  }

 private:
  std::function<CostScalar(const Vertex&)> calc_h_;
  std::function<CostScalar(const Vertex&, const Vertex&)> calc_d_;
  std::function<void(const Vertex&,
                     const std::function<void(const Vertex&, bool)>&)>
      get_n_;
  std::function<void(const Vertex&, bool)> on_each_n_;

  typename std::conditional<DYNAMIC_NODE_POOL, DynamicNodePool<Node>,
                            PreAllocNodePool<Node>>::type node_pool_;

  NodeMap open_list_;
  CostMap cost_list_;
  NodeMap close_list_;

  Node* cur_ = nullptr;
  bool reversed_ = false;
};

#endif
