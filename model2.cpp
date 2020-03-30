#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstring>
#include <deque>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <vector>

constexpr int64_t Ki = 1024;
constexpr int64_t Mi = Ki * Ki;
constexpr int64_t Gi = Mi * Ki;

inline double signum(double val) { return (0.0 < val) - (val < 0.0); }

inline double model_div(double n, double d) {
  double ret;
  if (n == 0.0 && d == 0.0) {
    ret = std::numeric_limits<double>::quiet_NaN();
  } else if (0.0 == std::abs(n)) {
    ret = n;
  } else if (0.0 == std::abs(d)) {
    ret = signum(n) * std::numeric_limits<double>::infinity();
  } else {
    ret = n / d;
  }
  return ret;
}

struct Task_s;
typedef std::shared_ptr<Task_s> Task;
struct PE_s;
typedef std::shared_ptr<PE_s> PE;
struct Link_s;
typedef std::shared_ptr<Link_s> Link;
struct Machine_s;
typedef std::shared_ptr<Machine_s> Machine;

enum class Op {
  Unknown,
  Send,
  WaitSends, // wait for all my outstanding sends
  Barrier,   // wait until everyone has arrived
};

// ideas for future semantics
#if 0
  WaitFor, // wait for another guy
  Signal,  // signal another Program to stop waiting

#endif

struct SendOperands {
  int64_t dst;
  int64_t size;
};

struct Instruction {
  Op code;
  union operands {
    SendOperands send;
  } ops;

  Instruction() : code(Op::Unknown) { std::memset(&ops, 0, sizeof(ops)); }

  static Instruction send(int64_t dst, int64_t size) {
    Instruction ret;
    ret.code = Op::Send;
    ret.ops.send.dst = dst;
    ret.ops.send.size = size;
    return ret;
  }

  static Instruction wait_sends() {
    Instruction ret;
    ret.code = Op::WaitSends;
    return ret;
  }

  static Instruction barrier() {
    Instruction ret;
    ret.code = Op::Barrier;
    return ret;
  };
};

class Program {
private:
  int64_t id_;
  std::vector<Instruction> instructions_;

public:
  Program(int64_t i) : id_(i) {}

  void push_back(const Instruction &i) { instructions_.push_back(i); }
  size_t size() const noexcept { return instructions_.size(); }
  const std::vector<Instruction> &instructions() const { return instructions_; }
};

typedef std::vector<Program> Application;

std::vector<Program> write_exchange() {
  std::vector<Program> programs;
  for (int64_t i = 0; i < 4; ++i) {
    programs.push_back(Program(i));
  }
  for (int64_t i = 0; i < 4; ++i) {
    programs[i].push_back(Instruction::send((i + 1) % 4, 10000));
    programs[i].push_back(Instruction::wait_sends());
    programs[i].push_back(Instruction::barrier());
    programs[i].push_back(Instruction::send((i + 2) % 4, 10000));
    programs[i].push_back(Instruction::wait_sends());
  }
  return programs;
}

struct PE_s {
  std::string name;
  int64_t memory;
  double injection;
};

struct Link_s {
  std::string name;
  double bandwidth;

  std::string str() const {
    return name + "[" + std::to_string(bandwidth) + "]";
  }
};

typedef std::vector<Link> Path;

std::string path_str(std::vector<Link> path) {
  (void)path;
  std::stringstream ss;
  return ss.str();
}

struct Machine_s {
  // the PEs
  std::vector<PE> pes;
  // the links connecting each pe
  std::map<PE, std::map<PE, Path>> topo;

  // add link l to the path between a and b
  void join(PE a, PE b, Link l) {
    insert(a);
    insert(b);
    topo[a][b].push_back(l);
  }

  void insert(PE p) {
    if (pes.end() == std::find(pes.begin(), pes.end(), p)) {
      pes.push_back(p);
    }
  }
};

inline PE make_V100(const std::string &name) {
  auto ret = std::make_shared<PE_s>();
  ret->name = name;
  ret->memory = 16 * Gi;
  ret->injection = std::numeric_limits<double>::infinity();
  return ret;
}

inline PE make_POWER9(const std::string &name) {
  auto ret = std::make_shared<PE_s>();
  ret->name = name;
  ret->memory = 0;
  ret->injection = std::numeric_limits<double>::infinity();
  return ret;
}

inline Link make_nvlink(const std::string &name) {
  auto ret = std::make_shared<Link_s>();
  ret->name = name;
  ret->bandwidth = 75 * Gi;
  return ret;
}

inline Link make_xbus(const std::string &name) {
  auto ret = std::make_shared<Link_s>();
  ret->name = name;
  ret->bandwidth = 32 * Gi;
  return ret;
}

inline Machine make_newell() {
  Machine ret = std::make_shared<Machine_s>();
  // PEs
  PE gpu0 = make_V100("gpu0");
  PE gpu1 = make_V100("gpu1");
  PE gpu2 = make_V100("gpu2");
  PE gpu3 = make_V100("gpu3");
  PE cpu0 = make_POWER9("cpu0");
  PE cpu1 = make_POWER9("cpu1");
  ret->pes.push_back(gpu0);
  ret->pes.push_back(gpu1);
  ret->pes.push_back(gpu2);
  ret->pes.push_back(gpu3);
  ret->pes.push_back(cpu0);
  ret->pes.push_back(cpu1);
  assert(ret->pes[0]);

  // Undirectional Links
  Link g0g1 = make_nvlink("g0g1");
  Link g0c0 = make_nvlink("g0c0");
  Link g1c0 = make_nvlink("g1c0");
  Link g1g0 = make_nvlink("g1g0");
  Link g2g3 = make_nvlink("g2g3");
  Link g2c1 = make_nvlink("g2c1");
  Link g3c1 = make_nvlink("g3c1");
  Link g3g2 = make_nvlink("g3g2");
  Link c0c1 = make_xbus("c0c1");
  Link c0g0 = make_nvlink("c0g0");
  Link c0g1 = make_nvlink("c0g1");
  Link c1c0 = make_xbus("c1c0");
  Link c1g2 = make_nvlink("c1g2");
  Link c1g3 = make_nvlink("c1g3");

  // 0 -> 1
  ret->topo[gpu0][gpu1].push_back(g0g1);

  // 0 -> 2
  ret->topo[gpu0][gpu2].push_back(g0c0);
  ret->topo[gpu0][gpu2].push_back(c0c1);
  ret->topo[gpu0][gpu2].push_back(c1g2);

  // 0 -> 3
  ret->topo[gpu0][gpu3].push_back(g0c0);
  ret->topo[gpu0][gpu3].push_back(c0c1);
  ret->topo[gpu0][gpu3].push_back(c1g3);

  // 1 -> 0
  ret->topo[gpu1][gpu0].push_back(g1g0);

  // 1 -> 2
  ret->topo[gpu1][gpu2].push_back(g1c0);
  ret->topo[gpu1][gpu2].push_back(c0c1);
  ret->topo[gpu1][gpu2].push_back(c1g2);

  // 1 -> 3
  ret->topo[gpu1][gpu3].push_back(g1c0);
  ret->topo[gpu1][gpu3].push_back(c0c1);
  ret->topo[gpu1][gpu3].push_back(c1g3);

  // 2 -> 0
  ret->topo[gpu2][gpu0].push_back(g2c1);
  ret->topo[gpu2][gpu0].push_back(c1c0);
  ret->topo[gpu2][gpu0].push_back(c0g0);

  // 2 -> 1
  ret->topo[gpu2][gpu1].push_back(g2c1);
  ret->topo[gpu2][gpu1].push_back(c1c0);
  ret->topo[gpu2][gpu1].push_back(c0g1);

  // 2 -> 3
  ret->topo[gpu2][gpu3].push_back(g2g3);

  // 3 -> 0
  ret->topo[gpu3][gpu0].push_back(g3c1);
  ret->topo[gpu3][gpu0].push_back(c1c0);
  ret->topo[gpu3][gpu0].push_back(c0g0);

  // 3 -> 1
  ret->topo[gpu3][gpu1].push_back(g3c1);
  ret->topo[gpu3][gpu1].push_back(c1c0);
  ret->topo[gpu3][gpu1].push_back(c0g1);

  // 3 -> 2
  ret->topo[gpu3][gpu2].push_back(g3g2);

  return ret;
};

struct EventData {
  double time;
};

struct Event {
  std::shared_ptr<EventData> data;
  bool operator<(const Event &rhs) const noexcept {
    return data->time > rhs.data->time;
  }
  bool operator==(const Event &rhs) const noexcept { return data == rhs.data; }

  Event(double time) {
    data = std::make_shared<EventData>();
    data->time = time;
  }
};

/* representing transmission of `size` at `start`
 */
struct tx_s {
  int64_t size;            // bytes left to transfer
  std::vector<Link> links; // the links involved
  Event event;             // the event representing the end of this transfer

  tx_s(int64_t size, const std::vector<Link> &links, Event event)
      : size(size), links(links), event(event) {}

  std::string str() const {
    return std::to_string((uintptr_t)this) + "[" + std::to_string(size) + "]";
  }
};
typedef std::shared_ptr<tx_s> Tx;

struct Executor {

  enum class State {
    Ready,
    WaitEvents,
    WaitBarrier,
    Done, // program is done
  };

  State state;
  int64_t id;
  Program program;
  int64_t pc; // program counter

  // data for various states
  std::vector<Event> events;

  std::vector<Tx> txs; // active transmissions

  Executor(int64_t i, const Program &program) : id(i), program(program), pc(0) {
    state = State::Ready;
  }

  /* fill i with next instruction.
  return true if success, return false if done
  */
  bool next(Instruction &i) {
    if (pc < int64_t(program.size())) {
      i = program.instructions()[pc++];
      return true;
    } else {
      return false;
    }
  }

  void notify(const Event &event) {
    if (State::WaitEvents == state) {
      auto &waitEvents = events;
      auto it = std::find(waitEvents.begin(), waitEvents.end(), event);
      if (it != waitEvents.end()) {
        std::cerr << "notify(): no longer waiting for event " << it->data
                  << "\n";
        waitEvents.erase(it);
      }
    }
  }
};

class Model {

private:
  /* event heap
   */
  std::vector<Event> events_;

  /* Program executors
   */
  std::vector<Executor> executors_;

  std::map<size_t, PE> map_;

  /* All active transmissions
   */
  std::set<Tx> activity;

  /* Which transmissions are active on each link
   */
  std::map<Link, std::set<Tx>> linkActivity;

  double time_;

  Machine machine_;

  /* add a new event
   */
  void add_event(const Event &e) {
    std::cerr << "add_event(): " << e.data << "\n";
    events_.push_back(e);
    std::push_heap(events_.begin(), events_.end());
  }

  /* pop the next event
   */
  Event pop_event() {
    std::pop_heap(events_.begin(), events_.end());
    Event ret = events_.back();
    events_.pop_back();
    return ret;
  }

public:
  Model(const Machine &m, const Application &a, std::map<size_t, PE> map)
      : map_(map), time_(0), machine_(m) {
    for (size_t i = 0; i < map.size(); ++i) {
      executors_.push_back(Executor(i, a[i]));
    }
  }

  /* Run the executors until they are stalled
   */
  void run_executors() {
    std::cerr << "Model::run_executors()\n";

    // Release barrier if everyone has entered
    bool releaseBarrier = true;
    for (size_t ei = 0; ei < executors_.size(); ++ei) {
      Executor &exec = executors_[ei];
      if (Executor::State::WaitBarrier != exec.state) {
        releaseBarrier = false;
        break;
      }
    }
    if (releaseBarrier) {
      std::cerr << "Model::run_executors(): release barrier\n";
      for (auto &exec : executors_) {
        exec.state = Executor::State::Ready;
      }
    }

    // check if anyone is waiting for sends that are already done
    for (size_t ei = 0; ei < executors_.size(); ++ei) {
      Executor &exec = executors_[ei];
      if (Executor::State::WaitEvents == exec.state) {
        if (exec.events.empty()) {
          std::cerr << "Model::run_executors(): " << ei
                    << " WaitEvents -> Ready\n";
          exec.state = Executor::State::Ready;
        }
      }
    }

    for (size_t ei = 0; ei < executors_.size(); ++ei) {
      Executor &exec = executors_[ei];

      while (Executor::State::Ready == exec.state) {
        Instruction instr;
        if (!exec.next(instr)) {
          exec.state = Executor::State::Done;
          continue;
        }

        switch (instr.code) {
        case Op::Send: {
          std::cerr << "Model::run_executors(): " << ei << " op=send\n";
          PE &srcPE = map_[ei];
          assert(srcPE);
          auto &dstPE = map_[instr.ops.send.dst];
          assert(dstPE);
          std::cerr << srcPE->name << " -> " << dstPE->name << "\n";
          Path links = machine_->topo[srcPE][dstPE];
          if (srcPE != dstPE) {
            assert(!links.empty());
          }
          std::cerr << srcPE->name << " -> " << dstPE->name
                    << "(links=" << links.size() << ")\n";
          Event event(0);
          Tx tx = std::make_shared<tx_s>(instr.ops.send.size, links, event);
          add_transmission(tx);
          exec.txs.push_back(tx);
          break;
        }
        case Op::WaitSends: {
          std::cerr << "Model::run_executors(): " << ei << " op=WaitSend\n";
          // wait on all active sends
          exec.state = Executor::State::WaitEvents;
          for (auto &tx : exec.txs) {
            std::cerr << "Model::run_executors(): " << ei << " wait on tx "
                      << tx->str() << "\n";
            exec.events.push_back(tx->event);
          }
          break;
        }

        case Op::Barrier:
          std::cerr << "Model::run_executors(): " << ei << " op=barrier\n";
          exec.state = Executor::State::WaitBarrier;
          break;
        default:
          assert(0);
          break;
        }
      }
    }

    std::cerr << "states\n";
    for (auto &exec : executors_) {
      std::cerr << int(exec.state) << "\n";
    }
  }

  /* real bandwidth for a transmission
   */
  double tx_bandwidth(Tx tx) {
    std::cerr << "Model::tx_bandwidth(): tx=" << tx->str() << " links\n";
    for (auto &link : tx->links) {
      std::cerr << link->str() << "\n";
    }

    double bw = std::numeric_limits<double>::infinity();
    for (auto &link : tx->links) {
      bw = std::min(bw, link->bandwidth / linkActivity[link].size());
    }
    return bw;
  }

  void add_transmission(Tx tx) {
    std::cerr << "Model::add_transmission()\n";
    add_event(tx->event);
    activity.insert(tx);
    for (auto &link : tx->links) {
      linkActivity[link].insert(tx);
    }
    std::cerr << "linkActivity.size() = " << linkActivity.size() << "\n";

    // update all transmission end times thanks to the new transmission
    update_transmissions();
  }

  void update_transmissions() {
    std::cerr << "Model::update_transmissions()\n";

    /* Go through each transmission, compute its bandwidth under contention, and
     * update
     */
    for (auto &tx : activity) {
      tx->event.data->time = time_ + tx->size / tx_bandwidth(tx);
    }

    // since times have changed, remake the heap
    std::make_heap(events_.begin(), events_.end());

    std::cerr << "Model::update_transmissions: events\n";
    for (auto &event : events_) {
      std::cerr << event.data << " " << event.data->time << "\n";
    }
  }

  /* Advance the simulation by a single event
   */
  void next_event() {
    std::cerr << "Model::next_event():\n";

    assert(!events_.empty());
    const Event event = pop_event();

    std::cerr << "Model::next_event(): @ " << event.data->time << "\n";

    // if event was a transmission, remove it from active transmissions
    for (auto &tx : activity) {
      if (tx->event == event) {
        for (auto &link : tx->links) {
          linkActivity[link].erase(tx);
        }
        activity.erase(tx);
        break;
      }
    }

    const double elapsed = event.data->time - time_;
    // update data remaining on all active transmissions
    for (auto &tx : activity) {
      const double bw = tx_bandwidth(tx);
      std::cerr << "Model::next_event(): tx->size=" << tx->size << " bw=" << bw
                << "\n";

      tx->size -= bw * elapsed;
      assert(tx->size >= 0);
    }

    // notify Executors of event
    for (auto &e : executors_) {
      e.notify(event);
    }

    // update simulation time
    time_ = event.data->time;
    std::cerr << "Model::next_event(): sim @" << time_ << "\n";
  }

  /* Advance simulation by finishing all events
   */
  void drain_events() {
    std::cerr << "Model::drain_events():\n";
    while (!events_.empty()) {
      next_event();
    }
  }

  /* When all programs are finished, return true
   */
  bool programs_done() const {
    for (auto &exec : executors_) {
      if (Executor::State::Done != exec.state) {
        return false;
      }
    }
    return true;
  }

  double run() {
    std::cerr << "Model::run()\n";

    while (!programs_done()) {
      std::cerr << "TIME=" << time_ << "\n";
      run_executors();

      /* executors may not create any events
       */
      if (!events_.empty()) {
        next_event();
      }
    }

    drain_events();

    return time_;
  }
};

int main() {
  Machine m = make_newell();

  Application a = write_exchange();

  // map tasks to PEs
  std::map<size_t, PE> map;
  assert(m->pes[0]);
  map[0] = m->pes[0];
  assert(m->pes[1]);
  map[1] = m->pes[1];
  assert(m->pes[2]);
  map[2] = m->pes[2];
  assert(m->pes[3]);
  map[3] = m->pes[3];

  Model model(m, a, map);

  double time = model.run();

  fprintf(stderr, "time = %e\n", time);

  return 0;
}