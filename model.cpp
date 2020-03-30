#include <cstdio>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <vector>
#include <sstream>
#include <algorithm>

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

struct Message {
  /* how fast the task can produce this message (B/s)
   */
  double injection;
  /* how large this message is (B)
   */
  int64_t size;

  Task dst; // destination of message
};

struct Task_s {
  int64_t size;              // task size in memory
  std::vector<Message> msgs; // messages this task wants to send
};

struct PE_s {
  std::string name;
  int64_t memory;
  double injection;
};

struct Link_s {
  double bandwidth;
};

struct Application {
  std::vector<Task> tasks;
};

typedef std::vector<Link> Path;

std::string path_str(std::vector<Link> path) {
  (void) path;
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

inline Message make_message(int64_t size, Task dst) {
  Message ret = {
    std::numeric_limits<double>::infinity(),
    size,
    dst,
  };
  return ret;
}

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

inline Link make_nvlink() {
  auto ret = std::make_shared<Link_s>();
  ret->bandwidth = 75 * Gi;
  return ret;
}

inline Link make_xbus() {
  auto ret = std::make_shared<Link_s>();
  ret->bandwidth = 32 * Gi;
  return ret;
}

inline Machine make_newell() {
  auto ret = std::make_shared<Machine_s>();
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


  // Undirectional Links
  Link g0g1 = make_nvlink();
  Link g0c0 = make_nvlink();
  Link g1c0 = make_nvlink();
  Link g1g0 = make_nvlink();
  Link g2g3 = make_nvlink();
  Link g2c1 = make_nvlink();
  Link g3c1 = make_nvlink();
  Link g3g2 = make_nvlink();
  Link c0c1 = make_xbus();
  Link c0g0 = make_nvlink();
  Link c0g1 = make_nvlink();
  Link c1c0 = make_xbus();
  Link c1g2 = make_nvlink();
  Link c1g3 = make_nvlink();




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
  ret->topo[gpu3][gpu2].push_back(g3c1);
  ret->topo[gpu3][gpu2].push_back(c1c0);
  ret->topo[gpu3][gpu2].push_back(c0g0);

  // 3 -> 1
  ret->topo[gpu3][gpu3].push_back(g3c1);
  ret->topo[gpu3][gpu3].push_back(c1c0);
  ret->topo[gpu3][gpu3].push_back(c0g1);

  // 3 -> 2
  ret->topo[gpu3][gpu2].push_back(g3g2);

  return ret;
};

/* provide the estimated time for machine m to run application a, with map[ti] =
pi mapping task ti to PE pi
 */
inline double model(const Machine &m, const Application &a,
                    std::map<Task, PE> map) {

  std::map<PE, int64_t> peComm;     // how many bytes each PE will send/recv
  std::map<Link, int64_t> linkComm; // bytes on each link

  /* Accumulate total bytes in/out of each PE
     Accumulate total bytes sent on each link
  */
  for (const auto &srcTask : a.tasks) {
    const auto &srcPE = map[srcTask];
    for (const auto &msg : srcTask->msgs) {
      const auto &dstTask = msg.dst;
      const auto &dstPE = map[dstTask];

      fprintf(stderr, "Task %p -> %p PE %s -> %s (%ld B) \n", srcTask.get(),
              dstTask.get(), srcPE->name.c_str(), dstPE->name.c_str(),
              msg.size);

      // add message size to PE
      auto ret = peComm.insert(std::make_pair(srcPE, msg.size));
      if (false == ret.second) {
        ret.first->second += msg.size;
      }
      ret = peComm.insert(std::make_pair(dstPE, msg.size));
      if (false == ret.second) {
        ret.first->second += msg.size;
      }

      // add message size to all involved links
      auto &links = m->topo[srcPE][dstPE];
      for (const auto &link : links) {
        auto ret = linkComm.insert(std::make_pair(link, msg.size));
        if (false == ret.second) {
          ret.first->second += msg.size;
        }
      }
    }
  }

  double time = -1;

  // limit the time according to the link bandwidth
  for (auto kv : linkComm) {
    const Link &link = kv.first;
    const int64_t bytes = kv.second;
    double limitTime = model_div(bytes, link->bandwidth);
    fprintf(stderr, "link has %ld bytes @ %e B/s = %es\n", bytes,
            link->bandwidth, limitTime);
    if (limitTime > time) {
      time = limitTime;
    }
  }

  // limit the time according to the message injections
  for (const auto &srcTask : a.tasks) {
    for (const auto &msg : srcTask->msgs) {
      double limitTime = model_div(msg.size, msg.injection);
      if (limitTime > time) {
        fprintf(stderr, "msg has %ld bytes @ %e B/s\n", msg.size,
                msg.injection);
        time = limitTime;
      }
    }
  }

  // limit the time according to the PE injection
  for (auto kv : peComm) {
    const PE &pe = kv.first;
    const int64_t bytes = kv.second;
    double limitTime = model_div(bytes, pe->injection);
    if (limitTime > time) {
      fprintf(stderr, "pe has %ld bytes @ %e B/s\n", bytes, pe->injection);
      time = limitTime;
    }
  }

  return time;
}

int main() {
  Machine m = make_newell();

  Application a;
  // four stencil tasks
  a.tasks.push_back(std::make_shared<Task_s>());
  a.tasks.push_back(std::make_shared<Task_s>());
  a.tasks.push_back(std::make_shared<Task_s>());
  a.tasks.push_back(std::make_shared<Task_s>());

  // three messages per task
  a.tasks[0]->msgs.push_back(make_message(512 * 512, a.tasks[1]));
  a.tasks[0]->msgs.push_back(make_message(512 * 512, a.tasks[2]));
  a.tasks[0]->msgs.push_back(make_message(512, a.tasks[3]));
  a.tasks[1]->msgs.push_back(make_message(512 * 512, a.tasks[0]));
  a.tasks[1]->msgs.push_back(make_message(512, a.tasks[2]));
  a.tasks[1]->msgs.push_back(make_message(512 * 512, a.tasks[3]));
  a.tasks[2]->msgs.push_back(make_message(512 * 512, a.tasks[0]));
  a.tasks[2]->msgs.push_back(make_message(512, a.tasks[1]));
  a.tasks[2]->msgs.push_back(make_message(512 * 512, a.tasks[3]));
  a.tasks[3]->msgs.push_back(make_message(512, a.tasks[0]));
  a.tasks[3]->msgs.push_back(make_message(512 * 512, a.tasks[1]));
  a.tasks[3]->msgs.push_back(make_message(512 * 512, a.tasks[2]));

  // map tasks to PEs
  std::map<Task, PE> map;
  map[a.tasks[0]] = m->pes[0];
  map[a.tasks[1]] = m->pes[1];
  map[a.tasks[2]] = m->pes[2];
  map[a.tasks[3]] = m->pes[3];

  double time = model(m, a, map);

  fprintf(stderr, "%e\n", time);

  return 0;
}