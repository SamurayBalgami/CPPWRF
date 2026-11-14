// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// The transport between the workers of a decomposed run (core/tiles.cpp): whole messages, in order per (sender,
// receiver) pair. Tiles: in-process mailboxes (tiles.cpp); MPI ranks: core/tiles_mpi.cpp.
#include <functional>
#include <stdexcept>
#include <vector>

#include "core/tiles.hpp"

namespace cppwrf::tile::detail {

struct Aborted : std::runtime_error {
  Aborted() : std::runtime_error("another worker stopped") {}
};

struct Comm {
  virtual ~Comm() = default;
  // never blocks on the receiver (the message is queued); the next message of this worker to `to`
  virtual void send(int to, std::vector<char>&& m) = 0;
  // waits for the next message from `from`; Aborted when another worker stopped
  virtual std::vector<char> recv(int from) = 0;
  // every queued message has left this worker's buffers (after each collective)
  virtual void flush() {}
};

// The calling thread as worker `id` of n = o.px * o.py over the transport c: body(id) with the decomposition's hooks.
int run_worker(Comm& c, int id, const Options& o, const std::function<int(int)>& body);

}  // namespace cppwrf::tile::detail
