/*
 * Regression test for OpenBucketHashtable::erase.
 *
 * Guards against erase() leaving entries unreachable: find() stopping at the
 * new hole before reaching them, after which at() inserts a duplicate key.
 * Checks every invariant against std::unordered_map after each operation.
 */
#include <cstdint>
#include <cstdio>
#include <iterator>
#include <random>
#include <unordered_map>

#include "../../open_bucket_hashtable.h"

template <typename GID>
static bool consistent(const OpenBucketHashtable<GID, uint32_t>& h, const std::unordered_map<GID, uint32_t>& ref) {
   if (h.size() != ref.size()) {
      return false;
   }
   for (const auto& [key, val] : ref) {
      const auto it = h.find(key);
      if (it == h.end() || it->second != val) {
         return false;
      }
   }
   size_t n = 0;
   for (auto it = h.begin(); it != h.end(); ++it, ++n) {
      if (ref.count(it->first) == 0) {
         return false;
      }
   }
   return n == ref.size();
}

// Home bucket of a key in a fresh (16-bucket) table.
static size_t homeOf(uint32_t k) {
   OpenBucketHashtable<uint32_t, uint32_t> t;
   t.at(k) = 0;
   return t.begin().getIndex();
}

// Smallest key >= start, not in used, whose home bucket is target.
static uint32_t keyWithHome(size_t target, uint32_t start, const std::unordered_map<uint32_t, uint32_t>& used) {
   uint32_t k = start;
   while (homeOf(k) != target || used.count(k)) {
      ++k;
   }
   return k;
}

// Two keys sharing a home bucket in a fresh table: erasing the first must keep the second findable.
static bool collidingPair() {
   const uint32_t a = 1;
   const uint32_t b = keyWithHome(homeOf(a), a + 1, {});
   OpenBucketHashtable<uint32_t, uint32_t> h;
   std::unordered_map<uint32_t, uint32_t> ref;
   h.at(a) = 10; ref[a] = 10;
   h.at(b) = 20; ref[b] = 20;
   h.erase(a); ref.erase(a);
   return consistent(h, ref);
}

// Nearly full table whose cluster spans the whole table: buckets 0..13 hold keys at home, bucket 14
// holds a second key with home 13. Erasing bucket 0 must not move that key back across the table
// into bucket 0, which would leave an empty bucket 14 between its home and it.
static bool wrappedCluster() {
   OpenBucketHashtable<uint32_t, uint32_t> h;
   std::unordered_map<uint32_t, uint32_t> ref;
   for (size_t bucket = 0; bucket <= 14; ++bucket) {
      const uint32_t k = keyWithHome(bucket == 14 ? 13 : bucket, 1, ref);
      h.at(k) = bucket;
      ref[k] = bucket;
   }
   if (h.bucket_count() != 16) {
      return false; // layout assumption broken, not the bug under test
   }
   const auto first = h.begin();
   ref.erase(first->first);
   h.erase(first);
   return consistent(h, ref);
}

// Random interleaved inserts and erases. Small maxKeys keeps 16-bucket tables nearly full,
// where clusters wrap around the table.
template <typename GID>
static long randomOps(uint32_t maxKeys, int runs, uint32_t seed) {
   std::mt19937 rng(seed);
   std::uniform_int_distribution<GID> gid(0, 4 * maxKeys + 50);
   long failingRuns = 0;
   for (int run = 0; run < runs; ++run) {
      OpenBucketHashtable<GID, uint32_t> h;
      std::unordered_map<GID, uint32_t> ref;
      for (uint32_t step = 0; step < 4 * maxKeys; ++step) {
         if ((rng() & 1) && ref.size() < maxKeys) {
            const GID k = gid(rng);
            h.at(k) = step;
            ref[k] = step;
         } else if (!ref.empty()) {
            auto it = ref.begin();
            std::advance(it, rng() % ref.size());
            h.erase(it->first);
            ref.erase(it);
         }
         if (!consistent(h, ref)) {
            ++failingRuns;
            break;
         }
      }
   }
   return failingRuns;
}

int main() {
   int failures = 0;
   if (!collidingPair()) {
      std::printf("FAIL: colliding pair\n");
      ++failures;
   }
   if (!wrappedCluster()) {
      std::printf("FAIL: wrapped cluster\n");
      ++failures;
   }
   for (const uint32_t maxKeys : {2u, 4u, 8u, 12u, 14u, 15u, 40u, 200u}) {
      const int runs = maxKeys > 40 ? 1000 : 5000;
      const long bad32 = randomOps<uint32_t>(maxKeys, runs, 777 + maxKeys);
      const long bad64 = randomOps<uint64_t>(maxKeys, runs, 999 + maxKeys); // fnv_1a hash path
      if (bad32 || bad64) {
         std::printf("FAIL: maxKeys=%u failing runs uint32=%ld uint64=%ld (of %d)\n", maxKeys, bad32, bad64, runs);
         ++failures;
      }
   }
   std::printf(failures ? "open_bucket_hashtable erase test FAILED\n" : "open_bucket_hashtable erase test passed\n");
   return failures ? 1 : 0;
}
