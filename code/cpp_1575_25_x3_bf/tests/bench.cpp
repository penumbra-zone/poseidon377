/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

extern "C" {
    #include "tools/bench_timing.h"
    #include "tools/bench_utils.h"
}

#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <HadesMIMC.h>
#include <random>
#include <array>
#include <utils/F_2_63_ops.h>

typedef struct {
  uint64_t instance_permute, instance_sbox, instance_lin, instance_ark, mul_reduce, reduce;
} timing_t;

static void print_timings(timing_t* timings, unsigned int iter) {
  timing_t total = {0,0,0,0};
  for (unsigned int i = 0; i < iter; i++) {
    printf("%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 "\n", timings[i].instance_permute,
           timings[i].instance_sbox, timings[i].instance_lin, timings[i].instance_ark, timings[i].mul_reduce, timings[i].reduce);

    total.instance_permute += timings[i].instance_permute;
    total.instance_sbox += timings[i].instance_sbox;
    total.instance_lin += timings[i].instance_lin;
    total.instance_ark += timings[i].instance_ark;
    total.mul_reduce += timings[i].mul_reduce;
    total.reduce += timings[i].reduce;
  }

  printf("%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 "\n", total.instance_permute,
         total.instance_sbox, total.instance_lin, total.instance_ark, total.mul_reduce, total.reduce);
}

static void bench_instance(const bench_options_t* options) {
  timing_t* timings = new timing_t[options->iter];


  timing_context_t ctx;
  if (!timing_init(&ctx)) {
    printf("Failed to initialize timing functionality.\n");
    return;
  }

  hashMIMC::HadesMIMCPermutation* perm;
  switch(options->instance) {
    case HASHMIMC_INSTANCES::INSTANCE_63_25_12_34:
      perm = new hashMIMC::HadesMIMCPermutation_63_25_12_34();
      break;
    default:
      printf("Unknown Instance.\n");
      return;
  }

  std::vector<uint64_t> arr(perm->get_num_sboxes());
  std::random_device rd1;
  std::uniform_int_distribution<uint64_t> dis;
  __m128i test;

  for (unsigned int i = 0; i != options->iter; ++i) {
    for(unsigned int j = 0; j < arr.size(); j++) {
      arr[j] = dis(rd1);
    }
    test = _mm_set1_epi64x(dis(rd1));

    timing_t* timing = &timings[i];

    uint64_t start_time = timing_read(&ctx);
    perm->permute(arr.data());
    uint64_t tmp_time = timing_read(&ctx);
    timing->instance_permute = tmp_time - start_time;

    start_time        = timing_read(&ctx);
    perm->SBOX(arr.data());
    tmp_time = timing_read(&ctx);
    timing->instance_sbox    = tmp_time - start_time;

    start_time        = timing_read(&ctx);
    perm->LIN(arr.data());
    tmp_time = timing_read(&ctx);
    timing->instance_lin    = tmp_time - start_time;

    start_time        = timing_read(&ctx);
    perm->ARK(arr.data(), 0);
    tmp_time = timing_read(&ctx);
    timing->instance_ark    = tmp_time - start_time;

    start_time        = timing_read(&ctx);
    hashMIMC::F_2_63::mulReduce(arr[0], arr[1]);
    tmp_time = timing_read(&ctx);
    timing->mul_reduce    = tmp_time - start_time;

    start_time        = timing_read(&ctx);
    hashMIMC::F_2_63::reduceOPT2(test);
    tmp_time = timing_read(&ctx);
    timing->reduce    = tmp_time - start_time;

  }

  timing_close(&ctx);
  print_timings(timings, options->iter);

  delete[] timings;
  delete perm;
}

int main(int argc, char** argv) {
  bench_options_t opts = {INSTANCE_INVALID, 0};
  int ret              = parse_args(&opts, argc, argv) ? 0 : -1;

  if (!ret) {
    bench_instance(&opts);
  }

  return ret;
}
