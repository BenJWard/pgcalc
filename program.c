// This file is a part of Julia. License is MIT: http://julialang.org/license

// Standard headers
#include <string.h>
#include <stdint.h>

// Julia headers (for initialization and gc commands)
#include "uv.h"
#include "julia.h"

// Declare C prototype of a function defined in Julia
extern void julia_main();

int main(int argc, char *argv[])
{
  intptr_t v;

  // Initialize Julia
  uv_setup_args(argc, argv); // no-op on Windows
  libsupport_init();
  jl_options.image_file = "libpgcalc";
  julia_init(JL_IMAGE_CWD);

  // build arguments array: `String[ unsafe_string(argv[i]) for i in 1:argc ]`
  jl_array_t *ARGS = jl_alloc_array_1d(jl_apply_array_type(jl_string_type, 1), 0);
  JL_GC_PUSH1(&ARGS);
  jl_array_grow_end(ARGS, argc - 1);
  for (int i = 1; i < argc; i++) {
      jl_value_t *s = (jl_value_t*)jl_cstr_to_string(argv[i]);
      jl_arrayset(ARGS, s, i - 1);
  }

  // Do some work
  julia_main(ARGS);

  // Cleanup and graceful exit
  jl_atexit_hook(0);
  return 0;
}
