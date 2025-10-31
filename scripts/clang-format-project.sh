#!/bin/bash
# set -x
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd -- "${script_dir}/.." && pwd)"

style="file:${repo_root}/.clang-format"

find "${repo_root}/src" -type f \
  \( -name '*.c' -o -name '*.cc' -o -name '*.cpp' -o -name '*.h' -o -name '*.hpp' \) -print0 |
  xargs -0 clang-format -i -style="${style}"