# Third-party notices

## nanoflann

`src/nanoflann.hpp` is a vendored copy of the `nanoflann` header-only C++
kd-tree library.

- Project: https://github.com/jlblancoc/nanoflann
- Vendored file: `include/nanoflann.hpp`
- Source URL: https://raw.githubusercontent.com/jlblancoc/nanoflann/ba47cfcb127c3597d69196d87f5aa9ca8811b0a9/include/nanoflann.hpp
- Pinned commit: `ba47cfcb127c3597d69196d87f5aa9ca8811b0a9`
- Retrieved for this package: 2026-05-15
- Vendored file SHA-256: `e91b910b6fcdeaefb0253d50438c519de675980cbf3374397da730dfc23bbfe7`
- License: BSD license, as stated in the header.

The vendored header was scanned for obvious unsafe behavior on 2026-05-15. No
process execution, network access, filesystem mutation, dynamic loading, or
inline assembly was found. The header includes stream serialization helpers and
a pooled allocator using `malloc`/`free`; these are normal library internals and
are not used for external IO by `seg`.
