# sa-is

Suffix Array Induced Sorting (SA-IS) algorithm based on the paper ["Linear Suffix Array Construction by Almost Pure Induced-Sorting"](https://code.google.com/archive/p/ge-nong/downloads) by G. Nong, S. Zhang and W. H. Chan and [the C++ implementation](https://github.com/chromium/chromium/blob/main/components/zucchini/suffix_array.h) in Chromium.

## Usage

```rust
use sa_is::make_suffix_array;

let answer = make_suffix_array(b"banana", char::MAX as usize);
assert_eq!(answer, [5, 3, 1, 0, 4, 2]);
```

## License

This project is licensed under both the BSD-3-Clause and MIT licenses:

1. The original code from the Chromium project is licensed under the BSD-3-Clause license.
2. Any new contributions made to this project are licensed under the MIT license.

See the [LICENSE](./LICENSE) file for details.
