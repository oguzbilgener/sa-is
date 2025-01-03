#![doc = include_str!("../README.md")]

use bitvec::{slice::BitSlice, vec::BitVec};
use num_traits::ToPrimitive;

const S_TYPE: bool = false;
const L_TYPE: bool = true;

/// Suffix Array Induced Sorting (SA-IS) algorithm.
///
/// # Examples
///
/// ```rust
/// use sa_is::make_suffix_array;
///
/// let answer = make_suffix_array(b"banana", char::MAX as usize);
/// assert_eq!(answer, [5, 3, 1, 0, 4, 2]);
/// ```
pub fn make_suffix_array<S>(input: &[S], key_bound: usize) -> Vec<usize>
where
    S: Ord + Clone + ToPrimitive,
{
    let mut suffix_array = vec![0; input.len()];
    suffix_sort_rec(input, key_bound, &mut suffix_array);
    suffix_array
}

fn suffix_sort_rec<S>(input: &[S], key_bound: usize, suffix_array: &mut [usize])
where
    S: Ord + Clone + ToPrimitive,
{
    if input.len() == 1 {
        suffix_array[0] = 0;
    }
    if input.len() < 2 {
        return;
    }

    let (sl_partition, lms_count) = build_sl_partition(input);
    let mut lms_indices = find_lms_suffixes(&sl_partition, lms_count);
    let buckets = make_bucket_count(input, key_bound);

    if lms_indices.len() > 1 {
        // Given |lms_indices| in the same order they appear in |str|, induce
        // LMS substrings relative order and write result to |suffix_array|.
        induced_sort(input, &sl_partition, &lms_indices, &buckets, suffix_array);

        // Given LMS substrings in relative order found in |suffix_array|,
        // map LMS substrings to unique labels to form a new string, |lms_str|.
        let (mut lms_str, label_count) =
            label_lms_substrings(input, &sl_partition, suffix_array, &mut lms_indices);

        if label_count < lms_str.len() {
            // Reorder |lms_str| to have LMS suffixes in the same order they
            // appear in |str|.
            for (i, &lms_index) in lms_indices.iter().enumerate() {
                suffix_array[lms_index] = lms_str[i];
            }

            let mut previous_type = S_TYPE;
            let mut j = 0;
            for i in 0..sl_partition.len() {
                let current_type = sl_partition[i];
                if current_type == S_TYPE && previous_type == L_TYPE {
                    lms_str[j] = suffix_array[i];
                    lms_indices[j] = i;
                    j += 1;
                }
                previous_type = current_type;
            }

            // Recursively apply SuffixSort on |lms_str|, which is formed from
            // labeled LMS suffixes in the same order they appear in |str|.
            // Note that |KeyType| will be size_type because |lms_str| contains
            // indices. |lms_str| is at most half the length of |str|.
            suffix_sort_rec(&lms_str, key_bound, suffix_array);

            // Map LMS labels back to indices in |str| and write result to
            // |lms_indices|. We're using |suffix_array| as a temporary buffer.
            for i in 0..lms_indices.len() {
                suffix_array[i] = lms_indices[suffix_array[i]];
            }

            let length = lms_indices.len();
            lms_indices[..length].copy_from_slice(&suffix_array[..length]);

            // At this point, |lms_indices| contains sorted LMS suffixes of |str|.
        }
    }
    // Given |lms_indices| where LMS suffixes are sorted, induce the full
    // order of suffixes in |str|.
    induced_sort(input, &sl_partition, &lms_indices, &buckets, suffix_array);
}

fn induced_sort<S>(
    input: &[S],
    sl_partition: &BitSlice,
    lms_indices: &[usize],
    buckets: &[usize],
    suffix_array: &mut [usize],
) where
    S: Ord + Clone + ToPrimitive,
{
    let length = input.len();
    // All indices are first marked as unset with the illegal value |length|.
    suffix_array[0..length].fill(length);

    debug_assert!(!buckets.is_empty());
    // Used to mark bucket boundaries (head or end) as indices in str.
    let mut bucket_bounds: Vec<usize> = vec![0; buckets.len()];

    // Step 1: Assign indices for LMS suffixes, populating the end of
    // respective buckets but keeping relative order.

    partial_sum(buckets.iter(), bucket_bounds.iter_mut());

    // Process each `lms_indices` in reverse and assign them to the end of their
    // respective buckets, so relative order is preserved.
    for &lms_index in lms_indices.iter().rev() {
        let key = item_to_key(&input[lms_index]);
        bucket_bounds[key] -= 1;
        suffix_array[bucket_bounds[key]] = lms_index;
    }

    // Step 2
    // Scan forward |suffix_array|; for each modified suf(S,i) for which
    // suf(S,SA(i) - 1) is L-type, place suf(S,SA(i) - 1) to the current
    // head of the corresponding bucket and forward the bucket head to the
    // right.

    // Find the head of each bucket and write it to |bucket_bounds|. Since
    // only LMS suffixes where inserted in |suffix_array| during Step 1,
    // |bucket_bounds| does not contains the head of each bucket and needs to
    // be updated.
    bucket_bounds[0] = 0;
    partial_sum(buckets.iter(), bucket_bounds.iter_mut().skip(1));

    // From Step 1, the sentinel $, which we treat implicitly, would have
    // been placed at the beginning of |suffix_array|, since $ is always
    // considered as the smallest character. We then have to deal with the
    // previous (last) suffix.
    if sl_partition[length - 1] == L_TYPE {
        let key = item_to_key(&input[length - 1]);
        suffix_array[bucket_bounds[key]] = length - 1;
        bucket_bounds[key] += 1;
    }

    for i in 0..length {
        let suffix_index = suffix_array[i];

        // While the original algorithm marks unset suffixes with -1,
        // we found that marking them with |length| is also possible and more
        // convenient because we are working with unsigned integers.
        if suffix_index != length && suffix_index > 0 {
            let suffix_index = suffix_index - 1;
            if sl_partition[suffix_index] == L_TYPE {
                let key = item_to_key(&input[suffix_index]);
                suffix_array[bucket_bounds[key]] = suffix_index;
                bucket_bounds[key] += 1;
            }
        }
    }

    // Step 3
    // Scan backward |suffix_array|; for each modified suf(S, i) for which
    // suf(S,SA(i) - 1) is S-type, place suf(S,SA(i) - 1) to the current
    // end of the corresponding bucket and forward the bucket head to the
    // left.

    // Find the end of each bucket and write it to |bucket_bounds|. Since
    // only L-type suffixes where inserted in |suffix_array| during Step 2,
    // |bucket_bounds| does not contain the end of each bucket and needs to
    // be updated.
    partial_sum(buckets.iter(), bucket_bounds.iter_mut());

    for i in (0..length).rev() {
        let suffix_index = suffix_array[i];
        if suffix_index != length && suffix_index > 0 {
            let suffix_index = suffix_index - 1;
            if sl_partition[suffix_index] == S_TYPE {
                let key = item_to_key(&input[suffix_index]);
                bucket_bounds[key] -= 1;
                suffix_array[bucket_bounds[key]] = suffix_index;
            }
        }
    }

    // Deals with the last suffix, because of the sentinel.
    if sl_partition[length - 1] == S_TYPE {
        let key = item_to_key(&input[length - 1]);
        bucket_bounds[key] -= 1;
        suffix_array[bucket_bounds[key]] = length - 1;
    }
}

/// Partition every suffix based on the SL-type. Return the number of LMS suffixes.
fn build_sl_partition<S>(input: &[S]) -> (BitVec, usize)
where
    S: Ord + Clone,
{
    let length = input.len();
    let mut lms_count = 0;
    let mut sl_partition = bitvec::bitvec![0; length];
    let mut previous_type = L_TYPE;
    let mut previous_key: Option<&S> = None;

    // We're traveling backwards to determine the partition,
    // as if we prepend one character ata time to the string, ex:
    // b$ is L-type because b > $.
    // ab$ is S-type because a < b, implying ab$ < b$.
    // bab$ is L-type because b > a, implying bab$ > ab$.
    // bbab$ is L-type, because bab$ was also L-type, implying bbab$ > bab$.
    for i in (0..length).rev() {
        let current_key = &input[i];

        if previous_key.is_none() || current_key > previous_key.unwrap() {
            // S[i] > S[i + 1] or S[i] is last character.
            if previous_type == S_TYPE {
                // suf(S,i) is L-type and suf(S,i + 1) is S-type, therefore,
                // suf(S,i+1) was a LMS suffix.
                lms_count += 1;
            }
            previous_type = L_TYPE;
        } else if current_key < previous_key.unwrap() {
            // S[i] < S[i + 1]
            previous_type = S_TYPE;
        }
        // Else, S[i] == S[i + 1]:
        // The next character that differs determines the SL-type,
        // so we reuse the last seen type.
        sl_partition.set(i, previous_type);
        previous_key = Some(current_key);
    }

    (sl_partition, lms_count)
}

fn label_lms_substrings<S>(
    input: &[S],
    sl_partition: &BitVec,
    suffix_array: &[usize],
    lms_indices: &mut [usize],
) -> (Vec<usize>, usize)
where
    S: Ord + Clone,
{
    let length = input.len();
    let mut lms_str = vec![0; lms_indices.len()];
    // Labelling starts at 0.
    let mut label = 0;
    // |previous_lms| is initialized to 0 to indicate it is unset.
    // Note that suf(S,0) is never a LMS suffix. Substrings will be visited in
    // lexicographical order.
    let mut previous_lms = 0;
    let mut j = 0;
    for &current_lms in &suffix_array[..length] {
        if current_lms > 0
            && sl_partition[current_lms] == S_TYPE
            && sl_partition[current_lms - 1] == L_TYPE
        {
            // suf(S, i) is an LMS suffix.

            if previous_lms != 0 {
                // There was a previous LMS suffix. Check if the current LMS
                // substring is equal to the previous one.
                let mut current_lms_type = S_TYPE;
                let mut previous_lms_type = S_TYPE;
                let mut k = 0;
                loop {
                    // |current_lms_end| and |previous_lms_end| denote whether we have
                    // reached the end of the current and previous LMS substring,
                    // respectively
                    let mut current_lms_end = false;
                    let mut previous_lms_end = false;

                    // Check for both previous and current substring ends.
                    // Note that it is more convenient to check if
                    // suf(S,current_lms + k) is an LMS suffix than to retrieve it
                    // from lms_indices.
                    if current_lms + k >= length
                        || (current_lms_type == L_TYPE && sl_partition[current_lms + k] == S_TYPE)
                    {
                        current_lms_end = true;
                    }
                    if previous_lms + k >= length
                        || (previous_lms_type == L_TYPE && sl_partition[previous_lms + k] == S_TYPE)
                    {
                        previous_lms_end = true;
                    }

                    if current_lms_end && previous_lms_end {
                        break; // Previous and current substrings are identical.
                    }
                    if current_lms_end != previous_lms_end
                        || input[current_lms + k] != input[previous_lms + k]
                    {
                        // Previous and current substrings differ, a new label is used.
                        label += 1;
                        break;
                    }

                    current_lms_type = sl_partition[current_lms + k];
                    previous_lms_type = sl_partition[previous_lms + k];

                    k += 1;
                }
            }

            lms_indices[j] = current_lms;
            lms_str[j] = label;
            j += 1;

            previous_lms = current_lms;
        }
    }

    (lms_str, label + 1)
}

fn find_lms_suffixes(sl_partition: &BitSlice, lms_count: usize) -> Vec<usize> {
    let mut previous_type = S_TYPE;
    let mut lms_indices = Vec::with_capacity(lms_count);

    for i in 0..sl_partition.len() {
        let current_type = sl_partition[i];
        if current_type == S_TYPE && previous_type == L_TYPE {
            lms_indices.push(i);
        }
        previous_type = current_type;
    }

    lms_indices
}

fn make_bucket_count<S>(input: &[S], key_bound: usize) -> Vec<usize>
where
    S: Ord + Clone + ToPrimitive,
{
    let mut buckets = vec![0usize; key_bound];
    for c in input {
        let key = item_to_key(c);
        buckets[key] += 1;
    }

    buckets
}

fn item_to_key<S>(item: &S) -> usize
where
    S: ToPrimitive,
{
    item.to_usize().expect("input is not convertible to usize")
}

fn partial_sum<'i, T>(
    iter_from: impl Iterator<Item = &'i T>,
    iter_to: impl Iterator<Item = &'i mut T>,
) where
    T: 'i + std::ops::AddAssign + Default + Copy,
{
    let mut sum: T = T::default();
    for (item_from, item_to) in iter_from.zip(iter_to) {
        sum += *item_from;
        *item_to = sum;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    const TEST_STRS_1: [&[u8]; 20] = [
        b"",
        b"a",
        b"aa",
        b"za",
        b"CACAO",
        b"aaaaa",
        b"banana",
        b"tobeornottobe",
        b"The quick brown fox jumps over the lazy dog.",
        b"elephantelephantelephantelephantelephant",
        b"walawalawashington",
        b"-------------------------",
        b"011010011001011010010110011010010",
        b"3141592653589793238462643383279502884197169399375105",
        b"\xFF\xFE\xFF\xFE\xFD\x80\x30\x31\x32\x80\x30\xFF\x01\xAB\xCD",
        b"abccbaabccbaabccbaabccbaabccbaabccbaabccbaabccba",
        b"0123456789876543210",
        b"9876543210123456789",
        b"aababcabcdabcdeabcdefabcdefg",
        b"asdhklgalksdjghalksdjghalksdjgh",
    ];

    #[rstest]
    #[case(&bits(&[]), &[], b"")]
    #[case(&bits(&[L_TYPE]), &[], b"a")]
    #[case(&bits(&[L_TYPE, L_TYPE]), &[], b"ba")]
    #[case(&bits(&[S_TYPE, L_TYPE]), &[], b"ab")]
    #[case(&bits(&[S_TYPE, S_TYPE, L_TYPE]), &[], b"aab")]
    #[case(&bits(&[L_TYPE, L_TYPE, L_TYPE]), &[], b"bba")]
    #[case(&bits(&[L_TYPE, S_TYPE, L_TYPE]), &[1], b"bab")]
    #[case(&bits(&[L_TYPE, S_TYPE, S_TYPE, L_TYPE]), &[1], b"baab")]
    #[case(&bits(&[
        L_TYPE,  // zucchini
        L_TYPE,  // ucchini
        S_TYPE,  // cchini
        S_TYPE,  // chini
        S_TYPE,  // hini
        S_TYPE,  // ini
        L_TYPE,  // ni
        L_TYPE,  // i
    ]), &[2], b"zucchini")]
    fn test_sl_partition(
        #[case] expected_sl_partition: &BitSlice,
        #[case] expected_lms_indices: &[usize],
        #[case] input: &[u8],
    ) {
        let (actual_sl_partition, lms_count) = build_sl_partition(input);
        assert_eq!(expected_lms_indices.len(), lms_count);
        assert_eq!(expected_sl_partition, actual_sl_partition.as_bitslice());

        let suffixes = find_lms_suffixes(expected_sl_partition, lms_count);
        assert_eq!(expected_lms_indices, suffixes.as_slice());
    }

    #[rstest]
    #[case(&[0, 0, 0, 0], &[], 4)]
    #[case(&[1, 0, 0, 0], &[0], 4)]
    #[case(&[0, 2, 0, 1], &[1, 1, 3], 4)]
    #[case(&[0, 0, 2, 1, 1], &[2, 3, 4, 2], 5)]
    #[case(&[0, 1, 3, 0, 0], &[2, 2, 1, 2], 5)]
    #[case(&[0, 0, 3, 1, 0], &[2, 2, 3, 2], 5)]
    fn test_bucket_count(
        #[case] expected: &[usize],
        #[case] input: &[u8],
        #[case] key_bound: usize,
    ) {
        let buckets = make_bucket_count(input, key_bound);
        assert_eq!(key_bound, buckets.len());
        assert_eq!(expected, buckets.as_slice());
    }

    fn induced_sort_substring(input: &str) -> Vec<usize> {
        let (sl_partition, lms_count) = build_sl_partition(input.as_bytes());
        let suffixes = find_lms_suffixes(&sl_partition, lms_count);
        let buckets = make_bucket_count(input.as_bytes(), 256);
        let mut suffix_array = vec![0; input.len()];
        induced_sort(
            input.as_bytes(),
            &sl_partition,
            &suffixes,
            &buckets,
            &mut suffix_array,
        );
        suffix_array
    }

    #[rstest]
    // L; a$
    #[case(&[0], "a")]
    // SL; ab$, b$
    #[case(&[0, 1], "ab")]
    // LL; a$, ba$
    #[case(&[1, 0], "ba")]
    // SLL; a$, aba$, ba$
    #[case(&[2, 0, 1], "aba")]
    // LSL; ab$, b$, ba
    #[case(&[1, 2, 0], "bab")]
    // SSL; aab$, ab$, b$
    #[case(&[0, 1, 2], "aab")]
    // LSSL; aab$, ab$, b$, ba
    #[case(&[1, 2, 3, 0], "baab")]
    fn test_induced_sort_substring(#[case] expected: &[usize], #[case] input: &str) {
        let actual = induced_sort_substring(input);
        assert_eq!(expected, actual.as_slice());
    }

    #[test]
    fn test_suffix_sort_1() {
        run_test_suffix_sort(&TEST_STRS_1);
    }

    #[test]
    fn test_suffix_sort_all_chars() {
        let mut input = Vec::with_capacity(256);
        for i in 0u8..=255 {
            input.push(i);
        }
        run_test_suffix_sort(&[input.as_slice()]);

        let mut input = Vec::with_capacity(256);
        for i in (0u8..=255).rev() {
            input.push(i);
        }

        run_test_suffix_sort(&[input.as_slice()]);
    }

    fn run_test_suffix_sort<const N: usize>(cases: &[&[u8]; N]) {
        for input in cases {
            let input_str = String::from_utf8_lossy(input);
            let suffix_array = make_suffix_array(input, 256);
            assert_eq!(input.len(), suffix_array.len(), "input: {input_str}");

            // Expect that I[] is a permutation of [0, len]
            let mut sorted_suffix = suffix_array.clone();
            sorted_suffix.sort_unstable();
            let expected = (0..input.len()).collect::<Vec<_>>();
            assert_eq!(expected, sorted_suffix, "input: {input_str}");

            // Expect that all suffixes are strictly ordered.
            for i in 1..suffix_array.len() {
                let suffix_a = &input[suffix_array[i - 1]..];
                let suffix_b = &input[suffix_array[i]..];
                assert!(suffix_a < suffix_b, "input: {input_str}");
            }
        }
    }

    #[rstest]
    #[case(&[], &[])]
    #[case(&[0], &[0])]
    #[case(&[0, 1], &[0, 1])]
    #[case(&[0, 1, 2], &[0, 1, 3])]
    #[case(&[0, 1, 2, 3], &[0, 1, 3, 6])]
    fn test_partial_sum(#[case] input: &[usize], #[case] expected: &[usize]) {
        let mut actual = vec![0; input.len()];
        partial_sum(input.iter(), actual.iter_mut());
        assert_eq!(expected, actual.as_slice());
    }

    #[test]
    fn test_partial_sum_different_lengths() {
        let left = &[0, 1, 2, 3];
        let dst = &mut [1, 1, 1];
        partial_sum(left.iter(), dst.iter_mut());
        assert_eq!(&[0, 1, 3], dst);
    }

    fn bits(input: &[bool]) -> BitVec {
        let mut bits = BitVec::with_capacity(input.len());
        for bit in input {
            bits.push(*bit);
        }
        bits
    }
}
