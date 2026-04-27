import sys

def validate_sequence(sequence, k):
    """Validate that a sequence can be processed for k-mer analysis.

    Parameters:
        sequence (str): DNA sequence string to validate.
        k (int): Requested k-mer length.

    Returns:
        bool: True when k is a positive integer, the sequence length is at
            least k, and every character in the sequence is a valid nucleotide
            (A, C, G, or T, case-insensitive); otherwise False. 
    """
    # Check if the sequence and k are the correct types
    if not isinstance(sequence, str):
        return False
    if not isinstance(k, int) or isinstance(k, bool):
        return False

    # A sequence must be at least k characters long to contain any k-mer.
    if k <= 0 or len(sequence) < k:
        return False

    # Scan each character and reject immediately if a digit is present.
    for nucleotide in sequence:
        if nucleotide not in 'ACGT':
            return False

    # Sequence passed all checks performed by this function.
    return True

def update_kmer_count(kmer_data, kmer, next_char):
    """Update k-mer count totals and next-character frequency data.

    Parameters:
        kmer_data (dict): Mapping of k-mer strings to count/context metadata.
        kmer (str): K-mer key to update in the mapping.
        next_char (str): Character observed immediately after the k-mer.

    Returns:
        dict: The updated mapping where the k-mer total count increases by one
            and the frequency for next_char under that k-mer increases by one.
    """
    # Initialize storage for a k-mer the first time we encounter it.
    if kmer not in kmer_data:
        kmer_data[kmer] = {'count': 1, 'next_chars': {}}
    # If the k-mer is already in the dictionary, increase the count.
    else:
        kmer_data[kmer]['count'] += 1
    
    # Initialize count bucket for this following character if it doesn't exist.
    if next_char not in kmer_data[kmer]['next_chars']:
        kmer_data[kmer]['next_chars'][next_char] = 1
    # If the next-character is already in the dictionary, increase the count.
    else:
        kmer_data[kmer]['next_chars'][next_char] += 1

    # Return the updated dictionary so callers can chain/use the result.
    return kmer_data

def count_kmers_with_context(sequence, k):
    """Count k-mers in a sequence with the frequency of following characters.

    Parameters:
        sequence (str): DNA sequence to analyze.
        k (int): Length of each k-mer window.

    Returns:
        dict: Mapping of each observed k-mer (that has at least one following
            character) to its total occurrence count and next-character
            frequencies across the sequence.
    """
    # Store aggregate counts for every observed k-mer and following character.
    kmer_data = {}
    
    # Stop at len(sequence) - k so each k-mer has a valid next character.
    for i in range(len(sequence) - k):
        # Slice the current k-mer window and capture the character after it.
        kmer = sequence[i:i+k]
        next_char = sequence[i+k]
        
        # Delegate counting logic to the helper that updates nested metadata.
        kmer_data = update_kmer_count(kmer_data, kmer, next_char)
    
    # Return all k-mer context frequency data for this sequence.
    return kmer_data


def write_results_to_file(kmer_data, output_filename):
    """Write k-mer context data to a text file in sorted order.

    Parameters:
        kmer_data (dict): Mapping of k-mers and their next-character frequencies.
        output_filename (str): Path to the output text file.

    Returns:
        None: Writes one line per k-mer with the k-mer followed by
            space-separated `character:frequency` entries, sorted by k-mer and
            then by character.
    """
    # Sort k-mers so output order is deterministic and easy to compare.
    sorted_kmers = sorted(kmer_data.keys())
    
    # Write one line per k-mer in the required text format.
    with open(output_filename, 'w') as f:
        for kmer in sorted_kmers:
            # Read total frequency and per-k-mer following-character counts.
            total_count = kmer_data[kmer]['count']
            # Read per-k-mer following-character counts.
            next_chars = kmer_data[kmer]['next_chars']
            
            # Build "char:count" tokens sorted by character, then join with spaces.
            next_char_str = " ".join(
                f"{char}:{freq}" 
                for char, freq in sorted(next_chars.items())
            )
            
            # Emit one output line: kmer, total count, then context frequencies.
            f.write(f"{kmer} {total_count} {next_char_str}\n")


def main():
    """Run the command-line workflow for k-mer context analysis.
    
    Parameters:
        None
    
    Returns:
        None: Reads sequence fragments from the input file, validates each
            sequence, computes k-mer context frequencies, and writes the final
            results to the specified output file.
    """
    # Parse required command-line arguments.
    sequence_file = sys.argv[1]
    k = int(sys.argv[2])
    output_file = sys.argv[3]

    # Basic progress message for the user.
    print(f"Reading sequences from {sequence_file}...")

    # Accumulate k-mer data across all sequences before writing.
    kmer_data = {}

    # Process each sequence fragment from the input file.
    with open(sequence_file, 'r') as f:
        for sequence in f:
            # Remove trailing whitespace/newline characters and convert to uppercase.
            sequence = sequence.strip().upper()
            
            # Skip malformed sequences and continue processing remaining lines.
            if not validate_sequence(sequence, k):
                print(f"  Warning: Skipping sequence")
                continue

            # Merge this sequence's k-mer counts into the running totals.
            sequence_kmer_data = count_kmers_with_context(sequence, k)
            for kmer, data in sequence_kmer_data.items():
                if kmer not in kmer_data:
                    kmer_data[kmer] = {'count': 0, 'next_chars': {}}

                # Add the total count contributed by this sequence.
                kmer_data[kmer]['count'] += data['count']

                # Add each next-character frequency into the aggregate map.
                for next_char, freq in data['next_chars'].items():
                    if next_char not in kmer_data[kmer]['next_chars']:
                        kmer_data[kmer]['next_chars'][next_char] = 0
                    kmer_data[kmer]['next_chars'][next_char] += freq

    # Write the final aggregated results once after all sequences are processed.
    write_results_to_file(kmer_data, output_file)

if __name__ == '__main__':
    main()
