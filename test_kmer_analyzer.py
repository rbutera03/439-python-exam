import pytest
from kmer_analyzer import validate_sequence, update_kmer_count, count_kmers_with_context, write_results_to_file

class TestValidateSequence:
    # --- Valid sequences ---
    def test_valid_dna_sequence(self):
        """Standard valid DNA sequence longer than k"""
        assert validate_sequence("ATGC", 2) == True

    def test_valid_all_same_nucleotide(self):
        """Sequence of repeated single nucleotide"""
        assert validate_sequence("AAAA", 2) == True

    def test_valid_sequence_length_equals_k(self):
        """Sequence exactly as long as k should be valid"""
        assert validate_sequence("ATG", 3) == True
        
    def test_sequence_with_lowercase(self):
        """Lowercase nucleotides should be valid"""
        assert validate_sequence("atgc", 2) == True

    def test_sequence_with_mixed_case(self):
        """Mixed case sequences should be valid"""
        assert validate_sequence("AtGc", 2) == True

    # --- Length checks ---
    def test_sequence_shorter_than_k(self):
        """Sequence shorter than k should be invalid"""
        assert validate_sequence("AT", 5) == False

    def test_empty_sequence(self):
        """Empty sequence should always be invalid"""
        assert validate_sequence("", 2) == False
        
    def test_sequence_with_zero_k(self):
        """k=0 is meaningless and should be invalid"""
        assert validate_sequence("ATGC", 0) == False
        
    def test_sequence_with_negative_k(self):
        """Negative k should be invalid"""
        assert validate_sequence("ATGC", -1) == False


    # --- Invalid characters ---
    def test_sequence_with_digit(self):
        """Digits should make a sequence invalid"""
        assert validate_sequence("ATG1C", 2) == False

    def test_sequence_with_special_character(self):
        """Special characters should make a sequence invalid"""
        assert validate_sequence("ATG!C", 2) == False
        
    def test_sequence_with_multiple_special_characters(self):
        """Multiple special characters should make a sequence invalid"""
        assert validate_sequence("ATG!@#$%^&*()_+-=[]{}|\\:;\"'<>,.?/~`", 2) == False
        
    def test_sequence_with_space(self):
        """Space should make a sequence invalid"""
        assert validate_sequence("ATG C", 2) == False
        
    def test_sequence_with_invalid_letter(self):
        """Letters outside ACGT should make a sequence invalid"""
        assert validate_sequence("ATGX", 2) == False
        
    def test_sequence_with_multiple_invalid_letters(self):
        """Multiple invalid letters should make a sequence invalid"""
        assert validate_sequence("ATGXZ", 2) == False
        
    def test_sequence_with_multiple_invalid_letters_and_special_characters(self):
        """Multiple invalid letters and special characters should make a sequence invalid"""
        assert validate_sequence("ATGX!@#$%^&*()_+-=[]{}|\\:;\"'<>,.?/~`", 2) == False
        
    def test_sequence_with_multiple_invalid_letters_and_spaces(self):
        """Multiple invalid letters and spaces should make a sequence invalid"""
        assert validate_sequence("ATGX C", 2) == False
        
    # --- Invalid types ---
    def test_sequence_is_not_string(self):
        """Non-string sequence should be invalid"""
        assert validate_sequence(123, 2) == False

    def test_k_is_not_int(self):
        """Non-integer k should be invalid"""
        assert validate_sequence("ATGC", 2.0) == False

    def test_k_is_bool(self):
        """Boolean k should be invalid"""
        assert validate_sequence("ATGC", True) == False

    def test_k_is_zero(self):
        """k=0 is meaningless and should be invalid"""
        assert validate_sequence("ATGC", 0) == False
        
class TestUpdateKmerCount:
    # --- New kmer ---
    def test_new_kmer_is_added(self):
        """A kmer not yet in kmer_data should be added"""
        kmer_data = {}
        result = update_kmer_count(kmer_data, "AT", "G")
        assert "AT" in result

    def test_new_kmer_count_is_one(self):
        """A brand new kmer should have a count of 1 — bug: will be 2"""
        kmer_data = {}
        result = update_kmer_count(kmer_data, "AT", "G")
        assert result["AT"]["count"] == 1

    def test_new_kmer_next_char_recorded(self):
        """Next character after a new kmer should be recorded with frequency 1"""
        kmer_data = {}
        result = update_kmer_count(kmer_data, "AT", "G")
        assert result["AT"]["next_chars"]["G"] == 1

    # --- Existing kmer, same next char ---
    def test_existing_kmer_count_increments(self):
        """Seeing the same kmer twice should give count of 2"""
        kmer_data = {}
        update_kmer_count(kmer_data, "AT", "G")
        result = update_kmer_count(kmer_data, "AT", "G")
        assert result["AT"]["count"] == 2

    def test_existing_kmer_next_char_increments(self):
        """Seeing the same kmer+next_char twice should give next_char frequency 2"""
        kmer_data = {}
        update_kmer_count(kmer_data, "AT", "G")
        result = update_kmer_count(kmer_data, "AT", "G")
        assert result["AT"]["next_chars"]["G"] == 2

    # --- Existing kmer, different next char ---
    def test_existing_kmer_new_next_char(self):
        """Same kmer followed by a different character should add a new next_char entry"""
        kmer_data = {}
        update_kmer_count(kmer_data, "AT", "G")
        result = update_kmer_count(kmer_data, "AT", "C")
        assert "C" in result["AT"]["next_chars"]

    def test_existing_kmer_different_next_chars_tracked_separately(self):
        """Two different next chars for the same kmer should each have frequency 1"""
        kmer_data = {}
        update_kmer_count(kmer_data, "AT", "G")
        result = update_kmer_count(kmer_data, "AT", "C")
        assert result["AT"]["next_chars"]["G"] == 1
        assert result["AT"]["next_chars"]["C"] == 1

    # --- Multiple kmers ---
    def test_multiple_different_kmers_tracked_independently(self):
        """Different kmers should not affect each other's counts"""
        kmer_data = {}
        update_kmer_count(kmer_data, "AT", "G")
        result = update_kmer_count(kmer_data, "TG", "C")
        assert "AT" in result
        assert "TG" in result
        assert result["AT"]["count"] == 1  
        assert result["TG"]["count"] == 1  

    # --- Return value ---
    def test_returns_kmer_data_dict(self):
        """Function should return the kmer_data dictionary"""
        kmer_data = {}
        result = update_kmer_count(kmer_data, "AT", "G")
        assert isinstance(result, dict)

class TestCountKmersWithContext:
    # --- Basic behavior ---
    def test_returns_dict(self):
        """Should return a dictionary"""
        result = count_kmers_with_context("ATGC", 2)
        assert isinstance(result, dict)

    def test_correct_kmers_extracted(self):
        """Should extract the correct set of kmers from a sequence"""
        result = count_kmers_with_context("ATGC", 2)
        assert set(result.keys()) == {"AT", "TG"}

    def test_kmer_count_single_occurrence(self):
        """A kmer appearing once should have count 1"""
        result = count_kmers_with_context("ATGC", 2)
        assert result["AT"]["count"] == 1

    def test_next_char_recorded_correctly(self):
        """The character following each kmer should be recorded"""
        result = count_kmers_with_context("ATGC", 2)
        assert result["AT"]["next_chars"] == {"G": 1}
        assert result["TG"]["next_chars"] == {"C": 1}

    def test_repeated_kmer_count(self):
        """A kmer appearing multiple times should have the correct count"""
        result = count_kmers_with_context("ATGAT", 2)
        assert result["AT"]["count"] == 2

    def test_repeated_kmer_same_next_char(self):
        """Repeated kmer with the same next char should accumulate frequency"""
        result = count_kmers_with_context("ATGATG", 2)
        assert result["AT"]["next_chars"]["G"] == 2

    def test_repeated_kmer_different_next_chars(self):
        """Same kmer followed by different chars should track both"""
        result = count_kmers_with_context("ATATC", 2)
        assert result["AT"]["next_chars"]["A"] == 1
        assert result["AT"]["next_chars"]["C"] == 1

    def test_last_kmer_excluded(self):
        """The final kmer in a sequence has no following character and should not appear"""
        result = count_kmers_with_context("ATGC", 2)
        assert "GC" not in result

    # --- k equals sequence length ---
    def test_k_equals_sequence_length(self):
        """When k equals sequence length there are no valid kmers with a following char"""
        result = count_kmers_with_context("ATG", 3)
        assert result == {}

    # --- k = 1 ---
    def test_k_equals_one(self):
        """k=1 should treat each single nucleotide as a kmer"""
        result = count_kmers_with_context("ATGC", 1)
        assert set(result.keys()) == {"A", "T", "G"}
        assert result["A"]["next_chars"] == {"T": 1}

class TestCountKmersWithContext:
    # --- Basic behavior ---
    def test_returns_dict(self):
        """Should return a dictionary"""
        result = count_kmers_with_context("ATGC", 2)
        assert isinstance(result, dict)

    def test_correct_kmers_extracted(self):
        """Should extract the correct set of kmers from a sequence"""
        result = count_kmers_with_context("ATGC", 2)
        assert set(result.keys()) == {"AT", "TG"}

    def test_kmer_count_single_occurrence(self):
        """A kmer appearing once should have count 1"""
        result = count_kmers_with_context("ATGC", 2)
        assert result["AT"]["count"] == 1

    def test_next_char_recorded_correctly(self):
        """The character following each kmer should be recorded"""
        result = count_kmers_with_context("ATGC", 2)
        assert result["AT"]["next_chars"] == {"G": 1}
        assert result["TG"]["next_chars"] == {"C": 1}

    def test_repeated_kmer_count(self):
        """A kmer appearing multiple times should have the correct count"""
        result = count_kmers_with_context("ATGATG", 2)
        assert result["AT"]["count"] == 2

    def test_repeated_kmer_same_next_char(self):
        """Repeated kmer with the same next char should accumulate frequency"""
        result = count_kmers_with_context("ATGATG", 2)
        assert result["AT"]["next_chars"]["G"] == 2

    def test_repeated_kmer_different_next_chars(self):
        """Same kmer followed by different chars should track both"""
        result = count_kmers_with_context("ATATC", 2)
        assert result["AT"]["next_chars"]["A"] == 1
        assert result["AT"]["next_chars"]["C"] == 1

    def test_last_kmer_excluded(self):
        """The final kmer in a sequence has no following character and should not appear"""
        result = count_kmers_with_context("ATGC", 2)
        assert "GC" not in result
        
    def test_kmer_with_no_following_char_not_counted(self):
        """A kmer that only appears at the end of the sequence with no following
        character should not be counted at all"""
        result = count_kmers_with_context("ATGAT", 2)
        assert result["AT"]["count"] == 1

    # --- k equals sequence length ---
    def test_k_equals_sequence_length(self):
        """When k equals sequence length there are no valid kmers with a following char"""
        result = count_kmers_with_context("ATG", 3)
        assert result == {}

    # --- k = 1 ---
    def test_k_equals_one(self):
        """k=1 should treat each single nucleotide as a kmer"""
        result = count_kmers_with_context("ATGC", 1)
        assert set(result.keys()) == {"A", "T", "G"}
        assert result["A"]["next_chars"] == {"T": 1}

class TestWriteResultsToFile:
    def test_output_file_exact_formatting(self, tmp_path):
        """
        Verifies the exact line format: KMER TOTAL_COUNT CHAR:FREQ CHAR:FREQ
        This ensures the 'Total Frequency' requirement from the rubric is met.
        """
        output_file = tmp_path / "output.txt"
        # Data representing: AT appears 3 times (2 followed by G, 1 followed by C)
        kmer_data = {
            "AT": {
                "count": 3, 
                "next_chars": {"G": 2, "C": 1}
            }
        }
        
        write_results_to_file(kmer_data, str(output_file))
        content = output_file.read_text().strip()
        
        # We expect: Kmer <space> TotalCount <space> SortedNextChars
        # Note: The next_chars should be sorted alphabetically (C before G)
        expected_line = "AT 3 C:1 G:2"
        assert content == expected_line

    def test_multiple_kmers_sorting_and_format(self, tmp_path):
        """
        Ensures multiple kmers are on separate lines and sorted alphabetically.
        """
        output_file = tmp_path / "output.txt"
        kmer_data = {
            "TG": {"count": 1, "next_chars": {"A": 1}},
            "AT": {"count": 2, "next_chars": {"G": 2}},
        }
        
        write_results_to_file(kmer_data, str(output_file))
        lines = output_file.read_text().strip().split("\n")
        
        assert len(lines) == 2
        assert lines[0].startswith("AT 2")  # AT comes before TG
        assert lines[1].startswith("TG 1")
        assert lines[0] == "AT 2 G:2"
        assert lines[1] == "TG 1 A:1"

    def test_empty_data_produces_empty_file(self, tmp_path):
        """Ensures no data results in a clean, empty file."""
        output_file = tmp_path / "output.txt"
        write_results_to_file({}, str(output_file))
        assert output_file.read_text() == ""