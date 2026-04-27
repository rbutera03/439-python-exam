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
        assert result["AT"]["count"] == 1  # will fail due to off-by-one bug
        assert result["TG"]["count"] == 1  # will fail due to off-by-one bug

    # --- Return value ---
    def test_returns_kmer_data_dict(self):
        """Function should return the kmer_data dictionary"""
        kmer_data = {}
        result = update_kmer_count(kmer_data, "AT", "G")
        assert isinstance(result, dict)


