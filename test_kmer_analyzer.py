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