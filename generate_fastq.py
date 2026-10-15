import random  # Added this missing import

# Function to generate quality scores matching sequence length
def generate_quality(length):
    # Using ASCII characters from ! (33) to I (73) for quality scores
    return ''.join(chr(random.randint(33, 73)) for _ in range(length))

# Function to create FASTQ file with manual sequences
def create_fastq_file(filename, sequences):
    with open(filename, 'w') as f:
        for i, sequence in enumerate(sequences, 1):
            # FASTQ format: @identifier, sequence, +, quality scores
            seq_id = f"@SEQ{i}"
            quality = generate_quality(len(sequence))
            
            # Write the four lines of FASTQ format
            f.write(f"{seq_id}\n")
            f.write(f"{sequence.upper()}\n")  # Convert to uppercase
            f.write("+\n")
            f.write(f"{quality}\n")

# Get 15 sequences from user input
sequences = []
print("Enter 15 DNA sequences (one per line):")
for i in range(1, 16):
    while True:
        seq = input(f"Sequence {i}: ").strip()
        # Basic validation: check if sequence contains only A, T, C, G
        if all(base in 'ATCGatcg' for base in seq) and len(seq) > 0:
            sequences.append(seq)
            break
        else:
            print("Invalid sequence! Please use only A, T, C, G characters.")

# Generate the file
create_fastq_file("manual_sequences.fastq", sequences)
print("FASTQ file 'manual_sequences.fastq' has been created with your 15 sequences.")