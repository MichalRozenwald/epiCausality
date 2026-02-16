import random

def find_cpgs(dna: str) -> list[int]:
    """Return the start positions of all CpG dinucleotides in a DNA string."""
    dna = dna.upper()
    return [i for i in range(len(dna) - 1) if dna[i] == "C" and dna[i + 1] == "G"]


def count_surrounding_cpgs(dna: str, n: int) -> dict[int, int]:
    """
    For each CpG in dna, count how many *other* CpGs have their start position
    within N bases (i.e. |pos_other - pos_self| <= N, excluding self).

    Returns a dict mapping each CpG start position -> surrounding CpG count.
    """
    cpg_positions = find_cpgs(dna)
    result = {}
    for pos in cpg_positions:
        count = sum(
            1
            for other in cpg_positions
            if other != pos and abs(other - pos) <= n
        )
        result[pos] = count
    return result


def generate_dna(length: int, cpg_rate: float = 0.05, seed: int = 42) -> str:
    """
    Generate a random DNA string of a given length.
    cpg_rate controls the approximate fraction of positions that start a CpG.
    """
    rng = random.Random(seed)
    bases = list("ACGT")
    seq = []
    i = 0
    while i < length:
        if i < length - 1 and rng.random() < cpg_rate:
            seq.append("C")
            seq.append("G")
            i += 2
        else:
            seq.append(rng.choice(bases))
            i += 1
    return "".join(seq[:length])


if __name__ == "__main__":
    # --- Generate a 100 bp test sequence ---
    SEQ_LEN = 100
    N = 10  # window size (bases on each side)

    dna = generate_dna(SEQ_LEN, cpg_rate=0.08, seed=42)
    print(f"DNA sequence ({SEQ_LEN} bp):")
    print(dna)
    print()

    cpg_positions = find_cpgs(dna)
    print(f"CpG positions (0-based): {cpg_positions}")
    print(f"Total CpGs found: {len(cpg_positions)}")
    print()

    surrounding = count_surrounding_cpgs(dna, n=N)
    print(f"Surrounding CpG counts within +/-{N} bases:")
    print(f"{'CpG pos':>10}  {'context':>20}  {'surrounding CpGs':>16}")
    print("-" * 52)
    for pos, count in surrounding.items():
        # Show a small context window around the CpG
        start = max(0, pos - 3)
        end = min(SEQ_LEN, pos + 5)
        context = dna[start:end]
        print(f"{pos:>10}  {context:>20}  {count:>16}")
