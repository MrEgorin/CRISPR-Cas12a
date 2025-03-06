import re
from Bio import Entrez, SeqIO
import pandas as pd

def fetch_sequence_from_ncbi(accession_number, email):
    """
    Завантажує послідовність із NCBI.
    
    :param accession_number: Номер доступу (NM_000142)
    :param email: Ваш email для NCBI
    :return: Послідовність у верхньому регістрі
    """
    Entrez.email = email
    handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    handle.close()
    return str(record.seq).upper()

def has_poly_repeats(sequence, max_repeat_length=3):
    """
    Перевіряє наявність полінуклеотидних повторів.
    
    :param sequence: Послідовність
    :param max_repeat_length: Максимальна довжина повторів
    :return: True, якщо є повтори, False — якщо немає
    """
    for i in range(len(sequence) - max_repeat_length):
        if sequence[i] == sequence[i + 1] == sequence[i + 2] == sequence[i + 3]:
            return True
    return False

def find_gRNA(dna_sequence, pam="TTTN", gRNA_length=20, gc_content_range=(40, 70), max_repeat_length=3):
    """
    Знаходить gRNA у послідовності для Cas12a.
    
    :param dna_sequence: Послідовність ДНК
    :param pam: PAM-послідовність (TTTN)
    :param gRNA_length: Довжина gRNA
    :param gc_content_range: Діапазон GC-вмісту
    :param max_repeat_length: Максимальна довжина повторів
    :return: Список gRNA із властивостями
    """
    pam_regex = pam.replace("N", "[ACGT]")
    potential_gRNAs = []
    for match in re.finditer(pam_regex, dna_sequence):
        pam_start = match.start()
        pam_end = match.end()
        if pam_end + gRNA_length <= len(dna_sequence):
            gRNA = dna_sequence[pam_end:pam_end + gRNA_length]
            pam_seq = dna_sequence[pam_start:pam_end]
            gc_count = gRNA.count("G") + gRNA.count("C")
            gc_content = (gc_count / gRNA_length) * 100
            if gc_content_range[0] <= gc_content <= gc_content_range[1]:
                if not has_poly_repeats(gRNA, max_repeat_length):
                    potential_gRNAs.append({
                        "gRNA": gRNA,
                        "PAM": pam_seq,
                        "Position": f"{pam_start + 1}-{pam_end + gRNA_length}",
                        "GC_Content": round(gc_content, 2)
                    })
    return potential_gRNAs

def select_top_gRNAs(potential_gRNAs, top_n=5):
    """
    Вибирає 5 найкращих gRNA за близькістю GC-вмісту до 50%.
    
    :param potential_gRNAs: Список gRNA
    :param top_n: Кількість gRNA для вибору
    :return: Список найкращих gRNA
    """
    for gRNA in potential_gRNAs:
        gRNA["GC_Distance"] = abs(gRNA["GC_Content"] - 50)
    sorted_gRNAs = sorted(potential_gRNAs, key=lambda x: x["GC_Distance"])
    return sorted_gRNAs[:top_n]

if __name__ == "__main__":
    accession_number = "NM_000142"
    email = "egorushkasuper12@gmail.com"  # Замініть на Ваш email
    try:
        dna_sequence = fetch_sequence_from_ncbi(accession_number, email)
        print(f"Послідовність FGFR3 ({accession_number}) завантажена. Довжина: {len(dna_sequence)} нуклеотидів.")
    except Exception as e:
        print(f"Помилка: {e}")
        dna_sequence = input("Введіть послідовність вручну (5' -> 3'): ").strip().upper()

    results = find_gRNA(dna_sequence)
    top_gRNAs = select_top_gRNAs(results)

    if top_gRNAs:
        df = pd.DataFrame(top_gRNAs)
        df = df[["gRNA", "PAM", "Position", "GC_Content"]]
        print("\n5 найкращих gRNA для FGFR3 (NM_000142):")
        print(df.to_string(index=False))
    else:
        print("gRNA не знайдено.")