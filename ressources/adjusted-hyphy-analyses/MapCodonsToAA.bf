/**
 * @name alignments.alignment.adjusted_MapCodonsToAA
 * Map in-frame nucleotides onto a protein alignment string

 * @param {String} codon_sequence - the codon sequence to map
 * @param {String} aa_sequence - the matching aligned a.a. sequence
 * @param {Number} no more than this many mismatches - the codon sequence to map
 * @param {Dict} mapping - code ["terms.code.mapping"]

 * @returns {String} the mapped sequence

 * @example
    GCAAAATCATTAGGGACTATGGAAAACAGA
    -AKSLGTMEN-R

    maps to

    ---GCAAAATCATTAGGGACTATGGAAAAC---AGA

 */

lfunction alignment.adjusted_MapCodonsToAA(codon_sequence, aa_sequence, this_many_mm, mapping) {

    seqLen = Abs(aa_sequence);
    translString = "";
    translString * (seqLen);
    seqLenN = Abs(codon_sequence);

    aaPos = 0;
    seqPos = 0;
    codon = codon_sequence[seqPos][seqPos + 2];
    currentAA = mapping[codon];

    mismatch_count = 0;

    for (aaPos = 0; aaPos < seqLen && seqPos < seqLenN; aaPos += 1) {
        advance = 1;
        copy_codon = 1;

        if (aa_sequence[aaPos] == "-") {
            translString * "---";
            advance = 0;
        } else {
            if (currentAA != 0) {
                mismatch_count += (aa_sequence[aaPos] != currentAA);
                if (this_many_mm == 1) {
                    if (mismatch_count == this_many_mm) {
                        translString * 0;
                        console.log (translString);
                        console.log (codon_sequence);
                    }
                    assert(mismatch_count < this_many_mm, "A mismatch between codon and protein sequences at position " + aaPos + " (codon `seqPos`) : codon '" + codon_sequence[seqPos][seqPos + 2] + "'(`currentAA`) a.a. '`aa_sequence[aaPos]`'");
                } else {
                    if (mismatch_count >= this_many_mm) {
                        translString * 0;
                        return None;
                    }
                }
                if (currentAA == "X") {
                    translString * "---";
                } else {
                    translString * codon;
                }
            } else {
                if (currentAA == "X") {
                    translString * codon;
                }
                // TODO HIER EINE FEHLERMELDUNG EINFÜGEN !!!
            }
            seqPos += 3;
            codon = codon_sequence[seqPos][seqPos + 2];
            currentAA = mapping[codon];
        }
    }

    for (; aaPos < seqLen; aaPos += 1) {
        translString * "---";
    }

    translString * 0;
    return translString;
}