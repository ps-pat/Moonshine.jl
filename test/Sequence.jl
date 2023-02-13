sequence_qc = Quickcheck("Tests Sequence Properties", serialize_fails = false)

const Seq = Sequence{UInt}

## Equality.
@add_predicate(sequence_qc,
               "Equality Reflexive",
               seq::Seq -> seq == seq)

## String
@add_predicate(sequence_qc,
               "To String and Back",
               seq::Seq -> convert(Seq, string(seq)) == seq)

## Bitwise operations.
@add_predicate(sequence_qc,
               "Negation Involutive",
               seq::Seq -> ~~seq == seq)

# @add_predicate(sequence_qc,
#                "Xor commutativity",
#                (seq1::Seq, seq2::Seq) -> seq1 ⊻ seq2 == seq2 ⊻ seq1)

# @add_predicate(sequence_qc,
#                "Xor associativity",
#                (seq1::Seq, seq2::Seq, seq3::Seq) -> seq1 ⊻ (seq2 ⊻ seq3) == (seq1 ⊻ seq2) ⊻ seq3)

@add_predicate(sequence_qc,
               "Xor identity element",
               seq::Seq -> seq ⊻ fillseq(false, length(seq)) == seq)

@add_predicate(sequence_qc,
               "Xor 2-idempotency",
               seq::Seq -> seq ⊻ seq == fillseq(false, length(seq)))

@add_predicate(sequence_qc,
               "Xor bit switch",
               seq::Seq -> seq ⊻ fillseq(true, length(seq)) == ~seq)

@add_predicate(sequence_qc,
               "& identity element",
               seq::Seq -> seq & fillseq(true, length(seq)) == seq)

@add_predicate(sequence_qc,
               "& absorbing element",
               seq::Seq -> seq & fillseq(false, length(seq)) ==
                   fillseq(false, length(seq)))

@add_predicate(sequence_qc,
               "| identity element",
               seq::Seq -> seq | fillseq(false, length(seq)) == seq)

@add_predicate(sequence_qc,
               "| absorbing element",
               seq::Seq -> seq | fillseq(true, length(seq)) ==
                   fillseq(true, length(seq)))

@testset "Sequence" begin
    @quickcheck sequence_qc
end
