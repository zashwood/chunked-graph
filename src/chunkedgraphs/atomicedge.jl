using NamedTuples

import Base: isless, isequal, ==, <, hash

struct AtomicEdge <: NamedTuples.NamedTuple
	u::Label
	v::Label
	affinity::Affinity
	AtomicEdge(u::Label, v::Label, affinity::Affinity) = u < v ? new(u, v, affinity) : new(v, u, affinity)
end
AtomicEdge(u, v) = AtomicEdge(Label(u), Label(v), Affinity(1))
AtomicEdge(u, v, affinity) = AtomicEdge(Label(u), Label(v), Affinity(affinity))

# Comparison
(==)(lhs::AtomicEdge, rhs::AtomicEdge) = isequal(lhs.u, rhs.u) && isequal(lhs.v, rhs.v)
isequal(lhs::AtomicEdge, rhs::AtomicEdge) = lhs == rhs

(<)(lhs::AtomicEdge, rhs::AtomicEdge) = isequal(lhs.u, rhs.u) ? lhs.v < rhs.v : lhs.u < rhs.u
isless(lhs::AtomicEdge, rhs::AtomicEdge) = lhs < rhs

hash(e::AtomicEdge, seed::UInt) = hash(e.u, hash(e.v, hash(:AtomicEdge, seed)))
