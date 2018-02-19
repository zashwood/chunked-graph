using NamedTuples

import Base: isless, isequal, ==, <, hash, convert

struct AtomicEdge <: NamedTuples.NamedTuple
	u::Label
	v::Label
	affinity::Affinity
	AtomicEdge(u::Label, v::Label, affinity::Affinity) = u < v ? new(u, v, affinity) : new(v, u, affinity)
end

AtomicEdge(u, v) = AtomicEdge(Label(u), Label(v), Affinity(1))
AtomicEdge(u, v, affinity) = AtomicEdge(Label(u), Label(v), Affinity(affinity))

# Comparison
(==)(lhs::AtomicEdge, rhs::AtomicEdge) = isequal(head(lhs), head(rhs)) && isequal(tail(lhs), tail(rhs))
isequal(lhs::AtomicEdge, rhs::AtomicEdge) = lhs == rhs

#TODO: Make this definition independent of bit layout!
#We mainly want to compare on chunkid first
(<)(lhs::AtomicEdge, rhs::AtomicEdge) = isequal(head(lhs), head(rhs)) ? tail(lhs) < tail(rhs) : head(lhs) < head(rhs)
isless(lhs::AtomicEdge, rhs::AtomicEdge) = lhs < rhs

#we don't have heterogeneous sets, so we don't need to hash in the type name
#hash(e::AtomicEdge, seed::UInt) = hash(e.u, hash(e.v, hash(:AtomicEdge, seed)))
hash(e::AtomicEdge, seed::UInt) = hash(e.u, hash(e.v, seed))
