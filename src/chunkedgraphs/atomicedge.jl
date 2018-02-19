import Base: isless, isequal, ==, <, hash

abstract type AbstractEdge end

struct AtomicEdge <: AbstractEdge
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

# TODO: Make this definition independent of bit layout!
# We mainly want to compare on chunkid first
(<)(lhs::AtomicEdge, rhs::AtomicEdge) = isequal(lhs.u, rhs.u) ? lhs.v < rhs.v : lhs.u < rhs.u
isless(lhs::AtomicEdge, rhs::AtomicEdge) = lhs < rhs

# We don't have heterogeneous sets, so we don't need to hash in the type name
# hash(e::AtomicEdge, seed::UInt) = hash(e.u, hash(e.v, hash(:AtomicEdge, seed)))
hash(e::AtomicEdge, seed::UInt) = hash(e.u, hash(e.v, seed))

head(e::AbstractEdge) = e.u
tail(e::AbstractEdge) = e.v
