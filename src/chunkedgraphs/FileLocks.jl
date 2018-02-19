__precompile__(true)
module FileLocks

import Base:lock, unlock, close
export FileLock

mutable struct FileLock
	f::IOStream
	s::String
	locked::Bool
	function FileLock(s)
		f = open(s, "w")
		return new(f, s, false)
	end
end

const LOCK_SH = 0x01		# shared file lock
const LOCK_EX = 0x02		# exclusive file lock
const LOCK_NB = 0x04		# don't block when locking
const LOCK_UN = 0x08		# unlock file

function lock(fl::FileLock)
	if !fl.locked
		x = ccall(:flock, Int32, (Int32, Int32), fd(fl.f), LOCK_EX)
		@assert x == 0
		info("acquired lock $(fl.s)")
		fl.locked = true
	end
end

function unlock(fl::FileLock)
	x = ccall(:flock, Int32, (Int32, Int32), fd(fl.f), LOCK_UN)
	@assert x == 0
	info("released lock $(fl.s)")
	fl.locked = false
end

function close(fl::FileLock)
	unlock(fl)
	close(fl.f)
end

end
