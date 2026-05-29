from dataclasses import dataclass

@dataclass
class Gene:
	strand: str  # '+' or '-'
	start: int   # 1-based inclusive
	end: int	 # 1-based inclusive

	@property
	def length_nt(self) -> int:
		return self.end - self.start + 1